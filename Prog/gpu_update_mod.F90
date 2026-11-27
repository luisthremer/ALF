!  Copyright (C) 2016 - 2025 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!--------------------------------------------------------------------
!> @brief
!> GPU port (P2.4 step 1c): WRAPGRUP / WRAPGRDO on a device-resident Green function.  The per-site
!> wraps and the rank-1 update run on the device, the small ratio block and the accept/reject decision on the host, with
!> the random numbers drawn in exactly ALF's order, so that for the same compiler and seeds the trajectory is the CPU
!> one until round-off first changes a decision.  Compiled only with -DALF_GPU; entered from Wrapgr_mod through guarded
!> early branches.
!>
!> Per site and flavour (nf):
!>   1. wrap: Op_Wrapup/Op_Wrapdo (N_Type 1) of operator n, one fused kernel launch for N = 1 operators (rows and columns
!>      of the one site together, rows before columns for the corner element as in ALF), else the 204 two-kernel routine.
!>      Non-diagonal / type-4 operators (P2.6.1): dense-block kernels of gpu_wrap_mod (nd_wrap), N_Type 1 and, after the update
!>      (WRAPGRUP) / before it (WRAPGRDO), N_Type 2 as in Wrapgr_mod; N_Type 2 does nothing for a diagonal operator.
!>   2. ratio (host, the code of Upgrade2): needs G(P(a),P(b),nf), a, b <= r = N_non_zero: one complex element D->H per
!>      flavour for r = 1, N = 1, else the r x r block (one gather kernel, one D->H).  Op_V%alpha and exp(Z1*E) are host data.
!>   3. decision (host): ranf_wrap() drawn at the same point as Upgrade2 (after the T0 test of WRAPGR).
!>   4. accept: per flavour 2 launches (copy the column G(:,P) and the row y = (-G(P,:), 1-G(P,P) at P) into scratch, then
!>      G += Zg * x y^T, the ZGERU of Upgrade2 for op_dim 1), then the field value on the host and one scalar H->D copy
!>      of nsigma_d(n,t) (the wrap kernels read the field on the device, so it must not go stale).
!>   Rank r > 1 or N > 1 (P2.6.1, port of gpu-port ticket 401-A): G <- G - (G(:,P) X) (Tri^T R), R = E_P^T - G(P,:), with
!>   the r x r factors X = x_v(P,:) and Tri (y_v = R^T Tri) of Upgrade2's recursion (host, rk_xtri), one H->D, two kernels.
!>   Upgrade2's grouping, not the Woodbury form G(:,P) W R with W = X Tri^T: that form loses ~30x in Precision Green for the
!>   type-4 J vertices of Spin_Peierls 12x12 (KNOWN_ISSUES G10; the CPU delayed/submatrix updates use it and show the same).
!>
!> G is device-authoritative while the sweep runs: WRAPGRUP/DO copy the host GR to G_d only when the device copy is not current
!> (g_in) and leave the host GR stale (g_out).  The stabilisation (gpu_stab_mod, step 1d) keeps G on the device; the sweep
!> (qmc_sweep_mod) fetches it (gpu_g_host) before each host reader of GR, and host code that writes GR calls gpu_g_host_wrote.
!>
!> Environment (debug runs): ALF_GPU_NOUPDATE=1 keeps the host WRAPGRUP/DO in a GPU build; ALF_GPU_FAULT=upd|nsig injects a
!> fault into the device update (negative controls).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_update_mod

  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use Fields_mod, only: fields
  Use Hamiltonian_main, only: calc_fl_map, ham, n_fl, n_fl_eff, n_sun, ndim, nsigma, op_v, reconstruction_needed
  Use Control, only: Control_upgrade, Control_upgrade_eff
  Use Random_Wrap, only: ranf_wrap
  Use gpu_state_mod, only: g_d, g_dev_ok, g_host_ok, gpu_active, gpu_env_is, info_line, nsigma_d, scr_vec, &
       & upd_cnt, v_exp, v_meta, v_type, v_nmax
  Use MyMats, only: Det
#if defined(ALF_LA_HIP)
  Use alf_hipk_mod, only: hk_enabled, hk_rk_gpp, hk_rk_accept
#endif
  Use gpu_wrap_mod, only: gpu_mmthr, gpu_mmthl_m1, gpu_mmthl, gpu_mmthr_m1, gpu_Vop_wrap
  Use gpu_lowrank_mod, only: gpu_lr_accept, gpu_lr_blk_enter, gpu_lr_blk_leave, gpu_lr_blockmode, gpu_lr_commit, &
       & gpu_lr_defer, gpu_lr_fast, gpu_lr_flush, gpu_lr_flush_all, gpu_lr_gpp, gpu_lr_mode, gpu_lr_site_ok, &
       & gpu_lr_step_down, gpu_lr_step_up, gpu_lr_wrap, gpu_lr_gen_ok, gpu_lr_gen_wrap, gpu_lr_gen_block, gpu_lr_gen_accept

  Implicit none

  Private
  Public :: gpu_upd_enabled, gpu_g_host, gpu_g_host_wrote, gpu_wrapgrup, gpu_wrapgrdo

  Logical, Save :: checked = .false., enabled = .false.
  Logical, Save :: fault_upd = .false., fault_nsig = .false.     ! negative controls ALF_GPU_FAULT=upd|nsig (debug runs only)
  Type (Fields), Save :: nsigma_new
  Logical, Save :: nsigma_new_made = .false.
  ! P2.6.1: device scratch of the rank-r update (block of G, W, sites, C = G(:,P), T = W R), allocated at the first use
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: rk_g(:), rk_w(:), rk_c(:), rk_t(:)
  Integer, Allocatable, Save :: rk_p(:)
  ! P2.6 global: .true. while G is moved by gpu_placegr or a global tau move runs: the wraps and the update take the standard
  ! device path (delayed / submatrix state flushed first, as the CPU calls Lowrank_flush before the global tau moves)
  Logical, Save :: glob_std = .false.

Contains

!--------------------------------------------------------------------
!> .true. iff the device update path replaces WRAPGRUP/WRAPGRDO: GPU path active and ALF_GPU_NOUPDATE /= 1.
  Logical function gpu_upd_enabled()
    Implicit none
    if (.not. checked) then
       enabled = gpu_active
       if (enabled) then
          if (gpu_env_is('ALF_GPU_NOUPDATE', '1')) then
             enabled = .false.
             call info_line('GPU update path: DISABLED by ALF_GPU_NOUPDATE=1, host WRAPGRUP/DO')
          else
             call info_line('GPU update path: ACTIVE (device wraps and local updates, host decision; G stays on the device)')
             fault_upd  = gpu_env_is('ALF_GPU_FAULT', 'upd')
             fault_nsig = gpu_env_is('ALF_GPU_FAULT', 'nsig')
          endif
       endif
       checked = .true.
    endif
    gpu_upd_enabled = enabled
  end function gpu_upd_enabled

!====================================================================
!  Interim copies of G (host GR <-> G_d)
!====================================================================
  !> G is device-authoritative during the sweep.  g_in copies the host GR to G_d only for a flavour whose device
  !> copy is not current (start of the run, after a host write); g_out leaves G_d current and marks the host GR stale.  A host reader
  !> calls gpu_g_host(GR) first; code that writes the host GR calls gpu_g_host_wrote().
  Subroutine g_in(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(:,:,:)
    Integer :: nf_eff, nf
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       if (g_dev_ok(nf)) cycle
       G_d(:,:,nf) = GR(:,:,nf)
       !$omp target update to(G_d(:,:,nf))
       g_dev_ok(nf) = .true.
       upd_cnt(7) = upd_cnt(7) + 1
       upd_cnt(9) = upd_cnt(9) + 16_int64*Ndim*Ndim
    enddo
  end subroutine g_in

  Subroutine g_out(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:)
    Integer :: nf_eff, nf
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       g_host_ok(nf) = .false.
    enddo
  end subroutine g_out

  !> D->H of G for every flavour whose host copy is stale (before a host reader of GR)
  Subroutine gpu_g_host(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:)
    Integer :: nf_eff, nf
    if (.not. allocated(g_host_ok)) return
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       if (g_host_ok(nf)) cycle
       !$omp target update from(G_d(:,:,nf))
       GR(:,:,nf) = G_d(:,:,nf)
       g_host_ok(nf) = .true.
       upd_cnt(8) = upd_cnt(8) + 1
       upd_cnt(9) = upd_cnt(9) + 16_int64*Ndim*Ndim
    enddo
  end subroutine gpu_g_host

  !> The host GR was written (host-only code): the device copy is stale
  Subroutine gpu_g_host_wrote()
    Implicit none
    if (.not. allocated(g_dev_ok)) return
    g_dev_ok = .false.
    g_host_ok = .true.
  end subroutine gpu_g_host_wrote

!====================================================================
!  Device kernels
!====================================================================
!> Op_Wrapup (dir = +1) / Op_Wrapdo (dir = -1), N_Type 1, of the single-site operator n at slice t on G, one launch: the
!> row r and the column r of G; element (r,r) gets the row factor, then the column factor (ALF's ZSCAL order).
  Subroutine k_site_wrap(G, ld, nd, n, t, dir, r, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nd, n, t, dir, r, eo, nt
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,nd)
    Integer :: j, sp
    Complex (Kind=Kind(0.d0)) :: za, zb
    !$omp target teams distribute parallel do private(sp, za, zb)
    do j = 1, nd
       sp = nint(real(nsigma_d(n,t)))
       za = V_exp(eo + 1 + (dir*sp + nt))
       zb = V_exp(eo + 1 + (-dir*sp + nt))
       if (j == r) then
          G(r,r) = zb * ( za * G(r,r) )
       else
          G(r,j) = za * G(r,j)
          G(j,r) = zb * G(j,r)
       endif
    enddo
  end subroutine k_site_wrap

!> x = G(:,p), y = -G(p,:) with y(p) = 1 - G(p,p)  (u, v, x_v, y_v of Upgrade2 for op_dim 1)
  Subroutine k_upd_prep(G, ld, nd, p, x, y)
    Implicit none
    Integer, intent(in) :: ld, nd, p
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,nd)
    Complex (Kind=Kind(0.d0)), intent(inout) :: x(nd), y(nd)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nd
       x(i) = G(i,p)
       y(i) = - G(p,i)
       if (i == p) y(i) = 1.d0 - G(p,p)
    enddo
  end subroutine k_upd_prep

!> G := G + Zg * x y^T (ZGERU: temp = alpha*y(j); G(i,j) = G(i,j) + x(i)*temp), Zg = (zr, zi).
!> The complex scalar is passed as two reals: with nvfortran 26.9 a complex scalar dummy (intent(in)) used in the target
!> region lost its imaginary part (measured: complex G, real-G models were not affected).
  Subroutine k_upd_geru(G, ld, nd, zr, zi, x, y)
    Implicit none
    Integer, intent(in) :: ld, nd
    Real (Kind=Kind(0.d0)), intent(in) :: zr, zi
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,nd)
    Complex (Kind=Kind(0.d0)), intent(in) :: x(nd), y(nd)
    Integer :: i, j
    Complex (Kind=Kind(0.d0)) :: tmp, Zg
    !$omp target teams distribute parallel do collapse(2) private(tmp, Zg)
    do j = 1, nd
       do i = 1, nd
          Zg = cmplx(zr, zi, kind(0.d0))
          tmp = Zg * y(j)
          G(i,j) = G(i,j) + x(i) * tmp
       enddo
    enddo
  end subroutine k_upd_geru

!> P2.6.1: buf(a + (b-1)*r) = G(P(a), P(b)), a, b = 1..r
  Subroutine k_rk_gpp(G, ld, P, r, buf)
    Implicit none
    Integer, intent(in) :: ld, r
    Integer, intent(in) :: P(V_NMAX)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: buf(V_NMAX*V_NMAX)
    Integer :: a, b
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_rk_gpp(G, ld, P, r, buf)
       return
    endif
#endif
    !$omp target teams distribute parallel do collapse(2)
    do b = 1, r
       do a = 1, r
          buf(a + (b-1)*r) = G(P(a), P(b))
       enddo
    enddo
  end subroutine k_rk_gpp

!> P2.6.1: C(i,a) = sum_c G(i,P(c)) X(c,a) (Upgrade2's xp_v) and T(a,i) = sum_b Tri(b,a) R(b,i) (y_v(i,a)), R(b,i) =
!> delta(i,P(b)) - G(P(b),i), i = 1..nd (before G changes).  XT(1:r*r) = X, XT(r*r+1:2*r*r) = Tri, column-major r x r.
  Subroutine k_rk_prep(G, ld, nd, P, r, XT, C, T)
    Implicit none
    Integer, intent(in) :: ld, nd, r
    Integer, intent(in) :: P(V_NMAX)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,nd), XT(2*V_NMAX*V_NMAX)
    Complex (Kind=Kind(0.d0)), intent(inout) :: C(nd*V_NMAX), T(V_NMAX*nd)
    Integer :: i, a, b
    Complex (Kind=Kind(0.d0)) :: acc, rb
    !$omp target teams distribute parallel do private(a, b, acc, rb)
    do i = 1, nd
       do a = 1, r
          acc = G(i, P(1)) * XT(1 + (a-1)*r)
          do b = 2, r
             acc = acc + G(i, P(b)) * XT(b + (a-1)*r)
          enddo
          C(i + (a-1)*nd) = acc
       enddo
       do a = 1, r
          acc = (0.d0, 0.d0)
          do b = 1, r
             rb = - G(P(b), i)
             if (i == P(b)) rb = 1.d0 - G(P(b), i)
             acc = acc + XT(r*r + b + (a-1)*r) * rb
          enddo
          T(a + (i-1)*r) = acc
       enddo
    enddo
  end subroutine k_rk_prep

!> P2.6.1: G(i,j) <- G(i,j) - sum_a C(i,a) T(a,j)
  Subroutine k_rk_upd(G, ld, nd, r, C, T)
    Implicit none
    Integer, intent(in) :: ld, nd, r
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,nd)
    Complex (Kind=Kind(0.d0)), intent(in) :: C(nd*V_NMAX), T(V_NMAX*nd)
    Integer :: i, j, a
    Complex (Kind=Kind(0.d0)) :: acc
    !$omp target teams distribute parallel do collapse(2) private(a, acc)
    do j = 1, nd
       do i = 1, nd
          acc = C(i) * T(1 + (j-1)*r)
          do a = 2, r
             acc = acc + C(i + (a-1)*nd) * T(a + (j-1)*r)
          enddo
          G(i,j) = G(i,j) - acc
       enddo
    enddo
  end subroutine k_rk_upd

!====================================================================
!  Host-side orchestration
!====================================================================
!> P2.6.1: device scratch of the rank-r update, allocated once (no transfer)
  Subroutine rk_alloc()
    Implicit none
    if (allocated(rk_g)) return
    allocate(rk_g(V_NMAX*V_NMAX), rk_w(2*V_NMAX*V_NMAX), rk_c(Ndim*V_NMAX), rk_t(V_NMAX*Ndim), rk_p(V_NMAX))
    rk_g = (0.d0,0.d0); rk_w = (0.d0,0.d0); rk_c = (0.d0,0.d0); rk_t = (0.d0,0.d0); rk_p = 1
    !$omp target enter data map(to: rk_g, rk_w, rk_c, rk_t, rk_p)
  end subroutine rk_alloc

!> P2.6.1: G(P(a),P(b)) of flavour nf, a, b = 1..r, to the host (gb, column-major r x r): one kernel, one D->H of r*r elements
  Subroutine rk_block(nf, n_op, r, gb)
    Implicit none
    Integer, intent(in) :: nf, n_op, r
    Complex (Kind=Kind(0.d0)), intent(out) :: gb(r,r)
    call rk_alloc()
    rk_p(1:r) = Op_V(n_op,nf)%P(1:r)
    !$omp target update to(rk_p(1:r))
    call k_rk_gpp(G_d(1,1,nf), Ndim, rk_p, r, rk_g)
    !$omp target update from(rk_g(1:r*r))
    gb = reshape(rk_g(1:r*r), (/r, r/))
    upd_cnt(5) = upd_cnt(5) + 1
  end subroutine rk_block

!> P2.6.1: the r x r factors of the accept G <- G - (G(:,P) X) (Tri^T R) (R(k,j) = delta(j,P(k)) - G(P(k),j)): Upgrade2's
!> recursion for x_v, y_v with y_v = V Tri (V(:,k) = v(:,k) = R(k,:)^T) restricted to the rows P (x_v vanishes outside P),
!> X = x_v(P,:).  The delayed general path uses W = X Tri^T.  Inputs: gb(a,b) = G(P(a),P(b)), Delta(1:r).
  Subroutine rk_xtri(r, gb, Delta, X, Tri)
    Implicit none
    Integer, intent(in) :: r
    Complex (Kind=Kind(0.d0)), intent(in) :: gb(r,r), Delta(r)
    Complex (Kind=Kind(0.d0)), intent(out) :: X(r,r), Tri(r,r)
    Complex (Kind=Kind(0.d0)) :: Vp(r,r), syu(r), sxv(r), Z
    Integer :: n, m, k
    ! Vp(a,k) = v(P(a),k) = delta(a,k) - G(P(k),P(a))
    do k = 1, r
       do n = 1, r
          Vp(n,k) = - gb(k,n)
       enddo
       Vp(k,k) = 1.d0 - gb(k,k)
    enddo
    X = (0.d0,0.d0); Tri = (0.d0,0.d0)
    X(1,1) = Delta(1)/(1.d0 + Vp(1,1)*Delta(1))
    Tri(1,1) = 1.d0
    do n = 2, r
       X(n,n) = Delta(n)
       Tri(n,n) = 1.d0
       Z = 1.d0 + Delta(n)*Vp(n,n)
       do m = 1, n-1
          ! syu(m) = - y_v(P(n),m) u(P(n),n),  y_v(P(n),m) = sum_k Vp(n,k) Tri(k,m);  sxv(m) = - sum_k x_v(P(k),m) v(P(k),n)
          syu(m) = - sum(Vp(n,1:n-1)*Tri(1:n-1,m)) * Delta(n)
          sxv(m) = - sum(X(1:n-1,m)*Vp(1:n-1,n))
       enddo
       do m = 1, n-1
          X(:,n) = X(:,n) + X(:,m)*syu(m)
          Tri(:,n) = Tri(:,n) + Tri(:,m)*sxv(m)
       enddo
       do m = 1, n-1
          Z = Z - syu(m)*sxv(m)
       enddo
       Z = 1.d0/Z                         ! as Upgrade2: reciprocal, then scale
       X(:,n) = X(:,n) * Z
    enddo
  end subroutine rk_xtri

!> P2.6.1: Op_Wrapup (dir = +1, after the update) / Op_Wrapdo (dir = -1, before it) with N_Type 2 of operator n: nothing for a
!> diagonal operator; otherwise G must be flushed (delayed / submatrix state) first.
  Subroutine site_wrap2(nf, n, t, dir)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    if (V_meta(6,n,nf) == 1) return
    if (gpu_lr_gen_ok(n) .and. .not. glob_std) then
       call gpu_lr_gen_wrap(nf, n, t, dir, 2)    ! delayed general vertex: pending X, Y follow the wrap
       return
    endif
    if (gpu_lr_mode() /= 0) call gpu_lr_flush(nf, 2)
    call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, t, dir, 2)
    upd_cnt(3) = upd_cnt(3) + 2
  end subroutine site_wrap2

  Subroutine site_wrap(nf, n, t, dir)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    if (gpu_lr_mode() /= 0 .and. .not. glob_std) then
       if (V_meta(1,n,nf) == 1) then
          call gpu_lr_wrap(nf, n, t, dir)       ! pending delayed / submatrix state follows the wrap
          return
       endif
       if (gpu_lr_gen_ok(n)) then
          call gpu_lr_gen_wrap(nf, n, t, dir, 1)   ! P2.6.1: delayed general vertex
          return
       endif
       call gpu_lr_flush(nf, 2)                 ! N > 1 vertex: the general wrap needs a flushed G
    endif
    if (V_meta(1,n,nf) == 1) then
       call k_site_wrap(G_d(1,1,nf), Ndim, Ndim, n, t, dir, Op_V(n,nf)%P(1), V_meta(3,n,nf), V_type(n))
       upd_cnt(3) = upd_cnt(3) + 1
    else
       call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, t, dir, 1)
       upd_cnt(3) = upd_cnt(3) + 2
    endif
  end subroutine site_wrap

!--------------------------------------------------------------------
!> Wrapgr_PlaceGR on the device: G(tau, m) -> G(tau, m1), the same wrap order (upward: N_Type 1 of all flavours, then
!> N_Type 2; downward: N_Type 2, then N_Type 1).  The delayed / submatrix state is flushed first; standard wraps only.
  Subroutine gpu_placegr(m, m1, nt)
    Implicit none
    Integer, intent(in) :: m, m1, nt
    Integer :: n, nf, nf_eff
    Logical :: std0
    if (m == m1) return
    if (.not. glob_std) call gpu_lr_flush_all(2)
    std0 = glob_std
    glob_std = .true.
    if (m1 > m) then
       do n = m+1, m1
          do nf_eff = 1, N_FL_eff
             nf = Calc_Fl_map(nf_eff)
             call site_wrap(nf, n, nt, 1)
          enddo
          do nf_eff = 1, N_FL_eff
             nf = Calc_Fl_map(nf_eff)
             call site_wrap2(nf, n, nt, 1)
          enddo
       enddo
    else
       do n = m, m1+1, -1
          do nf_eff = 1, N_FL_eff
             nf = Calc_Fl_map(nf_eff)
             call site_wrap2(nf, n, nt, -1)
          enddo
          do nf_eff = 1, N_FL_eff
             nf = Calc_Fl_map(nf_eff)
             call site_wrap(nf, n, nt, -1)
          enddo
       enddo
    endif
    glob_std = std0
  end subroutine gpu_placegr

!--------------------------------------------------------------------
!> Upgrade2 (mode FINAL, Prev_Ratiotot = 1) for the device-resident G: same ratio code, the decision draw at the same
!> point, the update on the device.  Op_V%N_non_zero <= 1 for every operator (gpu_supported).
  Subroutine site_update(N_op, NT, PHASE, Hs_new, S0_ratio, T0_proposal_ratio, toggle)
    Implicit none
    Integer, intent(in) :: N_op, NT
    Complex (Kind=Kind(0.d0)), intent(inout) :: PHASE
    Complex (Kind=Kind(0.d0)), intent(in) :: Hs_new
    Real (Kind=Kind(0.d0)), intent(in) :: S0_ratio, T0_proposal_ratio
    Logical, intent(inout) :: toggle

    Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Ratiotot, Z1, Prev_Ratiotot, Delta(N_FL), gpp(N_FL)
    Complex (Kind=Kind(0.d0)) :: Z, D_Mat, myexp, g_loc, Mat11, vii, xv, Zg, s1, s2
    Complex (Kind=Kind(0.d0)) :: gb(V_NMAX,V_NMAX,N_FL), Dl(V_NMAX,N_FL), Mat(V_NMAX,V_NMAX), W(V_NMAX*V_NMAX)
    Complex (Kind=Kind(0.d0)) :: Xm(V_NMAX,V_NMAX), Tm(V_NMAX,V_NMAX)
    Real (Kind=Kind(0.d0)) :: Weight
    Integer :: nf, nf_eff, op_dim_nf, p1, m, k, r2
    Logical :: lr, lrf, lrg, gen(N_FL)

    toggle = .false.
    lr = ( gpu_lr_mode() /= 0 )
    if (lr) lr = gpu_lr_site_ok(N_op)
    lrf = lr .and. gpu_lr_fast()
    lrg = (.not. lr) .and. gpu_lr_gen_ok(N_op)    ! P2.6.1: delayed general vertex (G = G0 - X Y^T kept)
    if (.not. (lr .or. lrg)) call gpu_lr_flush_all(2)   ! standard update of this site: G_d must be current
    if (.not. nsigma_new_made) then
       Call nsigma_new%make(1,1)
       nsigma_new_made = .true.
    endif
    Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
    Ratio = cmplx(1.d0,0.d0,kind(0.d0)); Delta = cmplx(0.d0,0.d0,kind(0.d0)); gpp = Delta

    ! Compute the ratio (Upgrade2)
    nf = 1
    nsigma_new%f(1,1)  = Hs_new
    nsigma_new%t(1)    = Op_V(n_op,nf)%Type
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       g_loc = Op_V(n_op,nf)%g
       if (op_v(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
       Z1 = g_loc * ( nsigma_new%Phi(1,1) -  nsigma%Phi(n_op,nt) )
       op_dim_nf = Op_V(n_op,nf)%N_non_zero
       ! P2.6.1: general vertex (rank > 1 or N > 1): the r x r block, Upgrade2's Mat and determinant branches
       gen(nf) = (op_dim_nf > 1 .or. (op_dim_nf == 1 .and. Op_V(n_op,nf)%N > 1))
       if (gen(nf)) then
          if (lrg) then
             call gpu_lr_gen_block(nf, n_op, op_dim_nf, gb(1:op_dim_nf,1:op_dim_nf,nf))
          else
             call rk_block(nf, n_op, op_dim_nf, gb(1:op_dim_nf,1:op_dim_nf,nf))
          endif
          do m = 1, op_dim_nf
             myexp = exp( Z1* Op_V(n_op,nf)%E(m) )
             Z = myexp - 1.d0
             Dl(m,nf) = Z
             do k = 1, op_dim_nf
                Mat(k,m) = - Z * gb(k,m,nf)
             enddo
             Mat(m,m) = myexp + Mat(m,m)
          enddo
          if (op_dim_nf == 1) then
             D_mat = Mat(1,1)
          elseif (op_dim_nf == 2) then
             s1 = Mat(1,1)*Mat(2,2)
             s2 = Mat(2,1)*Mat(1,2)
             If (Abs(s1) > Abs(s2)) then
                D_mat = s1*(1.D0 - s2/s1)
             else
                D_mat = s2*(s1/s2 - 1.D0)
             Endif
          else
             D_mat = Det(Mat(1:op_dim_nf,1:op_dim_nf), op_dim_nf)
          endif
       elseif (op_dim_nf == 0) then
          D_mat = 1.0d0
       else
          p1 = Op_V(n_op,nf)%P(1)
          ! The ratio needs one element of G on the host, G(P(1),P(1),nf): 16 bytes D->H per flavour and site (unavoidable
          ! while the decision is on the host, ADR-001 section 5).
          if (lr) then
             ! delayed / submatrix (gpu_lowrank_mod): G(p,p) of G0 + pending update, one device region, one D->H
             call gpu_lr_gpp(nf, nf_eff, p1, gpp(nf_eff))
          else
             !$omp target update from(G_d(p1,p1,nf))
             upd_cnt(5) = upd_cnt(5) + 1
             gpp(nf_eff) = G_d(p1,p1,nf)
          endif
          myexp = exp( Z1* Op_V(n_op,nf)%E(1) )
          Z = myexp - 1.d0
          Delta(nf_eff) = Z
          Mat11 = - Z * gpp(nf_eff)
          Mat11 = myexp + Mat11
          D_mat = Mat11
       endif
       Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
    Enddo

    ! Fill the non-calculated flavours, as Upgrade2
    if (reconstruction_needed) call ham%weight_reconstruction(Ratio)

    Ratiotot = Product(Ratio)
    nf = 1
    Ratiotot = (Ratiotot**dble(N_SUN)) * nsigma_new%Gama(1,1)/nsigma%Gama(n_op,nt)
    Ratiotot = Ratiotot * Prev_Ratiotot
    Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )

    upd_cnt(1) = upd_cnt(1) + 1
    if ( Weight > ranf_wrap() )  Then
       toggle = .true.
       upd_cnt(2) = upd_cnt(2) + 1
       Phase = Phase * Ratiotot/sqrt(Ratiotot*conjg(Ratiotot))
       Do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          if (gen(nf)) then
             op_dim_nf = Op_V(n_op,nf)%N_non_zero
             r2 = op_dim_nf**2
             call rk_xtri(op_dim_nf, gb(1:op_dim_nf,1:op_dim_nf,nf), Dl(1:op_dim_nf,nf), Xm(1:op_dim_nf,1:op_dim_nf), &
                  & Tm(1:op_dim_nf,1:op_dim_nf))
             ! W = X Tri^T: the HIP accept (GPU=amd, not yet in Upgrade2's grouping) take W
             W(1:r2) = reshape(matmul(Xm(1:op_dim_nf,1:op_dim_nf), transpose(Tm(1:op_dim_nf,1:op_dim_nf))), (/r2/))
             if (lrg) then
                call gpu_lr_gen_accept(nf, n_op, op_dim_nf, reshape(Xm(1:op_dim_nf,1:op_dim_nf), (/r2/)), &
                     & reshape(Tm(1:op_dim_nf,1:op_dim_nf), (/r2/)))
                cycle
             endif
             rk_p(1:op_dim_nf) = Op_V(n_op,nf)%P(1:op_dim_nf)
#if defined(ALF_LA_HIP)
             if (hk_enabled()) then
                rk_w(1:r2) = W(1:r2)
                !$omp target update to(rk_w(1:r2), rk_p(1:op_dim_nf))
                call hk_rk_accept(G_d(1,1,nf), Ndim, Ndim, rk_p, op_dim_nf, rk_w, rk_c, rk_t)
             else
#endif
             rk_w(1:r2) = reshape(Xm(1:op_dim_nf,1:op_dim_nf), (/r2/))
             rk_w(r2+1:2*r2) = reshape(Tm(1:op_dim_nf,1:op_dim_nf), (/r2/))
             !$omp target update to(rk_w(1:2*r2), rk_p(1:op_dim_nf))
             call k_rk_prep(G_d(1,1,nf), Ndim, Ndim, rk_p, op_dim_nf, rk_w, rk_c, rk_t)
             call k_rk_upd(G_d(1,1,nf), Ndim, Ndim, op_dim_nf, rk_c, rk_t)
#if defined(ALF_LA_HIP)
             endif
#endif
             upd_cnt(4) = upd_cnt(4) + 2
          elseif (Op_V(n_op,nf)%N_non_zero > 0) then
             p1 = Op_V(n_op,nf)%P(1)
             ! u = Delta, v(P) = 1 - G(P,P): x_v = u/(1 + v u), ZGERU alpha = -x_v(P)
             vii = 1.d0 - gpp(nf_eff)
             xv  = Delta(nf_eff)/(1.d0 + vii*Delta(nf_eff))
             Zg  = - xv
             if (lr) then
                call gpu_lr_accept(nf, nf_eff, p1, xv)
             else
                if (fault_upd) Zg = Zg * (1.d0 + 1.d-6)        ! negative control (ALF_GPU_FAULT=upd)
                call k_upd_prep(G_d(1,1,nf), Ndim, Ndim, p1, scr_vec(1,1), scr_vec(1,2))
                call k_upd_geru(G_d(1,1,nf), Ndim, Ndim, real(Zg), aimag(Zg), scr_vec(1,1), scr_vec(1,2))
                upd_cnt(4) = upd_cnt(4) + 2
             endif
          endif
       enddo
       if (lrf) call gpu_lr_commit()           ! delayed fast path: both flavours' accepts in one launch
       ! Flip the spin: host field, and the device copy that the wrap kernels read
       nsigma%f(n_op,nt) = nsigma_new%f(1,1)
       nsigma_d(n_op,nt) = nsigma%f(n_op,nt)
       if (lrf) then
          continue                              ! the device copy is written by the next launch (gpu_lr_defer)
       elseif (.not. fault_nsig) then
          ! one scalar H->D: the device copy of the field must follow the accepted flip (no stale nsigma_d reads)
          !$omp target update to(nsigma_d(n_op,nt))
          upd_cnt(6) = upd_cnt(6) + 1
       endif
    endif

    Call Control_upgrade(toggle)
    Call Control_upgrade_eff(toggle)
  End Subroutine site_update

!--------------------------------------------------------------------
!> WRAPGRUP of Wrapgr_mod on the device: G(tau) -> G(tau+1) with the sequential updates of slice tau+1.
!> No global moves (refused by gpu_supported).
  Subroutine gpu_wrapgrup(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:)
    Complex (Kind=Kind(0.d0)), intent(inout) :: PHASE
    Integer, intent(in) :: NTAU
    Logical, intent(in) :: Propose_S0
    Integer, intent(in) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

    Integer :: nf, nf_eff, NTAU1, n, m
    Complex (Kind=Kind(0.d0)) :: HS_New
    Real (Kind=Kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
    Logical :: Acc, toggle1, fast_lr, bm

    if (N_Global_tau > 0) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    upd_cnt(10) = upd_cnt(10) + 1
    fast_lr = gpu_lr_fast()
    bm = gpu_lr_blockmode()
    NTAU1 = NTAU + 1
    call g_in(GR)
    ! Wrap up with B(tau1)
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       call gpu_mmthr   (G_d(1,1,nf), Ndim, Ndim, nf, ntau1)
       call gpu_mmthl_m1(G_d(1,1,nf), Ndim, Ndim, nf, ntau1)
    Enddo
    ! G is at m = 0; the sequential range starts at m = Nt_sequential_start - 1 (KNOWN_ISSUES.md, H5)
    if (Nt_sequential_end >= Nt_sequential_start) call gpu_placegr(0, Nt_sequential_start-1, ntau1)
    Do n = Nt_sequential_start,Nt_sequential_end
       if (bm) then
          call gpu_lr_blk_enter(n, Nt_sequential_end, 1, ntau1, .true.)   ! submatrix block mode: host block, wrap on the host copy
       elseif (fast_lr) then
          call gpu_lr_step_up(n, ntau1)             ! wrap + gather / G(p,p) of all flavours, one launch (old field)
       else
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             call site_wrap(nf, n, ntau1, 1)           ! Op_Wrapup, N_Type 1 (reads nsigma_d(n,ntau1), the old field)
          enddo
       endif
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       Hs_New =   nsigma%flip(n,ntau1)
       S0_ratio          = ham%S0(n,ntau1, Hs_New )
       if ( Propose_S0 ) then
          If ( Op_V(n,nf)%type == 1)  then
             T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
             T0_Proposal_ratio = 1.d0/S0_ratio
          endif
       Endif
       Acc = .false.
       If ( T0_proposal > ranf_wrap() ) Then
          call site_update(n, ntau1, PHASE, HS_new, S0_ratio, T0_Proposal_ratio, Acc)
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
       if (bm) then
          call gpu_lr_blk_leave(n, ntau1, .false.)
       elseif (fast_lr) then
          call gpu_lr_defer(n, ntau1, Acc, .false.)
       endif
       ! Op_Wrapup with N_Type 2 (P2.6.1: non-diagonal operators; nothing for diagonal ones)
       if (.not. (bm .or. fast_lr)) then
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             call site_wrap2(nf, n, ntau1, 1)
          enddo
       endif
    Enddo
    call gpu_lr_flush_all(2)                    ! delayed / submatrix: G_d must be current at the end of the slice
    m = 0
    if (Nt_sequential_end >= Nt_sequential_start) m = Nt_sequential_end
    ! G leaves at m = Size(Op_V,1) (KNOWN_ISSUES.md, H5)
    call gpu_placegr(m, Size(Op_V,1), ntau1)
    call g_out(GR)
  end Subroutine gpu_wrapgrup

!--------------------------------------------------------------------
!> WRAPGRDO of Wrapgr_mod on the device: updates of slice tau, then G(tau) -> G(tau-1).
  Subroutine gpu_wrapgrdo(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:)
    Complex (Kind=Kind(0.d0)), intent(inout) :: PHASE
    Integer, intent(in) :: NTAU
    Logical, intent(in) :: Propose_S0
    Integer, intent(in) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

    Integer :: nf, nf_eff, n, m, m1
    Complex (Kind=Kind(0.d0)) :: HS_New
    Real (Kind=Kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
    Logical :: Acc, toggle1, fast_lr, bm

    if (N_Global_tau > 0) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    upd_cnt(10) = upd_cnt(10) + 1
    fast_lr = gpu_lr_fast()
    bm = gpu_lr_blockmode()
    call g_in(GR)
    m = Size(Op_V,1)
    m1 = 0
    if (Nt_sequential_end >= Nt_sequential_start) m1 = Nt_sequential_end
    call gpu_placegr(m, m1, ntau)
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       ! Op_Wrapdo with N_Type 2 (P2.6.1: non-diagonal operators; nothing for diagonal ones)
       if (.not. (bm .or. fast_lr)) then
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             call site_wrap2(nf, n, ntau, -1)
          enddo
       endif
       if (bm) then
          call gpu_lr_blk_enter(n, Nt_sequential_start, -1, ntau, .false.)
       elseif (fast_lr) then
          call gpu_lr_step_down(n, ntau)   ! deferred wrap of the previous vertex + gather / G(p,p), one launch
       endif
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       HS_new            = nsigma%flip(n,ntau)
       S0_ratio          = ham%S0(n,ntau,HS_new)
       if ( Propose_S0 ) then
          If ( Op_V(n,nf)%type == 1)  then
             T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
             T0_Proposal_ratio = 1.d0/S0_ratio
          endif
       Endif
       Acc = .false.
       If ( T0_proposal > ranf_wrap() ) Then
          call site_update(n, ntau, PHASE, HS_new, S0_ratio, T0_Proposal_ratio, Acc)
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
       ! the field has changed after the update: Op_Wrapdo N_Type 1 reads the new nsigma_d(n,ntau)
       if (bm) then
          call gpu_lr_blk_leave(n, ntau, .true.)      ! host wrap (final field); the block closes after its last vertex
       elseif (fast_lr) then
          call gpu_lr_defer(n, ntau, Acc, .true.)     ! the wrap follows in the next launch (or the end-of-slice sync)
       else
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             call site_wrap(nf, n, ntau, -1)
          enddo
       endif
    Enddo
    call gpu_lr_flush_all(2)                    ! delayed / submatrix: flush before the hopping wrap
    ! G is at m = Nt_sequential_start - 1; the hopping propagation needs m = 0 (KNOWN_ISSUES.md, H5)
    if (Nt_sequential_end >= Nt_sequential_start) call gpu_placegr(Nt_sequential_start-1, 0, ntau)
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       call gpu_mmthl   (G_d(1,1,nf), Ndim, Ndim, nf, ntau)
       call gpu_mmthr_m1(G_d(1,1,nf), Ndim, Ndim, nf, ntau)
    enddo
    call g_out(GR)
  end Subroutine gpu_wrapgrdo


end Module gpu_update_mod
#endif

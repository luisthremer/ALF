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
!> GPU port, ticket 205 (ADR-001 section 5, MVP): WRAPGRUP / WRAPGRDO on a device-resident Green function.  The per-site
!> wraps and the rank-1 update run on the device, the small ratio block and the accept/reject decision on the host, with
!> the random numbers drawn in exactly ALF's order, so that for the same compiler and seeds the trajectory is the CPU
!> one until round-off first changes a decision.  Compiled only with -DALF_GPU; entered from Wrapgr_mod through guarded
!> early branches.
!>
!> Per site and flavour (nf):
!>   1. wrap: Op_Wrapup/Op_Wrapdo (N_Type 1) of operator n, one fused kernel launch for N = 1 operators (rows and columns
!>      of the one site together, rows before columns for the corner element as in ALF), else the 204 two-kernel routine.
!>      N_Type 2 applies nothing for a diagonal operator (gpu_supported guarantees diagonal operators).
!>   2. ratio (host, the code of Upgrade2 verbatim for N_non_zero <= 1): needs G(P(1),P(1),nf): one complex element D->H
!>      per flavour (the only per-site device->host transfer).  Op_V%alpha and exp(Z1*E) are host data.
!>   3. decision (host): ranf_wrap() drawn at the same point as Upgrade2 (after the T0 test of WRAPGR).
!>   4. accept: per flavour 2 launches (copy the column G(:,P) and the row y = (-G(P,:), 1-G(P,P) at P) into scratch, then
!>      G += Zg * x y^T, the ZGERU of Upgrade2 for op_dim 1), then the field value on the host and one scalar H->D copy
!>      of nsigma_d(n,t) (the wrap kernels read the field on the device, so it must not go stale).
!>   Ticket 401-A (ADR-003 4.2/4.3): non-diagonal operators: Op_Wrapup/Op_Wrapdo N_Type 1 and 2 by the dense-block kernels of
!>   gpu_wrap_mod (gpu_Vop_wrap); N_Type 2 is applied after the update (WRAPGRUP) / before it (WRAPGRDO) as in Wrapgr_mod.
!>   Rank r = N_non_zero >= 2: the r x r block G(P,P) is gathered by one kernel, one D->H of r*r elements per flavour; the ratio is
!>   the code of Upgrade2 (determinant branches r = 2 and LAPACK Det for r > 2) on the host; accept: G <- G - G(:,P) W (E_P^T - G(P,:))
!>   with W = Delta Mat^-1 (Woodbury; checked against the x_v/y_v recursion of Upgrade2 by the rank self-test, ALF_GPU_UPDTEST=1),
!>   W computed on the host (Gauss-Jordan, partial pivoting), one H->D of r*r elements, two kernels.  r = 1 keeps the path above.
!>
!> Ticket 208-1: G is device-authoritative (G_d).  WRAPGRUP/DO copy the host GR to G_d only when the device copy is not current
!> and never copy back; host readers call gpu_g_host(GR) (D->H of the stale flavours), host writers gpu_g_host_wrote().
!>
!> Environment (testing): ALF_GPU_NOUPDATE=1 keeps the host WRAPGRUP/DO in a GPU build; ALF_GPU_TRACE=<file> writes, per
!> WRAPGRUP/DO call (both paths), the flipped sites of the slice and the Phase (first-divergence tests);
!> ALF_GPU_UPDTEST=1 runs the one-slice test (Wrapgr_mod, at the first WRAPGRUP call), ALF_GPU_FAULT=upd|nsig injects a
!> fault into the device update (negative controls), ALF_GPU_EXIT=1 ends the run after the test.
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_update_mod

  Use omp_lib
  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Operator_mod
  Use Fields_mod
  Use Hamiltonian_main
  Use QMC_runtime_var
  Use Control
  Use Random_Wrap
  Use gpu_state_mod
  Use gpu_wrap_mod
  Use gpu_lowrank_mod
  Use MyMats, only: Det

  Implicit none

  Private
  Public :: gpu_upd_host_only, gpu_upd_enabled, gpu_upd_override, gpu_upd_resync_fields, gpu_upd_check_fields, gpu_upd_selftest_due, gpu_upd_fault_set, gpu_upd_fault_clear
  Public :: gpu_g_host, gpu_g_host_wrote
  Public :: gpu_wrapgrup, gpu_wrapgrdo, gpu_trace_on, gpu_trace_begin, gpu_trace_end
#ifdef TIMING
  Public :: gpu_upd_time, gpu_tacc, gpu_tget
  Real (Kind=Kind(0.d0)), Save :: t_acc(24) = 0.d0
  Real (Kind=Kind(0.d0)), Save :: t_up = 0.d0
#endif

  Logical, Save :: checked = .false., enabled = .false., forced = .false., forced_val = .false.
  Logical, Save :: fault_upd = .false., fault_nsig = .false., test_host_only = .false.
  Logical, Save :: sel_done = .false., sel_checked = .false.
  Logical, Save :: tr_checked = .false., tr_on = .false.
  Integer, Save :: tr_unit = -1
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: tr_f(:)
  Type (Fields), Save :: nsigma_new
  Logical, Save :: nsigma_new_made = .false.
  ! ticket 401-A: device buffers of the rank-r update (allocated at the first rank-r site): G(P,P) block, W, C = G(:,P), T = W R
  Integer, Parameter :: RK = V_NMAX
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: rk_g(:), rk_w(:), rk_c(:), rk_t(:)

Contains

!--------------------------------------------------------------------
!> .true. iff the device update path replaces WRAPGRUP/WRAPGRDO: GPU path active and ALF_GPU_NOUPDATE /= 1.
  Logical function gpu_upd_enabled()
    Implicit none
    if (forced) then
       gpu_upd_enabled = forced_val .and. gpu_active
       return
    endif
    if (.not. checked) then
       enabled = gpu_active
       if (enabled) then
          if (gpu_env_is('ALF_GPU_NOUPDATE', '1')) then
             enabled = .false.
             call info_line('GPU update path (ticket 205): DISABLED by ALF_GPU_NOUPDATE=1, host WRAPGRUP/DO')
          else
             call info_line('GPU update path (ticket 205): ACTIVE (device wraps + rank-1 update, host decision; INTERIM per-slice G copies)')
          endif
       endif
       checked = .true.
    endif
    gpu_upd_enabled = enabled
  end function gpu_upd_enabled

!> self-test helper: force the path on/off temporarily (override = .true.), or release it
  Subroutine gpu_upd_override(on, val)
    Implicit none
    Logical, intent(in) :: on, val
    forced = on; forced_val = val
  end subroutine gpu_upd_override

!> self-test helper: nsigma_d (host shadow and device) := nsigma%f, after the test restored nsigma%f
  Subroutine gpu_upd_resync_fields()
    Implicit none
    nsigma_d = nsigma%f
    !$omp target update to(nsigma_d)     ! self-test only: the test rewinds nsigma%f on the host
  end subroutine gpu_upd_resync_fields

!> self-test helper: number of entries where the device copy of the fields differs from nsigma%f (stale nsigma_d)
  Integer function gpu_upd_check_fields()
    Implicit none
    !$omp target update from(nsigma_d)     ! self-test only: read back the device copy
    gpu_upd_check_fields = count( nsigma_d /= nsigma%f )
  end function gpu_upd_check_fields

  Logical function gpu_upd_host_only()
    Implicit none
    gpu_upd_host_only = test_host_only
  end function gpu_upd_host_only

  Subroutine gpu_upd_fault_set(name)
    Implicit none
    Character (len=*), intent(in) :: name
    select case (trim(name))
    case ('upd');  fault_upd = .true.
    case ('nsig'); fault_nsig = .true.
    case ('hostonly'); test_host_only = .true.      ! harness check: the 'device' run of the test is a second host run
    case default                ! '' or a fault of another self-test (ticket 204)
    end select
  end subroutine gpu_upd_fault_set

!> switches the update-test faults off again after the self-test
  Subroutine gpu_upd_fault_clear()
    Implicit none
    fault_upd = .false.; fault_nsig = .false.; test_host_only = .false.
  end subroutine gpu_upd_fault_clear

!> .true. once, at the first call with ALF_GPU_UPDTEST=1 (the update path must be active)
  Logical function gpu_upd_selftest_due()
    Implicit none
    gpu_upd_selftest_due = .false.
    if (sel_done) return
    if (.not. sel_checked) then
       sel_checked = .true.
       sel_done = .not. gpu_env_is('ALF_GPU_UPDTEST', '1')
    endif
    if (.not. sel_done) then
       sel_done = .true.
       gpu_upd_selftest_due = .true.
       call rank_selftest()          ! ticket 401-A (gate K2): rank-r ratio and accept, before the one-slice test of Wrapgr_mod
    endif
  end function gpu_upd_selftest_due

!--------------------------------------------------------------------
!> Ticket 401-A, gate K2 (acceptance/401.json): rank-r ratio and accept, r = 1, 2, 3, 4 (and 8 if Ndim >= 8), random complex G,
!> random sites P and eigenvalues, against the code of Upgrade2 (Mat, determinant branches, the x_v / y_v recursion and the
!> ZGEMM update, copied in ref_upgrade2).  Three numbers per r: ratio rel diff (device gather + rk_ratio vs reference), G after the
!> accept (device rk_accept with W = Delta Mat^-1 vs the recursion; max|dG|/max|G|), and the derivation check on the host only
!> (G - C W R in host arithmetic vs the recursion).  Threshold 1e-10.  The model's own vertices are covered by the one-slice test.
  Subroutine rank_selftest()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-10
    Integer, parameter :: rlist(5) = (/ 1, 2, 3, 4, 8 /)
    Integer(int64) :: sd
    Integer :: ir, r, it, i, j, a, b, pool
    Integer, allocatable :: PP(:)
    Complex (Kind=Kind(0.d0)), allocatable :: G0(:,:), Gr(:,:), Gh(:,:), Dt(:,:)
    Complex (Kind=Kind(0.d0)) :: Z1, Dref, Ddev, gb(RK*RK), dl(RK), mk(RK,RK), w(RK*RK), dl2(RK)
    Real (Kind=Kind(0.d0)) :: E(RK), dr, dg, dh, mr, mg, mh
    Character (len=256) :: line
    sd = 4012026_int64
    call rk_alloc()
    allocate(G0(Ndim,Ndim), Gr(Ndim,Ndim), Gh(Ndim,Ndim), Dt(Ndim,Ndim), PP(RK))
    !$omp target enter data map(alloc: Dt, PP)
    mr = 0.d0; mg = 0.d0; mh = 0.d0
    do ir = 1, size(rlist)
       r = rlist(ir)
       if (r > Ndim) cycle
       dr = 0.d0; dg = 0.d0; dh = 0.d0
       do it = 1, 3
          call rand_mat(G0, sd)
          G0 = 0.3d0*G0
          do i = 1, Ndim
             G0(i,i) = G0(i,i) + 0.5d0
          enddo
          ! r distinct sites (partial Fisher-Yates over 1..Ndim)
          pool = Ndim
          block
            Integer :: perm(Ndim), k, tmp
            do i = 1, Ndim; perm(i) = i; enddo
            do i = 1, r
               k = i + int((rnd(sd) + 0.5d0)*real(pool - i + 1, kind(0.d0)))
               if (k > pool) k = pool
               tmp = perm(i); perm(i) = perm(k); perm(k) = tmp
               PP(i) = perm(i)
            enddo
          end block
          do i = 1, r
             E(i) = rnd(sd)*2.d0
          enddo
          Z1 = cmplx(0.4d0 + 0.2d0*rnd(sd), 0.3d0*rnd(sd), kind(0.d0))
          ! reference: Upgrade2's code
          Gr = G0
          call ref_upgrade2(Gr, PP, r, Z1, E, Dref)
          ! device: gather, host ratio, W, accept
          Dt = G0
          !$omp target update to(Dt, PP)     ! self-test only: random input
          call k_rk_gpp(Dt, Ndim, PP, r, rk_g)
          !$omp target update from(rk_g(1:r*r))
          gb(1:r*r) = rk_g(1:r*r)
          call rk_ratio(r, gb, Z1, E, Ddev, dl, mk)
          call rk_wmat(r, mk, dl, w)
          call rk_accept(Dt, Ndim, PP, r, w)
          !$omp target update from(Dt)       ! self-test only: result for the comparison
          dr = max(dr, abs(Ddev - Dref)/abs(Dref))
          dg = max(dg, maxval(abs(Dt - Gr))/maxval(abs(Gr)))
          ! derivation check, host only: G - C W R with C = G(:,P), R = E_P^T - G(P,:)
          do b = 1, r
             do a = 1, r
                gb(a + (b-1)*r) = G0(PP(a), PP(b))
             enddo
          enddo
          call rk_ratio(r, gb, Z1, E, Ddev, dl2, mk)
          call rk_wmat(r, mk, dl2, w)
          Gh = G0
          do j = 1, Ndim
             do i = 1, Ndim
                do a = 1, r
                   do b = 1, r
                      Gh(i,j) = Gh(i,j) - G0(i,PP(a)) * w(a + (b-1)*r) * (merge(1.d0, 0.d0, j == PP(b)) - G0(PP(b),j))
                   enddo
                enddo
             enddo
          enddo
          dh = max(dh, maxval(abs(Gh - Gr))/maxval(abs(Gr)))
       enddo
       write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1,A)') 'GPU rank-r update self-test K2 (401-A) r = ', r, &
            & ': ratio rel diff ', dr, ', G after accept (device) ', dg, ', host W = Delta Mat^-1 formula vs Upgrade2 recursion ', dh, &
            & ' (tol ', tol, ')'
       call info_line(line); write(output_unit,'(1X,A)') trim(line)
       if (.not. (dr <= mr)) mr = dr
       if (.not. (dg <= mg)) mg = dg
       if (.not. (dh <= mh)) mh = dh
    enddo
    !$omp target exit data map(delete: Dt, PP)
    deallocate(G0, Gr, Gh, Dt, PP)
    if (.not. (mr <= tol .and. mg <= tol .and. mh <= tol)) call selftest_fail('K2: rank-r ratio / accept vs Upgrade2 (401-A)')
    call info_line('GPU rank-r update self-test K2 PASSED'); write(output_unit,'(A)') ' GPU rank-r update self-test K2 PASSED'
  end subroutine rank_selftest

!> Reference for rank_selftest: the ratio (alpha = 0) and accept of Upgrade2 (upgrade_mod.F90) for one flavour, sites PP(1:r) in
!> the eigenbasis, eigenvalues E(1:r), Z1; code copied from Upgrade2 (ratio lines 177-203, accept lines 248-297)
  Subroutine ref_upgrade2(GR, PP, op_dim_nf, Z1, E, D_mat)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:)
    Integer, intent(in) :: PP(:), op_dim_nf
    Complex (Kind=Kind(0.d0)), intent(in) :: Z1
    Real (Kind=Kind(0.d0)), intent(in) :: E(:)
    Complex (Kind=Kind(0.d0)), intent(out) :: D_mat
    Complex (Kind=Kind(0.d0)) :: Mat(op_dim_nf,op_dim_nf), Delta(op_dim_nf), myexp, Z, s1, s2, alpha, beta
    Complex (Kind=Kind(0.d0)), allocatable :: u(:,:), v(:,:), x_v(:,:), y_v(:,:), xp_v(:,:), syu(:), sxv(:), zarr(:,:), grarr(:,:)
    Integer :: n, m, i, nd
    nd = size(GR,1)
    Do m = 1,op_dim_nf
       myexp = exp( Z1* E(m) )
       Z = myexp - 1.d0
       Delta(m) = Z
       do n = 1,op_dim_nf
          Mat(n,m) = - Z * GR( PP(n), PP(m) )
       Enddo
       Mat(m,m) = myexp + Mat(m,m)
    Enddo
    If (op_dim_nf == 1 ) then
       D_mat = Mat(1,1)
    elseif (op_dim_nf == 2 ) then
       s1 = Mat(1,1)*Mat(2,2)
       s2 = Mat(2,1)*Mat(1,2)
       If (Abs(s1) > Abs(s2)) then
          D_mat = s1*(1.D0 - s2/s1)
       else
          D_mat = s2*(s1/s2 - 1.D0)
       Endif
    else
       D_mat = Det(Mat,op_dim_nf)
    endif
    allocate(u(nd,op_dim_nf), v(nd,op_dim_nf), x_v(nd,op_dim_nf), y_v(nd,op_dim_nf), xp_v(nd,op_dim_nf))
    u = 0.d0; v = 0.d0
    do n = 1,op_dim_nf
       u( PP(n), n) = Delta(n)
       do i = 1,nd
          v(i,n) = - GR( PP(n), i )
       enddo
       v(PP(n), n)  = 1.d0 - GR( PP(n),  PP(n))
    enddo
    x_v = 0.d0; y_v = 0.d0
    i = PP(1)
    x_v(i, 1) = u(i, 1)/(1.d0 + v(i,1)*u(i,1) )
    call zcopy(nd, v(:, 1), 1, y_v(:, 1), 1)
    beta = 0.d0
    do n = 2,op_dim_nf
       call zcopy(nd, u(:, n), 1, x_v(:, n), 1)
       call zcopy(nd, v(:, n), 1, y_v(:, n), 1)
       Z = 1.d0 + u( PP(n), n)*v(PP(n),n)
       alpha = -1.D0
       Allocate(syu(n), sxv(n))
       call zgemv('T', nd, n-1, alpha, y_v, nd, u(1,n), 1, beta , syu, 1)
       call zgemv('T', nd, n-1, alpha, x_v, nd, v(1,n), 1, beta , sxv, 1)
       alpha = 1.D0
       call zgemv('N', nd, n-1, alpha, x_v, nd, syu, 1, alpha, x_v(1, n), 1)
       call zgemv('N', nd, n-1, alpha, y_v, nd, sxv, 1, alpha, y_v(1, n), 1)
       do m = 1,n-1
          Z = Z - syu(m)*sxv(m)
       enddo
       Z = 1.D0/Z
       call zscal(nd, Z, x_v(1, n), 1)
       Deallocate(syu, sxv)
    enddo
    IF (op_dim_nf == 1) THEN
       CALL ZCOPY(nd, gr(1, PP(1)), 1, xp_v(1, 1), 1)
       Z = -x_v(PP(1), 1)
       CALL ZGERU(nd, nd, Z, xp_v(1,1), 1, y_v(1, 1), 1, gr(1,1), nd)
    ELSE
       Allocate (zarr(op_dim_nf, op_dim_nf), grarr(nd, op_dim_nf))
       Zarr = x_v(PP(1:op_dim_nf), :)
       grarr = gr(:, PP(1:op_dim_nf))
       beta  = 0.d0
       alpha = 1.D0
       CALL ZGEMM('N', 'N', nd, op_dim_nf, op_dim_nf, alpha, grarr, nd, Zarr, op_dim_nf, beta, xp_v, nd)
       Deallocate(Zarr, grarr)
       beta  = cmplx ( 1.0d0, 0.0d0, kind(0.D0))
       alpha = -1.D0
       CALL ZGEMM('N','T',nd,nd,op_dim_nf,alpha,xp_v, nd, y_v, nd,beta,gr(1,1), nd)
    ENDIF
    deallocate(u, v, x_v, y_v, xp_v)
  end subroutine ref_upgrade2

!====================================================================
!  Trace (first-divergence tests)
!====================================================================
  Logical function gpu_trace_on()
    Implicit none
    Character (len=256) :: env
    Integer :: st
    if (.not. tr_checked) then
       tr_checked = .true.
       if (gpu_env_get('ALF_GPU_TRACE', env)) then
          open(newunit=tr_unit, file=trim(env), status='unknown', action='write')
          tr_on = .true.
       endif
    endif
    gpu_trace_on = tr_on
  end function gpu_trace_on

  Subroutine gpu_trace_begin(ntau)
    Implicit none
    Integer, intent(in) :: ntau
    if (.not. allocated(tr_f)) allocate(tr_f(size(Op_V,1)))
    tr_f = nsigma%f(:,ntau)
  end subroutine gpu_trace_begin

  ! one line per call: direction, slice, flipped sites of the slice as 0/1 string, Phase
  Subroutine gpu_trace_end(dir, ntau, Phase)
    Implicit none
    Character (len=1), intent(in) :: dir
    Integer, intent(in) :: ntau
    Complex (Kind=Kind(0.d0)), intent(in) :: Phase
    Character (len=:), Allocatable :: str
    Integer :: n
    allocate(character(len=size(Op_V,1)) :: str)
    do n = 1, size(Op_V,1)
       if (nsigma%f(n,ntau) /= tr_f(n)) then
          str(n:n) = '1'
       else
          str(n:n) = '0'
       endif
    enddo
    write(tr_unit,'(A1,1X,I5,1X,A,1X,ES25.17,1X,ES25.17)') dir, ntau, str, real(Phase), aimag(Phase)
    flush(tr_unit)
  end subroutine gpu_trace_end

!====================================================================
!  Interim copies of G (host GR <-> G_d)
!====================================================================
  !> Ticket 208-1: G is device-authoritative during the sweep.  g_in copies the host GR to G_d only for a flavour whose device
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
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
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: x(*), y(*)
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
!> region lost its imaginary part (measured: complex G, ticket 205; real-G models were not affected).
  Subroutine k_upd_geru(G, ld, nd, zr, zi, x, y)
    Implicit none
    Integer, intent(in) :: ld, nd
    Real (Kind=Kind(0.d0)), intent(in) :: zr, zi
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: x(*), y(*)
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

!> Ticket 401-A: buf(a + (b-1)*r) = G(PP(a), PP(b)), a, b = 1..r (one team)
  Subroutine k_rk_gpp(G, ld, PP, r, buf)
    Implicit none
    Integer, intent(in) :: ld, r
    Integer, intent(in) :: PP(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: buf(*)
    Integer :: a, b
    !$omp target teams distribute parallel do collapse(2) num_teams(1)
    do b = 1, r
       do a = 1, r
#ifdef ALF_PLANT_GPPT
          buf(a + (b-1)*r) = G(PP(b), PP(a))
#else
          buf(a + (b-1)*r) = G(PP(a), PP(b))
#endif
       enddo
    enddo
  end subroutine k_rk_gpp

!> Ticket 401-A: C(i,a) = G(i,PP(a)) and T(a,i) = sum_b W(a,b) R(b,i), R(b,i) = delta(i,PP(b)) - G(PP(b),i), for i = 1..nd
!> (before G changes; thread i writes row i of C and column i of T)
  Subroutine k_rk_prep(G, ld, nd, PP, r, W, C, T)
    Implicit none
    Integer, intent(in) :: ld, nd, r
    Integer, intent(in) :: PP(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*), W(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: C(*), T(*)
    Integer :: i, a, b
    Complex (Kind=Kind(0.d0)) :: acc, rb
    !$omp target teams distribute parallel do private(a, b, acc, rb)
    do i = 1, nd
       do a = 1, r
          C(i + (a-1)*nd) = G(i, PP(a))
       enddo
       do a = 1, r
          acc = (0.d0, 0.d0)
          do b = 1, r
             rb = - G(PP(b), i)
             if (i == PP(b)) rb = 1.d0 - G(PP(b), i)
             acc = acc + W(a + (b-1)*r) * rb
          enddo
          T(a + (i-1)*r) = acc
       enddo
    enddo
  end subroutine k_rk_prep

!> Ticket 401-A: G(i,j) <- G(i,j) - sum_a C(i,a) T(a,j)
  Subroutine k_rk_upd(G, ld, nd, r, C, T)
    Implicit none
    Integer, intent(in) :: ld, nd, r
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: C(*), T(*)
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

  Subroutine rk_alloc()
    Implicit none
    if (allocated(rk_g)) return
    allocate(rk_g(RK*RK), rk_w(RK*RK), rk_c(Ndim*RK), rk_t(RK*Ndim))
    rk_g = (0.d0,0.d0); rk_w = (0.d0,0.d0); rk_c = (0.d0,0.d0); rk_t = (0.d0,0.d0)
    !$omp target enter data map(alloc: rk_g, rk_w, rk_c, rk_t)     ! device scratch of the rank-r update, no transfer
  end subroutine rk_alloc

!> Ticket 401-A: Mat and determinant of Upgrade2 for one flavour, rank r >= 1 (the same lines: Mat(n,m) = -Z_m G(P_n,P_m),
!> Mat(m,m) += exp(Z1 E_m), branches r = 1, 2, else Det).  gb(a,b) = G(P_a,P_b) (column-major r x r), Matc = Mat before Det.
  Subroutine rk_ratio(r, gb, Z1, E, D_mat, Delta, Matc)
    Implicit none
    Integer, intent(in) :: r
    Complex (Kind=Kind(0.d0)), intent(in) :: gb(*), Z1
    Real (Kind=Kind(0.d0)), intent(in) :: E(*)
    Complex (Kind=Kind(0.d0)), intent(out) :: D_mat, Delta(*), Matc(RK,RK)
    Complex (Kind=Kind(0.d0)) :: Mat(r,r), myexp, Z, s1, s2
    Integer :: n, m
    Do m = 1, r
       myexp = exp( Z1* E(m) )
       Z = myexp - 1.d0
       Delta(m) = Z
       do n = 1, r
          Mat(n,m) = - Z * gb(n + (m-1)*r)
       Enddo
       Mat(m,m) = myexp + Mat(m,m)
    Enddo
    Matc(1:r,1:r) = Mat
    If (r == 1) then
       D_mat = Mat(1,1)
    elseif (r == 2) then
       s1 = Mat(1,1)*Mat(2,2)
       s2 = Mat(2,1)*Mat(1,2)
       If (Abs(s1) > Abs(s2)) then
          D_mat = s1*(1.D0 - s2/s1)
       else
          D_mat = s2*(s1/s2 - 1.D0)
       Endif
    else
       D_mat = Det(Mat, r)
    endif
  end subroutine rk_ratio

!> Ticket 401-A: W = Delta Mat^-1 (r x r, column-major into w), Gauss-Jordan with partial pivoting on a copy of Mat
  Subroutine rk_wmat(r, Matc, Delta, w)
    Implicit none
    Integer, intent(in) :: r
    Complex (Kind=Kind(0.d0)), intent(in) :: Matc(RK,RK), Delta(*)
    Complex (Kind=Kind(0.d0)), intent(out) :: w(*)
    Complex (Kind=Kind(0.d0)) :: A(r,r), X(r,r), tmp(r), f
    Integer :: i, j, k, piv
    A = Matc(1:r,1:r)
    X = (0.d0,0.d0)
    do i = 1, r
       X(i,i) = (1.d0,0.d0)
    enddo
    do k = 1, r
       piv = k
       do i = k+1, r
          if (abs(A(i,k)) > abs(A(piv,k))) piv = i
       enddo
       if (piv /= k) then
          tmp = A(k,:); A(k,:) = A(piv,:); A(piv,:) = tmp
          tmp = X(k,:); X(k,:) = X(piv,:); X(piv,:) = tmp
       endif
       f = 1.d0 / A(k,k)
       A(k,:) = A(k,:) * f; X(k,:) = X(k,:) * f
       do i = 1, r
          if (i /= k) then
             f = A(i,k)
             A(i,:) = A(i,:) - f*A(k,:); X(i,:) = X(i,:) - f*X(k,:)
          endif
       enddo
    enddo
    do j = 1, r
       do i = 1, r
#ifdef ALF_PLANT_WT
          w(j + (i-1)*r) = Delta(i) * X(i,j)
#else
          w(i + (j-1)*r) = Delta(i) * X(i,j)
#endif
       enddo
    enddo
  end subroutine rk_wmat

!> Ticket 401-A: device rank-r accept on G (ld x Ndim) for the site list PP(1..r) with W (host, r x r): upload W, two kernels
  Subroutine rk_accept(G, ld, PP, r, w)
    Implicit none
    Integer, intent(in) :: ld, r
    Integer, intent(in) :: PP(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: w(*)
    rk_w(1:r*r) = w(1:r*r)
    !$omp target update to(rk_w(1:r*r))    ! r*r complex H->D per flavour and accepted rank-r flip: W is computed on the host (ADR-003 4.3)
    call k_rk_prep(G, ld, Ndim, PP, r, rk_w, rk_c, rk_t)
    call k_rk_upd(G, ld, Ndim, r, rk_c, rk_t)
  end subroutine rk_accept

!====================================================================
!  Host-side orchestration
!====================================================================
  Subroutine site_wrap(nf, n, t, dir)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    if (gpu_lr_mode() /= 0) then
       if (V_meta(1,n,nf) == 1) then
          call gpu_lr_wrap(nf, n, t, dir)       ! pending delayed / submatrix state follows the wrap
          return
       endif
       call gpu_lr_flush(nf, 2)                 ! N > 1 vertex: the general wrap needs a flushed G
    endif
    if (V_meta(1,n,nf) == 1) then
       call k_site_wrap(G_d(1,1,nf), Ndim, Ndim, n, t, dir, Op_V(n,nf)%P(1), V_meta(3,n,nf), V_type(n))
       upd_cnt(3) = upd_cnt(3) + 1
    else
       call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, t, dir, 1)     ! N > 1: diagonal (row / column scalings) or dense block (401-A)
       upd_cnt(3) = upd_cnt(3) + 2
    endif
  end subroutine site_wrap

!> Ticket 401-A: Op_Wrapup (dir = +1, after the update) / Op_Wrapdo (dir = -1, before it) with N_Type 2: G <- U G U^H resp.
!> U^H G U on P; nothing for a diagonal operator (as ALF)
  Subroutine site_wrap2(nf, n, t, dir)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    if (V_meta(6,n,nf) == 1) return
#ifdef ALF_PLANT_NT2
    return
#endif
    if (gpu_lr_mode() /= 0) call gpu_lr_flush(nf, 2)       ! never pending here (non-diagonal vertices flush first); kept for safety
    call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, t, dir, 2)
    upd_cnt(3) = upd_cnt(3) + 2
  end subroutine site_wrap2

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
    Complex (Kind=Kind(0.d0)) :: Z, D_Mat, myexp, g_loc, Mat11, vii, xv, Zg
    Real (Kind=Kind(0.d0)) :: Weight
    Integer :: nf, nf_eff, op_dim_nf, p1, po
    Logical :: lr, lrf
    ! ticket 401-A, rank r >= 2: per computed flavour the block G(P,P), Delta, Mat (Upgrade2), W
    Complex (Kind=Kind(0.d0)) :: gbk(RK*RK, N_FL), dlt(RK, N_FL), matk(RK, RK, N_FL), wk(RK*RK)

    toggle = .false.
    lr = ( gpu_lr_mode() /= 0 )
    if (lr) lr = gpu_lr_site_ok(N_op)
    lrf = lr .and. gpu_lr_fast()
    if (.not. lr) call gpu_lr_flush_all(2)       ! standard rank-1 update of this site: G_d must be current
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
       if (op_dim_nf == 0) then
          D_mat = 1.0d0
       elseif (op_dim_nf >= 2) then
          ! ticket 401-A: rank r = op_dim_nf >= 2 (lr is .false.: gpu_lr_site_ok needs N = 1).  The ratio needs the r x r block
          ! G(P,P) on the host: one gather kernel and one D->H of r*r elements per flavour (unavoidable while the decision is on the host)
          call rk_alloc()
          po = V_meta(2,n_op,nf)
          call k_rk_gpp(G_d(1,1,nf), Ndim, V_P(po+1), op_dim_nf, rk_g)
          !$omp target update from(rk_g(1:op_dim_nf*op_dim_nf))
          upd_cnt(4) = upd_cnt(4) + 1
          upd_cnt(5) = upd_cnt(5) + op_dim_nf*op_dim_nf
          gbk(1:op_dim_nf*op_dim_nf, nf_eff) = rk_g(1:op_dim_nf*op_dim_nf)
          call rk_ratio(op_dim_nf, gbk(1,nf_eff), Z1, Op_V(n_op,nf)%E, D_mat, dlt(1,nf_eff), matk(1,1,nf_eff))
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

    ! ticket 401-B: fill the ratio of the flavours that are not calculated (as upgrade_mod: before Product(Ratio))
    if (reconstruction_needed) call ham%weight_reconstruction(Ratio)

    Ratiotot = Product(Ratio)
    nf = 1
    Ratiotot = (Ratiotot**dble(N_SUN)) * nsigma_new%Gama(1,1)/nsigma%Gama(n_op,nt)
    Ratiotot = Ratiotot * Prev_Ratiotot
    Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )

    upd_cnt(1) = upd_cnt(1) + 1
    if (upd_cnt(1) < 40) write(error_unit,'(A,3I4,6ES12.4)') 'DBG401 ', n_op, nt, Op_V(n_op,1)%N_non_zero, Ratiotot, Weight, gbk(1,1), Phase  ! DBG401
    if ( Weight > ranf_wrap() )  Then
       toggle = .true.
       upd_cnt(2) = upd_cnt(2) + 1
       Phase = Phase * Ratiotot/sqrt(Ratiotot*conjg(Ratiotot))
       Do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          if (Op_V(n_op,nf)%N_non_zero >= 2) then
             ! ticket 401-A: G <- G - G(:,P) W (E_P^T - G(P,:)), W = Delta Mat^-1 (r x r, host)
             op_dim_nf = Op_V(n_op,nf)%N_non_zero
             call rk_wmat(op_dim_nf, matk(1,1,nf_eff), dlt(1,nf_eff), wk)
             if (fault_upd) wk(1) = wk(1) * (1.d0 + 1.d-6)        ! negative control (ALF_GPU_FAULT=upd)
             po = V_meta(2,n_op,nf)
             call rk_accept(G_d(1,1,nf), Ndim, V_P(po+1), op_dim_nf, wk)
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
#ifdef ALF_PLANT_STALE
       elseif (.not. fault_nsig .and. Op_V(n_op,1)%N_non_zero < 2) then
#else
       elseif (.not. fault_nsig) then
#endif
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

    Integer :: nf, nf_eff, NTAU1, n
    Complex (Kind=Kind(0.d0)) :: HS_New
    Real (Kind=Kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
    Logical :: Acc, toggle1, fast_lr, bm
#ifdef TIMING
    Real (Kind=Kind(0.d0)) :: t0w
#endif

    if (N_Global_tau > 0) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    upd_cnt(10) = upd_cnt(10) + 1
#ifdef TIMING
    t0w = omp_get_wtime()
#endif
    fast_lr = gpu_lr_fast()
    bm = gpu_lr_blockmode()
    NTAU1 = NTAU + 1
    call g_in(GR)
    ! Wrap up with B(tau1)
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       call gpu_mmthr   (G_d(1,1,nf), Ndim, Ndim, nf)
       call gpu_mmthl_m1(G_d(1,1,nf), Ndim, Ndim, nf)
    Enddo
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
       ! Op_Wrapup with N_Type 2: nothing for diagonal operators (block / fast modes have diagonal N = 1 vertices only)
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          call site_wrap2(nf, n, ntau1, 1)
       enddo
    Enddo
    call gpu_lr_flush_all(2)                    ! delayed / submatrix: G_d must be current at the end of the slice
    call g_out(GR)
#ifdef TIMING
    t_up = t_up + (omp_get_wtime() - t0w)
#endif
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

    Integer :: nf, nf_eff, n
    Complex (Kind=Kind(0.d0)) :: HS_New
    Real (Kind=Kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
    Logical :: Acc, toggle1, fast_lr, bm
#ifdef TIMING
    Real (Kind=Kind(0.d0)) :: t0w
#endif

    if (N_Global_tau > 0) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    upd_cnt(10) = upd_cnt(10) + 1
#ifdef TIMING
    t0w = omp_get_wtime()
#endif
    fast_lr = gpu_lr_fast()
    bm = gpu_lr_blockmode()
    call g_in(GR)
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       ! Op_Wrapdo with N_Type 2: nothing for diagonal operators (block / fast modes have diagonal N = 1 vertices only)
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          call site_wrap2(nf, n, ntau, -1)
       enddo
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
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       call gpu_mmthl   (G_d(1,1,nf), Ndim, Ndim, nf)
       call gpu_mmthr_m1(G_d(1,1,nf), Ndim, Ndim, nf)
    enddo
    call g_out(GR)
#ifdef TIMING
    t_up = t_up + (omp_get_wtime() - t0w)
#endif
  end Subroutine gpu_wrapgrdo

#ifdef TIMING
  Subroutine gpu_tacc(i, sgn)       ! stopwatch: call with -1 before, +1 after
    Implicit none
    Integer, intent(in) :: i
    Real (Kind=Kind(0.d0)), intent(in) :: sgn
    t_acc(i) = t_acc(i) + sgn*omp_get_wtime()
  end subroutine gpu_tacc
  Real (Kind=Kind(0.d0)) function gpu_tget(i)
    Implicit none
    Integer, intent(in) :: i
    gpu_tget = t_acc(i)
  end function gpu_tget

  Real (Kind=Kind(0.d0)) function gpu_upd_time()
    Implicit none
    gpu_upd_time = t_up
  end function gpu_upd_time
#endif

end Module gpu_update_mod
#endif

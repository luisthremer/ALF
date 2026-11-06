!  Copyright (C) 2016 - 2026 The ALF project
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
!> GPU version of the delayed and submatrix local updates (test branch gpu-port-delayed-submatrix).
!> The algebra is that of Prog/Lowrank_update_mod.F90 (docs/algorithm-notes.md sections 3, 4, 12) restricted to what the
!> GPU path supports (gpu_supported): rank-1 vertices (N_non_zero <= 1), diagonal operators.  Selected by the parameters
!> Sequential_scheme = "delayed" | "submatrix" and N_delay, the same keys as on the CPU.
!>
!> Delayed (per flavour nf), G = G0 - X Y^T with X, Y (Ndim x N_delay) device-resident, G0 = G_d:
!>   proposal : G(p,p) = G0(p,p) - X(p,:) Y(p,:)^T          one single-thread region, 2 reals D->H (replaces the 1-element copy)
!>   accept   : c = G0(:,p) - X Y(p,:)^T, r = G0(p,:) - X(p,:) Y^T, X(:,j+1) = (Delta/d) c, Y(:,j+1) = e_p - r: ONE kernel
!>              launch per flavour (the rank-1 update needs 2 launches and touches the whole N x N matrix; here O(N j))
!>   wrap     : G0 <- L G0 R as before, X(r,:) *= za, Y(r,:) *= zb (rows of the wrapped site), one kernel
!>   flush    : G0 <- G0 - X Y^T, one ZGEMM (la_gemm), when j = N_delay, before an N > 1 vertex, and at the end of
!>              WRAPGRUP / WRAPGRDO (the hopping wraps, stabilization, measurements all read a flushed G_d)
!> Submatrix, M-form with the pending sites S (distinct, k <= N_delay), M (k x k) kept on the HOST:
!>   G = G0 - G0(:,S) M (E_S^T - G0(S,:))
!>   proposal : gather G0(p,p), G0(S,p), G0(p,S) (2k+1 complex, one region, one D->H), host O(k^2) -> G(p,p)
!>   accept   : host only: M <- M + (Delta/d) a b^T (padding S with p if needed); NO device work at all
!>   wrap     : G0 as before (one kernel); host M(m,:) *= za, M(:,m) *= zb for a wrapped site in S
!>   flush    : Qt = E_S - G0(S,:)^T (gather), Y = Qt M^T (ZGEMM), G0 <- G0 - G0(:,S) Y^T (ZGEMM), as the delayed flush
!> Neither scheme changes the random-number order: the same ratio code, the same draws.  Numerical agreement is not the
!> goal of this test branch (see the ticket text); the CPU schemes are the reference.
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_lowrank_mod

  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Files_mod, only: str_to_upper
  Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map, Op_V
  Use QMC_runtime_var, only: get_Sequential_scheme, get_N_delay
  Use alf_la_mod, only: la_gemm, LA_N, LA_T
  Use gpu_state_mod

  Implicit none

  Private
  Public :: gpu_lr_mode, gpu_lr_init, gpu_lr_finalize, gpu_lr_site_ok, gpu_lr_pending, gpu_lr_flush_all, gpu_lr_flush
  Public :: gpu_lr_wrap, gpu_lr_gpp, gpu_lr_accept

  Integer, parameter :: LR_STD = 0, LR_DELAYED = 1, LR_SUBMATRIX = 2
  Complex (Kind=Kind(0.d0)), parameter :: czero_c = (0.d0,0.d0)

  Logical, Save :: inited = .false.
  Integer, Save :: mode = LR_STD, nd = 0, ndim_ = 0, nfl_ = 0
  Integer, Save :: n_flush(3) = 0       ! flushes: capacity, forced (N>1 vertex or slice end), total
  Integer(int64), Save :: n_acc = 0, n_prop = 0

  ! per flavour (index nf of G_d)
  Integer, Allocatable, Save :: jc(:)                          ! pending columns (delayed) / pending sites k (submatrix)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: X_d(:,:,:), Y_d(:,:,:)       ! (Ndim, nd, N_FL) device
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: gbuf(:)      ! (2*nd+1) gather buffer, device + host
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: res2(:)      ! (2) single-thread result
  Integer, Allocatable, Save :: S_h(:,:), slot_h(:,:)          ! (nd, N_FL), (Ndim, N_FL) host
  Integer, Allocatable, Save :: S_d(:)                          ! (nd) device copy of S of one flavour
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: M_h(:,:,:)   ! (nd, nd, N_FL) host
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: M_d(:,:)     ! (nd, nd) device copy of M of one flavour (flush)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: W_h(:,:), grow_h(:,:)   ! cached proposal data (nd, N_FL) host, per flavour
  Integer, Allocatable, Save :: prop_p(:), prop_k(:)          ! site and pending count of the cached proposal, per flavour

Contains

!--------------------------------------------------------------------
  Integer function gpu_lr_mode()
    Implicit none
    if (.not. inited) call gpu_lr_init()
    gpu_lr_mode = mode
  end function gpu_lr_mode

!> static rule of Lowrank_static_N_delay (docs/algorithm-notes.md), integer arithmetic
  Integer function static_nd(sub, Ndim, max_pend)
    Implicit none
    Logical, intent(in) :: sub
    Integer, intent(in) :: Ndim, max_pend
    Integer :: raw, U, n, p
    Integer(int64) :: nd8
    nd8 = int(Ndim, int64)
    raw = 1
    if (.not. sub) then
       U = 128
       do while (raw < U .and. nd8*nd8 >= 2_int64*(20_int64*raw)**2)
          raw = 2*raw
       enddo
    else
       U = 256
       do while (raw < U .and. nd8 >= 12_int64*(2_int64*raw))
          raw = 2*raw
       enddo
    endif
    n = max(16, min(U, raw))
    n = min(n, max(max_pend,1))
    if (sub) n = min(n, Ndim)
    p = 1
    do while (2*p <= n)
       p = 2*p
    enddo
    static_nd = p
  end function static_nd

!--------------------------------------------------------------------
  Subroutine gpu_lr_init()
    Implicit none
    Character (len=:), Allocatable :: name
    Character (len=200) :: line
    Integer :: nfl, nop, n_nz, nreq, nf, n
    if (inited) return
    inited = .true.
    name = str_to_upper(trim(adjustl(get_Sequential_scheme())))
    select case (name)
    case ('STANDARD')
       mode = LR_STD; return
    case ('DELAYED')
       mode = LR_DELAYED
    case ('SUBMATRIX')
       mode = LR_SUBMATRIX
    case default
       write(error_unit,*) 'gpu_lr_init: unknown Sequential_scheme "', trim(get_Sequential_scheme()), '"'
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    end select
    ndim_ = Ndim; nfl_ = N_FL
    nop = size(Op_V,1)
    n_nz = 0
    do n = 1, nop
       do nf = 1, N_FL
          if (Op_V(n,nf)%N_non_zero > 0) then
             n_nz = n_nz + 1
             exit
          endif
       enddo
    enddo
    nreq = get_N_delay()
    if (nreq >= 1) then
       nd = nreq
    else
       nd = static_nd(mode == LR_SUBMATRIX, Ndim, n_nz)
    endif
    nd = max(1, min(nd, Ndim))
    allocate( jc(N_FL), S_h(nd,N_FL), slot_h(Ndim,N_FL), M_h(nd,nd,N_FL), W_h(nd,N_FL), grow_h(nd,N_FL), prop_p(N_FL), prop_k(N_FL) )
    prop_p = 0; prop_k = -1
    jc = 0; S_h = 0; slot_h = 0; M_h = cmplx(0.d0,0.d0,kind(0.d0))
    allocate( gbuf(2*nd+1), res2(2), S_d(nd), M_d(nd,nd) )
    gbuf = cmplx(0.d0,0.d0,kind(0.d0)); res2 = 0.d0; S_d = 0; M_d = cmplx(0.d0,0.d0,kind(0.d0))
    if (mode == LR_DELAYED) then
       allocate( X_d(Ndim,nd,N_FL), Y_d(Ndim,nd,N_FL) )
       X_d = cmplx(0.d0,0.d0,kind(0.d0)); Y_d = X_d
       !$omp target enter data map(alloc: X_d, Y_d)
    else
       allocate( X_d(Ndim,nd,N_FL), Y_d(1,1,1) )     ! X_d holds G0(:,S) in the flush
       X_d = cmplx(0.d0,0.d0,kind(0.d0)); Y_d = X_d(1,1,1)
       !$omp target enter data map(alloc: X_d)
    endif
    !$omp target enter data map(alloc: gbuf, res2, S_d, M_d)
    if (mode == LR_DELAYED) then
       write(line,'(A,I0,A,I0,A)') 'GPU delayed update (test branch): N_delay = ', nd, ', Ndim = ', Ndim, &
            & ', G = G0 - X Y^T on the device, one ZGEMM per flush'
    else
       write(line,'(A,I0,A,I0,A)') 'GPU submatrix update (test branch): N_delay = ', nd, ', Ndim = ', Ndim, &
            & ', M on the host, G0 untouched between flushes'
    endif
    call info_line(trim(line))
  end subroutine gpu_lr_init

!--------------------------------------------------------------------
  Subroutine gpu_lr_finalize()
    Implicit none
    Character (len=200) :: line
    if (.not. inited) return
    if (mode /= LR_STD) then
       write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'GPU low-rank update counts: proposals ', n_prop, ', accepted ', n_acc, &
            & ', flushes ', n_flush(3), ' (capacity ', n_flush(1), ', forced ', n_flush(2)
       call info_line(trim(line)//')')
       !$omp target exit data map(delete: gbuf, res2, S_d, M_d)
       !$omp target exit data map(delete: X_d)
       if (mode == LR_DELAYED) then
          !$omp target exit data map(delete: Y_d)
       endif
       deallocate( jc, S_h, slot_h, M_h, W_h, grow_h, prop_p, prop_k, gbuf, res2, S_d, M_d, X_d, Y_d )
    endif
    inited = .false.; mode = LR_STD
  end subroutine gpu_lr_finalize

!> .true. iff vertex n can use the low-rank path: every computed flavour has a single-site operator
  Logical function gpu_lr_site_ok(n)
    Implicit none
    Integer, intent(in) :: n
    Integer :: nf_eff, nf
    gpu_lr_site_ok = .true.
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       if (V_meta(1,n,nf) /= 1) then
          gpu_lr_site_ok = .false.; return
       endif
    enddo
  end function gpu_lr_site_ok

  Integer function gpu_lr_pending()
    Implicit none
    gpu_lr_pending = 0
    if (inited .and. mode /= LR_STD) gpu_lr_pending = maxval(jc)
  end function gpu_lr_pending

!====================================================================
!  Device kernels
!====================================================================
!> delayed wrap: G0 as k_site_wrap (row r, column r), plus X(r,1:jc) *= za, Y(r,1:jc) *= zb
  Subroutine k_dly_wrap(G, X, Y, ld, nd_, jcnt, n, t, dir, r, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nd_, jcnt, n, t, dir, r, eo, nt
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*), X(ld,*), Y(ld,*)
    Integer :: j, sp, nloop
    Complex (Kind=Kind(0.d0)) :: za, zb
    nloop = max(nd_, jcnt)
    !$omp target teams distribute parallel do private(sp, za, zb)
    do j = 1, nloop
       sp = nint(real(nsigma_d(n,t)))
       za = V_exp(eo + 1 + (dir*sp + nt))
       zb = V_exp(eo + 1 + (-dir*sp + nt))
       if (j <= nd_) then
          if (j == r) then
             G(r,r) = zb * ( za * G(r,r) )
          else
             G(r,j) = za * G(r,j)
             G(j,r) = zb * G(j,r)
          endif
       endif
       if (j <= jcnt) then
          X(r,j) = za * X(r,j)
          Y(r,j) = zb * Y(r,j)
       endif
    enddo
  end subroutine k_dly_wrap

!> res2 = G(p,p) - X(p,1:jcnt) . Y(p,1:jcnt)   (one thread, result to the host)
  Subroutine k_dly_gpp(G, X, Y, ld, jcnt, p, res)
    Implicit none
    Integer, intent(in) :: ld, jcnt, p
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*), X(ld,*), Y(ld,*)
    Real (Kind=Kind(0.d0)), intent(inout) :: res(*)
    Integer :: k
    Complex (Kind=Kind(0.d0)) :: s
    !$omp target private(s)
    s = G(p,p)
    do k = 1, jcnt
       s = s - X(p,k) * Y(p,k)
    enddo
    res(1) = real(s, kind(0.d0))
    res(2) = aimag(s)
    !$omp end target
  end subroutine k_dly_gpp

!> delayed accept: column j1 = jcnt+1 of X and Y.  c = G(:,p) - X Y(p,:)^T, r = G(p,:) - X(p,:) Y^T,
!> X(:,j1) = z c, Y(:,j1) = e_p - r.  z = (zr, zi) passed as two reals (nvfortran complex-scalar dummy issue, ticket 205).
  Subroutine k_dly_accept(G, X, Y, ld, nrow, jcnt, p, zr, zi)
    Implicit none
    Integer, intent(in) :: ld, nrow, jcnt, p
    Real (Kind=Kind(0.d0)), intent(in) :: zr, zi
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: X(ld,*), Y(ld,*)
    Integer :: i, k
    Complex (Kind=Kind(0.d0)) :: c, rr, z
    !$omp target teams distribute parallel do private(k, c, rr, z)
    do i = 1, nrow
       z = cmplx(zr, zi, kind(0.d0))
       c  = G(i,p)
       rr = G(p,i)
       do k = 1, jcnt
          c  = c  - X(i,k) * Y(p,k)
          rr = rr - X(p,k) * Y(i,k)
       enddo
       X(i,jcnt+1) = z * c
       Y(i,jcnt+1) = - rr
       if (i == p) Y(i,jcnt+1) = Y(i,jcnt+1) + 1.d0
    enddo
  end subroutine k_dly_accept

!> submatrix gather for a proposal at p: buf(1) = G(p,p), buf(1+m) = G(S(m),p), buf(1+k+m) = G(p,S(m))
  Subroutine k_sub_gather(G, ld, S, k, p, buf)
    Implicit none
    Integer, intent(in) :: ld, k, p
    Integer, intent(in) :: S(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: buf(*)
    Integer :: m
    !$omp target teams distribute parallel do
    do m = 0, k
       if (m == 0) then
          buf(1) = G(p,p)
       else
          buf(1+m)   = G(S(m),p)
          buf(1+k+m) = G(p,S(m))
       endif
    enddo
  end subroutine k_sub_gather

!> flush gathers: C(:,m) = G(:,S(m)),  Qt(:,m) = e_{S(m)} - G(S(m),:)^T
  Subroutine k_sub_flush_gather(G, ld, nrow, S, k, C, Qt)
    Implicit none
    Integer, intent(in) :: ld, nrow, k
    Integer, intent(in) :: S(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: C(ld,*), Qt(ld,*)
    Integer :: i, m
    !$omp target teams distribute parallel do collapse(2)
    do m = 1, k
       do i = 1, nrow
          C(i,m)  = G(i,S(m))
          Qt(i,m) = - G(S(m),i)
          if (i == S(m)) Qt(i,m) = Qt(i,m) + 1.d0
       enddo
    enddo
  end subroutine k_sub_flush_gather

!====================================================================
!  Host orchestration
!====================================================================
  !> Flush the pending update of flavour nf into G_d(:,:,nf)
  Subroutine gpu_lr_flush(nf, reason)
    Implicit none
    Integer, intent(in) :: nf
    Integer, intent(in) :: reason          ! 1 capacity, 2 forced
    Integer :: k, m
    Complex (Kind=Kind(0.d0)), parameter :: cone = (1.d0,0.d0), cmone = (-1.d0,0.d0), czero = (0.d0,0.d0)
    if (.not. inited .or. mode == LR_STD) return
    k = jc(nf)
    if (k == 0) return
    n_flush(reason) = n_flush(reason) + 1
    n_flush(3) = n_flush(3) + 1
    if (mode == LR_DELAYED) then
       call la_gemm(LA_N, LA_T, Ndim, Ndim, k, cmone, X_d(1,1,nf), Ndim, Y_d(1,1,nf), Ndim, cone, G_d(1,1,nf), Ndim)
       upd_cnt(4) = upd_cnt(4) + 1
    else
       ! S, M to the device (k ints, k*k complex), gather C = G0(:,S), Qt = E_S - G0(S,:)^T, Y = Qt M^T, G0 -= C Y^T
       S_d(1:k) = S_h(1:k,nf)
       M_d(1:k,1:k) = M_h(1:k,1:k,nf)
       !$omp target update to(S_d(1:k))
       do m = 1, k
          !$omp target update to(M_d(1:k,m))
       enddo
       call k_sub_flush_gather(G_d(1,1,nf), Ndim, Ndim, S_d, k, X_d(1,1,nf), scr_c(1,1,1))
       call la_gemm(LA_N, LA_T, Ndim, k, k, cone, scr_c(1,1,1), Ndim, M_d, nd, czero, scr_c(1,1,2), Ndim)
       call la_gemm(LA_N, LA_T, Ndim, Ndim, k, cmone, X_d(1,1,nf), Ndim, scr_c(1,1,2), Ndim, cone, G_d(1,1,nf), Ndim)
       upd_cnt(4) = upd_cnt(4) + 3
       do m = 1, k
          slot_h(S_h(m,nf),nf) = 0
          S_h(m,nf) = 0
       enddo
       M_h(1:k,1:k,nf) = czero
    endif
    jc(nf) = 0
  end subroutine gpu_lr_flush

  Subroutine gpu_lr_flush_all(reason)
    Implicit none
    Integer, intent(in) :: reason
    Integer :: nf_eff
    if (.not. inited .or. mode == LR_STD) return
    do nf_eff = 1, N_FL_eff
       call gpu_lr_flush(Calc_Fl_map(nf_eff), reason)
    enddo
  end subroutine gpu_lr_flush_all

!--------------------------------------------------------------------
!> Wrap of single-site vertex n on flavour nf (Op_Wrapup, dir = +1 / Op_Wrapdo, dir = -1), pending state included.
  Subroutine gpu_lr_wrap(nf, n, t, dir)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    Integer :: r, m, sp, eo, nt
    Complex (Kind=Kind(0.d0)) :: za, zb
    r  = Op_V(n,nf)%P(1)
    eo = V_meta(3,n,nf)
    nt = V_type(n)
    if (mode == LR_DELAYED .and. jc(nf) > 0) then
       call k_dly_wrap(G_d(1,1,nf), X_d(1,1,nf), Y_d(1,1,nf), Ndim, Ndim, jc(nf), n, t, dir, r, eo, nt)
       upd_cnt(3) = upd_cnt(3) + 1
       return
    endif
    call k_wrap_g0(G_d(1,1,nf), Ndim, n, t, dir, r, eo, nt)
    upd_cnt(3) = upd_cnt(3) + 1
    if (mode == LR_SUBMATRIX .and. jc(nf) > 0) then
       m = slot_h(r,nf)
       if (m /= 0) then
          ! host copy of the field (nsigma_d is kept in step with the device copy): same factors as the kernel
          sp = nint(real(nsigma_d(n,t)))
          za = V_exp(eo + 1 + (dir*sp + nt))
          zb = V_exp(eo + 1 + (-dir*sp + nt))
          M_h(m,1:jc(nf),nf) = za * M_h(m,1:jc(nf),nf)
          M_h(1:jc(nf),m,nf) = zb * M_h(1:jc(nf),m,nf)
       endif
    endif
  end subroutine gpu_lr_wrap

!> G0 wrap kernel, same as gpu_update_mod's k_site_wrap
  Subroutine k_wrap_g0(G, ld, n, t, dir, r, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, n, t, dir, r, eo, nt
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Integer :: j, sp
    Complex (Kind=Kind(0.d0)) :: za, zb
    !$omp target teams distribute parallel do private(sp, za, zb)
    do j = 1, ld
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
  end subroutine k_wrap_g0

!--------------------------------------------------------------------
!> G(p,p) of the current (G0 + pending) Green function of flavour nf, for the ratio of a proposal at site p
  Subroutine gpu_lr_gpp(nf, p, gpp)
    Implicit none
    Integer, intent(in) :: nf, p
    Complex (Kind=Kind(0.d0)), intent(out) :: gpp
    Integer :: k, a, m, i
    Complex (Kind=Kind(0.d0)) :: s, wa
    n_prop = n_prop + 1
    if (mode == LR_DELAYED) then
       call k_dly_gpp(G_d(1,1,nf), X_d(1,1,nf), Y_d(1,1,nf), Ndim, jc(nf), p, res2)
       !$omp target update from(res2)
       upd_cnt(5) = upd_cnt(5) + 1
       gpp = cmplx(res2(1), res2(2), kind(0.d0))
       return
    endif
    ! submatrix
    k = jc(nf)
    S_d(1:k) = S_h(1:k,nf)
    if (k > 0) then
       !$omp target update to(S_d(1:k))
    endif
    call k_sub_gather(G_d(1,1,nf), Ndim, S_d, k, p, gbuf)
    upd_cnt(5) = upd_cnt(5) + 1
    !$omp target update from(gbuf(1:2*k+1))
    gpp = gbuf(1)
    prop_p(nf) = p; prop_k(nf) = k
    if (k == 0) return
    ! W = M (J - G0(S,p)),  G = G0(p,p) - G0(p,S) W
    i = slot_h(p,nf)
    do a = 1, k
       s = czero_c
       do m = 1, k
          s = s - M_h(a,m,nf) * gbuf(1+m)
       enddo
       if (i /= 0) s = s + M_h(a,i,nf)
       W_h(a,nf) = s
       grow_h(a,nf) = gbuf(1+k+a)
    enddo
    s = czero_c
    do a = 1, k
       s = s + grow_h(a,nf) * W_h(a,nf)
    enddo
    gpp = gpp - s
  end subroutine gpu_lr_gpp

!--------------------------------------------------------------------
!> Accept at site p in flavour nf with z = Delta/d.  Must follow gpu_lr_gpp for the same (nf, p).
  Subroutine gpu_lr_accept(nf, p, z)
    Implicit none
    Integer, intent(in) :: nf, p
    Complex (Kind=Kind(0.d0)), intent(in) :: z
    Integer :: j, k, a, b, i, kp
    Complex (Kind=Kind(0.d0)) :: av(nd+1), bt(nd+1)
    if (nf == 1) n_acc = n_acc + 1
    if (mode == LR_DELAYED) then
       if (jc(nf) + 1 > nd) call gpu_lr_flush(nf, 1)
       j = jc(nf)
       call k_dly_accept(G_d(1,1,nf), X_d(1,1,nf), Y_d(1,1,nf), Ndim, Ndim, j, p, real(z), aimag(z))
       upd_cnt(4) = upd_cnt(4) + 1
       jc(nf) = j + 1
       if (jc(nf) >= nd) call gpu_lr_flush(nf, 1)
       return
    endif
    ! submatrix, host only
    if (slot_h(p,nf) == 0 .and. jc(nf) + 1 > nd) call gpu_lr_flush(nf, 1)
    k = jc(nf)
    i = slot_h(p,nf)
    kp = k
    if (i == 0) then
       kp = k + 1
       S_h(kp,nf) = p
       slot_h(p,nf) = kp
       M_h(1:kp,kp,nf) = czero_c
       M_h(kp,1:kp,nf) = czero_c
       i = kp
    endif
    ! a = J - [W; 0],  b^T = J^T + [G0(p,S) M, 0]  (k = 0 after a flush: a = b = e_i; W, G0(p,S) are then not used)
    av(1:kp) = czero_c
    bt(1:kp) = czero_c
    if (k > 0 .and. prop_k(nf) == k .and. prop_p(nf) == p) then
       do a = 1, k
          av(a) = - W_h(a,nf)
       enddo
       do b = 1, k
          do a = 1, k
             bt(b) = bt(b) + grow_h(a,nf) * M_h(a,b,nf)
          enddo
       enddo
    endif
    av(i) = av(i) + 1.d0
    bt(i) = bt(i) + 1.d0
    do b = 1, kp
       do a = 1, kp
          M_h(a,b,nf) = M_h(a,b,nf) + z * av(a) * bt(b)
       enddo
    enddo
    jc(nf) = kp
    if (kp >= nd) call gpu_lr_flush(nf, 1)
  end subroutine gpu_lr_accept

end Module gpu_lowrank_mod
#endif

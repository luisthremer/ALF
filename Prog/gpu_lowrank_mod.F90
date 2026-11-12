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
  Use omp_lib, only: omp_get_wtime
  Use runtime_error_mod
  Use Files_mod, only: str_to_upper
  Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map, Op_V, nsigma
  Use QMC_runtime_var, only: get_Sequential_scheme, get_N_delay
  Use alf_la_mod, only: la_gemm, la_sync, LA_N, LA_T
#if defined(ALF_LA_HIP)
  Use alf_hipk_mod
#endif
  Use gpu_state_mod
#if defined(ALF_LA_ONEMKL)
  Use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_intptr_t, c_double, c_loc, c_f_pointer, c_associated, c_null_ptr
#endif

  Implicit none

  Private
  Public :: gpu_lr_mode, gpu_lr_init, gpu_lr_finalize, gpu_lr_site_ok, gpu_lr_pending, gpu_lr_flush_all, gpu_lr_flush
  Public :: gpu_lr_wrap, gpu_lr_gpp, gpu_lr_accept
  Public :: gpu_lr_blockmode, gpu_lr_blk_enter, gpu_lr_blk_leave
#ifdef ALF_TEAM_TIMERS
  Public :: p3t, p3n, p3_now
#endif
  Public :: gpu_lr_fast, gpu_lr_step_up, gpu_lr_step_down, gpu_lr_defer, gpu_lr_commit, gpu_lr_sync

  Integer, parameter :: LR_STD = 0, LR_DELAYED = 1, LR_SUBMATRIX = 2
  Complex (Kind=Kind(0.d0)), parameter :: czero_c = (0.d0,0.d0)

#ifdef ALF_TEAM_TIMERS
  ! Phase 3b timers: 1 blk_begin total, 2 its device part (to bsd, gather, from gbd), 3 blk_end total, 4 its la_sync, 5 its device part
  ! (uploads, scale, flushes), 6 slice level (g_in, hopping wraps, final flush, g_out); p3n(1) blocks
  Real (Kind=Kind(0.d0)), Save :: p3t(0:7) = 0.d0
  Integer(int64), Save :: p3n(0:7) = 0
#endif
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
  ! ---- fast path (all flavours in one launch, no host->device transfer per accept) ----
  Logical, Save :: fast = .false., stepped = .false.
  Integer, Save :: gs = 0                                      ! gbuf stride per flavour (2 nd + 1)
  Integer, Allocatable, Save :: slot_d(:,:)                    ! (Ndim, N_FL) device: site -> position in S (kept by the kernels)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: M_pk(:)      ! (nd*nd) device/host packed M of one flavour (flush)
  ! deferred work, applied by the next launch (arguments of the kernel, not transfers)
  Logical, Save :: dfr_wrap = .false.
  Integer, Save :: dfr_n = 0, dfr_t = 0, dfr_rw(2) = 0, dfr_pp(2) = 0, dfr_pk(2) = 0
  Real (Kind=Kind(0.d0)), Save :: dfr_vr = 0.d0, dfr_vi = 0.d0
  Complex (Kind=Kind(0.d0)), Save :: dfr_za(2), dfr_zb(2)
  Integer, Save :: dfr_pn(2) = 0                                ! gather site of the step per flavour
  Integer, Save :: acc_p(2) = 0, acc_flag(2) = 0                ! delayed: accepts collected for gpu_lr_commit
  Complex (Kind=Kind(0.d0)), Save :: acc_z(2)
  Integer(int64), Save :: n_fused = 0, n_blocks = 0
  ! block mode (submatrix): host copy of the principal block G0(B,B) per flavour
  Logical, Save :: blk_open = .false.
  Integer, Save :: blk_v1 = 0, blk_dv = 1, blk_b = 0, blk_cnt = 0, blk_a = 0
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: Gb(:,:,:), zab_h(:,:,:), gbd(:,:,:), cbuf(:,:)
  Integer, Allocatable, Save :: bidx(:,:), bs_h(:,:), bsd(:,:), bS_d(:,:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: W_h(:,:), grow_h(:,:)   ! cached proposal data (nd, N_FL) host, per flavour
  Integer, Allocatable, Save :: prop_p(:), prop_k(:)          ! site and pending count of the cached proposal, per flavour
#if defined(ALF_LA_ONEMKL)
  ! Asynchronous block path (env ALFSH_ASYNCBLK=1, read once in ab_setup): the device sequence of a block as SYCL kernels +
  ! oneMKL GEMM on the interop queue (libalfsb, Libraries/Modules/alf_la_sycl_blk.cpp), no wait except the D2H of the block gather
  Logical, Save :: ab_on = .false., ab_pending = .false.
  Type (c_ptr), Save :: ab_pG, ab_pX, ab_pS, ab_pBs, ab_pBS2, ab_pCb, ab_pSl, ab_pNs, ab_pGb, ab_ph
  Complex (Kind=Kind(0.d0)), Pointer, Save :: gbh(:,:,:) => null()      ! host-USM landing buffer of the block gather
  Interface
    Type (c_ptr) Function alfsh_get_queue() Bind(C, name='alfsh_get_queue')
      Import :: c_ptr
    End Function
    Integer (c_int) Function alfsb_init(q) Bind(C, name='alfsb_init')
      Import :: c_ptr, c_int
      Type (c_ptr), Value :: q
    End Function
    Type (c_ptr) Function alfsb_alloc_host(nbytes) Bind(C, name='alfsb_alloc_host')
      Import :: c_ptr, c_int64_t
      Integer (c_int64_t), Value :: nbytes
    End Function
    Integer (c_int) Function alfsb_wait() Bind(C, name='alfsb_wait')
      Import :: c_int
    End Function
    Integer (c_int) Function alfsb_copy(dst, src, nbytes) Bind(C, name='alfsb_copy')
      Import :: c_ptr, c_int, c_int64_t
      Type (c_ptr), Value :: dst, src
      Integer (c_int64_t), Value :: nbytes
    End Function
    Integer (c_int) Function alfsb_gather(nfl, nfa, nfb, ld, b, nd, gbk, bsk, Gk) Bind(C, name='alfsb_gather')
      Import :: c_ptr, c_int
      Integer (c_int), Value :: nfl, nfa, nfb, ld, b, nd
      Type (c_ptr), Value :: gbk, bsk, Gk
    End Function
    Integer (c_int) Function alfsb_scale(nfl, nfa, nfb, ld, b, nd, ldc, Gk, bsk, cbk) Bind(C, name='alfsb_scale')
      Import :: c_ptr, c_int
      Integer (c_int), Value :: nfl, nfa, nfb, ld, b, nd, ldc
      Type (c_ptr), Value :: Gk, bsk, cbk
    End Function
    Integer (c_int) Function alfsb_fgather(G, ld, nrow, S, k, nf, C, Qt, slot) Bind(C, name='alfsb_fgather')
      Import :: c_ptr, c_int
      Integer (c_int), Value :: ld, nrow, k, nf
      Type (c_ptr), Value :: G, S, C, Qt, slot
    End Function
    Integer (c_int) Function alfsb_gemm(ta, tb, m, n, k, ar, ai, A, lda, B, ldb, br, bi, C, ldc) Bind(C, name='alfsb_gemm')
      Import :: c_ptr, c_int, c_double
      Integer (c_int), Value :: ta, tb, m, n, k, lda, ldb, ldc
      Real (c_double), Value :: ar, ai, br, bi
      Type (c_ptr), Value :: A, B, C
    End Function
  End Interface
#endif

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
       if (mode == LR_SUBMATRIX) nd = max(nd, min(32, n_nz))     ! GPU: larger blocks amortise the device round trips
    endif
    nd = max(1, min(nd, Ndim))
    allocate( jc(N_FL), S_h(nd,N_FL), slot_h(Ndim,N_FL), M_h(nd,nd,N_FL), W_h(nd,N_FL), grow_h(nd,N_FL), prop_p(N_FL), prop_k(N_FL) )
    prop_p = 0; prop_k = -1
    jc = 0; S_h = 0; slot_h = 0; M_h = cmplx(0.d0,0.d0,kind(0.d0))
    gs = 2*nd + 1
    allocate( gbuf(N_FL*gs), res2(2), S_d(nd), M_d(nd,nd), slot_d(Ndim,N_FL), M_pk(nd*nd) )
    slot_d = 0; M_pk = cmplx(0.d0,0.d0,kind(0.d0))
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
    !$omp target enter data map(alloc: gbuf, res2, S_d, M_d, M_pk)
    !$omp target enter data map(to: slot_d)
    call fast_check()
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
#if defined(ALF_LA_ONEMKL)
    call ab_drain()
#endif
    if (mode /= LR_STD) then
       write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') 'GPU low-rank update counts: proposals ', n_prop, ', accepted ', n_acc, &
            & ', fused steps ', n_fused, ', blocks ', n_blocks, ', flushes ', n_flush(3), ' (capacity ', n_flush(1), ', forced ', n_flush(2)
       call info_line(trim(line)//')')
       !$omp target exit data map(delete: gbuf, res2, S_d, M_d, M_pk, slot_d)
       !$omp target exit data map(delete: X_d)
       if (mode == LR_DELAYED) then
          !$omp target exit data map(delete: Y_d)
       endif
       if (allocated(Gb)) then
          !$omp target exit data map(delete: bsd, bS_d, gbd, cbuf)
          deallocate( Gb, bidx, bs_h, zab_h, bsd, bS_d, gbd, cbuf )
       endif
       deallocate( jc, S_h, slot_h, M_h, W_h, grow_h, prop_p, prop_k, gbuf, res2, S_d, M_d, X_d, Y_d, slot_d, M_pk )
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,1), X(ld,1), Y(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,1), X(ld,1), Y(ld,1)
    Real (Kind=Kind(0.d0)), intent(inout) :: res(2)
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
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: X(ld,1), Y(ld,1)
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

!> submatrix gather for a proposal at p: buf(1) = G(p,p), buf(1+m) = G(S(m),p), buf(1+nd+m) = G(p,S(m))
  Subroutine k_sub_gather(G, ld, S, k, nd_, p, buf)
    Implicit none
    Integer, intent(in) :: ld, k, nd_, p
    Integer, intent(in) :: S(nd_)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,ld)
    Complex (Kind=Kind(0.d0)), intent(inout) :: buf(2*nd_+1)
    Integer :: m
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do map(present,alloc: G, S, buf)
#else
    !$omp target teams distribute parallel do map(alloc: G, S, buf)
#endif
    do m = 0, k
       if (m == 0) then
          buf(1) = G(p,p)
       else
          buf(1+m)     = G(S(m),p)
          buf(1+nd_+m) = G(p,S(m))
       endif
    enddo
  end subroutine k_sub_gather

!> flush gathers: C(:,m) = G(:,S(m)),  Qt(:,m) = e_{S(m)} - G(S(m),:)^T;  slot_d of the flushed sites is cleared
  Subroutine k_sub_flush_gather(G, ld, nrow, S, k, nf, C, Qt, slot, nfg)
    Implicit none
    Integer, intent(in) :: ld, nrow, k, nf, nfg
    Integer, intent(in) :: S(k)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,ld)
    Complex (Kind=Kind(0.d0)), intent(inout) :: C(ld,k), Qt(ld,k)
    Integer, intent(inout) :: slot(ld,nfg)
    Integer :: i, m
#if defined(ALF_LA_HIP)
    call hk_flush_gather(ld, nrow, k, nf, size(slot,1), G(1,1), S, C(1,1), Qt(1,1), slot(1,1))
#else
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do collapse(2) map(present,alloc: G, S, C, Qt, slot)
#else
    !$omp target teams distribute parallel do collapse(2) map(alloc: G, S, C, Qt, slot)
#endif
    do m = 1, k
       do i = 1, nrow
          C(i,m)  = G(i,S(m))
          Qt(i,m) = - G(S(m),i)
          if (i == S(m)) then
             Qt(i,m) = Qt(i,m) + 1.d0
             slot(i,nf) = 0
          endif
       enddo
    enddo
#endif
  end subroutine k_sub_flush_gather

!====================================================================
!  Host orchestration
!====================================================================
  !> Flush the pending update of flavour nf into G_d(:,:,nf)
  Subroutine gpu_lr_flush(nf, reason)
    Implicit none
    Integer, intent(in) :: nf
    Integer, intent(in) :: reason          ! 1 capacity, 2 forced
    Integer :: k, m, b
    Complex (Kind=Kind(0.d0)), parameter :: cone = (1.d0,0.d0), cmone = (-1.d0,0.d0), czero = (0.d0,0.d0)
    if (.not. inited .or. mode == LR_STD) return
#if defined(ALF_LA_ONEMKL)
    call ab_drain()
#endif
    k = jc(nf)
    if (k == 0) return
    if (fast) call gpu_lr_sync()                    ! deferred wrap / field / slot work first: G0 must be current
    n_flush(reason) = n_flush(reason) + 1
    n_flush(3) = n_flush(3) + 1
    if (mode == LR_DELAYED) then
       call la_gemm(LA_N, LA_T, Ndim, Ndim, k, cmone, X_d(1,1,nf), Ndim, Y_d(1,1,nf), Ndim, cone, G_d(1,1,nf), Ndim)
       upd_cnt(4) = upd_cnt(4) + 1
    else
       ! S, M (packed k x k) to the device, gather C = G0(:,S), Qt = E_S - G0(S,:)^T, Y = Qt M^T, G0 -= C Y^T
       call la_sync()                               ! asynchronous library calls of a previous flush still read S_d / M_pk
       S_d(1:k) = S_h(1:k,nf)
       !$omp target update to(S_d(1:k))
       do m = 1, k
          M_pk((m-1)*k+1:m*k) = M_h(1:k,m,nf)
       enddo
       !$omp target update to(M_pk(1:k*k))
       call k_sub_flush_gather(G_d(1,1,nf), Ndim, Ndim, S_d, k, nf, X_d(1,1,nf), scr_c(1,1,1), slot_d, size(slot_d,2))
       call la_gemm(LA_N, LA_T, Ndim, k, k, cone, scr_c(1,1,1), Ndim, M_pk, k, czero, scr_c(1,1,2), Ndim)
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
#if defined(ALF_LA_ONEMKL)
    call ab_drain()
#endif
    if (fast) call gpu_lr_sync()
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,1)
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
  Subroutine gpu_lr_gpp(nf, nf_eff, p, gpp)
    Implicit none
    Integer, intent(in) :: nf, nf_eff, p
    Complex (Kind=Kind(0.d0)), intent(out) :: gpp
    Integer :: k, a, m, i, base
    Complex (Kind=Kind(0.d0)) :: s
    n_prop = n_prop + 1
    base = 0
    if (blk_open) then
       call blk_gpp(nf, nf_eff, p, gpp)
       return
    endif
    if (fast .and. stepped) then
       base = (nf_eff-1)*gs
       if (mode == LR_DELAYED) then
          gpp = gbuf(base+1)                 ! G(p,p) of G0 + pending, from the fused step kernel
          return
       endif
    else
       if (mode == LR_DELAYED) then
          call k_dly_gpp(G_d(1,1,nf), X_d(1,1,nf), Y_d(1,1,nf), Ndim, jc(nf), p, res2)
          !$omp target update from(res2)
          upd_cnt(5) = upd_cnt(5) + 1
          gpp = cmplx(res2(1), res2(2), kind(0.d0))
          return
       endif
       k = jc(nf)
       S_d(1:k) = S_h(1:k,nf)
       if (k > 0) then
          !$omp target update to(S_d(1:k))
       endif
       call k_sub_gather(G_d(1,1,nf), Ndim, S_d, k, nd, p, gbuf)
       upd_cnt(5) = upd_cnt(5) + 1
       !$omp target update from(gbuf(1:2*nd+1))
    endif
    ! submatrix: G(p,p) = G0(p,p) - G0(p,S) W,  W = M (J - G0(S,p))   (host, O(k^2))
    k = jc(nf)
    gpp = gbuf(base+1)
    prop_p(nf) = p; prop_k(nf) = k
    if (k == 0) return
    i = slot_h(p,nf)
    do a = 1, k
       s = czero_c
       do m = 1, k
          s = s - M_h(a,m,nf) * gbuf(base+1+m)
       enddo
       if (i /= 0) s = s + M_h(a,i,nf)
       W_h(a,nf) = s
       grow_h(a,nf) = gbuf(base+1+nd+a)
    enddo
    s = czero_c
    do a = 1, k
       s = s + grow_h(a,nf) * W_h(a,nf)
    enddo
    gpp = gpp - s
  end subroutine gpu_lr_gpp

!--------------------------------------------------------------------
!> Accept at site p in flavour nf with z = Delta/d.  Must follow gpu_lr_gpp for the same (nf, p).
  Subroutine gpu_lr_accept(nf, nf_eff, p, z)
    Implicit none
    Integer, intent(in) :: nf, nf_eff, p
    Complex (Kind=Kind(0.d0)), intent(in) :: z
    Integer :: j, k, a, b, i, kp
    Logical :: padded
    Complex (Kind=Kind(0.d0)) :: av(nd+1), bt(nd+1)
    if (nf_eff == 1) n_acc = n_acc + 1
    if (mode == LR_DELAYED) then
       if (fast) then
          acc_p(nf_eff) = p; acc_z(nf_eff) = z; acc_flag(nf_eff) = 1     ! launched together by gpu_lr_commit
          return
       endif
       if (jc(nf) + 1 > nd) call gpu_lr_flush(nf, 1)
       j = jc(nf)
       call k_dly_accept(G_d(1,1,nf), X_d(1,1,nf), Y_d(1,1,nf), Ndim, Ndim, j, p, real(z), aimag(z))
       upd_cnt(4) = upd_cnt(4) + 1
       jc(nf) = j + 1
       if (jc(nf) >= nd) call gpu_lr_flush(nf, 1)
       return
    endif
    ! submatrix, host only
    if (blk_open .and. slot_h(p,nf) == 0 .and. jc(nf) + 1 > nd) then
       write(error_unit,*) 'gpu_lr_accept: pending sites exceed N_delay inside a block'
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    endif
    if (slot_h(p,nf) == 0 .and. jc(nf) + 1 > nd) call gpu_lr_flush(nf, 1)
    k = jc(nf)
    i = slot_h(p,nf)
    kp = k
    padded = .false.
    if (i == 0) then
       kp = k + 1
       S_h(kp,nf) = p
       slot_h(p,nf) = kp
       M_h(1:kp,kp,nf) = czero_c
       M_h(kp,1:kp,nf) = czero_c
       i = kp
       padded = .true.
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
    if (blk_open) then
       continue                          ! block mode: flushed by blk_end
    elseif (kp >= nd) then
       call gpu_lr_flush(nf, 1)          ! (fast: syncs the deferred work first; this accept's pad is dropped, the flush clears slot_d)
    elseif (fast .and. padded .and. .not. blk_open) then
       dfr_pp(nf_eff) = p; dfr_pk(nf_eff) = kp           ! slot_d(p) := kp, applied by the next launch
    endif
  end subroutine gpu_lr_accept


!====================================================================
!  Fast path: all flavours in ONE launch per site, the gathered data copied back by the same region (one
!  synchronisation per site), no host->device transfer per accept.  Work that does not need a decision is deferred and
!  passed as scalar arguments of the next launch: the field write nsigma_d(n,t), the slot_d entry of a newly pending site,
!  and (WRAPGRDO) the wrap of the site just updated.
!  Conditions (fast_check): N_FL_eff <= 2, every vertex a single-site operator in every computed flavour, and the site
!  P(1) of different vertices distinct (so a deferred wrap never touches the next gather site).
!====================================================================
  Subroutine fast_check()
    Implicit none
    Logical, Allocatable :: seen(:,:)
    Integer :: n, fl, nf
    fast = (N_FL_eff <= 2)
    if (.not. fast) return
    allocate( seen(Ndim,N_FL) )
    seen = .false.
    do n = 1, size(Op_V,1)
       do fl = 1, N_FL_eff
          nf = Calc_Fl_map(fl)
          if (V_meta(1,n,nf) /= 1) then
             fast = .false.
          else
             if (seen(Op_V(n,nf)%P(1),nf)) fast = .false.
             seen(Op_V(n,nf)%P(1),nf) = .true.
          endif
       enddo
    enddo
    deallocate(seen)
  end subroutine fast_check

  Logical function gpu_lr_fast()
    Implicit none
    if (.not. inited) call gpu_lr_init()
    gpu_lr_fast = (mode /= LR_STD) .and. fast
  end function gpu_lr_fast

!--------------------------------------------------------------------
!> UP, submatrix: wrap of site r (old field), gather of G(r,r), G(S,r), G(r,S) of the wrapped G, deferred field / slot work.
  Subroutine k_sub_up(nfl, nfa, nfb, ld, nd_, gs_, n, t, nt, ra, rb, eoa, eob, dpa, dka, dpb, dkb, dn, dt, dvr, dvi, &
       &              nfg, ngb, nsn1, nsn2, nve, Gk, gb, slot, nsig, vexp)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, nd_, gs_, n, t, nt, ra, rb, eoa, eob, dpa, dka, dpb, dkb, dn, dt
    Integer, intent(in) :: nfg, ngb, nsn1, nsn2, nve
    Real (Kind=Kind(0.d0)), intent(in) :: dvr, dvi
    Complex (Kind=Kind(0.d0)), intent(inout) :: Gk(ld,ld,nfg), gb(ngb), nsig(nsn1,nsn2)
    Complex (Kind=Kind(0.d0)), intent(in) :: vexp(nve)
    Integer, intent(inout) :: slot(ld,nfg)
    Integer :: fl, j, nf, r, eo, sl, base, dp, dk, sp
    Complex (Kind=Kind(0.d0)) :: za, zb, gv
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do collapse(2) private(nf, r, eo, sl, base, dp, dk, sp, za, zb, gv) &
    !$omp&   map(present,alloc: Gk, gb, slot, nsig, vexp)
#else
    !$omp target teams distribute parallel do collapse(2) private(nf, r, eo, sl, base, dp, dk, sp, za, zb, gv) &
    !$omp&   map(alloc: Gk, gb, slot, nsig, vexp)
#endif
    do fl = 1, nfl
       do j = 1, ld
          if (fl == 1) then
             nf = nfa; r = ra; eo = eoa; dp = dpa; dk = dka
          else
             nf = nfb; r = rb; eo = eob; dp = dpb; dk = dkb
          endif
          base = (fl-1)*gs_
          sp = nint(real(nsig(n,t)))
          za = vexp(eo + 1 + (sp + nt))
          zb = vexp(eo + 1 + (-sp + nt))
          sl = slot(j,nf)
          if (dk > 0 .and. j == dp) then
             sl = dk
             slot(j,nf) = dk
          endif
          if (j == r) then
             gv = zb * ( za * Gk(r,r,nf) )
             Gk(r,r,nf) = gv
             gb(base+1) = gv
             if (sl > 0) then
                gb(base+1+sl) = gv
                gb(base+1+nd_+sl) = gv
             endif
          else
             Gk(r,j,nf) = za * Gk(r,j,nf)
             Gk(j,r,nf) = zb * Gk(j,r,nf)
             if (sl > 0) then
                gb(base+1+sl) = Gk(j,r,nf)
                gb(base+1+nd_+sl) = Gk(r,j,nf)
             endif
          endif
          if (fl == 1 .and. j == 1 .and. dn > 0) nsig(dn,dt) = cmplx(dvr, dvi, kind(0.d0))
       enddo
    enddo
    !$omp target update from(gb(1:nfl*gs_))
  end subroutine k_sub_up

!> DOWN, submatrix: deferred wrap of site rw (factors za_w, zb_w), gather at site pn of the wrapped G, deferred field / slot work.
  Subroutine k_sub_down(nfl, nfa, nfb, ld, nd_, gs_, wv, pna, pnb, rwa, rwb, &
       &                zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dpa, dka, dpb, dkb, dn, dt, dvr, dvi, &
       &                nfg, ngb, nsn1, nsn2, Gk, gb, slot, nsig)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, nd_, gs_, wv, pna, pnb, rwa, rwb, dpa, dka, dpb, dkb, dn, dt
    Integer, intent(in) :: nfg, ngb, nsn1, nsn2
    Complex (Kind=Kind(0.d0)), intent(inout) :: Gk(ld,ld,nfg), gb(ngb), nsig(nsn1,nsn2)
    Integer, intent(inout) :: slot(ld,nfg)
    Real (Kind=Kind(0.d0)), intent(in) :: zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dvr, dvi
    Integer :: fl, j, nf, pn, rw, sl, slrw, base, dp, dk
    Complex (Kind=Kind(0.d0)) :: za, zb
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do collapse(2) private(nf, pn, rw, sl, slrw, base, dp, dk, za, zb) &
    !$omp&   map(present,alloc: Gk, gb, slot, nsig)
#else
    !$omp target teams distribute parallel do collapse(2) private(nf, pn, rw, sl, slrw, base, dp, dk, za, zb) &
    !$omp&   map(alloc: Gk, gb, slot, nsig)
#endif
    do fl = 1, nfl
       do j = 1, ld
          if (fl == 1) then
             nf = nfa; pn = pna; rw = rwa; dp = dpa; dk = dka
             za = cmplx(zar, zai, kind(0.d0)); zb = cmplx(zbr, zbi, kind(0.d0))
          else
             nf = nfb; pn = pnb; rw = rwb; dp = dpb; dk = dkb
             za = cmplx(zcr, zci, kind(0.d0)); zb = cmplx(zdr, zdi, kind(0.d0))
          endif
          base = (fl-1)*gs_
          sl = slot(j,nf)
          if (dk > 0 .and. j == dp) then
             sl = dk
             slot(j,nf) = dk
          endif
          if (wv /= 0) then
             if (j == rw) then
                Gk(rw,rw,nf) = zb * ( za * Gk(rw,rw,nf) )
             else
                Gk(rw,j,nf) = za * Gk(rw,j,nf)
                Gk(j,rw,nf) = zb * Gk(j,rw,nf)
             endif
          endif
          if (j == pn) then
             gb(base+1) = Gk(pn,pn,nf)
             if (sl > 0) then
                gb(base+1+sl) = Gk(pn,pn,nf)
                gb(base+1+nd_+sl) = Gk(pn,pn,nf)
             endif
             if (wv /= 0) then
                ! the entries (rw,pn), (pn,rw) are written by this thread: gather them here (no read of data another thread writes)
                if (dk > 0 .and. rw == dp) then
                   slrw = dk
                else
                   slrw = slot(rw,nf)
                endif
                if (slrw > 0) then
                   gb(base+1+slrw) = Gk(rw,pn,nf)
                   gb(base+1+nd_+slrw) = Gk(pn,rw,nf)
                endif
             endif
          else
             if (sl > 0) then
                if (.not. (wv /= 0 .and. j == rw)) then
                   gb(base+1+sl) = Gk(j,pn,nf)
                   gb(base+1+nd_+sl) = Gk(pn,j,nf)
                endif
             endif
          endif
          if (fl == 1 .and. j == 1 .and. dn > 0) nsig(dn,dt) = cmplx(dvr, dvi, kind(0.d0))
       enddo
    enddo
    !$omp target update from(gb(1:nfl*gs_))
  end subroutine k_sub_down

!> apply the deferred work only (sync): wrap of rw, slot pad, field write
  Subroutine k_apply(nfl, nfa, nfb, ld, wv, rwa, rwb, zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dpa, dka, dpb, dkb, dn, dt, dvr, dvi, &
       &               nfg, nsn1, nsn2, Gk, slot, nsig)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, wv, rwa, rwb, dpa, dka, dpb, dkb, dn, dt
    Integer, intent(in) :: nfg, nsn1, nsn2
    Complex (Kind=Kind(0.d0)), intent(inout) :: Gk(ld,ld,nfg), nsig(nsn1,nsn2)
    Integer, intent(inout) :: slot(ld,nfg)
    Real (Kind=Kind(0.d0)), intent(in) :: zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dvr, dvi
    Integer :: fl, j, nf, rw, dp, dk
    Complex (Kind=Kind(0.d0)) :: za, zb
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do collapse(2) private(nf, rw, dp, dk, za, zb) &
    !$omp&   map(present,alloc: Gk, slot, nsig)
#else
    !$omp target teams distribute parallel do collapse(2) private(nf, rw, dp, dk, za, zb) &
    !$omp&   map(alloc: Gk, slot, nsig)
#endif
    do fl = 1, nfl
       do j = 1, ld
          if (fl == 1) then
             nf = nfa; rw = rwa; dp = dpa; dk = dka
             za = cmplx(zar, zai, kind(0.d0)); zb = cmplx(zbr, zbi, kind(0.d0))
          else
             nf = nfb; rw = rwb; dp = dpb; dk = dkb
             za = cmplx(zcr, zci, kind(0.d0)); zb = cmplx(zdr, zdi, kind(0.d0))
          endif
          if (dk > 0 .and. j == dp) slot(j,nf) = dk
          if (wv /= 0) then
             if (j == rw) then
                Gk(rw,rw,nf) = zb * ( za * Gk(rw,rw,nf) )
             else
                Gk(rw,j,nf) = za * Gk(rw,j,nf)
                Gk(j,rw,nf) = zb * Gk(j,rw,nf)
             endif
          endif
          if (fl == 1 .and. j == 1 .and. dn > 0) nsig(dn,dt) = cmplx(dvr, dvi, kind(0.d0))
       enddo
    enddo
  end subroutine k_apply

!> UP, delayed: wrap of site r (G0 row/column, rows r of X and Y), G(r,r) of G0 + pending, deferred field write
  Subroutine k_dly_up(nfl, nfa, nfb, ld, gs_, n, t, nt, ra, rb, eoa, eob, jca, jcb, dn, dt, dvr, dvi)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, gs_, n, t, nt, ra, rb, eoa, eob, jca, jcb, dn, dt
    Real (Kind=Kind(0.d0)), intent(in) :: dvr, dvi
    Integer :: fl, j, nf, r, eo, jcn, base, sp, k
    Complex (Kind=Kind(0.d0)) :: za, zb, gv, s, xk, yk
    !$omp target teams distribute parallel do collapse(2) private(nf, r, eo, jcn, base, sp, k, za, zb, gv, s, xk, yk)
    do fl = 1, nfl
       do j = 1, ld
          if (fl == 1) then
             nf = nfa; r = ra; eo = eoa; jcn = jca
          else
             nf = nfb; r = rb; eo = eob; jcn = jcb
          endif
          base = (fl-1)*gs_
          sp = nint(real(nsigma_d(n,t)))
          za = V_exp(eo + 1 + (sp + nt))
          zb = V_exp(eo + 1 + (-sp + nt))
          if (j == r) then
             gv = zb * ( za * G_d(r,r,nf) )
             G_d(r,r,nf) = gv
             s = gv
             do k = 1, jcn
                xk = za * X_d(r,k,nf)
                yk = zb * Y_d(r,k,nf)
                X_d(r,k,nf) = xk
                Y_d(r,k,nf) = yk
                s = s - xk * yk
             enddo
             gbuf(base+1) = s
          else
             G_d(r,j,nf) = za * G_d(r,j,nf)
             G_d(j,r,nf) = zb * G_d(j,r,nf)
          endif
          if (fl == 1 .and. j == 1 .and. dn > 0) nsigma_d(dn,dt) = cmplx(dvr, dvi, kind(0.d0))
       enddo
    enddo
    !$omp target update from(gbuf(1:(nfl-1)*gs_+1))
  end subroutine k_dly_up

!> DOWN, delayed: deferred wrap of rw, G(pn,pn) of G0 + pending, deferred field write
  Subroutine k_dly_down(nfl, nfa, nfb, ld, gs_, wv, pna, pnb, rwa, rwb, jca, jcb, &
       &                zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dn, dt, dvr, dvi)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, gs_, wv, pna, pnb, rwa, rwb, jca, jcb, dn, dt
    Real (Kind=Kind(0.d0)), intent(in) :: zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dvr, dvi
    Integer :: fl, j, nf, pn, rw, jcn, base, k
    Complex (Kind=Kind(0.d0)) :: za, zb, s
    !$omp target teams distribute parallel do collapse(2) private(nf, pn, rw, jcn, base, k, za, zb, s)
    do fl = 1, nfl
       do j = 1, ld
          if (fl == 1) then
             nf = nfa; pn = pna; rw = rwa; jcn = jca
             za = cmplx(zar, zai, kind(0.d0)); zb = cmplx(zbr, zbi, kind(0.d0))
          else
             nf = nfb; pn = pnb; rw = rwb; jcn = jcb
             za = cmplx(zcr, zci, kind(0.d0)); zb = cmplx(zdr, zdi, kind(0.d0))
          endif
          base = (fl-1)*gs_
          if (wv /= 0) then
             if (j == rw) then
                G_d(rw,rw,nf) = zb * ( za * G_d(rw,rw,nf) )
                do k = 1, jcn
                   X_d(rw,k,nf) = za * X_d(rw,k,nf)
                   Y_d(rw,k,nf) = zb * Y_d(rw,k,nf)
                enddo
             else
                G_d(rw,j,nf) = za * G_d(rw,j,nf)
                G_d(j,rw,nf) = zb * G_d(j,rw,nf)
             endif
          endif
          if (j == pn) then
             s = G_d(pn,pn,nf)
             do k = 1, jcn
                s = s - X_d(pn,k,nf) * Y_d(pn,k,nf)
             enddo
             gbuf(base+1) = s
          endif
          if (fl == 1 .and. j == 1 .and. dn > 0) nsigma_d(dn,dt) = cmplx(dvr, dvi, kind(0.d0))
       enddo
    enddo
    !$omp target update from(gbuf(1:(nfl-1)*gs_+1))
  end subroutine k_dly_down

!> apply deferred work, delayed: wrap of rw (G0, X, Y rows), field write
  Subroutine k_dly_apply(nfl, nfa, nfb, ld, wv, rwa, rwb, jca, jcb, zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dn, dt, dvr, dvi)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, wv, rwa, rwb, jca, jcb, dn, dt
    Real (Kind=Kind(0.d0)), intent(in) :: zar, zai, zbr, zbi, zcr, zci, zdr, zdi, dvr, dvi
    Integer :: fl, j, nf, rw, jcn, k
    Complex (Kind=Kind(0.d0)) :: za, zb
    !$omp target teams distribute parallel do collapse(2) private(nf, rw, jcn, k, za, zb)
    do fl = 1, nfl
       do j = 1, ld
          if (fl == 1) then
             nf = nfa; rw = rwa; jcn = jca
             za = cmplx(zar, zai, kind(0.d0)); zb = cmplx(zbr, zbi, kind(0.d0))
          else
             nf = nfb; rw = rwb; jcn = jcb
             za = cmplx(zcr, zci, kind(0.d0)); zb = cmplx(zdr, zdi, kind(0.d0))
          endif
          if (wv /= 0) then
             if (j == rw) then
                G_d(rw,rw,nf) = zb * ( za * G_d(rw,rw,nf) )
                do k = 1, jcn
                   X_d(rw,k,nf) = za * X_d(rw,k,nf)
                   Y_d(rw,k,nf) = zb * Y_d(rw,k,nf)
                enddo
             else
                G_d(rw,j,nf) = za * G_d(rw,j,nf)
                G_d(j,rw,nf) = zb * G_d(j,rw,nf)
             endif
          endif
          if (fl == 1 .and. j == 1 .and. dn > 0) nsigma_d(dn,dt) = cmplx(dvr, dvi, kind(0.d0))
       enddo
    enddo
  end subroutine k_dly_apply

!> accepts of both flavours in one launch (delayed): column jcn+1 of X and Y, as k_dly_accept
  Subroutine k_dly_accept2(nfl, nfa, nfb, ld, pa, pb, ja, jb, fa, fb, zar, zai, zbr, zbi)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, pa, pb, ja, jb, fa, fb
    Real (Kind=Kind(0.d0)), intent(in) :: zar, zai, zbr, zbi
    Integer :: fl, i, k, nf, p, jn, fg
    Complex (Kind=Kind(0.d0)) :: c, rr, z
    !$omp target teams distribute parallel do collapse(2) private(nf, p, jn, fg, k, c, rr, z)
    do fl = 1, nfl
       do i = 1, ld
          if (fl == 1) then
             nf = nfa; p = pa; jn = ja; fg = fa; z = cmplx(zar, zai, kind(0.d0))
          else
             nf = nfb; p = pb; jn = jb; fg = fb; z = cmplx(zbr, zbi, kind(0.d0))
          endif
          if (fg /= 0) then
             c  = G_d(i,p,nf)
             rr = G_d(p,i,nf)
             do k = 1, jn
                c  = c  - X_d(i,k,nf) * Y_d(p,k,nf)
                rr = rr - X_d(p,k,nf) * Y_d(i,k,nf)
             enddo
             X_d(i,jn+1,nf) = z * c
             Y_d(i,jn+1,nf) = - rr
             if (i == p) Y_d(i,jn+1,nf) = Y_d(i,jn+1,nf) + 1.d0
          endif
       enddo
    enddo
  end subroutine k_dly_accept2

!--------------------------------------------------------------------
!> wrap factors of vertex n, flavour nf, at slice t with the field currently in nsigma_d (host copy)
  Subroutine wrap_factors(nf, n, t, dir, za, zb)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    Complex (Kind=Kind(0.d0)), intent(out) :: za, zb
    Integer :: sp, eo, nt
    sp = nint(real(nsigma_d(n,t)))
    eo = V_meta(3,n,nf)
    nt = V_type(n)
    za = V_exp(eo + 1 + (dir*sp + nt))
    zb = V_exp(eo + 1 + (-dir*sp + nt))
  end subroutine wrap_factors

  Subroutine scale_M(nf, r, za, zb)
    Implicit none
    Integer, intent(in) :: nf, r
    Complex (Kind=Kind(0.d0)), intent(in) :: za, zb
    Integer :: m
    if (mode /= LR_SUBMATRIX) return
    if (jc(nf) == 0) return
    m = slot_h(r,nf)
    if (m == 0) return
    M_h(m,1:jc(nf),nf) = za * M_h(m,1:jc(nf),nf)
    M_h(1:jc(nf),m,nf) = zb * M_h(1:jc(nf),m,nf)
  end subroutine scale_M

  Subroutine clear_deferred()
    Implicit none
    dfr_wrap = .false.; dfr_n = 0; dfr_pp = 0; dfr_pk = 0; dfr_rw = 0
  end subroutine clear_deferred

!--------------------------------------------------------------------
!> WRAPGRUP step at vertex n, slice t: one launch for all flavours (wrap, gather / G(p,p)), one synchronisation
  Subroutine gpu_lr_step_up(n, t)
    Implicit none
    Integer, intent(in) :: n, t
    Integer :: fl, nf, nfm(2), r(2), eo(2), jcv(2)
    Complex (Kind=Kind(0.d0)) :: za, zb
    nfm = 0; r = 1; eo = 0; jcv = 0
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       nfm(fl) = nf; r(fl) = Op_V(n,nf)%P(1); eo(fl) = V_meta(3,n,nf); jcv(fl) = jc(nf)
    enddo
    if (N_FL_eff == 1) then
       nfm(2) = nfm(1); r(2) = r(1); eo(2) = eo(1); jcv(2) = jcv(1)
    endif
    if (mode == LR_SUBMATRIX) then
       call k_sub_up(N_FL_eff, nfm(1), nfm(2), Ndim, nd, gs, n, t, V_type(n), r(1), r(2), eo(1), eo(2), &
            & dfr_pp(1), dfr_pk(1), dfr_pp(2), dfr_pk(2), dfr_n, dfr_t, dfr_vr, dfr_vi, &
            & size(G_d,3), size(gbuf), size(nsigma_d,1), size(nsigma_d,2), size(V_exp), G_d, gbuf, slot_d, nsigma_d, V_exp)
    else
       call k_dly_up(N_FL_eff, nfm(1), nfm(2), Ndim, gs, n, t, V_type(n), r(1), r(2), eo(1), eo(2), jcv(1), jcv(2), &
            & dfr_n, dfr_t, dfr_vr, dfr_vi)
    endif
    upd_cnt(3) = upd_cnt(3) + 1
    n_fused = n_fused + 1
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       call wrap_factors(nf, n, t, 1, za, zb)
       call scale_M(nf, r(fl), za, zb)          ! submatrix: M <- L_S M R_S for the wrapped site
    enddo
    call clear_deferred()
    stepped = .true.
  end subroutine gpu_lr_step_up

!> WRAPGRDO step at vertex n, slice t: the deferred wrap of the previous vertex, then gather / G(p,p) at n, one launch
  Subroutine gpu_lr_step_down(n, t)
    Implicit none
    Integer, intent(in) :: n, t
    Integer :: fl, nf, nfm(2), pn(2), rw(2), jcv(2), wv
    Real (Kind=Kind(0.d0)) :: z(8)
    nfm = 0; pn = 1; rw = 1; jcv = 0; z = 0.d0
    wv = 0
    if (dfr_wrap) wv = 1
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       nfm(fl) = nf; pn(fl) = Op_V(n,nf)%P(1); jcv(fl) = jc(nf)
       rw(fl) = max(dfr_rw(fl), 1)
       if (wv /= 0) then
          z(4*fl-3) = real(dfr_za(fl)); z(4*fl-2) = aimag(dfr_za(fl))
          z(4*fl-1) = real(dfr_zb(fl)); z(4*fl)   = aimag(dfr_zb(fl))
       endif
    enddo
    if (N_FL_eff == 1) then
       nfm(2) = nfm(1); pn(2) = pn(1); rw(2) = rw(1); jcv(2) = jcv(1)
    endif
    if (mode == LR_SUBMATRIX) then
       call k_sub_down(N_FL_eff, nfm(1), nfm(2), Ndim, nd, gs, wv, pn(1), pn(2), rw(1), rw(2), &
            & z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), dfr_pp(1), dfr_pk(1), dfr_pp(2), dfr_pk(2), dfr_n, dfr_t, dfr_vr, dfr_vi, &
            & size(G_d,3), size(gbuf), size(nsigma_d,1), size(nsigma_d,2), G_d, gbuf, slot_d, nsigma_d)
    else
       call k_dly_down(N_FL_eff, nfm(1), nfm(2), Ndim, gs, wv, pn(1), pn(2), rw(1), rw(2), jcv(1), jcv(2), &
            & z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), dfr_n, dfr_t, dfr_vr, dfr_vi)
    endif
    upd_cnt(3) = upd_cnt(3) + 1
    n_fused = n_fused + 1
    call clear_deferred()
    stepped = .true.
  end subroutine gpu_lr_step_down

!> End of the site: record the deferred work.  accepted: the field of (n,t) changed (nsigma_d host copy is current);
!> down: WRAPGRDO, the wrap of vertex n (final field) follows and is deferred to the next launch.
  Subroutine gpu_lr_defer(n, t, accepted, down)
    Implicit none
    Integer, intent(in) :: n, t
    Logical, intent(in) :: accepted, down
    Integer :: fl, nf
    Complex (Kind=Kind(0.d0)) :: za, zb
    stepped = .false.
    if (accepted) then
       dfr_n = n; dfr_t = t
       dfr_vr = real(nsigma_d(n,t)); dfr_vi = aimag(nsigma_d(n,t))
    endif
    if (down) then
       dfr_wrap = .true.
       do fl = 1, N_FL_eff
          nf = Calc_Fl_map(fl)
          dfr_rw(fl) = Op_V(n,nf)%P(1)
          call wrap_factors(nf, n, t, -1, za, zb)
          dfr_za(fl) = za; dfr_zb(fl) = zb
          call scale_M(nf, dfr_rw(fl), za, zb)
       enddo
       if (N_FL_eff == 1) then
          dfr_rw(2) = dfr_rw(1); dfr_za(2) = dfr_za(1); dfr_zb(2) = dfr_zb(1)
       endif
    endif
  end subroutine gpu_lr_defer

!> delayed: launch the accepts collected by gpu_lr_accept (both flavours in one kernel), flush at capacity
  Subroutine gpu_lr_commit()
    Implicit none
    Integer :: fl, nf, nfm(2), pv(2), jv(2), fg(2)
    Complex (Kind=Kind(0.d0)) :: zz(2)
    if (.not. (fast .and. mode == LR_DELAYED)) return
    if (acc_flag(1) == 0 .and. acc_flag(2) == 0) return
    nfm = 1; pv = 1; jv = 0; fg = 0; zz = (0.d0,0.d0)
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       if (acc_flag(fl) /= 0 .and. jc(nf) + 1 > nd) call gpu_lr_flush(nf, 1)     ! capacity
    enddo
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       nfm(fl) = nf; pv(fl) = acc_p(fl); jv(fl) = jc(nf); fg(fl) = acc_flag(fl); zz(fl) = acc_z(fl)
    enddo
    if (N_FL_eff == 1) then
       nfm(2) = nfm(1); pv(2) = pv(1); jv(2) = jv(1); fg(2) = 0
    endif
    call k_dly_accept2(N_FL_eff, nfm(1), nfm(2), Ndim, pv(1), pv(2), jv(1), jv(2), fg(1), fg(2), &
         & real(zz(1)), aimag(zz(1)), real(zz(2)), aimag(zz(2)))
    upd_cnt(4) = upd_cnt(4) + 1
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       if (acc_flag(fl) /= 0) jc(nf) = jc(nf) + 1
    enddo
    acc_flag = 0
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       if (jc(nf) >= nd) call gpu_lr_flush(nf, 1)
    enddo
  end subroutine gpu_lr_commit

!> apply the deferred work now (before a flush, an end-of-slice hop, or any host reader of G_d)
  Subroutine gpu_lr_sync()
    Implicit none
    Integer :: fl, nf, nfm(2), rw(2), jcv(2), wv
    Real (Kind=Kind(0.d0)) :: z(8)
    if (.not. (dfr_wrap .or. dfr_n > 0 .or. any(dfr_pk > 0))) return
    nfm = 0; rw = 1; jcv = 0; z = 0.d0
    wv = 0
    if (dfr_wrap) wv = 1
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       nfm(fl) = nf; jcv(fl) = jc(nf); rw(fl) = max(dfr_rw(fl), 1)
       if (wv /= 0) then
          z(4*fl-3) = real(dfr_za(fl)); z(4*fl-2) = aimag(dfr_za(fl))
          z(4*fl-1) = real(dfr_zb(fl)); z(4*fl)   = aimag(dfr_zb(fl))
       endif
    enddo
    if (N_FL_eff == 1) then
       nfm(2) = nfm(1); rw(2) = rw(1); jcv(2) = jcv(1)
    endif
    if (mode == LR_SUBMATRIX) then
       call k_apply(N_FL_eff, nfm(1), nfm(2), Ndim, wv, rw(1), rw(2), z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), &
            & dfr_pp(1), dfr_pk(1), dfr_pp(2), dfr_pk(2), dfr_n, dfr_t, dfr_vr, dfr_vi, &
            & size(G_d,3), size(nsigma_d,1), size(nsigma_d,2), G_d, slot_d, nsigma_d)
    else
       call k_dly_apply(N_FL_eff, nfm(1), nfm(2), Ndim, wv, rw(1), rw(2), jcv(1), jcv(2), z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), &
            & dfr_n, dfr_t, dfr_vr, dfr_vi)
    endif
    upd_cnt(4) = upd_cnt(4) + 1
    call clear_deferred()
    stepped = .false.
  end subroutine gpu_lr_sync


!====================================================================
!  Block mode (submatrix, fast): a block of b <= N_delay consecutive vertices is handled entirely on the host.
!  Between flushes G0 changes only by the diagonal wrap scalings, and a proposal at p needs G0(p,p), G0(p,S), G0(S,p) with
!  S inside the block: the principal block Gb = G0(B,B) is copied ONCE (b x b per flavour), scaled on the host by every wrap,
!  and used for all b proposals.  At the end of the block the device applies the scalings to G0 (rows, then columns of the
!  block sites) and flushes M (S, M uploaded once).  Device round trips: ~12 per block instead of ~2 per site; k <= b <= N_delay,
!  so no capacity flush can happen inside a block.
!====================================================================
  Logical function gpu_lr_blockmode()
    Implicit none
    if (.not. inited) call gpu_lr_init()
    gpu_lr_blockmode = (mode == LR_SUBMATRIX) .and. fast
  end function gpu_lr_blockmode

#ifdef ALF_TEAM_TIMERS
  Real (Kind=Kind(0.d0)) function p3_now()
    Integer(int64) :: c, r
    call system_clock(c, r)
    p3_now = real(c,8)/real(r,8)
  end function p3_now
#endif

#if defined(ALF_LA_ONEMKL)
  Subroutine ab_check(rc, what)
    Integer (c_int), intent(in) :: rc
    Character(len=*), intent(in) :: what
    if (rc /= 0) then
       write(error_unit,*) 'gpu_lowrank: asynchronous block path (ALFSH_ASYNCBLK) failed: ', what
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    endif
  end subroutine ab_check

  Function devz(A) result(p)
    Complex (Kind=Kind(0.d0)), Target :: A(*)
    Type (c_ptr) :: p
    !$omp target data use_device_addr(A)
    p = c_loc(A(1))
    !$omp end target data
  end function devz

  Function devi(A) result(p)
    Integer, Target :: A(*)
    Type (c_ptr) :: p
    !$omp target data use_device_addr(A)
    p = c_loc(A(1))
    !$omp end target data
  end function devi

  Function offp(p, nbytes) result(q)
    Type (c_ptr), intent(in) :: p
    Integer (int64), intent(in) :: nbytes
    Type (c_ptr) :: q
    q = transfer(transfer(p, 0_c_intptr_t) + int(nbytes, c_intptr_t), q)
  end function offp

  Subroutine ab_up_z(dev, h, n)       ! host -> device copy, queued
    Type (c_ptr), intent(in) :: dev
    Complex (Kind=Kind(0.d0)), Target :: h(*)
    Integer, intent(in) :: n
    call ab_check(alfsb_copy(dev, c_loc(h(1)), 16_c_int64_t*n), 'copy z')
  end subroutine ab_up_z

  Subroutine ab_up_i(dev, h, n)
    Type (c_ptr), intent(in) :: dev
    Integer, Target :: h(*)
    Integer, intent(in) :: n
    call ab_check(alfsb_copy(dev, c_loc(h(1)), 4_c_int64_t*n), 'copy i')
  end subroutine ab_up_i

  Subroutine ab_drain()               ! every queued block operation done (before any other device work / host read of G_d)
    if (.not. ab_pending) return
    call ab_check(alfsb_wait(), 'wait')
    ab_pending = .false.
  end subroutine ab_drain

  Subroutine ab_setup()
    Character(len=8) :: v
    Integer :: st
    Type (c_ptr) :: qh
    ab_on = .false.
    v = ''
    call get_environment_variable('ALFSH_ASYNCBLK', v, status=st)
    if (st == 0 .and. trim(v) == '1') ab_on = .true.
    if (.not. ab_on) return
    qh = alfsh_get_queue()
    if (.not. c_associated(qh)) call ab_check(-1_c_int, 'no queue (alf_la not initialised)')
    st = alfsb_init(qh)
    if (st /= 0) call ab_check(-1_c_int, 'the interop queue is not in-order')
    ab_pG = devz(G_d); ab_pX = devz(X_d); ab_pS = devz(scr_c); ab_pSl = devi(slot_d); ab_pNs = devz(nsigma_d)
    ab_pBs = devi(bsd); ab_pBS2 = devi(bS_d); ab_pCb = devz(cbuf); ab_pGb = devz(gbd)
    ab_ph = alfsb_alloc_host(16_c_int64_t*nd*nd*2)
    if (.not. c_associated(ab_ph)) call ab_check(-1_c_int, 'host USM allocation')
    call c_f_pointer(ab_ph, gbh, [nd, nd, 2])
    write(output_unit,*) ' gpu_lowrank: ALFSH_ASYNCBLK=1, submatrix blocks run asynchronously on the SYCL queue'
  end subroutine ab_setup
#endif

  Subroutine blk_alloc()
    Implicit none
    if (allocated(Gb)) return
    allocate( Gb(nd,nd,2), bidx(Ndim,N_FL), bs_h(nd,2), zab_h(nd,2,2) )
    allocate( bsd(nd,2), bS_d(nd,2), gbd(nd,nd,2), cbuf(2*nd+nd*nd,2) )
    bidx = 0; bs_h = 0; Gb = czero_c; zab_h = czero_c
    bsd = 0; bS_d = 0; gbd = czero_c; cbuf = czero_c
    !$omp target enter data map(alloc: bsd, bS_d, gbd, cbuf)
#if defined(ALF_LA_ONEMKL)
    call ab_setup()
#endif
  end subroutine blk_alloc

!> k_blk_gather: gbd(a,c,fl) = G(s(a),s(c)) for the block sites.
!> All device arrays are explicit-shape dummy arguments (map(present,alloc:); nvfortran 26.9 rejects the OpenMP 5.1 present
!> modifier, so the CUDA path uses map(alloc:), which moves no data for arrays already mapped by enter data): referencing module allocatables inside the target
!> region makes ifx allocate and copy their descriptors as first-private data at every launch (measured 16.0 vs 10.7 us/launch).
  Subroutine k_blk_gather(nfl, nfa, nfb, ld, b, nd_, nfg, gbk, bsk, Gk)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, b, nd_, nfg
    Complex (Kind=Kind(0.d0)), intent(inout) :: gbk(nd_,nd_,2), Gk(ld,ld,nfg)
    Integer, intent(in) :: bsk(nd_,2)
    Integer :: fl, a, c, nf
#if defined(ALF_LA_HIP)
    call hk_blk_gather(nfl, nfa, nfb, b, nd_, ld, Gk(1,1,1), bsk(1,1), gbk(1,1,1))
#else
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do collapse(3) private(nf) map(present,alloc: gbk, bsk, Gk)
#else
    !$omp target teams distribute parallel do collapse(3) private(nf) map(alloc: gbk, bsk, Gk)
#endif
    do fl = 1, nfl
       do c = 1, b
          do a = 1, b
             nf = nfa
             if (fl == 2) nf = nfb
             gbk(a,c,fl) = Gk(bsk(a,fl),bsk(c,fl),nf)
          enddo
       enddo
    enddo
#endif
  end subroutine k_blk_gather

!> End of a block: rows, then columns of G(:,:,nf) scaled by the block's wrap factors (cbuf: za(1:b), zb(1:b) of each flavour),
!> both flavours in ONE launch.  Every element is written by exactly one iteration and gets exactly the factors, in exactly the
!> order, of the two-kernel version (rows kernel, then columns kernel): an element (s(a),j) gets cbuf(a)*G, and if j = s(c) is
!> itself a block site then cbuf(nd+c)*(cbuf(a)*G); an element (j,s(a)) with j not a block site gets cbuf(nd+a)*G.  The block
!> sites of one flavour are distinct (gpu_lr fast_check), so "j is a block site" is found by a search over bsk(1:b,fl).
  Subroutine k_blk_scale(nfl, nfa, nfb, ld, b, nd_, nfg, ldc, Gk, bsk, cbk)
    Implicit none
    Integer, intent(in) :: nfl, nfa, nfb, ld, b, nd_, nfg, ldc
    Complex (Kind=Kind(0.d0)), intent(inout) :: Gk(ld,ld,nfg)
    Integer, intent(in) :: bsk(nd_,2)
    Complex (Kind=Kind(0.d0)), intent(in) :: cbk(ldc,2)
    Integer :: fl, a, j, c, cj, nf, sa
    Complex (Kind=Kind(0.d0)) :: v
#if defined(ALF_LA_ONEMKL)
    !$omp target teams distribute parallel do collapse(3) private(nf, sa, c, cj, v) map(present,alloc: Gk, bsk, cbk)
#else
    !$omp target teams distribute parallel do collapse(3) private(nf, sa, c, cj, v) map(alloc: Gk, bsk, cbk)
#endif
    do fl = 1, nfl
       do a = 1, b
          do j = 1, ld
             nf = nfa
             if (fl == 2) nf = nfb
             sa = bsk(a,fl)
             cj = 0
             do c = 1, b
                if (bsk(c,fl) == j) cj = c
             enddo
             v = cbk(a,fl) * Gk(sa,j,nf)                  ! row factor
             if (cj > 0) v = cbk(nd_+cj,fl) * v           ! column factor of the same element (j is a block site)
             Gk(sa,j,nf) = v
             if (cj == 0) Gk(j,sa,nf) = cbk(nd_+a,fl) * Gk(j,sa,nf)
          enddo
       enddo
    enddo
  end subroutine k_blk_scale

!> open a block at vertex v1 (v1, v1 + dv, ...), b vertices
  Subroutine blk_begin(v1, dv, b)
    Implicit none
    Integer, intent(in) :: v1, dv, b
    Integer :: fl, nf, a, v, nfm(2)
#ifdef ALF_TEAM_TIMERS
    Real (Kind=Kind(0.d0)) :: ta, tb, tc
    ta = p3_now()
#endif
    call blk_alloc()
#if defined(ALF_LA_ONEMKL)
    if (.not. ab_on) &
#endif
    call la_sync()          ! previous flush (asynchronous library calls) before the new uploads
    blk_open = .true.; blk_v1 = v1; blk_dv = dv; blk_b = b; blk_cnt = 0
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl); nfm(fl) = nf
       bsd(1:nd,fl) = 0
       do a = 1, b
          v = v1 + dv*(a-1)
          bs_h(a,fl) = Op_V(v,nf)%P(1)
          bsd(a,fl) = bs_h(a,fl)
          bidx(bs_h(a,fl),nf) = a
       enddo
       zab_h(1:nd,1:2,fl) = (1.d0,0.d0)
    enddo
    if (N_FL_eff == 1) nfm(2) = nfm(1)
#ifdef ALF_TEAM_TIMERS
    tb = p3_now()
#endif
#if defined(ALF_LA_ONEMKL)
    if (ab_on) then
       call ab_up_i(ab_pBs, bsd, size(bsd))
       call ab_check(alfsb_gather(N_FL_eff, nfm(1), nfm(2), Ndim, b, nd, ab_pGb, ab_pBs, ab_pG), 'gather')
       call ab_check(alfsb_copy(ab_ph, ab_pGb, 16_c_int64_t*nd*nd*2), 'copy gbd')
       ab_pending = .true.
       call ab_drain()      ! the one unavoidable wait of a block: the host needs G0(B,B)
    else
#endif
    !$omp target update to(bsd)
    call k_blk_gather(N_FL_eff, nfm(1), nfm(2), Ndim, b, nd, size(G_d,3), gbd, bsd, G_d)
    !$omp target update from(gbd)
#if defined(ALF_LA_ONEMKL)
    endif
#endif
#ifdef ALF_TEAM_TIMERS
    tc = p3_now(); p3t(2) = p3t(2) + (tc - tb)
#endif
    do fl = 1, N_FL_eff
#if defined(ALF_LA_ONEMKL)
       if (ab_on) then
          Gb(1:b,1:b,fl) = gbh(1:b,1:b,fl)
       else
#endif
       Gb(1:b,1:b,fl) = gbd(1:b,1:b,fl)
#if defined(ALF_LA_ONEMKL)
       endif
#endif
    enddo
    n_blocks = n_blocks + 1
    upd_cnt(3) = upd_cnt(3) + 1
#ifdef ALF_TEAM_TIMERS
    p3t(1) = p3t(1) + (p3_now() - ta); p3n(1) = p3n(1) + 1
#endif
  end subroutine blk_begin

!> block end: scalings of G0, flush of M, field values to the device
!> NOTE: nvfortran 26.9 at -O2/-O3 mis-compiles this routine (wrong G after the first stabilisation; -O1 and -O0 are correct,
!> found by bisecting the optimisation level per routine).  Kept at -O1 until the cause is understood.
!pgi$r opt=1
  Subroutine blk_end(t)
    Implicit none
    Integer, intent(in) :: t
    Integer :: fl, nf, a, k, b, m, v, vmin, vmax, nfm(2), kk
#ifdef ALF_TEAM_TIMERS
    Real (Kind=Kind(0.d0)) :: ta, tb, tc
    ta = p3_now()
#endif
    b = blk_b
#if !defined(ALF_LA_CUDA)
#if defined(ALF_LA_ONEMKL)
    if (.not. ab_on) &
#endif
    call la_sync()          ! the previous block's flush (cuBLAS) is asynchronous: it must be done before cbuf / bS_d are overwritten
#else
    ! CUDA (speedup S2): no la_sync here.  (1) The previous block's flush (the only library work reading the device cbuf / bS_d)
    ! was waited for by the la_sync at the start of this block's blk_begin, and blk_begin ends with a blocking D->H of gbd.
    ! (2) Between blk_begin and here only host work, synchronous target regions/updates and (non-block sites only)
    ! gpu_lr_flush can run; gpu_lr_flush does not touch cbuf / bS_d and syncs itself before reusing S_d / M_pk.  (3) cbuf / bS_d
    ! below are written on the HOST copy (no unified memory in the CUDA build), the device copy changes only by the target update,
    ! and cuBLAS/cuSOLVER run on the stream of the synchronous target regions (alf_la_mod, la_init), so that update and
    ! k_blk_scale are ordered after any library call still in flight (G_d writes included).  The removed sync waited on an idle
    ! stream (its check_info_slots found islot = 0, reset by blk_begin's la_sync).
#endif
#ifdef ALF_TEAM_TIMERS
    p3t(4) = p3t(4) + (p3_now() - ta)
#endif
    ! cbuf: za, zb of every block site (1 if not wrapped), M packed after them; bS_d: pending sites
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       k = jc(nf)
       cbuf(1:b,fl) = zab_h(1:b,1,fl)
       cbuf(nd+1:nd+b,fl) = zab_h(1:b,2,fl)
       bS_d(1:nd,fl) = 0
       if (k > 0) then
          bS_d(1:k,fl) = S_h(1:k,nf)
          do m = 1, k
             cbuf(2*nd+(m-1)*k+1:2*nd+m*k,fl) = M_h(1:k,m,nf)
          enddo
          if (k == 1 .and. nd >= 2) then          ! k = 1: flush as k = 2 with a zero row / column (see blk_flush)
             bS_d(2,fl) = S_h(1,nf)
             cbuf(2*nd+1:2*nd+4,fl) = czero_c
             cbuf(2*nd+1,fl) = M_h(1,1,nf)
          endif
       endif
    enddo
#ifdef ALF_TEAM_TIMERS
    tb = p3_now()
#endif
#if defined(ALF_LA_ONEMKL)
    if (ab_on) then
       call ab_up_z(ab_pCb, cbuf, size(cbuf)); call ab_up_i(ab_pBS2, bS_d, size(bS_d))
    else
#endif
    !$omp target update to(cbuf, bS_d)
#if defined(ALF_LA_ONEMKL)
    endif
#endif
    do fl = 1, N_FL_eff
       nfm(fl) = Calc_Fl_map(fl)
    enddo
    if (N_FL_eff == 1) nfm(2) = nfm(1)
#if defined(ALF_LA_ONEMKL)
    if (ab_on) then
       call ab_check(alfsb_scale(N_FL_eff, nfm(1), nfm(2), Ndim, b, nd, size(cbuf,1), ab_pG, ab_pBs, ab_pCb), 'scale')
    else
#endif
#if !defined(ALF_LA_HIP)
    ! HIP scales inside the fused hk_blk_end kernel below
    call k_blk_scale(N_FL_eff, nfm(1), nfm(2), Ndim, b, nd, size(G_d,3), size(cbuf,1), G_d, bsd, cbuf)
#endif
#if defined(ALF_LA_ONEMKL)
    endif
#endif
#if !defined(ALF_LA_HIP)
    upd_cnt(4) = upd_cnt(4) + 1
#endif
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
#if defined(ALF_LA_HIP)
       ! one fused kernel: row and column scalings plus the flush gathers (AMD only: launch count)
       kk = jc(nf)
       if (kk == 1 .and. nd >= 2) kk = 2
       call hk_blk_end(nf, fl, Ndim, b, nd, size(cbuf,1), kk, size(slot_d,1), G_d(1,1,1), bsd(1,1), cbuf(1,1), &
            &          bS_d(1,fl), X_d(1,1,nf), scr_c(1,1,1), slot_d(1,1))
       upd_cnt(4) = upd_cnt(4) + 1
       if (jc(nf) > 0) call blk_flush(nf, fl, .true.)
#else
       if (jc(nf) > 0) then
          call blk_flush(nf, fl, .false.)
       endif
#endif
       do a = 1, b
          bidx(bs_h(a,fl),nf) = 0
       enddo
    enddo
    vmin = min(blk_v1, blk_v1 + blk_dv*(b-1)); vmax = max(blk_v1, blk_v1 + blk_dv*(b-1))
#if defined(ALF_LA_ONEMKL)
    if (ab_on) then
       call ab_up_z(offp(ab_pNs, 16_int64*((vmin-1) + int(size(nsigma_d,1),int64)*(t-1))), nsigma_d(vmin,t), vmax-vmin+1)
       ab_pending = .true.          ! drained by the next blk_begin, or by the flush / sync entry points at the end of the slice
    else
#endif
    !$omp target update to(nsigma_d(vmin:vmax,t))
#if defined(ALF_LA_ONEMKL)
    endif
#endif
    upd_cnt(6) = upd_cnt(6) + 1
    blk_open = .false.; blk_cnt = 0
#ifdef ALF_TEAM_TIMERS
    tc = p3_now(); p3t(5) = p3t(5) + (tc - tb); p3t(3) = p3t(3) + (tc - ta)
#endif
  end subroutine blk_end

#if defined(ALF_LA_ONEMKL)
!> the device sequence of blk_flush on the SYCL queue: same kernel, same two GEMMs (shapes, ops, scalars, leading dimensions)
  Subroutine ab_flush_dev(nf, fl, kk)
    Integer, intent(in) :: nf, fl, kk
    Integer (int64) :: oG, oX, oS2, oC, oB
    oG  = 16_int64*size(G_d,1)*size(G_d,2)*(nf-1)
    oX  = 16_int64*size(X_d,1)*size(X_d,2)*(nf-1)
    oS2 = 16_int64*size(scr_c,1)*size(scr_c,2)
    oC  = 16_int64*(2*nd + int(size(cbuf,1),int64)*(fl-1))
    oB  = 4_int64*size(bS_d,1)*(fl-1)
    call ab_check(alfsb_fgather(offp(ab_pG,oG), Ndim, Ndim, offp(ab_pBS2,oB), kk, nf, offp(ab_pX,oX), ab_pS, ab_pSl), 'fgather')
    call ab_check(alfsb_gemm(0, 1, Ndim, kk, kk, 1.d0, 0.d0, ab_pS, Ndim, offp(ab_pCb,oC), kk, 0.d0, 0.d0, offp(ab_pS,oS2), Ndim), 'gemm1')
    call ab_check(alfsb_gemm(0, 1, Ndim, Ndim, kk, -1.d0, 0.d0, offp(ab_pX,oX), Ndim, offp(ab_pS,oS2), Ndim, 1.d0, 0.d0, &
         offp(ab_pG,oG), Ndim), 'gemm2')
  end subroutine ab_flush_dev
#endif

!> flush from the already uploaded S and packed M (bS_d, cbuf): same device sequence as gpu_lr_flush, no further uploads
  Subroutine blk_flush(nf, fl, gathered)
    Implicit none
    Integer, intent(in) :: nf, fl
    Logical, intent(in) :: gathered       ! gathers (X, scr 1) already done by hk_blk_end
    Integer :: k, m, kk
    Complex (Kind=Kind(0.d0)), parameter :: cone = (1.d0,0.d0), cmone = (-1.d0,0.d0), czero = (0.d0,0.d0)
    k = jc(nf)
    kk = k
    if (k == 1 .and. nd >= 2) kk = 2
    n_flush(3) = n_flush(3) + 1; n_flush(2) = n_flush(2) + 1
#if defined(ALF_LA_ONEMKL)
    if (ab_on) then
       call ab_flush_dev(nf, fl, kk)
    else
#endif
    if (.not. gathered) call k_sub_flush_gather(G_d(1,1,nf), Ndim, Ndim, bS_d(1,fl), kk, nf, X_d(1,1,nf), scr_c(1,1,1), slot_d, size(slot_d,2))
    call la_gemm(LA_N, LA_T, Ndim, kk, kk, cone, scr_c(1,1,1), Ndim, cbuf(2*nd+1,fl), kk, czero, scr_c(1,1,2), Ndim)
    call la_gemm(LA_N, LA_T, Ndim, Ndim, kk, cmone, X_d(1,1,nf), Ndim, scr_c(1,1,2), Ndim, cone, G_d(1,1,nf), Ndim)
#if defined(ALF_LA_ONEMKL)
    endif
#endif
    upd_cnt(4) = upd_cnt(4) + 3
    do m = 1, k
       slot_h(S_h(m,nf),nf) = 0
       S_h(m,nf) = 0
    enddo
    M_h(1:k,1:k,nf) = czero
    jc(nf) = 0
  end subroutine blk_flush

!> vertex n of the sequence n, n + dv, ... up to nlast, slice t: opens a block if none is open; dv = +1 (WRAPGRUP):
!> the wrap of vertex n (old field) is applied to the host copy before the proposal
  Subroutine gpu_lr_blk_enter(n, nlast, dv, t, wrap_first)
    Implicit none
    Integer, intent(in) :: n, nlast, dv, t
    Logical, intent(in) :: wrap_first
    Integer :: b, fl, nf
    if (.not. blk_open) then
       b = min(nd, abs(nlast - n) + 1)
       call blk_begin(n, dv, b)
    endif
    blk_a = abs(n - blk_v1) + 1
    if (wrap_first) call blk_wrap_vertex(n, t, 1)
  end subroutine gpu_lr_blk_enter

!> end of the vertex: WRAPGRDO applies the wrap after the update (final field); closes the block when it is complete
  Subroutine gpu_lr_blk_leave(n, t, wrap_after)
    Implicit none
    Integer, intent(in) :: n, t
    Logical, intent(in) :: wrap_after
    if (wrap_after) call blk_wrap_vertex(n, t, -1)
    blk_cnt = blk_cnt + 1
    if (blk_cnt >= blk_b) call blk_end(t)
  end subroutine gpu_lr_blk_leave

  Subroutine blk_wrap_vertex(n, t, dir)
    Implicit none
    Integer, intent(in) :: n, t, dir
    Integer :: fl, nf, a
    Complex (Kind=Kind(0.d0)) :: za, zb
    a = blk_a
    do fl = 1, N_FL_eff
       nf = Calc_Fl_map(fl)
       call wrap_factors(nf, n, t, dir, za, zb)
       Gb(a,1:blk_b,fl) = za * Gb(a,1:blk_b,fl)
       Gb(1:blk_b,a,fl) = zb * Gb(1:blk_b,a,fl)
       zab_h(a,1,fl) = za; zab_h(a,2,fl) = zb
       call scale_M(nf, bs_h(a,fl), za, zb)
    enddo
  end subroutine blk_wrap_vertex

!> G(p,p) of the current G from the host block (submatrix, block mode): p = bs_h(blk_a, fl)
  Subroutine blk_gpp(nf, nf_eff, p, gpp)
    Implicit none
    Integer, intent(in) :: nf, nf_eff, p
    Complex (Kind=Kind(0.d0)), intent(out) :: gpp
    Integer :: k, a, m, i, ap, bm
    Complex (Kind=Kind(0.d0)) :: s
    k = jc(nf)
    ap = blk_a
    gpp = Gb(ap,ap,nf_eff)
    prop_p(nf) = p; prop_k(nf) = k
    if (k == 0) return
    i = slot_h(p,nf)
    do a = 1, k
       s = czero_c
       do m = 1, k
          bm = bidx(S_h(m,nf),nf)
          s = s - M_h(a,m,nf) * Gb(bm,ap,nf_eff)
       enddo
       if (i /= 0) s = s + M_h(a,i,nf)
       W_h(a,nf) = s
       grow_h(a,nf) = Gb(ap,bidx(S_h(a,nf),nf),nf_eff)
    enddo
    s = czero_c
    do a = 1, k
       s = s + grow_h(a,nf) * W_h(a,nf)
    enddo
    gpp = gpp - s
  end subroutine blk_gpp

end Module gpu_lowrank_mod
#endif

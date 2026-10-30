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
!> GPU port, ticket 206 (ADR-001 section 4, ADR-002): stabilization on the device.  Compiled only with -DALF_GPU.
!>
!> What runs on the device (finite temperature, default and STABLOG stabilization; STAB1/2/3 and the projector code are
!> refused by gpu_supported):
!>   * gpu_stab_wrapur / gpu_stab_wrapul: the slice products of WRAPUR / WRAPUL (ticket 204 composites) on the U matrices of
!>     the device states, then the UDV decomposition (decompose_UDV_state) and the scale test (testscale);
!>   * the decomposition: default = column-norm pre-pivot + unpivoted QR (alf_la la_qrpp, owner-approved, ticket 104), D
!>     separated by a kernel, phase and permutation sign, V = R * P^-1 V by one GEMM on the triangular factor, Q = reflectors applied
!>     to the unit matrix by la_ormqr (NOT la_ungqr: see dev_decompose, accuracy of the sweep's Precision Green).
!>     STABLOG = the CPU algorithm itself: log-norm selection sort (the norms D->H, the swap permutation H->D, columns
!>     gathered on the device), unpivoted QR (la_geqrf), the exp(L_j - L_i) row scaling and the log update;
!>   * gpu_stab_cgr: CGR for NVAR = 1 and 2 (cgr1_mod CGR, branches default and STABLOG).  Default: la_qrpp.  STABLOG: host
!>     ZGEQP3 through la_geqp3_host (ADR-001 section 4: upstream numerics; the N x N matrix goes D->H and the factors H->D).
!>     D separation, phase (permutation sign, R diagonal, reflector determinants, det of U_R^dagger U_L through la_getrf),
!>     ZUNMQR, ZTRSM, ZLAPMR/T and the final GEMM on the device;
!>   * the precision check of Control_PrecisionG as a device reduction (max |A-B|, sum |A-B|, NaN counts): one D->H of four
!>     numbers, handed to Control_PrecisionG_vals (control_mod).
!>
!> Device states (gpu_state_mod): udw_*(.,.,1,nfe) = udvl, udw_*(.,.,2,nfe) = udvr, udv_*(.,.,NST,nfe) = udvst.  Host UDV_State
!> objects of main are NOT updated by the device path, except (interim) udvst(NST) in the backward sweep when Ltau = 1,
!> because TAU_M (host, ticket 207) reads it.  TAU_M calls the host WRAPUR/CGR2_2 on its own objects, so WRAPUR/WRAPUL
!> themselves stay CPU code; the device branches are in main.
!>
!> Transfers (every one is marked in the code; per stabilization of one flavour, N = Ndim):
!>   decompose, default : 8 N B D->H (column norms), 4 N B H->D (permutation), 24 B D->H (scale test) [inside la_qrpp: the
!>                        norm/permutation transfers]
!>   decompose, STABLOG : 8 N B D->H (log norms), 4 N B H->D, 0 B for the test
!>   CGR default        : 8 N + 4 N B (la_qrpp), 16 B D->H (phase)
!>   CGR STABLOG        : 16 N^2 B D->H + 16 N^2 + 16 N + 4 N B H->D (la_geqp3_host), 16 B D->H (phase)
!>   INTERIM 206 (removed by 207/208): G D->H 16 N^2 B (GR is host-authoritative), wrapped G H->D 16 N^2 B (for the precision
!>   reduction), 32 B D->H (precision numbers); with Ltau = 1 in the backward sweep: udvst(NST) D->H, 32 N^2 + 16 N B.
!>
!> Self-test ALF_GPU_SELFTEST=3 (oracle = the upstream host routines WRAPUR/WRAPUL/CGR/decompose on host UDV_State objects;
!> the default path needs a build with -DQR_PREPIVOT, the CPU variant of the same algorithm), timing ALF_GPU_STABTIME=1,
!> ALF_GPU_EXIT=1 ends the run after the hook.  ALF_GPU_FAULT_STAB=perm|phase|nan injects faults (negative controls).
!> Debug: ALF_GPU_STABCHECK=1|2 repeats every CGR (1) and every WRAPUR/WRAPUL block (2) with the host routines on the same device
!> states and reports the deviations in info; ALF_GPU_STABMIX=1|2 uses the host CGR / host block result instead (isolation).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_stab_mod

  Use omp_lib
  Use iso_c_binding, only: c_int
  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Hamiltonian_main
  Use QMC_runtime_var
  Use MyMats, only: COMPARE
  Use UDV_State_mod
  Use QDRP_mod, only: Pivot_phase
  Use wrapur_mod, only: WRAPUR
  Use wrapul_mod, only: WRAPUL
  Use cgr1_mod, only: CGR
  Use Control, only: Control_PrecisionG_vals
  Use alf_la_mod
  Use gpu_state_mod
  Use gpu_update_mod, only: gpu_upd_enabled
  Use gpu_wrap_mod, only: gpu_wrapur_slices, gpu_wrapul_slices

  Implicit none

  Private
  Public :: gpu_stab_on, gpu_stab_setup, gpu_stab_finalize, gpu_stab_hook
  Public :: gpu_stab_wrapur, gpu_stab_wrapul, gpu_stab_exchange, gpu_stab_store_l, gpu_stab_cgr
  Public :: gpu_stab_reset_l, gpu_stab_reset_r, gpu_stab_reset_st

  Logical, Save :: gpu_stab_on = .false.

  Integer, Save :: s_N = 0, s_nfe = 0, s_nstm = 0
  Integer, Allocatable, Save :: s_nt(:)                     ! Stab_nt(0:NSTM)
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

  ! ---- own device buffers (allocated once in gpu_stab_setup) ----
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: sA(:,:), sB(:,:), sC(:,:)   ! N x N work matrices
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_Gn(:,:), st_Gt(:,:)       ! new G, wrapped G (device + host shadow)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_tau(:), st_Dl(:)          ! QR reflector scalars, local D of the QR
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_c(:)                      ! device scalars: 1 phase, 2 1/phase, 3 CGR phase
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_r(:)                      ! device scalars: scale test, precision numbers
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_nrm(:), st_Lt(:)          ! column norms / log norms; L work vector
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_pw(:)                     ! per-column partial sums of k_prec (ticket 212)
  Integer (c_int), Allocatable, Save :: st_ipvt(:), st_perm(:), st_ipiv(:)      ! pivot (ZGEQP3 semantics), sort permutation, LU pivots
  Logical, Save :: warn_once = .true.
  Logical, Save :: ns_on = .true.    ! Newton-Schulz step on Q in decompose (ALF_GPU_NS=0 switches it off; see DISCREPANCIES.md, cuSOLVER Q orthogonality)
  Logical, Save :: in_selftest = .false.                                        ! inside gpu_stab_selftest (ticket 212)
  Integer, Save :: fault = 0                                                    ! 1 perm, 2 phase, 3 nan (negative controls)
  ! ---- DEBUG AIDS, not used in production: ALF_GPU_STABCHECK=1|2: every CGR (1) and every WRAPUR/WRAPUL block (2) of the run is repeated on the host
  ! with the upstream routines from the same device state, and the deviations are summarised in info (host work, slow)
  Integer, Save :: chk_level = 0, chk_ncgr = 0, chk_nblk = 0, chk_nbad = 0, chk_mix = 0   ! chk_mix 1: use the host CGR result, 2: the host WRAPUR/UL result (isolation experiments)
  Real (Kind=Kind(0.d0)), Save :: chk_g(2) = 0.d0, chk_ph = 0.d0, chk_u = 0.d0, chk_v = 0.d0, chk_s = 0.d0
  Class (UDV_State), Allocatable, Save :: ck_l(:), ck_r(:), ck_d(:)

Contains

!====================================================================
!  Setup
!====================================================================
  Subroutine gpu_stab_setup(Stab_nt, NSTM)
    Implicit none
    Integer, intent(in) :: NSTM
    Integer, intent(in) :: Stab_nt(0:NSTM)
    Integer :: n, nfe, lw
    Character (len=256) :: line

    gpu_stab_on = .false.
    if (.not. gpu_active) return
#if defined(ALF_LA_ONEMKL)
    call info_line('GPU stabilization (206): INACTIVE, the alf_la backend of GPU=intel is not implemented yet (ticket 302)')
    return
#else
    s_N = Ndim; s_nfe = N_FL_eff; s_nstm = NSTM
    allocate(s_nt(0:NSTM)); s_nt = Stab_nt
    n = s_N
    allocate(sA(n,n), sB(n,n), sC(n,n), st_Gn(n,n), st_Gt(n,n), st_tau(n), st_Dl(n), st_c(8), st_r(16), st_nrm(n), st_Lt(n), st_pw(n))
    allocate(st_ipvt(n), st_perm(n), st_ipiv(n))
    st_ipvt = 1; st_perm = 1; st_ipiv = 1
    ! device allocation once, no transfer; the first use of each buffer is a device write
    !$omp target enter data map(alloc: sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_pw, st_ipvt, st_perm, st_ipiv)
    call la_init()
    if (gpu_env_is('ALF_GPU_NS', '0')) ns_on = .false.      ! (ticket 212: the value 0 is listed in the debug-env line)
    ! pre-grow the cuSOLVER work space (no device allocation inside the sweep)
    lw = la_geqrf_lwork(n, n, sA, n)
    lw = max(lw, la_ormqr_lwork('L', LA_C, n, n, n, sA, n, st_tau, sB, n))
    lw = max(lw, la_ormqr_lwork('R', LA_N, n, n, n, sA, n, st_tau, sB, n))
    lw = max(lw, la_getrf_lwork(n, n, sA, n))
    call la_reserve_work(lw)
    do nfe = 1, s_nfe
       call gpu_stab_reset_l(nfe)
       call gpu_stab_reset_r(nfe)
       call gpu_stab_reset_st(s_nstm, nfe)
    enddo
    gpu_stab_on = .true.
    write(line,'(A,I0,A,I0,A)') 'GPU stabilization (206): ACTIVE (decompose, CGR and precision check on the device), N = ', n, &
         & ', Nstm = ', s_nstm, ', pre-pivoted QR'
#if defined(STABLOG)
    line = trim(line)//' in decompose, host ZGEQP3 in CGR (STABLOG)'
#else
    line = trim(line)//' in decompose and CGR (default)'
#endif
    call info_line(line)
#endif
  end subroutine gpu_stab_setup

  Subroutine gpu_stab_finalize()
    Implicit none
    if (.not. gpu_stab_on) return
    call check_report()
    !$omp target exit data map(delete: sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_pw, st_ipvt, st_perm, st_ipiv)
    deallocate(sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_pw, st_ipvt, st_perm, st_ipiv, s_nt)
    call la_finalize()
    gpu_stab_on = .false.
  end subroutine gpu_stab_finalize

!====================================================================
!  Slot management (device to device)
!====================================================================
  Subroutine k_ident(n, U, ld)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (i == j) then
             U(i,j) = z1
          else
             U(i,j) = z0
          endif
       enddo
    enddo
  end subroutine k_ident

  Subroutine k_setD(n, D)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       D(i) = z1
    enddo
  end subroutine k_setD

  Subroutine k_setL(n, L)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(inout) :: L(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       L(i) = 0.d0
    enddo
  end subroutine k_setL

  Subroutine k_copy_c(nel, a, b)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)), intent(in) :: a(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: b(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nel
       b(i) = a(i)
    enddo
  end subroutine k_copy_c

  Subroutine k_copy_r(nel, a, b)
    Implicit none
    Integer, intent(in) :: nel
    Real (Kind=Kind(0.d0)), intent(in) :: a(*)
    Real (Kind=Kind(0.d0)), intent(inout) :: b(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nel
       b(i) = a(i)
    enddo
  end subroutine k_copy_r

  !> U = V = 1, D = 1 / L = 0 for slot of udw (1 = udvl, 2 = udvr)
  Subroutine reset_w(slot, nfe)
    Implicit none
    Integer, intent(in) :: slot, nfe
    call k_ident(s_N, udw_U(1,1,slot,nfe), s_N)
    call k_ident(s_N, udw_V(1,1,slot,nfe), s_N)
#if defined(STABLOG)
    call k_setL(s_N, udw_L(1,slot,nfe))
#else
    call k_setD(s_N, udw_D(1,slot,nfe))
#endif
  end subroutine reset_w

  Subroutine gpu_stab_reset_l(nfe)
    Implicit none
    Integer, intent(in) :: nfe
    call reset_w(1, nfe)
  end subroutine gpu_stab_reset_l

  Subroutine gpu_stab_reset_r(nfe)
    Implicit none
    Integer, intent(in) :: nfe
    call reset_w(2, nfe)
  end subroutine gpu_stab_reset_r

  !> udvst(NST, nfe) = identity state
  Subroutine gpu_stab_reset_st(nst, nfe)
    Implicit none
    Integer, intent(in) :: nst, nfe
    call k_ident(s_N, udv_U(1,1,nst,nfe), s_N)
    call k_ident(s_N, udv_V(1,1,nst,nfe), s_N)
#if defined(STABLOG)
    call k_setL(s_N, udv_L(1,nst,nfe))
#else
    call k_setD(s_N, udv_D(1,nst,nfe))
#endif
  end subroutine gpu_stab_reset_st

  !> dst = src slot copies: wslot (1/2) <- udvst(nst), or udvst(nst) <- wslot
  Subroutine copy_w_to_st(wslot, nfe, nst)
    Implicit none
    Integer, intent(in) :: wslot, nfe, nst
    call k_copy_c(s_N*s_N, udw_U(1,1,wslot,nfe), udv_U(1,1,nst,nfe))
    call k_copy_c(s_N*s_N, udw_V(1,1,wslot,nfe), udv_V(1,1,nst,nfe))
#if defined(STABLOG)
    call k_copy_r(s_N, udw_L(1,wslot,nfe), udv_L(1,nst,nfe))
#else
    call k_copy_c(s_N, udw_D(1,wslot,nfe), udv_D(1,nst,nfe))
#endif
  end subroutine copy_w_to_st

  Subroutine copy_st_to_w(nst, nfe, wslot)
    Implicit none
    Integer, intent(in) :: wslot, nfe, nst
    call k_copy_c(s_N*s_N, udv_U(1,1,nst,nfe), udw_U(1,1,wslot,nfe))
    call k_copy_c(s_N*s_N, udv_V(1,1,nst,nfe), udw_V(1,1,wslot,nfe))
#if defined(STABLOG)
    call k_copy_r(s_N, udv_L(1,nst,nfe), udw_L(1,wslot,nfe))
#else
    call k_copy_c(s_N, udv_D(1,nst,nfe), udw_D(1,wslot,nfe))
#endif
  end subroutine copy_st_to_w

  !> main: udvst(NST,nfe) = udvl(nfe)  (initial pass; forward = .false. semantics of the store)
  Subroutine gpu_stab_store_l(nst, nfe)
    Implicit none
    Integer, intent(in) :: nst, nfe
    call copy_w_to_st(1, nfe, nst)
  end subroutine gpu_stab_store_l

  !> fwd = .true.  (forward sweep):  udvl = udvst(NST); udvst(NST) = udvr
  !> fwd = .false. (backward sweep): udvr = udvst(NST); udvst(NST) = udvl; with Ltau = 1 the new udvst(NST) is copied D->H
  !>                                 into the host object ust, because TAU_M (host) reads it  (INTERIM 206: removed by 207)
  Subroutine gpu_stab_exchange(nst, nfe, fwd, ust)
    Implicit none
    Integer, intent(in) :: nst, nfe
    Logical, intent(in) :: fwd
    Class (UDV_State), intent(inout) :: ust
    if (fwd) then
       call copy_st_to_w(nst, nfe, 1)
       call copy_w_to_st(2, nfe, nst)
    else
       call copy_st_to_w(nst, nfe, 2)
       call copy_w_to_st(1, nfe, nst)
       if (get_Ltau() == 1) then
          ! INTERIM 206: removed by 207/208: D->H of the stored left propagation, TAU_M (host) reads udvst
          !$omp target update from(udv_U(:,:,nst,nfe), udv_V(:,:,nst,nfe))
          ust%U = udv_U(:,:,nst,nfe)
          ust%V = udv_V(:,:,nst,nfe)
          ust%side = 'l'
#if defined(STABLOG)
          !$omp target update from(udv_L(:,nst,nfe))
          ust%L = udv_L(:,nst,nfe)
#else
          !$omp target update from(udv_D(:,nst,nfe))
          ust%D = udv_D(:,nst,nfe)
#endif
       endif
    endif
  end subroutine gpu_stab_exchange

!====================================================================
!  WRAPUR / WRAPUL
!====================================================================
  !> WRAPUR(NTAU, NTAU1, udvr) for all flavours: slices NTAU+1..NTAU1 on U of udvr, decompose, scale test.
  Subroutine gpu_stab_wrapur(ntau, ntau1)
    Implicit none
    Integer, intent(in) :: ntau, ntau1
    Integer :: nfe, nf
    call sync_fields(ntau+1, ntau1)
    if (chk_level >= 2) call check_block_in(2)
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe)
       call gpu_wrapur_slices(udw_U(1,1,2,nfe), s_N, s_N, nf, ntau, ntau1)
       call decompose_slot('R', 2, nfe)
    enddo
    if (chk_level >= 2) call check_block_out(2, ntau, ntau1)
  end subroutine gpu_stab_wrapur

  !> WRAPUL(NTAU1, NTAU, udvl) for all flavours: slices NTAU1..NTAU+1 (descending) on U of udvl, decompose, scale test.
  Subroutine gpu_stab_wrapul(ntau1, ntau)
    Implicit none
    Integer, intent(in) :: ntau1, ntau
    Integer :: nfe, nf
    call sync_fields(ntau+1, ntau1)
    if (chk_level >= 2) call check_block_in(1)
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe)
       call gpu_wrapul_slices(udw_U(1,1,1,nfe), s_N, s_N, nf, ntau, ntau1)
       call decompose_slot('L', 1, nfe)
    enddo
    if (chk_level >= 2) call check_block_out(1, ntau, ntau1)
  end subroutine gpu_stab_wrapul

  !> INTERIM 206: only for the host-update path (ALF_GPU_NOUPDATE=1); with device updates (205) nsigma_d is current and this is skipped.
  !> The HS fields of the slices t1..t2 are changed by the host updates, so they go H->D before the slice products (N_op * (t2-t1+1) * 16 B).
  Subroutine sync_fields(t1, t2)
    Implicit none
    Integer, intent(in) :: t1, t2
    if (t2 < t1) return
    if (gpu_upd_enabled()) return      ! device update path: nsigma_d is kept current by the flips (ticket 205)
    nsigma_d(:,t1:t2) = nsigma%f(:,t1:t2)
    ! INTERIM 206: removed by 207/208: H->D of the block's HS fields
    !$omp target update to(nsigma_d(:,t1:t2))
  end subroutine sync_fields

  Subroutine decompose_slot(side, slot, nfe)
    Implicit none
    Character (len=1), intent(in) :: side
    Integer, intent(in) :: slot, nfe
#if defined(STABLOG)
    call dev_decompose(side, udw_U(1,1,slot,nfe), udw_V(1,1,slot,nfe), udw_L(1,slot,nfe))
#else
    call dev_decompose(side, udw_U(1,1,slot,nfe), udw_V(1,1,slot,nfe), udw_D(1,slot,nfe))
    call dev_testscale(udw_D(1,slot,nfe))
#endif
  end subroutine decompose_slot

!--------------------------------------------------------------------
!> decompose_UDV_state on the device.  U, V: N x N device arrays; S: D (default, complex) or L (STABLOG, real).
!> On input U holds B * U_old (propagated), D/L and V the previous factors; side R: A = U D V, side L: A = U D V^dagger.
!--------------------------------------------------------------------
#if defined(STABLOG)
  Subroutine dev_decompose(side, U, V, S)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: S(*)
#else
  Subroutine dev_decompose(side, U, V, S)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: S(*)
#endif
    Character (len=1), intent(in) :: side
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(*), V(*)
    Integer :: n, i, nswap, sgn
    Integer(c_int) :: kk
    Logical :: rside
    Complex (Kind=Kind(0.d0)) :: ph1

    n = s_N
    rside = (side == 'R' .or. side == 'r')
#if !defined(STABLOG)
    ! U(:,i) = U(:,i) * D(i)   (V is allocated for finite T)
    call k_scale_cols(n, U, n, S)
    ! column-norm pre-pivot + unpivoted QR: U holds the reflectors and R, st_ipvt the permutation (ZGEQP3 semantics)
    call la_qrpp(n, n, U, n, st_ipvt, st_tau)
    call fault_perm()
    ! separate D (|R_ii|) from R: U(i,j) = U(i,j)/|R_ii|, j >= i
    call k_diagabs(n, U, n, S)
    call k_divrows(n, U, n, S)
    ph1 = z1
    call Pivot_phase(ph1, st_ipvt, n)                     ! permutation sign, host (st_ipvt is valid on the host)
    sgn = nint(real(ph1))
#else
    ! log-norm pre-sort: tmpnorm(i) = log(|U(:,i)|) + L(i), then the selection sort with column swaps of the CPU code
    call k_lognorm(n, U, n, S, st_nrm)
    call la_sync()
    ! D->H: n log norms for the host selection sort (STABLOG decompose, ADR-001 section 4)
    !$omp target update from(st_nrm(1:n))
    call sort_stablog(n, st_nrm, st_perm, nswap)
    if (fault == 1) then                                  ! negative control: wrong permutation
       kk = st_perm(1); st_perm(1) = st_perm(2); st_perm(2) = kk
    endif
    ! H->D: the n-entry permutation
    !$omp target update to(st_perm(1:n))
    sgn = merge(-1, 1, mod(nswap,2) == 1)
    call k_perm_cols(n, U, n, st_perm, sA, n)             ! sA(:,j) = U(:,perm(j))
    call k_copy_c(n*n, sA, U)
    if (rside) then
       call k_perm_rows(n, V, n, st_perm, sB, n)          ! swaps of the rows of V (side R)
    else
       call k_perm_cols(n, V, n, st_perm, sB, n)          ! swaps of the columns of V (side L)
    endif
    call k_perm_vec(n, S, st_perm, st_Lt)                 ! L(i) = L(perm(i))
    call k_copy_r(n, st_Lt, S)
    call la_geqrf(n, n, U, n, st_tau)                     ! unpivoted QR (CPU: ZGEQP3 with all columns fixed)
    call k_diagabs(n, U, n, st_Dl)
    call k_divrows(n, U, n, st_Dl)
#endif
    call k_phase(n, U, n, sgn, .not. rside, st_c)         ! Phase = prod U(i,i), permutation sign, conj for side L; beta = 1/Phase
    call k_scale_row1(n, U, n, st_c)                      ! scale the first row of R with 1/Phase
#if defined(STABLOG)
    call k_stablog_scale(n, U, n, S, st_Dl)               ! U(i,j) *= exp(L(j) - L(i)), j > i;  L(i) = log(D(i)) + L(i)
#else
    if (rside) then
       call k_perm_rows(n, V, n, st_ipvt, sB, n)          ! ZLAPMR forward
    else
       call k_perm_cols(n, V, n, st_ipvt, sB, n)          ! ZLAPMT forward
    endif
#endif
    ! V = R * V (side R) or V * R^dagger (side L); the triangular factor is multiplied as a full GEMM (R in sA, zeros below)
    call k_triu(n, U, n, sA, n)
    if (rside) then
       call la_gemm(LA_N, LA_N, n, n, n, z1, sA, n, sB, n, z0, V, n)
    else
       call la_gemm(LA_N, LA_C, n, n, n, z1, sB, n, sA, n, z0, V, n)
    endif
    ! U = Q = H_1 ... H_n applied to the unit matrix by ZUNMQR (sC is free here).  Not la_ungqr: with it the sweep's Precision Green
    ! was 4x (mean) / 11x (max) worse than on the CPU at N = 36 (isolation study: replacing only this step by host ZUNGQR or by
    ! ZUNMQR on 1 restores the CPU level at N = 36; at N = 144 the Newton-Schulz step below is needed in addition).
    call k_ident(n, sC, n)
    call la_ormqr('L', LA_N, n, n, n, U, n, st_tau, sC, n)
    call k_copy_c(n*n, sC, U)
    ! One Newton-Schulz step Q <- Q (3 - Q^H Q)/2 restores the orthogonality to O(defect^2).  Measured (ticket 206, 12x12, beta 4,
    ! same seeds): without it Precision Green mean 9.7e-11 / max 1.2e-6 (CPU oracle 4.1e-11 / 7.2e-7), with it 5.4e-11 / 3.8e-7;
    ! cuSOLVER's Q has an orthogonality defect about 5x that of LAPACK (1.6e-14 vs 3.3e-15 at N = 400) and the CGR solves amplify it
    ! (see docs/gpu/DISCREPANCIES.md).  ALF_GPU_NS=0 switches the step off (read once in gpu_stab_setup).
    if (ns_on) then
       call la_gemm(LA_C, LA_N, n, n, n, z1, U, n, U, n, z0, sA, n)
       call k_ns(n, sA)
       call la_gemm(LA_N, LA_N, n, n, n, z1, U, n, sA, n, z0, sB, n)
       call k_copy_c(n*n, sB, U)
    endif
    call k_scale_col1(n, U, n, st_c)                      ! scale the first column of U with Phase
#if defined(STABLOG)
    if (fault == 3) then                                  ! negative control: NaN in the new L(1)
       st_Lt(1) = transfer(9221120237041090560_int64, 0.d0)
       !$omp target update to(st_Lt(1:1))
       call k_copy_r(1, st_Lt, S)
    endif
#endif
  end subroutine dev_decompose

  !> Fault injection (negative control): a wrong pivot permutation, host and device copy (ALF_GPU_FAULT_STAB=perm)
  Subroutine fault_perm()
    Implicit none
    Integer(c_int) :: t
    if (fault /= 1) return
    t = st_ipvt(1); st_ipvt(1) = st_ipvt(2); st_ipvt(2) = t
    !$omp target update to(st_ipvt(1:s_N))
  end subroutine fault_perm

  !> Stable selection sort of the STABLOG CPU decompose: perm = composite of the column swaps, nswap = number of swaps.
  Subroutine sort_stablog(n, tn, perm, nswap)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(inout) :: tn(n)
    Integer(c_int), intent(out) :: perm(n)
    Integer, intent(out) :: nswap
    Integer :: i, j, pvt
    Integer(c_int) :: t
    do i = 1, n
       perm(i) = i
    enddo
    nswap = 0
    do i = 1, n
       pvt = i
       do j = i+1, n
          if (tn(j) > tn(pvt)) pvt = j
       enddo
       if (pvt /= i) then
          t = perm(i); perm(i) = perm(pvt); perm(pvt) = t
          tn(pvt) = tn(i)
          nswap = nswap + 1
       endif
    enddo
  end subroutine sort_stablog

  !> testscale_UDV_state for the default stabilization: NaN scale -> error, extreme scales -> warning (once).
  !> One D->H of three numbers.
  Subroutine dev_testscale(D)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Integer :: n
    n = s_N
    call k_testscale(n, D, st_r)
    call la_sync()
    ! D->H: NaN count and the extreme scales (testscale of UDV_State, 24 B)
    !$omp target update from(st_r(1:3))
    if (st_r(1) > 0.5d0 .or. fault == 3) then
       ! ticket 212: inside the self-test the NaN scale is a detected failure: print the FAILED line before the abort
       if (in_selftest) call selftest_fail('stab: NaN scale in decompose')
       write(error_unit,*)
       write(error_unit,*) "Error: At least one scale is NaN."
       write(error_unit,*) "       Switch to LOG is required."
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    endif
    if (st_r(2) > 0.1d0*huge(1.d0) .and. warn_once) then
       write(error_unit,*)
       write(error_unit,*) "Warning: Largest scale is approaching the largest representable value."
       write(error_unit,*) "         Consider switching to LOG."
       warn_once = .false.
    endif
    if (st_r(3) < 10.d0*tiny(1.d0) .and. warn_once) then
       write(error_unit,*)
       write(error_unit,*) "Warning: Smallest scale is approaching the smallest representable value."
       write(error_unit,*) "         Consider switching to LOG."
       warn_once = .false.
    endif
  end subroutine dev_testscale

!====================================================================
!  CGR
!====================================================================
!> CGR(Phase, NVAR, GRUP, udvr, udvl) of cgr1_mod for the device states of flavour nfe (udvl = slot 1, udvr = slot 2).
!> G(:,:) (host) receives the new Green function: INTERIM 206, GR is host-authoritative until 208.  With_prec: the
!> precision check of Control_PrecisionG against the incoming G (the wrapped one), as a device reduction.
  Subroutine gpu_stab_cgr(nfe, nf, nvar, G, phase, with_prec)
    Implicit none
    Integer, intent(in) :: nfe, nf, nvar
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(Ndim,Ndim,*)
    Complex (Kind=Kind(0.d0)), intent(out) :: phase
    Logical, intent(in) :: with_prec
    Integer :: n
    Real (Kind=Kind(0.d0)) :: xmax, xsum
    Integer :: nna, nnb

    n = s_N
    if (with_prec) then
       st_Gt = G(:,:,nf)
       ! INTERIM 206: removed by 207/208: H->D of the wrapped G (host-authoritative), 16 N^2 B, for the precision reduction
       !$omp target update to(st_Gt)
    endif
    call dev_cgr(nfe, nvar)
    call la_sync()
    ! D->H: the phase (16 B)
    !$omp target update from(st_c(3:3))
    phase = st_c(3)
    ! INTERIM 206: removed by 207/208: D->H of the new G, 16 N^2 B (GR is host-authoritative: observables on the host)
    !$omp target update from(st_Gn)
    G(:,:,nf) = st_Gn
    if (chk_level >= 1) call check_cgr(nfe, nvar, phase, G(:,:,nf))
    if (with_prec) then
       call dev_prec(n, st_Gn, st_Gt, xmax, xsum, nna, nnb)
       call Control_PrecisionG_vals(xmax, xsum, nna, nnb, n)
    endif
  end subroutine gpu_stab_cgr

  !> Control_PrecisionG as a reduction: xmax = max |A-B|, xsum = sum |A-B| (xmean = xsum/N^2), NaN entries of A and B.
  !> NaN entries do not enter xmax/xsum.  One D->H of four numbers (INTERIM 206 hook of control_mod).
  Subroutine dev_prec(n, A, B, xmax, xsum, nna, nnb)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: A(n,*), B(n,*)
    Real (Kind=Kind(0.d0)), intent(out) :: xmax, xsum
    Integer, intent(out) :: nna, nnb
    call k_prec(n, A, B, st_r)
    call la_sync()
    ! D->H: the four reduced numbers of the precision check (32 B)
    !$omp target update from(st_r(4:7))
    xmax = st_r(4); xsum = st_r(5); nna = nint(st_r(6)); nnb = nint(st_r(7))
  end subroutine dev_prec

!--------------------------------------------------------------------
  !> ALF_GPU_STABCHECK: host CGR of the same device states vs the device G (st_Gn, valid on the host) and phase
  Subroutine check_cgr(nfe, nvar, phase, Gout)
    Implicit none
    Integer, intent(in) :: nfe, nvar
    Complex (Kind=Kind(0.d0)), intent(inout) :: phase
    Complex (Kind=Kind(0.d0)), intent(inout) :: Gout(:,:)
    Complex (Kind=Kind(0.d0)), Allocatable :: Gh(:,:)
    Complex (Kind=Kind(0.d0)) :: ph
    Real (Kind=Kind(0.d0)) :: d
    Character (len=256) :: line
    allocate(Gh(s_N,s_N))
    call dl_state(1, 1, nfe, ck_l(nfe)); call dl_state(1, 2, nfe, ck_r(nfe))
    call CGR(ph, nvar, Gh, ck_r(nfe), ck_l(nfe))
    d = rel_diff(st_Gn, Gh)
    chk_ncgr = chk_ncgr + 1
    call mx(chk_g(nvar), d); call mx(chk_ph, abs(phase - ph))
    if (chk_mix == 1) then
       Gout = Gh; phase = ph        ! experiment: the host CGR result replaces the device one
    endif
    if (.not. (d <= 1.d-8)) then
       chk_nbad = chk_nbad + 1
       if (chk_nbad <= 20) then
          write(line,'(A,I0,A,I0,A,I0,A,ES10.3)') 'GPU stab check: CGR call ', chk_ncgr, ' flavour ', nfe, ' NVAR ', nvar, &
               & ': device G vs host CGR rel diff ', d
          call info_line(line)
       endif
    endif
    deallocate(Gh)
  end subroutine check_cgr

  !> max |U^dagger U - 1| of the device slot U (self-test diagnostic)
  Real (Kind=Kind(0.d0)) function dev_orth(slot, nfe)
    Implicit none
    Integer, intent(in) :: slot, nfe
    call la_gemm(LA_C, LA_N, s_N, s_N, s_N, z1, udw_U(1,1,slot,nfe), s_N, udw_U(1,1,slot,nfe), s_N, z0, sC, s_N)
    call k_orth(s_N, sC, st_r)
    call la_sync()
    !$omp target update from(st_r(8:8))
    dev_orth = st_r(8)
  end function dev_orth

  Real (Kind=Kind(0.d0)) function host_orth(U)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: U(:,:)
    Complex (Kind=Kind(0.d0)), Allocatable :: P(:,:)
    Integer :: i
    allocate(P(size(U,2),size(U,2)))
    P = matmul(conjg(transpose(U)), U)
    do i = 1, size(P,1)
       P(i,i) = P(i,i) - 1.d0
    enddo
    host_orth = maxval(abs(P))
  end function host_orth

  !> ALF_GPU_STABCHECK=2: compare a device block (slot, after the call) with the host WRAPUR / WRAPUL on the saved input
  Subroutine check_block_in(slot)
    Implicit none
    Integer, intent(in) :: slot
    Integer :: nfe
    do nfe = 1, s_nfe
       if (slot == 2) then
          call dl_state(1, 2, nfe, ck_r(nfe))
       else
          call dl_state(1, 1, nfe, ck_l(nfe))
       endif
    enddo
  end subroutine check_block_in

  Subroutine check_block_out(slot, ntau, ntau1)
    Implicit none
    Integer, intent(in) :: slot, ntau, ntau1
    Integer :: nfe
    Real (Kind=Kind(0.d0)) :: dU, dV, dS
    if (slot == 2) then
       call WRAPUR(ntau, ntau1, ck_r)
    else
       call WRAPUL(ntau1, ntau, ck_l)
    endif
    do nfe = 1, s_nfe
       call dl_state(1, slot, nfe, ck_d(nfe))
       if (slot == 2) then
          call state_diff(ck_d(nfe), ck_r(nfe), dU, dV, dS)
       else
          call state_diff(ck_d(nfe), ck_l(nfe), dU, dV, dS)
       endif
       call mx(chk_u, dU); call mx(chk_v, dV); call mx(chk_s, dS)
       if (chk_mix == 2) then
          if (slot == 2) then
             call ul_state(2, nfe, ck_r(nfe))      ! experiment: the host block result replaces the device one
          else
             call ul_state(1, nfe, ck_l(nfe))
          endif
       endif
    enddo
    chk_nblk = chk_nblk + 1
  end subroutine check_block_out

  Subroutine check_report()
    Implicit none
    Character (len=512) :: line
    if (chk_level == 0) return
    write(line,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A,I0)') 'GPU stab check (device vs host upstream routines on the same device states): CGR calls ', &
         & chk_ncgr, ', blocks ', chk_nblk, ': max rel diff G NVAR=1 ', chk_g(1), ' NVAR=2 ', chk_g(2), ' phase ', chk_ph, &
         & '; CGR calls with diff > 1e-8: ', chk_nbad, ' of ', chk_ncgr
    call say(line)
    if (chk_nblk > 0) then
       write(line,'(A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU stab check blocks: max rel diff U ', chk_u, ' V ', chk_v, ' scales ', chk_s
       call say(line)
    endif
  end subroutine check_report

!--------------------------------------------------------------------
  Subroutine dev_cgr(nfe, nvar)
    Implicit none
    Integer, intent(in) :: nfe, nvar
    Integer :: n

    n = s_N
    ! RHS = U_R^dagger U_L  (sA),  T = V_R V_L  (sB)
    call la_gemm(LA_C, LA_N, n, n, n, z1, udw_U(1,1,2,nfe), n, udw_U(1,1,1,nfe), n, z0, sA, n)
    call la_gemm(LA_N, LA_N, n, n, n, z1, udw_V(1,1,2,nfe), n, udw_V(1,1,1,nfe), n, z0, sB, n)
#if defined(STABLOG)
    call k_tpup_log(n, sB, sA, udw_L(1,2,nfe), udw_L(1,1,nfe))
#else
    call k_tpup(n, sB, sA, udw_D(1,2,nfe), udw_D(1,1,nfe))
#endif
    ! NOTE: upstream's one-time "Large number encountered" warning of CGR is omitted here (the scale test of decompose warns instead).
    ! det(U_R^dagger U_L): LU of RHS; Phase = conj(det)/|det|  (st_c(3))
    call la_getrf(n, n, sA, n, st_ipiv)
    call k_det_phase(n, sA, n, st_ipiv, st_c)
    if (nvar == 1) then
       call cgr_solve(nfe, nvar, sB)
    else
       call k_ct(n, sB, n, sC, n)                          ! T = T^dagger
       call cgr_solve(nfe, nvar, sC)
    endif
  end subroutine dev_cgr

  !> The part of CGR after the matrix T = TPUP is built: QR, phase, solve.  T is destroyed.  Result in st_Gn, phase st_c(3).
  Subroutine cgr_solve(nfe, nvar, T)
    Implicit none
    Integer, intent(in) :: nfe, nvar
    Complex (Kind=Kind(0.d0)), intent(inout) :: T(s_N,*)
    Integer :: n, sgn
    Complex (Kind=Kind(0.d0)) :: ph1

    n = s_N
#if defined(STABLOG)
    ! STABLOG: host ZGEQP3 (ADR-001 section 4): D->H of T, H->D of the factors, tau, ipvt (inside la_geqp3_host)
    call la_geqp3_host(n, n, T, n, st_ipvt, st_tau)
    call fault_perm()
#else
    call la_qrpp(n, n, T, n, st_ipvt, st_tau)
    call fault_perm()
#endif
    call k_diagabs(n, T, n, st_Dl)
    call k_divrows(n, T, n, st_Dl)
    ph1 = z1
    call Pivot_phase(ph1, st_ipvt, n)
    sgn = nint(real(ph1))
    call k_cgr_phase(n, T, n, st_tau, st_c, nvar, sgn, merge(1, 0, fault == 2))
    if (nvar == 1) then
       call k_ct(n, udw_U(1,1,2,nfe), n, sA, n)             ! RHS = U_R^dagger
#if defined(STABLOG)
       call k_rowscale_log(n, sA, n, udw_L(1,2,nfe))
#endif
       call la_ormqr('L', LA_C, n, n, n, T, n, st_tau, sA, n)
       call k_rowdiv(n, sA, n, st_Dl)
       call la_trsm('L', 'U', LA_N, 'N', n, n, z1, T, n, sA, n)
       call k_iperm_rows(n, sA, n, st_ipvt, sC, n)           ! ZLAPMR backward  (sC is free: T is sB)
#if defined(STABLOG)
       call k_rowscale_log(n, sC, n, udw_L(1,1,nfe))
#endif
       call la_gemm(LA_N, LA_N, n, n, n, z1, udw_U(1,1,1,nfe), n, sC, n, z0, st_Gn, n)
    else
       call k_copy_c(n*n, udw_U(1,1,1,nfe), sA)              ! RHS = U_L
#if defined(STABLOG)
       call k_colscale_log(n, sA, n, udw_L(1,1,nfe))
#endif
       call la_ormqr('R', LA_N, n, n, n, T, n, st_tau, sA, n)
       call k_colinv(n, sA, n, st_Dl)
       call la_trsm('R', 'U', LA_C, 'N', n, n, z1, T, n, sA, n)
       call k_iperm_cols(n, sA, n, st_ipvt, sB, n)           ! ZLAPMT backward  (sB is free: T is sC)
#if defined(STABLOG)
       call k_colscale_log(n, sB, n, udw_L(1,2,nfe))
#endif
       call la_gemm(LA_N, LA_C, n, n, n, z1, sB, n, udw_U(1,1,2,nfe), n, z0, st_Gn, n)
    endif
  end subroutine cgr_solve

!====================================================================
!  Device kernels (plain arrays; synchronous target regions)
!====================================================================
  ! U(:,j) = U(:,j) * D(j)
  Subroutine k_scale_cols(n, U, ld, D)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          U(i,j) = U(i,j)*D(j)
       enddo
    enddo
  end subroutine k_scale_cols

  ! D(i) = |U(i,i)| (complex, zero imaginary part)
  Subroutine k_diagabs(n, U, ld, D)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       D(i) = cmplx(abs(U(i,i)), 0.d0, kind(0.d0))
    enddo
  end subroutine k_diagabs

  ! U(i,j) = U(i,j) / D(i), j >= i   (D real valued)
  Subroutine k_divrows(n, U, ld, D)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: x
    !$omp target teams distribute parallel do collapse(2) private(x)
    do j = 1, n
       do i = 1, n
          if (j >= i) then
             x = real(D(i))
             U(i,j) = U(i,j) / x
          endif
       enddo
    enddo
  end subroutine k_divrows

  ! Phase = prod_i U(i,i), times the permutation sign, conjugated for side L; cs(1) = Phase, cs(2) = 1/Phase.
  ! One serial loop on the device (N complex products, no transfer).
  Subroutine k_phase(n, U, ld, sgn, cj, cs)
    Implicit none
    Integer, intent(in) :: n, ld, sgn
    Logical, intent(in) :: cj
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(*)
    Complex (Kind=Kind(0.d0)) :: ph
    Integer :: i
    !$omp target
    ph = cmplx(1.d0, 0.d0, kind(0.d0))
    do i = 1, n
       ph = ph*U(i,i)
    enddo
    if (sgn < 0) ph = -ph
    if (cj) ph = conjg(ph)
    cs(1) = ph
    cs(2) = cmplx(1.d0, 0.d0, kind(0.d0)) / ph
    !$omp end target
  end subroutine k_phase

  ! U(1,j) = cs(2) * U(1,j)   (ZSCAL of the first row of R with 1/Phase)
  Subroutine k_scale_row1(n, U, ld, cs)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: cs(*)
    Integer :: j
    !$omp target teams distribute parallel do
    do j = 1, n
       U(1,j) = cs(2)*U(1,j)
    enddo
  end subroutine k_scale_row1

  ! U(i,1) = cs(1) * U(i,1)   (ZSCAL of the first column of Q with Phase)
  Subroutine k_scale_col1(n, U, ld, cs)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: cs(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       U(i,1) = cs(1)*U(i,1)
    enddo
  end subroutine k_scale_col1

  ! R(i,j) = U(i,j) for i <= j, else 0
  Subroutine k_triu(n, U, ld, R, ldr)
    Implicit none
    Integer, intent(in) :: n, ld, ldr
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: R(ldr,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (i <= j) then
             R(i,j) = U(i,j)
          else
             R(i,j) = z0
          endif
       enddo
    enddo
  end subroutine k_triu

  ! W(i,j) = X(perm(i),j)        (ZLAPMR forward)
  Subroutine k_perm_rows(n, X, ldx, perm, W, ldw)
    Implicit none
    Integer, intent(in) :: n, ldx, ldw
    Integer(c_int), intent(in) :: perm(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          W(i,j) = X(perm(i),j)
       enddo
    enddo
  end subroutine k_perm_rows

  ! W(i,j) = X(i,perm(j))        (ZLAPMT forward)
  Subroutine k_perm_cols(n, X, ldx, perm, W, ldw)
    Implicit none
    Integer, intent(in) :: n, ldx, ldw
    Integer(c_int), intent(in) :: perm(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          W(i,j) = X(i,perm(j))
       enddo
    enddo
  end subroutine k_perm_cols

  ! W(perm(i),j) = X(i,j)        (ZLAPMR backward)
  Subroutine k_iperm_rows(n, X, ldx, perm, W, ldw)
    Implicit none
    Integer, intent(in) :: n, ldx, ldw
    Integer(c_int), intent(in) :: perm(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          W(perm(i),j) = X(i,j)
       enddo
    enddo
  end subroutine k_iperm_rows

  ! W(i,perm(j)) = X(i,j)        (ZLAPMT backward)
  Subroutine k_iperm_cols(n, X, ldx, perm, W, ldw)
    Implicit none
    Integer, intent(in) :: n, ldx, ldw
    Integer(c_int), intent(in) :: perm(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          W(i,perm(j)) = X(i,j)
       enddo
    enddo
  end subroutine k_iperm_cols

  ! W(i) = L(perm(i))
  Subroutine k_perm_vec(n, L, perm, W)
    Implicit none
    Integer, intent(in) :: n
    Integer(c_int), intent(in) :: perm(*)
    Real (Kind=Kind(0.d0)), intent(in) :: L(*)
    Real (Kind=Kind(0.d0)), intent(inout) :: W(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       W(i) = L(perm(i))
    enddo
  end subroutine k_perm_vec

  ! B(i,j) = conjg(A(j,i))
  Subroutine k_ct(n, A, lda, B, ldb)
    Implicit none
    Integer, intent(in) :: n, lda, ldb
    Complex (Kind=Kind(0.d0)), intent(in) :: A(lda,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: B(ldb,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          B(i,j) = conjg(A(j,i))
       enddo
    enddo
  end subroutine k_ct

  ! A(i,j) = A(i,j) / D(i)   (complex division, as the CPU code)
  Subroutine k_rowdiv(n, A, lda, D)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          A(i,j) = A(i,j) / D(i)
       enddo
    enddo
  end subroutine k_rowdiv

  ! A(i,j) = A(i,j) * (1/Re D(j))
  Subroutine k_colinv(n, A, lda, D)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: sv
    !$omp target teams distribute parallel do collapse(2) private(sv)
    do j = 1, n
       do i = 1, n
          sv = 1.d0 / real(D(j))
          A(i,j) = A(i,j) * sv
       enddo
    enddo
  end subroutine k_colinv

  ! default CGR: T(i,j) = Dr(i) * T(i,j) * Dl(j) + RHS(i,j)
  Subroutine k_tpup(n, T, R, Dr, Dl)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: T(n,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: R(n,*), Dr(*), Dl(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          T(i,j) = Dr(i) * T(i,j) * Dl(j)
          T(i,j) = T(i,j) + R(i,j)
       enddo
    enddo
  end subroutine k_tpup

  ! STABLOG CGR: the four cases of cgr1_mod (exponent signs chosen so that no exp overflows)
  Subroutine k_tpup_log(n, T, R, Lr, Ll)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: T(n,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: R(n,*)
    Real (Kind=Kind(0.d0)), intent(in) :: Lr(*), Ll(*)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: lr_i, ll_j
    Complex (Kind=Kind(0.d0)) :: dup, dlj
    !$omp target teams distribute parallel do collapse(2) private(lr_i, ll_j, dup, dlj)
    do j = 1, n
       do i = 1, n
          lr_i = Lr(i); ll_j = Ll(j)
          if (lr_i <= 0.d0) then
             dup = cmplx(exp(lr_i), 0.d0, kind(0.d0))
          else
             dup = cmplx(exp(-lr_i), 0.d0, kind(0.d0))
          endif
          if (ll_j <= 0.d0) then
             dlj = cmplx(exp(ll_j), 0.d0, kind(0.d0))
             if (lr_i <= 0.d0) then
                T(i,j) = R(i,j) + cmplx(exp(lr_i+ll_j), 0.d0, kind(0.d0))*T(i,j)
             else
                T(i,j) = dup*R(i,j) + dlj*T(i,j)
             endif
          else
             dlj = cmplx(exp(-ll_j), 0.d0, kind(0.d0))
             if (lr_i <= 0.d0) then
                T(i,j) = dlj*R(i,j) + dup*T(i,j)
             else
                T(i,j) = cmplx(exp(-lr_i-ll_j), 0.d0, kind(0.d0))*R(i,j) + T(i,j)
             endif
          endif
       enddo
    enddo
  end subroutine k_tpup_log

  ! rows: A(j,:) *= exp(-L(j)) if L(j) > 0
  Subroutine k_rowscale_log(n, A, lda, L)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,*)
    Real (Kind=Kind(0.d0)), intent(in) :: L(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (L(i) > 0.d0) A(i,j) = cmplx(exp(-L(i)), 0.d0, kind(0.d0)) * A(i,j)
       enddo
    enddo
  end subroutine k_rowscale_log

  ! columns: A(:,j) *= exp(-L(j)) if L(j) > 0
  Subroutine k_colscale_log(n, A, lda, L)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,*)
    Real (Kind=Kind(0.d0)), intent(in) :: L(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (L(j) > 0.d0) A(i,j) = cmplx(exp(-L(j)), 0.d0, kind(0.d0)) * A(i,j)
       enddo
    enddo
  end subroutine k_colscale_log

  ! cs(3) = conj(det A)/|det A| for the LU factors A, ipiv (DET_C of mymats_mod): serial loop
  Subroutine k_det_phase(n, A, lda, ipiv, cs)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(in) :: A(lda,*)
    Integer(c_int), intent(in) :: ipiv(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(*)
    Complex (Kind=Kind(0.d0)) :: det
    Integer :: i, sg
    !$omp target
    det = cmplx(1.d0, 0.d0, kind(0.d0))
    do i = 1, n
       det = det*A(i,i)
    enddo
    sg = 1
    do i = 1, n
       if (ipiv(i) /= i) sg = -sg
    enddo
    if (sg == -1) det = -det
    det = conjg(det)
    cs(3) = det / abs(det)
    !$omp end target
  end subroutine k_det_phase

  ! cs(3) *= permutation sign * prod_i (R_ii phase) * reflector determinants   (cgr1_mod, after the QR); serial loop.
  ! bad = 1: fault injection (negative control), phase sign flipped
  Subroutine k_cgr_phase(n, T, ld, tau, cs, nvar, sgn, bad)
    Implicit none
    Integer, intent(in) :: n, ld, nvar, sgn, bad
    Complex (Kind=Kind(0.d0)), intent(in) :: T(ld,*), tau(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(*)
    Complex (Kind=Kind(0.d0)) :: ph, z
    Real (Kind=Kind(0.d0)) :: x
    Integer :: i
    !$omp target
    ph = cs(3)
    if (sgn < 0) ph = -ph
    do i = 1, n
       z = tau(i)
       if (nvar == 1) then
          ph = ph * T(i,i)/abs(T(i,i))
       else
          ph = ph * conjg(T(i,i))/abs(T(i,i))
          z = conjg(z)
       endif
       if (z /= cmplx(0.d0, 0.d0, kind(0.d0))) then
          x = abs(z)
          z = 1.d0 - 2.d0 * (z/x) * (real(z)/x)
          ph = ph * z/abs(z)
       endif
    enddo
    if (bad == 1) ph = -ph
    cs(3) = ph
    !$omp end target
  end subroutine k_cgr_phase

  ! A <- 3/2 - A/2  (A = Q^H Q on input)
  Subroutine k_ns(n, A)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(n,*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          A(i,j) = -0.5d0*A(i,j)
          if (i == j) A(i,j) = A(i,j) + 1.5d0
       enddo
    enddo
  end subroutine k_ns

  ! r(8) = max |A - 1|
  Subroutine k_orth(n, A, r)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: A(n,*)
    Real (Kind=Kind(0.d0)), intent(inout) :: r(*)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: xm, d
    xm = 0.d0
    !$omp target teams distribute parallel do collapse(2) reduction(max:xm) private(d)
    do j = 1, n
       do i = 1, n
          if (i == j) then
             d = abs(A(i,j) - cmplx(1.d0, 0.d0, kind(0.d0)))
          else
             d = abs(A(i,j))
          endif
          xm = max(xm, d)
       enddo
    enddo
    !$omp target
    r(8) = xm
    !$omp end target
  end subroutine k_orth

  ! scale test: r(1) = number of NaN scales, r(2) = Re D(1), r(3) = Re D(n)
  Subroutine k_testscale(n, D, r)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Real (Kind=Kind(0.d0)), intent(inout) :: r(*)
    Integer :: i, nn
    Real (Kind=Kind(0.d0)) :: a, b
    nn = 0
    !$omp target teams distribute parallel do reduction(+:nn) private(a, b)
    do i = 1, n
       a = real(D(i)); b = aimag(D(i))
       if (isnan_bits(a) .or. isnan_bits(b)) nn = nn + 1
    enddo
    !$omp target
    r(1) = real(nn, kind(0.d0))
    r(2) = real(D(1))
    r(3) = real(D(n))
    !$omp end target
  end subroutine k_testscale

  ! precision check: r(4) = max |A-B|, r(5) = sum |A-B|, r(6)/r(7) = number of NaN entries in A / B (not in max/sum)
  Subroutine k_prec(n, A, B, r)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: A(n,*), B(n,*)
    Real (Kind=Kind(0.d0)), intent(inout) :: r(*)
    Integer :: i, j, na, nb
    Real (Kind=Kind(0.d0)) :: xm, xs, d, cs
    Logical :: fa, fb
    xm = 0.d0; xs = 0.d0; na = 0; nb = 0
    ! ticket 212: the sum has a fixed order: one thread per column sums serially into st_pw(j), then one thread adds the
    ! n column sums in order (max and the integer counts are order independent and stay reductions)
    !$omp target teams distribute parallel do reduction(max:xm) reduction(+:na, nb) private(i, d, fa, fb, cs)
    do j = 1, n
       cs = 0.d0
       do i = 1, n
          fa = isnan_bits(real(A(i,j))) .or. isnan_bits(aimag(A(i,j)))
          fb = isnan_bits(real(B(i,j))) .or. isnan_bits(aimag(B(i,j)))
          if (fa) na = na + 1
          if (fb) nb = nb + 1
          if (.not. (fa .or. fb)) then
             d = abs(A(i,j) - B(i,j))
             cs = cs + d
             xm = max(xm, d)
          endif
       enddo
       st_pw(j) = cs
    enddo
    !$omp target
    xs = 0.d0
    do j = 1, n
       xs = xs + st_pw(j)
    enddo
    r(5) = xs
    !$omp end target
    !$omp target
    r(4) = xm
    r(6) = real(na, kind(0.d0))
    r(7) = real(nb, kind(0.d0))
    !$omp end target
  end subroutine k_prec

  Logical function isnan_bits(x)
    Implicit none
    !$omp declare target
    Real (Kind=Kind(0.d0)), intent(in) :: x
    Integer(int64) :: b
    b = transfer(x, 0_int64)
    isnan_bits = (iand(ishft(b, -52), 2047_int64) == 2047_int64) .and. (iand(b, 4503599627370495_int64) /= 0_int64)
  end function isnan_bits

#if defined(STABLOG)
  ! out(j) = log(|U(:,j)|_2) + L(j)   (the CPU's tmpnorm); scaled norm as in alf_la
  Subroutine k_lognorm(n, U, ld, L, out)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,*)
    Real (Kind=Kind(0.d0)), intent(in) :: L(*)
    Real (Kind=Kind(0.d0)), intent(inout) :: out(*)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: mx, s, v
    ! ticket 212: one thread per column, serial sums (fixed order; a device reduction combines in an unspecified order)
    !$omp target teams distribute parallel do private(i, mx, s, v)
    do j = 1, n
       mx = 0.d0
       do i = 1, n
          mx = max(mx, abs(real(U(i,j))), abs(aimag(U(i,j))))
       enddo
       s = 0.d0
       if (mx > 0.d0) then
          do i = 1, n
             v = real(U(i,j)) / mx
             s = s + v*v
             v = aimag(U(i,j)) / mx
             s = s + v*v
          enddo
          out(j) = log(mx*sqrt(s)) + L(j)
       else
          out(j) = -huge(1.d0)           ! CPU: log(0) = -Inf
       endif
    enddo
  end subroutine k_lognorm

  ! after the QR: U(i,j) *= exp(L(j) - L(i)) for j > i, then L(i) = log(Re D(i)) + L(i)
  Subroutine k_stablog_scale(n, U, ld, L, D)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Real (Kind=Kind(0.d0)), intent(inout) :: L(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (j > i) U(i,j) = U(i,j) * cmplx(exp(L(j)-L(i)), 0.d0, kind(0.d0))
       enddo
    enddo
    !$omp target teams distribute parallel do
    do i = 1, n
       L(i) = log(real(D(i))) + L(i)
    enddo
  end subroutine k_stablog_scale
#endif

!====================================================================
!  Hook: self-test, timing
!====================================================================
  ! ticket 212: all ALF_GPU_* variables are read once by gpu_state_mod (gpu_env_get / gpu_env_is)

  Subroutine gpu_stab_hook()
    Implicit none
    Character (len=64) :: f
    Integer :: st, ib
    Logical :: did
    did = .false.
    if (.not. gpu_stab_on) return
    if (gpu_env_is('ALF_GPU_STABCHECK', '1')) chk_level = 1
    if (gpu_env_is('ALF_GPU_STABCHECK', '2')) chk_level = 2
    if (gpu_env_is('ALF_GPU_STABMIX', '1')) then
       chk_level = max(chk_level, 1); chk_mix = 1
    endif
    if (gpu_env_is('ALF_GPU_STABMIX', '2')) then
       chk_level = 2; chk_mix = 2
    endif
    if (chk_level > 0) then
       allocate(ck_l(s_nfe), ck_r(s_nfe), ck_d(s_nfe))
       do st = 1, s_nfe
          call ck_l(st)%init(s_N, 'l'); call ck_r(st)%init(s_N, 'r'); call ck_d(st)%init(s_N, 'r')
       enddo
    endif
    ! values were validated (and the armed line printed) by gpu_env_check; eps:<x> acts in the SELFTEST=3 comparison only
    if (gpu_env_get('ALF_GPU_FAULT_STAB', f)) then
       select case (trim(f))
       case ('perm');  fault = 1
       case ('phase'); fault = 2
       case ('nan');   fault = 3
       case default                  ! eps:<x>
       end select
    endif
    if (gpu_env_is('ALF_GPU_SELFTEST', '3')) then
       call gpu_stab_selftest()
       did = .true.
    endif
    if (gpu_env_is('ALF_GPU_STABTIME', '1')) then
       call gpu_stab_timing()
       did = .true.
    endif
    if (did) then
       ! the tests leave the device slots dirty: restore the state before the hook (identity states, as gpu_stab_setup)
       do st = 1, s_nfe
          call gpu_stab_reset_l(st); call gpu_stab_reset_r(st)
          do ib = 1, s_nstm
             call gpu_stab_reset_st(ib, st)
          enddo
       enddo
    endif
    if (did .and. gpu_env_is('ALF_GPU_EXIT', '1')) then
       call gpu_stab_finalize()
       call gpu_state_finalize()
       write(output_unit,'(A)') ' GPU stabilization hook: ALF_GPU_EXIT=1, stopping after the test'
       stop
    endif
  end subroutine gpu_stab_hook

  Subroutine say(line)
    Implicit none
    Character (len=*), intent(in) :: line
    call info_line(line)
    write(output_unit,'(1X,A)') trim(line)
  end subroutine say

  ! download / upload the state of a device slot to / from a host UDV_State object (self-test only)
  Subroutine dl_state(kind_, slot, nfe, ust)
    Implicit none
    Integer, intent(in) :: kind_, slot, nfe     ! kind_ 1 = udw slot, 2 = udv slot
    Class (UDV_State), intent(inout) :: ust
    if (kind_ == 1) then
       !$omp target update from(udw_U(:,:,slot,nfe), udw_V(:,:,slot,nfe))
       ust%U = udw_U(:,:,slot,nfe); ust%V = udw_V(:,:,slot,nfe)
#if defined(STABLOG)
       !$omp target update from(udw_L(:,slot,nfe))
       ust%L = udw_L(:,slot,nfe)
#else
       !$omp target update from(udw_D(:,slot,nfe))
       ust%D = udw_D(:,slot,nfe)
#endif
    else
       !$omp target update from(udv_U(:,:,slot,nfe), udv_V(:,:,slot,nfe))
       ust%U = udv_U(:,:,slot,nfe); ust%V = udv_V(:,:,slot,nfe)
#if defined(STABLOG)
       !$omp target update from(udv_L(:,slot,nfe))
       ust%L = udv_L(:,slot,nfe)
#else
       !$omp target update from(udv_D(:,slot,nfe))
       ust%D = udv_D(:,slot,nfe)
#endif
    endif
  end subroutine dl_state

  Subroutine ul_state(slot, nfe, ust)
    Implicit none
    Integer, intent(in) :: slot, nfe
    Class (UDV_State), intent(in) :: ust
    udw_U(:,:,slot,nfe) = ust%U; udw_V(:,:,slot,nfe) = ust%V
    !$omp target update to(udw_U(:,:,slot,nfe), udw_V(:,:,slot,nfe))
#if defined(STABLOG)
    udw_L(:,slot,nfe) = ust%L
    !$omp target update to(udw_L(:,slot,nfe))
#else
    udw_D(:,slot,nfe) = ust%D
    !$omp target update to(udw_D(:,slot,nfe))
#endif
  end subroutine ul_state

  ! relative difference of two states: U, V (max |a-b| / max |b|), scales (componentwise relative, D; L absolute)
  Subroutine state_diff(a, b, dU, dV, dS)
    Implicit none
    Class (UDV_State), intent(in) :: a, b
    Real (Kind=Kind(0.d0)), intent(out) :: dU, dV, dS
    Integer :: i
    Real (Kind=Kind(0.d0)) :: d
    dU = rel_diff(a%U, b%U)
    dV = rel_diff(a%V, b%V)
    dS = 0.d0
#if defined(STABLOG)
    do i = 1, size(a%L)
       if (.not. (is_finite(a%L(i)) .and. is_finite(b%L(i)))) then
          dS = huge(1.d0)
       else
          d = abs(a%L(i) - b%L(i)); if (d > dS) dS = d
       endif
    enddo
#else
    do i = 1, size(a%D)
       if (.not. (finite_c(a%D(i)) .and. finite_c(b%D(i)))) then
          dS = huge(1.d0)
       else
          d = abs(a%D(i) - b%D(i)) / max(abs(b%D(i)), tiny(1.d0)); if (d > dS) dS = d
       endif
    enddo
#endif
  end subroutine state_diff

  Subroutine mx(m, d)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: m
    Real (Kind=Kind(0.d0)), intent(in) :: d
    if (.not. is_finite(d)) then
       m = huge(1.d0)
    else if (d > m) then
       m = d
    endif
  end subroutine mx

!--------------------------------------------------------------------
!> Self-test (ALF_GPU_SELFTEST=3): the fixed initial fields; one pass of the sweep.
!>   left chain  WRAPUL + decompose for NST = NSTM-1 .. 1, right chain WRAPUR + decompose for NST = 1 .. NSTM-1, and at every
!>   stabilization point CGR (NVAR = 1 and 2) from the states exchanged as in the sweep, each against the host routines
!>   (WRAPUL/WRAPUR/decompose/CGR of upstream, run on host UDV_State objects that hold the same input) and the precision
!>   reduction against COMPARE.  Every block starts from the same input on both sides.
!--------------------------------------------------------------------
  Subroutine gpu_stab_selftest()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-10, tolr = 1.d-12
    Class (UDV_State), Allocatable :: hl(:), hr(:), dl_(:), dr_(:), hp(:)
    Integer :: nfe, nf, nst, nvar, ib, N
    Integer :: NT1, NT
    Real (Kind=Kind(0.d0)) :: dU, dV, dS, mu(2), mv(2), ms(2), mg(2), mph(2), d, xmax_h, xmean_h, xmax_d, xsum_d
    Real (Kind=Kind(0.d0)) :: dpr(4), sens_h(2), sens_d(2), orth_d, orth_h
    Complex (Kind=Kind(0.d0)), Allocatable :: Gh(:,:), Gd(:,:,:), Gp(:,:), Gd2(:,:,:)
    Complex (Kind=Kind(0.d0)) :: ph_h, ph_d, ph_p
    Integer(int64) :: seed
    Integer :: nna, nnb, nblk(2)
    Character (len=512) :: line

#if !defined(STABLOG) && !defined(QR_PREPIVOT)
    write(error_unit,*) ' ALF_GPU_SELFTEST=3 with the default stabilization needs a build with -DQR_PREPIVOT (ALF_FLAGS_EXT=-DQR_PREPIVOT before configure.sh): the host oracle'
    write(error_unit,*) ' must be the same pre-pivoted algorithm (ZGEQP3 gives other, equally valid, pivots).'
    call info_line('GPU stabilization self-test: needs a build with ALF_FLAGS_EXT=-DQR_PREPIVOT for the default stabilization')
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
    N = s_N
    seed = 20251101_int64
    in_selftest = .true.
    call say('GPU stabilization self-test (ticket 206) start')
    allocate(hl(s_nfe), hr(s_nfe), dl_(s_nfe), dr_(s_nfe), hp(s_nfe))
    allocate(Gh(N,N), Gp(N,N), Gd(N,N,N_FL), Gd2(N,N,N_FL))
    sens_h = 0.d0; sens_d = 0.d0; orth_d = 0.d0; orth_h = 0.d0
    mu = 0.d0; mv = 0.d0; ms = 0.d0; mg = 0.d0; mph = 0.d0; dpr = 0.d0; nblk = 0
    do nfe = 1, s_nfe
       call hl(nfe)%init(N, 'l'); call hr(nfe)%init(N, 'r')
       call dl_(nfe)%init(N, 'l'); call dr_(nfe)%init(N, 'r'); call hp(nfe)%init(N, 'r')
       call gpu_stab_reset_l(nfe); call gpu_stab_reset_r(nfe)
       call gpu_stab_reset_st(s_nstm, nfe)
    enddo
    ! ---- left chain
    do nst = s_nstm-1, 1, -1
       NT1 = s_nt(nst+1); NT = s_nt(nst)
       do nfe = 1, s_nfe
          call dl_state(1, 1, nfe, hl(nfe))               ! self-test: the host oracle starts from the device input
       enddo
       call WRAPUL(NT1, NT, hl)
       call gpu_stab_wrapul(NT1, NT)
       do nfe = 1, s_nfe
          call dl_state(1, 1, nfe, dl_(nfe))               ! self-test: device result for the comparison
          call state_diff(dl_(nfe), hl(nfe), dU, dV, dS)
          call mx(mu(1), dU); call mx(mv(1), dV); call mx(ms(1), dS)
          call mx(orth_d, dev_orth(1, nfe)); call mx(orth_h, host_orth(hl(nfe)%U))
          call gpu_stab_store_l(nst, nfe)
       enddo
       nblk(1) = nblk(1) + 1
    enddo
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test WRAPUL+decompose over ', nblk(1), &
         & ' blocks: max rel diff U ', mu(1), ', V ', mv(1), ', scales ', ms(1), ' (tol ', tol
    call say(line)
    ! ---- right chain with CGR at every stabilization point
    do nst = 1, s_nstm
       NT1 = s_nt(nst-1); NT = s_nt(nst)
       if (.true.) then
          do nfe = 1, s_nfe
             call dl_state(1, 2, nfe, hr(nfe))
          enddo
          call WRAPUR(NT1, NT, hr)
          call gpu_stab_wrapur(NT1, NT)
          do nfe = 1, s_nfe
             call dl_state(1, 2, nfe, dr_(nfe))
             call state_diff(dr_(nfe), hr(nfe), dU, dV, dS)
             call mx(mu(2), dU); call mx(mv(2), dV); call mx(ms(2), dS)
             call mx(orth_d, dev_orth(2, nfe)); call mx(orth_h, host_orth(hr(nfe)%U))
          enddo
          nblk(2) = nblk(2) + 1
       endif
       ! CGR as in the sweep: udvl = udvst(NST), udvst(NST) = udvr (device), then compare with the host CGR on the same states
       do nfe = 1, s_nfe
          nf = Calc_Fl_map(nfe)
          call gpu_stab_exchange(nst, nfe, .true., hr(nfe))   ! (the host object is only used in the backward direction)
          call dl_state(1, 1, nfe, hl(nfe))
          call dl_state(1, 2, nfe, hr(nfe))
          do nvar = 1, 2
             call CGR(ph_h, nvar, Gh, hr(nfe), hl(nfe))
             call gpu_stab_cgr(nfe, nf, nvar, Gd, ph_d, .false.)
             call fault_eps_apply(Gd(:,:,nf))
             call state_gdiff(Gd(:,:,nf), Gh, d)
             call mx(mg(nvar), d)
             call mx(mph(nvar), abs(ph_d - ph_h))
             ! intrinsic sensitivity (reported, not a criterion): U_R perturbed by 2e-16 relative per entry, host and device CGR
             call perturb_state(hr(nfe), hp(nfe), seed)
             call CGR(ph_p, nvar, Gp, hp(nfe), hl(nfe))
             call mx(sens_h(nvar), rel_diff(Gp, Gh))
             call ul_state(2, nfe, hp(nfe))                      ! self-test: upload the perturbed udvr
             call gpu_stab_cgr(nfe, nf, nvar, Gd2, ph_p, .false.)
             call mx(sens_d(nvar), rel_diff(Gd2(:,:,nf), Gd(:,:,nf)))
             call ul_state(2, nfe, hr(nfe))                      ! self-test: restore udvr
             ! precision reduction vs COMPARE: (new G, reference G) and (new G, G + noise of order 1e-3)
             call rand_perturb(Gh, Gp, seed)
             call dev_check_prec(Gd(:,:,nf), Gh, xmax_h, xmean_h, xmax_d, xsum_d, nna, nnb)
             call mx(dpr(1), abs(xmax_d - xmax_h) / max(xmax_h, tiny(1.d0)))
             call mx(dpr(2), abs(xsum_d/dble(N*N) - xmean_h) / max(xmean_h, tiny(1.d0)))
             call dev_check_prec(Gd(:,:,nf), Gp, xmax_h, xmean_h, xmax_d, xsum_d, nna, nnb)
             call mx(dpr(3), abs(xmax_d - xmax_h) / max(xmax_h, tiny(1.d0)))
             call mx(dpr(4), abs(xsum_d/dble(N*N) - xmean_h) / max(xmean_h, tiny(1.d0)))
             if (nna /= 0 .or. nnb /= 0) call selftest_fail('stab: NaN count of the precision reduction')
          enddo
       enddo
    enddo
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test WRAPUR+decompose over ', nblk(2), &
         & ' blocks: max rel diff U ', mu(2), ', V ', mv(2), ', scales ', ms(2), ' (tol ', tol
    call say(line)
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test CGR at ', nblk(2), &
         & ' stabilization points x flavours: max rel diff G NVAR=1 ', mg(1), ', NVAR=2 ', mg(2), '; phase NVAR=1 ', mph(1), &
         & ', NVAR=2 ', mph(2), ' (tol ', tol
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3)') 'GPU stab self-test (reported) orthogonality max|U^dagger U - 1| of the decomposed U (all blocks): device ', orth_d, &
         & ', host ', orth_h
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU stab self-test (reported) intrinsic sensitivity of G to a 2e-16 relative (1 ulp) perturbation of U_R: host CGR NVAR=1 ', &
         & sens_h(1), ', NVAR=2 ', sens_h(2), '; device CGR NVAR=1 ', sens_d(1), ', NVAR=2 ', sens_d(2)
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test precision reduction vs COMPARE: rel diff max/mean ', &
         & dpr(1), ' / ', dpr(2), ' (G vs G_ref), ', max(dpr(3),dpr(4)), ' (G vs perturbed) (tol ', tolr
    call say(line)
    if (.not. (mu(1) <= tol .and. mv(1) <= tol .and. ms(1) <= tol .and. mu(2) <= tol .and. mv(2) <= tol .and. ms(2) <= tol)) &
         & call selftest_fail('stab: WRAPUL/WRAPUR + decompose')
    if (.not. (mg(1) <= tol .and. mg(2) <= tol .and. mph(1) <= tol .and. mph(2) <= tol)) call selftest_fail('stab: CGR')
    if (.not. (dpr(1) <= tolr .and. dpr(2) <= tolr .and. dpr(3) <= tolr .and. dpr(4) <= tolr)) &
         & call selftest_fail('stab: precision reduction')
    do nfe = 1, s_nfe
       call hl(nfe)%dealloc; call hr(nfe)%dealloc; call dl_(nfe)%dealloc; call dr_(nfe)%dealloc; call hp(nfe)%dealloc
    enddo
    deallocate(hl, hr, dl_, dr_, hp, Gh, Gp, Gd, Gd2)
    in_selftest = .false.
    call say('GPU stab self-test PASSED')
  end subroutine gpu_stab_selftest

  ! graded fault (ALF_GPU_FAULT_STAB=eps:<x>): the element of largest |G| of the device G (host copy) times (1 + x)
  Subroutine fault_eps_apply(G)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(:,:)
    Integer :: i, j, im, jm
    Real (Kind=Kind(0.d0)) :: x, a, amax
    x = gpu_env_fault_stab_eps()
    if (x <= 0.d0) return
    amax = -1.d0; im = 1; jm = 1
    do j = 1, size(G,2)
       do i = 1, size(G,1)
          a = abs(G(i,j))
          if (a > amax) then
             amax = a; im = i; jm = j
          endif
       enddo
    enddo
    G(im,jm) = G(im,jm) * (1.d0 + x)
  end subroutine fault_eps_apply

  ! hp = hr with every entry of U multiplied by (1 + 4e-16 * r), r in (-1/2, 1/2): +-2e-16, about one ulp
  Subroutine perturb_state(a, p, seed)
    Implicit none
    Class (UDV_State), intent(in) :: a
    Class (UDV_State), intent(inout) :: p
    Integer(int64), intent(inout) :: seed
    Integer :: i, j
    p%U = a%U; p%V = a%V
#if defined(STABLOG)
    p%L = a%L
#else
    p%D = a%D
#endif
    do j = 1, size(p%U,2)
       do i = 1, size(p%U,1)
          p%U(i,j) = p%U(i,j) * (1.d0 + 4.d-16*rnd(seed))
       enddo
    enddo
  end subroutine perturb_state

  Subroutine state_gdiff(A, B, d)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:), B(:,:)
    Real (Kind=Kind(0.d0)), intent(out) :: d
    d = rel_diff(A, B)
  end subroutine state_gdiff

  Subroutine rand_perturb(A, B, seed)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:)
    Complex (Kind=Kind(0.d0)), intent(out) :: B(:,:)
    Integer(int64), intent(inout) :: seed
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: re, im
    do j = 1, size(A,2)
       do i = 1, size(A,1)
          re = rnd(seed); im = rnd(seed)
          B(i,j) = A(i,j) + 1.d-3*cmplx(re, im, kind(0.d0))
       enddo
    enddo
  end subroutine rand_perturb

  ! device reduction of (A, B) next to the host COMPARE of the same matrices
  Subroutine dev_check_prec(A, B, xmax_h, xmean_h, xmax_d, xsum_d, nna, nnb)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:), B(:,:)
    Real (Kind=Kind(0.d0)), intent(out) :: xmax_h, xmean_h, xmax_d, xsum_d
    Integer, intent(out) :: nna, nnb
    call COMPARE(A, B, xmax_h, xmean_h)
    st_Gn = A; st_Gt = B
    !$omp target update to(st_Gn, st_Gt)       ! self-test: upload the two matrices of the reduction test
    call dev_prec(size(A,1), st_Gn, st_Gt, xmax_d, xsum_d, nna, nnb)
  end subroutine dev_check_prec

!====================================================================
!  Timing (ALF_GPU_STABTIME=1): time per stabilization step, device vs host, same states
!====================================================================
  Subroutine gpu_stab_timing()
    Implicit none
    Class (UDV_State), Allocatable :: hl(:), hr(:)
    Integer :: nfe, nst, nrep, k, N, nf, NT1, NT
    Integer (int64) :: c0, c1, crate
    Real (Kind=Kind(0.d0)) :: tg(5), tc(5), tgc(5), tcc(5), tgc0(5)
    Complex (Kind=Kind(0.d0)) :: ph
    Complex (Kind=Kind(0.d0)), Allocatable :: G(:,:,:)
    Character (len=512) :: line
    Character (len=8) :: mode

    N = s_N
    nrep = 5
    allocate(hl(s_nfe), hr(s_nfe), G(N,N,N_FL))
    G = z0
    do nfe = 1, s_nfe
       call hl(nfe)%init(N, 'l'); call hr(nfe)%init(N, 'r')
       call gpu_stab_reset_l(nfe); call gpu_stab_reset_r(nfe); call gpu_stab_reset_st(s_nstm, nfe)
    enddo
    ! build the left states for all blocks on the device (as the initial pass), then take a mid-chain right state
    do nst = s_nstm-1, 1, -1
       call gpu_stab_wrapul(s_nt(nst+1), s_nt(nst))
       do nfe = 1, s_nfe
          call gpu_stab_store_l(nst, nfe)
       enddo
    enddo
    nst = max(1, s_nstm/2)
    if (nst > 1) call gpu_stab_wrapur(s_nt(0), s_nt(nst-1))
    NT1 = s_nt(nst-1); NT = s_nt(nst)
#if defined(STABLOG)
    mode = 'LOG'
#else
    mode = 'default'
#endif
    call la_sync()
    do k = 1, nrep
       ! one stabilization step = WRAPUR block (Nwrap slices) + decompose, all flavours: device
       do nfe = 1, s_nfe
          call dl_state(1, 2, nfe, hr(nfe))
       enddo
       call gpu_stab_wrapur(NT1, NT); call la_sync()                      ! warm-up (the first call allocates library work space)
       do nfe = 1, s_nfe
          call ul_state(2, nfe, hr(nfe))                                  ! reset the input (timing: not included)
       enddo
       call system_clock(c0, crate)
       call gpu_stab_wrapur(NT1, NT); call la_sync()
       call system_clock(c1)
       tg(k) = dble(c1-c0)/dble(crate)
       ! host
       call system_clock(c0)
       call WRAPUR(NT1, NT, hr)
       call system_clock(c1)
       tc(k) = dble(c1-c0)/dble(crate)
       ! CGR: device (with the interim copies of the sweep), device only, host
       do nfe = 1, s_nfe
          nf = Calc_Fl_map(nfe)
          call gpu_stab_exchange(nst, nfe, .true., hr(nfe))
       enddo
       do nfe = 1, s_nfe
          call dl_state(1, 1, nfe, hl(nfe)); call dl_state(1, 2, nfe, hr(nfe))
       enddo
       do nfe = 1, s_nfe
          call gpu_stab_cgr(nfe, Calc_Fl_map(nfe), 1, G, ph, .false.)    ! G = the Green function of these states (precision check: ~0)
       enddo
       call system_clock(c0)
       do nfe = 1, s_nfe
          nf = Calc_Fl_map(nfe)
          call gpu_stab_cgr(nfe, nf, 1, G, ph, .true.)
       enddo
       call system_clock(c1)
       tgc(k) = dble(c1-c0)/dble(crate)
       call system_clock(c0)
       do nfe = 1, s_nfe
          call dev_cgr(nfe, 1)
       enddo
       call la_sync()
       call system_clock(c1)
       tgc0(k) = dble(c1-c0)/dble(crate)
       call system_clock(c0)
       do nfe = 1, s_nfe
          call CGR(ph, 1, G(:,:,Calc_Fl_map(nfe)), hr(nfe), hl(nfe))
       enddo
       call system_clock(c1)
       tcc(k) = dble(c1-c0)/dble(crate)
    enddo
    call sort5(tg); call sort5(tc); call sort5(tgc); call sort5(tgc0); call sort5(tcc)
    write(line,'(A,I0,A,I0,A,A,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU stab timing N = ', N, ', Nwrap block = ', &
         & NT-NT1, ', ', trim(mode), ', flavours ', s_nfe, ' (median of 5, s per step for all flavours): WRAPUR+decompose GPU ', tg(3), &
         & ' CPU ', tc(3), ' | CGR GPU incl. interim copies ', tgc(3), ' GPU only ', tgc0(3), ' CPU ', tcc(3)
    call say(line)
    do nfe = 1, s_nfe
       call hl(nfe)%dealloc; call hr(nfe)%dealloc
    enddo
    deallocate(hl, hr, G)
  end subroutine gpu_stab_timing

  Subroutine sort5(t)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: t(5)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: x
    do i = 2, 5
       x = t(i); j = i-1
       do while (j >= 1)
          if (t(j) <= x) exit
          t(j+1) = t(j); j = j-1
       enddo
       t(j+1) = x
    enddo
  end subroutine sort5

End Module gpu_stab_mod
#endif

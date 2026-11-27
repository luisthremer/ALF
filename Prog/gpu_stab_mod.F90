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
!> GPU port (P2.4 step 1d): stabilisation on the device-resident UDV states.  Compiled only with -DALF_GPU.
!>
!> @details
!> What runs on the device (finite temperature and projector, default and STABLOG stabilisation; STAB1/2/3 are refused by
!> gpu_supported; the projector part is at the end of the module: Ndim x N_part states, dev_decompose_p, dev_cgrp):
!>   * gpu_stab_wrapur / gpu_stab_wrapul: the slice products of WRAPUR / WRAPUL on the U matrices of the device states
!>     (gpu_wrap_mod), then the UDV decomposition (decompose_UDV_state) and the scale test (testscale);
!>   * the decomposition: default = column-norm pre-pivot + unpivoted QR (la_zqrpp; the CPU analogue is a build with
!>     -DQR_PREPIVOT), D separated by a kernel, phase and permutation sign, V = R * P^-1 V by one GEMM on the triangular
!>     factor, Q = reflectors applied to the unit matrix by ZUNMQR (not ZUNGQR: accuracy of the Precision Green).
!>     STABLOG = the CPU algorithm itself: log-norm selection sort (norms D->H, swap permutation H->D, columns gathered on the
!>     device), unpivoted QR (ZGEQRF), the exp(L_j - L_i) row scaling and the log update;
!>   * gpu_stab_cgr: CGR for NVAR = 1 and 2 (cgr1_mod, branches default and STABLOG).  Default: la_zqrpp.  STABLOG: ZGEQP3 on
!>     the host (la_zgeqp3, upstream numerics).  D separation, phase (permutation sign, R diagonal, reflector determinants),
!>     ZUNMQR, ZTRSM, ZLAPMR/T and the final GEMM on the device;
!>   * the determinant phase of U_R^dagger U_L, from one of two sources (ALF_GPU_CGR_LU, "det tracking" below).  LU: ZGETRF at
!>     every CGR, the kernel k_det_phase returns the normalised phase and the raw determinant (32 B D->H), the host picks between
!>     them with the rule of LU_PHASE_SAFE (gpu_lu_phase_select, tested bitwise against it); the raw product alone underflows for
!>     N > ~2950 (fix of fix/det-phase-underflow, KNOWN_ISSUES "Determinant phase underflow").  Tracked: det(U)/|det U| of every
!>     device slot U is kept with the slot (no LU, no transfer).  CGR_PHASE_DEFAULT names the source of a run without
!>     ALF_GPU_CGR_LU: 'tracked' (step 1f, speedup S1; one commit, "git revert" of it makes 'lu' the default again).
!>     Why no LU: the LU (ZGETRF of the N x N matrix, a serial determinant kernel and a 32 B D->H at every CGR) computes a number
!>     the decompositions already know, since the reflectors of the two QRs fix det U_L and det U_R exactly; measured on an H100
!>     (Hubbard, ham_h0 = 0): L = 32 8.57 -> 7.49 s per sweep (-12.6 %).  The tracked phase never forms a raw product and needs
!>     only tau of the QR, which every vendor backend has;
!>   * the precision check of Control_PrecisionG as a device reduction (max |A-B|, sum |A-B|, non-finite counts): one D->H of
!>     four numbers, handed to Control_PrecisionG_vals (control_mod).
!>
!> Device states (gpu_state_mod): udw_*(.,.,1,nfe) = udvl, udw_*(.,.,2,nfe) = udvr, udv_*(.,.,NST,nfe) = udvst.  The host
!> UDV_State objects of the walker are not updated by this path, except udvst(NST) in the backward sweep when Ltau = 1,
!> because the host TAU_M reads it.  G is device-authoritative (G_d): the precision check reads it device to device and the
!> host copy is made on demand by the sweep (gpu_g_host).
!>
!> The la_* calls (linalg_mod) are the places a vendor backend replaces (P2.4 step 2); with GPU=host they run on the host.
!> Transfers (N = Ndim): decompose default 24 B D->H (scale test) + 8 N + 4 N B inside la_zqrpp; STABLOG 8 N B D->H (log norms)
!> and 4 N B H->D (permutation); CGR 32 B D->H (LU phase; none with the tracked phase) + 16 B D->H (phase) + 32 B D->H (precision
!> numbers); STABLOG CGR
!> additionally 16 N^2 B each way for the host ZGEQP3; backward sweep with Ltau = 1: udvst(NST) D->H, 32 N^2 + 16 N B.
!>
!> Det tracking (speedup S1).  ZGEQRF writes U = H_1 ... H_n R, H_i = I - tau_i v_i v_i^H with v_i(i) = 1 and |v_i|^2 = 2 Re(tau_i) /
!> |tau_i|^2, so det H_i = -tau_i / conj(tau_i) (1 for tau_i = 0): modulus 1.  The new U of a decomposition is Q diag(Phase,1,..,1) with
!> Q = H_1 ... H_n and Phase = cs(1) of k_phase (permutation sign and phases of the R_ii included), so det U_new = prod det(H_i) * Phase:
!> known exactly at the cost of n complex products (k_phase, cs(6)).  det(U_R^dagger U_L) = conj(det U_R) det U_L, and the CGR phase is
!> conj(det)/|det| of it, the number the LU gives (k_det_phase_u, cs(7)).  Only phases are formed, never a raw product (N > ~2950).
!> Invariant: every writer of a device U sets its det.  decompose_slot rebuilds U and det from the same reflectors (every wrap U <- B U
!> is followed by decompose_slot in the same routine: keep it so); resets set the identity and det 1; copies and the exchange copy the
!> det with U; gpu_stab_put_udv is the only upload and takes the det from a QR of each uploaded U (k_det_qr).  The flags w_det_ok /
!> st_det_ok record it; a CGR in a mode that uses the tracked phase stops with "untracked U in CGR" instead of a silent LU fallback.
!> The device tau pass (gpu_tau_mod) only reads U.  A wrong det is silent in a model with a sign problem (the phase sign flips): hence
!> the check mode and the fault control below.
!> Check mode (debug runs, needs ALF_GPU_DEBUG=1): ALF_GPU_CGR_LU = 0 tracked, 1 LU, 2 both (the tracked phase is used and compared
!> with the LU phase at every CGR; the number of comparisons and the largest difference are written to the info file at the end, a
!> difference above 1e-8 stops the run).  The self-test runs in mode 2 whatever the default is.
!>
!> Environment (debug runs): ALF_GPU_FAULT_STAB=perm|phase|nan|det injects a fault (negative controls of the self-test; det: sign of
!> the tracked phase flipped, seen by mode 2).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_stab_mod

  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, N_SUN, Calc_Fl_map, Nsigma, Op_V, Projector, WF_L, WF_R
  Use Operator_mod, only: Op_phase
  Use QDRP_mod, only: Pivot_phase
  Use Control, only: Control_PrecisionG_vals, Control_PrecisionP, Control_Check_Phase, Control_Is_Finite
  Use MyMats, only: COMPARE
  Use wrapur_mod, only: WRAPUR
  Use wrapul_mod, only: WRAPUL
  Use cgr1_mod, only: CGR, CGRP
  Use linalg_mod, only: la_zgemm, la_ztrsm, la_zunmqr, la_zgeqrf, la_zgetrf, la_zqrpp
#if defined(ALF_LA_CUDA)
  Use linalg_cuda_mod, only: la_cuda_sync
#endif
#if defined(ALF_LA_HIP)
  Use alf_hipk_mod, only: hk_enabled, hk_ident, hk_setd, hk_setl, hk_copyz, hk_scale_cols, hk_diagabs, hk_divrows, hk_phase, &
       &                  hk_scale_row1, hk_scale_col1, hk_triu, hk_perm, hk_perm_vec, hk_ct, hk_rowdiv, hk_colinv, hk_tpup, &
       &                  hk_det_phase, hk_cgr_phase, hk_testscale, hk_prec
#endif
#if defined(STABLOG)
  Use linalg_mod, only: la_zgeqp3
#endif
  Use QMC_runtime_var, only: get_Ltau
  Use qmc_walker_mod, only: qmc_walker_t
  Use UDV_State_mod, only: UDV_State
  Use gpu_state_mod, only: gpu_put_udv, gpu_active, G_d, g_dev_ok, g_host_ok, udv_U, udv_V, udv_D, udv_L, udw_U, udw_V, udw_D, udw_L, &
       &                   info_line, gpu_env_get, gpu_env_is, selftest_fail, rnd, rel_diff, finite_c, is_finite
  Use gpu_update_mod, only: gpu_upd_enabled
  Use gpu_wrap_mod, only: gpu_wrapur_slices, gpu_wrapul_slices

  Implicit none

  Private
  Public :: gpu_stab_on, gpu_stab_setup, gpu_stab_finalize
  Public :: gpu_stab_wrapur, gpu_stab_wrapul, gpu_stab_exchange, gpu_stab_store_l, gpu_stab_cgr, gpu_stab_recompute
  Public :: gpu_stab_reset_l, gpu_stab_reset_r, gpu_stab_reset_st, gpu_stab_put_udv, gpu_stab_untrack, gpu_stab_export_r
  Public :: gpu_lu_phase_select, k_det_phase, gpu_stab_hook, gpu_stab_selftest
  Public :: k_det_qr, k_det_phase_u, k_set_one      ! det tracking (tests)
#if defined(ALF_LA_HIP)
  Public :: k_phase, k_cgr_phase, k_scale_cols, k_rowdiv, k_colinv, k_perm_rows, k_iperm_cols, k_triu, k_testscale  ! 56-hipk-stab
#if !defined(STABLOG)
  Public :: k_tpup
#endif
#endif
  Public :: gpu_stab_tau_dev, dev_prec, k_ident, k_copy_c, k_ct, k_diagabs, k_divrows, say     ! used by gpu_tau_mod

  Logical, Save :: gpu_stab_on = .false.      ! the device stabilisation replaces WRAPUR/WRAPUL/CGR in the sweeps
  Logical, Save :: gpu_stab_tau_dev = .false.  ! the device tau pass reads udvst on the device: no D->H export of udvst (Ltau = 1)

  Integer, Save :: s_N = 0, s_nfe = 0, s_nstm = 0, s_lw = 0
  Integer, Allocatable, Save :: s_nt(:)                     ! Stab_nt(0:NSTM)
  ! projector code: U of every slot is Ndim x N_part (first s_nc(nfe) columns of the Ndim x Ndim slot), no V, no tracked det
  Logical, Save :: s_proj = .false.
  Integer, Allocatable, Save :: s_nc(:)                     ! columns of U per flavour: Ndim, or N_part (projector)
  Integer, Save :: s_npm = 0                                ! max N_part (projector)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: p_wf(:,:,:,:)   ! (Ndim,s_npm,2,nfe): 1 WF_L%P, 2 WF_R%P (device)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: p_S(:,:), p_R(:,:)  ! CGRP: S (np x np, LU factors), R (np x Ndim) (device)
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

  ! ---- own device buffers (allocated once in gpu_stab_setup) ----
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: sA(:,:), sB(:,:), sC(:,:)   ! N x N work matrices
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_Gn(:,:), st_Gt(:,:)       ! new G, wrapped G
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_tau(:), st_Dl(:)          ! QR reflector scalars, local D of the QR
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_c(:)                      ! scalars: 1 phase, 2 1/phase, 3 CGR phase, 4 LU phase, 5 raw det,
                                                                                ! 6 det(U) of the new U (k_phase), 7 tracked CGR determinant phase
  ! det tracking: det(U)/|det U| of every device slot U (udw: 1 udvl, 2 udvr; udv: udvst(nst)); the flag says that the value is valid
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: udw_det(:,:), udv_det(:,:)    ! (2,nfe), (nstm,nfe), device
  Logical, Allocatable, Save :: w_det_ok(:,:), st_det_ok(:,:)                    ! host flags
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: st_work(:)                   ! LAPACK work space (queried once)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_r(:)                      ! scalars: 1-3 scale test, 4-7 precision numbers
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_nrm(:), st_Lt(:)          ! log norms (STABLOG); L work vector
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_pw(:)                     ! per-column partial sums of k_prec
  Integer, Allocatable, Save :: st_ipvt(:), st_perm(:), st_ipiv(:)             ! pivot (ZGEQP3 semantics), sort permutation, LU pivots
#if defined(STABLOG)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: st_rw(:)                     ! RWORK of ZGEQP3
#endif
  Logical, Save :: warn_once = .true.
  Integer, Save :: fault = 0                                                  ! 1 perm, 2 phase, 3 nan, 4 det (negative controls)
  ! source of the CGR determinant phase: 0 tracked det(U) (no LU), 1 LU of the RHS, 2 both (check: tracked used, compared with the LU).
  ! CGR_PHASE_DEFAULT is the mode without ALF_GPU_CGR_LU (which needs ALF_GPU_DEBUG=1)
  Character (len=*), Parameter :: CGR_PHASE_DEFAULT = 'tracked'               ! 'lu' or 'tracked'
  Integer, Save :: cgr_mode = 1
  Integer, Save :: ph_ncmp = 0                                                 ! mode 2: CGR calls compared
  Real (Kind=Kind(0.d0)), Save :: ph_maxdiff = 0.d0                            ! mode 2: max |phase(tracked) - phase(LU)|
  Real (Kind=Kind(0.d0)), Parameter :: PH_ABORT = 1.d-8                        ! mode 2: abort above this difference
  Logical, Save :: in_selftest = .false.                                      ! inside gpu_stab_selftest

Contains

!====================================================================
!  Setup
!====================================================================
!--------------------------------------------------------------------
!> @brief
!> Allocates the work buffers, queries the LAPACK work space, sets all device states to the identity.
!> @details Called once from Setup_gpu after gpu_state_setup.  The stabilisation is on iff the GPU path and the device
!> update path are active.  Needs Walker_layout (Stab_nt).
!> @param[in] walker walker whose layout (NSTM, Stab_nt) is used
!--------------------------------------------------------------------
  Subroutine gpu_stab_setup(walker)
    Implicit none
    Type (qmc_walker_t), intent(in) :: walker
    Integer :: n, nfe, info, lw
    Complex (Kind=Kind(0.d0)) :: wq(1)
    Character (len=256) :: line
    Character (len=64) :: v

    gpu_stab_on = .false.
    if (.not. gpu_active) return
    if (.not. gpu_upd_enabled()) return
    cgr_mode = merge(0, 1, CGR_PHASE_DEFAULT == 'tracked')
    if (gpu_env_get('ALF_GPU_CGR_LU', v)) then                  ! values 0|1|2 are checked by gpu_env_check
       select case (trim(v))
       case ('0'); cgr_mode = 0
       case ('1'); cgr_mode = 1
       case ('2'); cgr_mode = 2
       end select
    endif
    ph_ncmp = 0; ph_maxdiff = 0.d0
    s_N = Ndim; s_nfe = N_FL_eff; s_nstm = walker%NSTM
    allocate(s_nt(0:s_nstm)); s_nt = walker%Stab_nt
    n = s_N
    allocate(s_nc(s_nfe)); s_nc = n
    s_proj = Projector
    if (s_proj) call proj_setup()
    allocate(sA(n,n), sB(n,n), sC(n,n), st_Gn(n,n), st_Gt(n,n), st_tau(n), st_Dl(n), st_c(8), st_r(16), st_nrm(n), st_Lt(n), st_pw(n))
    allocate(st_ipvt(n), st_perm(n), st_ipiv(n))
    allocate(udw_det(2,s_nfe), udv_det(s_nstm,s_nfe), w_det_ok(2,s_nfe), st_det_ok(s_nstm,s_nfe))
    w_det_ok = .false.; st_det_ok = .false.
    st_ipvt = 1; st_perm = 1; st_ipiv = 1
    ! device allocation once, no transfer; the first use of each buffer is a device write
    !$omp target enter data map(alloc: sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_pw, st_ipvt, st_perm, st_ipiv, udw_det, udv_det)
    ! LAPACK work space: the largest of the queries of the calls below
    call la_zqrpp(n, n, sA, n, st_ipvt, st_tau, wq, -1, info)
    lw = int(real(wq(1)))
    call la_zunmqr('L', 'C', n, n, n, sA, n, st_tau, sB, n, wq, -1, info)
    lw = max(lw, int(real(wq(1))))
    call la_zunmqr('L', 'N', n, n, n, sA, n, st_tau, sB, n, wq, -1, info)
    lw = max(lw, int(real(wq(1))))
    call la_zunmqr('R', 'N', n, n, n, sA, n, st_tau, sB, n, wq, -1, info)
    lw = max(lw, int(real(wq(1))))
#if defined(STABLOG)
    allocate(st_rw(2*n))
    call la_zgeqp3(n, n, sA, n, st_ipvt, st_tau, wq, -1, st_rw, info)
    lw = max(lw, int(real(wq(1))))
#endif
    s_lw = max(lw, 1)
    allocate(st_work(s_lw))
    do nfe = 1, s_nfe
       call gpu_stab_reset_l(nfe)
       call gpu_stab_reset_r(nfe)
       call gpu_stab_reset_st(s_nstm, nfe)
    enddo
    gpu_stab_on = .true.
    write(line,'(A,I0,A,I0)') 'GPU stabilization (step 1d): ACTIVE (decompose, CGR and precision check on the device), N = ', n, &
         & ', Nstm = ', s_nstm
#if defined(STABLOG)
    line = trim(line)//', STABLOG: log-norm pivoting in decompose, host ZGEQP3 in CGR'
#else
    line = trim(line)//', default: column-norm pre-pivoted QR in decompose and CGR'
#endif
    select case (cgr_mode)
    case (0); line = trim(line)//'; CGR phase: tracked reflector determinants (no LU)'
    case (1); line = trim(line)//'; CGR phase: LU (ZGETRF)'
    case (2); line = trim(line)//'; CGR phase: tracked and LU, compared at every CGR (ALF_GPU_CGR_LU=2)'
    end select
    call info_line(line)
  end subroutine gpu_stab_setup

!--------------------------------------------------------------------
!> @brief Releases the buffers of gpu_stab_setup.
!--------------------------------------------------------------------
  Subroutine gpu_stab_finalize()
    Implicit none
    Character (len=256) :: line
    if (.not. gpu_stab_on) return
    if (cgr_mode == 2) then
       write(line,'(A,I0,A,ES10.3)') 'GPU CGR phase check (ALF_GPU_CGR_LU=2): CGR calls compared ', ph_ncmp, &
            & ', max |phase(tracked) - phase(LU)| ', ph_maxdiff
       call info_line(line)
    endif
    !$omp target exit data map(delete: sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_pw, st_ipvt, st_perm, st_ipiv, udw_det, udv_det)
    deallocate(sA, sB, sC, st_Gn, st_Gt, st_tau, st_Dl, st_c, st_r, st_nrm, st_Lt, st_pw, st_ipvt, st_perm, st_ipiv, st_work, s_nt)
    deallocate(udw_det, udv_det, w_det_ok, st_det_ok, s_nc)
#if defined(STABLOG)
    deallocate(st_rw)
#endif
    if (s_proj) then
       !$omp target exit data map(delete: p_wf, p_S, p_R)
       deallocate(p_wf, p_S, p_R)
    endif
    gpu_stab_on = .false.
  end subroutine gpu_stab_finalize

!====================================================================
!  Slot management (device to device)
!====================================================================
  Subroutine k_ident(n, U, ld)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_ident(n, ld, U)
       return
    endif
#endif
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(n)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_setd(n, D)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do i = 1, n
       D(i) = z1
    enddo
  end subroutine k_setD

  Subroutine k_setL(n, L)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(inout) :: L(n)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_setl(n, L)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do i = 1, n
       L(i) = 0.d0
    enddo
  end subroutine k_setL

  Subroutine k_copy_c(nel, a, b)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)), intent(in) :: a(nel)
    Complex (Kind=Kind(0.d0)), intent(inout) :: b(nel)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_copyz(nel, a, b)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do i = 1, nel
       b(i) = a(i)
    enddo
  end subroutine k_copy_c

  Subroutine k_copy_r(nel, a, b)
    Implicit none
    Integer, intent(in) :: nel
    Real (Kind=Kind(0.d0)), intent(in) :: a(nel)
    Real (Kind=Kind(0.d0)), intent(inout) :: b(nel)
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
    if (s_proj) then                                      ! projector: U = WF_L%P (slot 1) or WF_R%P (slot 2)
       call k_copy_c(s_N*s_nc(nfe), p_wf(1,1,slot,nfe), udw_U(1,1,slot,nfe)); return
    endif
    call k_ident(s_N, udw_U(1,1,slot,nfe), s_N)
    call k_ident(s_N, udw_V(1,1,slot,nfe), s_N)
    call k_set_one(1, udw_det(slot,nfe)); w_det_ok(slot,nfe) = .true.     ! det(1) = 1
#if defined(STABLOG)
    call k_setL(s_N, udw_L(1,slot,nfe))
#else
    call k_setD(s_N, udw_D(1,slot,nfe))
#endif
  end subroutine reset_w

!> @brief udvl of flavour nfe = identity state (device)
  Subroutine gpu_stab_reset_l(nfe)
    Implicit none
    Integer, intent(in) :: nfe
    call reset_w(1, nfe)
  end subroutine gpu_stab_reset_l

!> @brief udvr of flavour nfe = identity state (device)
  Subroutine gpu_stab_reset_r(nfe)
    Implicit none
    Integer, intent(in) :: nfe
    call reset_w(2, nfe)
  end subroutine gpu_stab_reset_r

!> @brief udvst(nst, nfe) = identity state (device)
  Subroutine gpu_stab_reset_st(nst, nfe)
    Implicit none
    Integer, intent(in) :: nst, nfe
    if (s_proj) then                                      ! projector: left state at Ltrot = WF_L%P
       call k_copy_c(s_N*s_nc(nfe), p_wf(1,1,1,nfe), udv_U(1,1,nst,nfe)); return
    endif
    call k_ident(s_N, udv_U(1,1,nst,nfe), s_N)
    call k_ident(s_N, udv_V(1,1,nst,nfe), s_N)
    call k_set_one(1, udv_det(nst,nfe)); st_det_ok(nst,nfe) = .true.      ! det(1) = 1
#if defined(STABLOG)
    call k_setL(s_N, udv_L(1,nst,nfe))
#else
    call k_setD(s_N, udv_D(1,nst,nfe))
#endif
  end subroutine gpu_stab_reset_st

  !> dst = src slot copies: wslot (1/2) -> udvst(nst)
  Subroutine copy_w_to_st(wslot, nfe, nst)
    Implicit none
    Integer, intent(in) :: wslot, nfe, nst
    call k_copy_c(s_N*s_N, udw_U(1,1,wslot,nfe), udv_U(1,1,nst,nfe))
    call k_copy_c(s_N*s_N, udw_V(1,1,wslot,nfe), udv_V(1,1,nst,nfe))
    call k_copy_c(1, udw_det(wslot,nfe), udv_det(nst,nfe)); st_det_ok(nst,nfe) = w_det_ok(wslot,nfe)
#if defined(STABLOG)
    call k_copy_r(s_N, udw_L(1,wslot,nfe), udv_L(1,nst,nfe))
#else
    call k_copy_c(s_N, udw_D(1,wslot,nfe), udv_D(1,nst,nfe))
#endif
  end subroutine copy_w_to_st

  !> udvst(nst) -> wslot (1/2)
  Subroutine copy_st_to_w(nst, nfe, wslot)
    Implicit none
    Integer, intent(in) :: wslot, nfe, nst
    call k_copy_c(s_N*s_N, udv_U(1,1,nst,nfe), udw_U(1,1,wslot,nfe))
    call k_copy_c(s_N*s_N, udv_V(1,1,nst,nfe), udw_V(1,1,wslot,nfe))
    call k_copy_c(1, udv_det(nst,nfe), udw_det(wslot,nfe)); w_det_ok(wslot,nfe) = st_det_ok(nst,nfe)
#if defined(STABLOG)
    call k_copy_r(s_N, udv_L(1,nst,nfe), udw_L(1,wslot,nfe))
#else
    call k_copy_c(s_N, udv_D(1,nst,nfe), udw_D(1,wslot,nfe))
#endif
  end subroutine copy_st_to_w

!> @brief udvst(nst, nfe) = udvl(nfe) (initial pass, device to device)
  Subroutine gpu_stab_store_l(nst, nfe)
    Implicit none
    Integer, intent(in) :: nst, nfe
    call copy_w_to_st(1, nfe, nst)
  end subroutine gpu_stab_store_l

!--------------------------------------------------------------------
!> @brief
!> Exchange of the current and the stored propagation at a stabilisation point (device to device).
!> @details
!> fwd = .true.  (forward sweep):  udvl = udvst(nst); udvst(nst) = udvr.
!> fwd = .false. (backward sweep): udvr = udvst(nst); udvst(nst) = udvl; with Ltau = 1 the new udvst(nst) is also copied
!> D->H into the host object ust, because the host TAU_M reads it.
!> @param[in] nst, nfe stabilisation index, flavour index
!> @param[in] fwd direction
!> @param[inout] ust host object walker%udvst(nst,nfe) (written only in the backward direction with Ltau = 1)
!--------------------------------------------------------------------
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
       if (get_Ltau() == 1 .and. s_proj) then
          ! projector: Tau_p (host) reads U of udvst (N_part columns; no V, D unused by CGRP)
          !$omp target update from(udv_U(:,1:s_nc(nfe),nst,nfe))
          ust%U = udv_U(:,1:s_nc(nfe),nst,nfe)
          ust%side = 'l'
       else if (get_Ltau() == 1 .and. .not. gpu_stab_tau_dev) then
          ! D->H of the stored left propagation: the host TAU_M reads walker%udvst (unavoidable until the device tau pass, step 1e)
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
!> @brief WRAPUR(ntau, ntau1, udvr) for all flavours: slices ntau+1..ntau1 on U of udvr, decompose, scale test.
  Subroutine gpu_stab_wrapur(ntau, ntau1)
    Implicit none
    Integer, intent(in) :: ntau, ntau1
    Integer :: nfe, nf
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe)
       call gpu_wrapur_slices(udw_U(1,1,2,nfe), s_N, s_nc(nfe), nf, ntau, ntau1)
       call decompose_slot('R', 2, nfe)
    enddo
  end subroutine gpu_stab_wrapur

!> @brief WRAPUL(ntau1, ntau, udvl) for all flavours: slices ntau1..ntau+1 (descending) on U of udvl, decompose, scale test.
  Subroutine gpu_stab_wrapul(ntau1, ntau)
    Implicit none
    Integer, intent(in) :: ntau1, ntau
    Integer :: nfe, nf
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe)
       call gpu_wrapul_slices(udw_U(1,1,1,nfe), s_N, s_nc(nfe), nf, ntau, ntau1)
       call decompose_slot('L', 1, nfe)
    enddo
  end subroutine gpu_stab_wrapul

  Subroutine decompose_slot(side, slot, nfe)
    Implicit none
    Character (len=1), intent(in) :: side
    Integer, intent(in) :: slot, nfe
    if (s_proj) then
       call dev_decompose_p(side, s_nc(nfe), udw_U(1,1,slot,nfe)); return
    endif
#if defined(STABLOG)
    call dev_decompose(side, udw_U(1,1,slot,nfe), udw_V(1,1,slot,nfe), udw_L(1,slot,nfe))
#else
    call dev_decompose(side, udw_U(1,1,slot,nfe), udw_V(1,1,slot,nfe), udw_D(1,slot,nfe))
    call dev_testscale(udw_D(1,slot,nfe))
#endif
    ! det(U) of the new U (k_phase), kept with the slot: U is rebuilt from the reflectors here, so its det is known exactly
    call k_copy_c(1, st_c(6), udw_det(slot,nfe)); w_det_ok(slot,nfe) = .true.
  end subroutine decompose_slot

!--------------------------------------------------------------------
!> decompose_UDV_state on the device.  U, V: N x N device arrays; S: D (default, complex) or L (STABLOG, real).
!> On input U holds B * U_old (propagated), D/L and V the previous factors; side R: A = U D V, side L: A = U D V^dagger.
!--------------------------------------------------------------------
#if defined(STABLOG)
  Subroutine dev_decompose(side, U, V, S)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: S(s_N)
#else
  Subroutine dev_decompose(side, U, V, S)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: S(s_N)
#endif
    Character (len=1), intent(in) :: side
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(s_N,s_N), V(s_N,s_N)
    Integer :: n, nswap, sgn, info
    Logical :: rside
    Complex (Kind=Kind(0.d0)) :: ph1
#if defined(STABLOG)
    Integer :: kk
#endif

    n = s_N
    rside = (side == 'R' .or. side == 'r')
#if !defined(STABLOG)
    ! U(:,i) = U(:,i) * D(i)   (V is allocated for finite T)
    call k_scale_cols(n, U, n, S)
    ! column-norm pre-pivot + unpivoted QR: U holds the reflectors and R, st_ipvt the permutation (ZGEQP3 semantics)
    call la_zqrpp(n, n, U, n, st_ipvt, st_tau, st_work, s_lw, info)
    if (info /= 0) call lapack_fail('decompose: la_zqrpp (-101: NaN or Inf in the matrix to decompose)', info)
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
    ! D->H: n log norms for the host selection sort (STABLOG decompose: the sort is sequential)
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
    call la_zgeqrf(n, n, U, n, st_tau, st_work, s_lw, info)   ! unpivoted QR (CPU: ZGEQP3 with all columns fixed)
    if (info /= 0) call lapack_fail('decompose: zgeqrf', info)
    call k_diagabs(n, U, n, st_Dl)
    call k_divrows(n, U, n, st_Dl)
#endif
    call k_phase(n, U, n, sgn, .not. rside, st_c, st_tau)        ! Phase = prod U(i,i), permutation sign, conj for side L; beta = 1/Phase
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
       call la_zgemm('N', 'N', n, n, n, z1, sA, n, sB, n, z0, V, n)
    else
       call la_zgemm('N', 'C', n, n, n, z1, sB, n, sA, n, z0, V, n)
    endif
    ! U = Q = H_1 ... H_n applied to the unit matrix by ZUNMQR (sC is free here).  Not ZUNGQR: the Precision Green of the sweep
    ! is better with ZUNMQR on the unit matrix (measured on the GPU line, ticket 206).
    call k_ident(n, sC, n)
    call la_zunmqr('L', 'N', n, n, n, U, n, st_tau, sC, n, st_work, s_lw, info)
    if (info /= 0) call lapack_fail('decompose: zunmqr', info)
    call k_copy_c(n*n, sC, U)
    call k_scale_col1(n, U, n, st_c)                      ! scale the first column of U with Phase
#if defined(STABLOG)
    if (fault == 3) then                                  ! negative control: NaN in the new L(1)
       st_Lt(1) = transfer(9221120237041090560_int64, 0.d0)
       !$omp target update to(st_Lt(1:1))
       call k_copy_r(1, st_Lt, S)
    endif
#endif
  end subroutine dev_decompose

  !> Stops the run after a failed LAPACK call of the stabilisation
  Subroutine lapack_fail(what, info)
    Implicit none
    Character (len=*), intent(in) :: what
    Integer, intent(in) :: info
    write(error_unit,*)
    write(error_unit,*) 'Error in the GPU stabilisation: ', what, ', info = ', info
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine lapack_fail

  !> Fault injection (negative control): a wrong pivot permutation (ALF_GPU_FAULT_STAB=perm)
  Subroutine fault_perm()
    Implicit none
    Integer :: t
    if (fault /= 1) return
    t = st_ipvt(1); st_ipvt(1) = st_ipvt(2); st_ipvt(2) = t
    !$omp target update to(st_ipvt(1:s_N))
  end subroutine fault_perm

#if defined(STABLOG)
  !> Stable selection sort of the STABLOG CPU decompose: perm = composite of the column swaps, nswap = number of swaps.
  Subroutine sort_stablog(n, tn, perm, nswap)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(inout) :: tn(n)
    Integer, intent(out) :: perm(n)
    Integer, intent(out) :: nswap
    Integer :: i, j, pvt, t
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
#endif

  !> testscale_UDV_state for the default stabilisation: NaN or Inf scale -> error, extreme scales -> warning (once).
  !> One D->H of three numbers.
  Subroutine dev_testscale(D)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: D(s_N)
    call k_testscale(s_N, D, st_r)
    ! D->H: NaN/Inf count and the extreme scales (testscale of UDV_State, 24 B)
    !$omp target update from(st_r(1:3))
    if (st_r(1) > 0.5d0) then
       if (in_selftest) call selftest_fail('stab: NaN or Inf scale in decompose')
       write(error_unit,*)
       write(error_unit,*) "Error: At least one scale is NaN or Inf."
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
!--------------------------------------------------------------------
!> @brief
!> CGR(Phase, NVAR, GRUP, udvr, udvl) of cgr1_mod for the device states of flavour nfe (udvl = slot 1, udvr = slot 2).
!> @details
!> lazy = .true. (the sweep): G is device-authoritative.  The wrapped G for the precision check is G_d(:,:,nf) (device to
!> device; uploaded from G once if the device copy is stale), the new G goes to G_d(:,:,nf) and the host G is NOT written
!> (g_host_ok = .false.; host readers call gpu_g_host).  Otherwise (lazy absent or .false., self-test) G(:,:,nf) on the
!> host receives the new Green function and G_d is not touched.
!> with_prec: the precision check of Control_PrecisionG against the incoming G (the wrapped one), as a device reduction.
!> The returned phase is checked by Control_Check_Phase.
!> @param[in] nfe, nf, nvar flavour index (stored), flavour, CGR variant
!> @param[inout] G Green function (host)
!> @param[out] phase phase of the determinant, as CGR
!> @param[in] with_prec precision check on/off
!> @param[in] lazy device-authoritative G
!--------------------------------------------------------------------
  Subroutine gpu_stab_cgr(nfe, nf, nvar, G, phase, with_prec, lazy)
    Implicit none
    Integer, intent(in) :: nfe, nf, nvar
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(Ndim,Ndim,N_FL)
    Complex (Kind=Kind(0.d0)), intent(out) :: phase
    Logical, intent(in) :: with_prec
    Logical, intent(in), optional :: lazy
    Integer :: n, i, j
    Real (Kind=Kind(0.d0)) :: xmax, xsum
    Integer :: nna, nnb
    Logical :: sess

    n = s_N
    sess = .false.
    if (present(lazy)) sess = lazy
    if (with_prec) then
       if (sess) then
          if (.not. g_dev_ok(nf)) then
             ! the device copy is stale (host write since): H->D of the wrapped G once (16 N^2 B, unavoidable: host GR is the newer copy)
             G_d(:,:,nf) = G(:,:,nf)
             !$omp target update to(G_d(:,:,nf))
             g_dev_ok(nf) = .true.
          endif
          !$omp target teams distribute parallel do collapse(2)
          do j = 1, n
             do i = 1, n
                st_Gt(i,j) = G_d(i,j,nf)
             enddo
          enddo
       else
          st_Gt = G(:,:,nf)
          ! H->D of the wrapped G (host-authoritative), 16 N^2 B, for the precision reduction (self-test path)
          !$omp target update to(st_Gt)
       endif
    endif
    if (s_proj) then
       call dev_cgrp(nfe)
    else
       call dev_cgr(nfe, nvar)
    endif
#if defined(ALF_LA_CUDA)
    call la_cuda_sync()                            ! cuSOLVER info of this stabilisation, before its results are used
#endif
    ! D->H: the phase (16 B)
    !$omp target update from(st_c(3:3))
    phase = st_c(3)
    call Control_Check_Phase(phase)
    if (sess) then
       ! the new G becomes the device Green function (device to device)
       !$omp target teams distribute parallel do collapse(2)
       do j = 1, n
          do i = 1, n
             G_d(i,j,nf) = st_Gn(i,j)
          enddo
       enddo
       g_dev_ok(nf) = .true.
       g_host_ok(nf) = .false.
    else
       ! D->H of the new G, 16 N^2 B (self-test path: the host receives G)
       !$omp target update from(st_Gn)
       G(:,:,nf) = st_Gn
    endif
    if (with_prec) then
       call dev_prec(n, st_Gn, st_Gt, xmax, xsum, nna, nnb)
       call Control_PrecisionG_vals(xmax, xsum, nna, nnb, n)
    endif
  end subroutine gpu_stab_cgr

!--------------------------------------------------------------------
!> @brief
!> Stabilisation point of the sweep on the device: Recompute_GR_phase of qmc_walker_mod for the device states.
!> @details
!> Per computed flavour the device CGR with the precision check, then Op_phase, the product of the flavour phases,
!> Control_PrecisionP and walker%Phase on the host.  walker%GR is not written (G stays on the device).
!> @param[inout] walker the walker (Phase, Phase_array, GR as host copy)
!> @param[in] NVAR CGR variant (1 or 2)
!--------------------------------------------------------------------
  Subroutine gpu_stab_recompute(walker, NVAR)
    Implicit none
    Type (qmc_walker_t), intent(inout) :: walker
    Integer, intent(in) :: NVAR
    Integer :: nf, nf_eff
    Complex (Kind=Kind(0.d0)) :: Z, Z1

    walker%Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       call gpu_stab_cgr(nf_eff, nf, NVAR, walker%GR, Z1, .true., .true.)
       call Op_phase(Z1, OP_V, Nsigma, nf)
       walker%Phase_array(nf) = Z1
    enddo
    ! (gpu_supported refuses models with reconstruction_needed, so ham%weight_reconstruction is not needed here)
    Z = product(walker%Phase_array)
    Z = Z**N_SUN
    Call Control_PrecisionP(Z, walker%Phase)
    walker%Phase = Z
  end subroutine gpu_stab_recompute

  !> Control_PrecisionG as a reduction: xmax = max |A-B|, xsum = sum |A-B| (xmean = xsum/N^2), NaN or Inf entries of A and B.
  !> Such entries do not enter xmax/xsum.  One D->H of four numbers.
  Subroutine dev_prec(n, A, B, xmax, xsum, nna, nnb)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: A(n,n), B(n,n)
    Real (Kind=Kind(0.d0)), intent(out) :: xmax, xsum
    Integer, intent(out) :: nna, nnb
    call k_prec(n, A, B, st_r)
    ! D->H: the four reduced numbers of the precision check (32 B)
    !$omp target update from(st_r(4:7))
    xmax = st_r(4); xsum = st_r(5); nna = nint(st_r(6)); nnb = nint(st_r(7))
  end subroutine dev_prec

!--------------------------------------------------------------------
  Subroutine dev_cgr(nfe, nvar)
    Implicit none
    Integer, intent(in) :: nfe, nvar
    Integer :: n, info
    Complex (Kind=Kind(0.d0)) :: ph0
    Logical :: use_lu, use_tr

    n = s_N
    ! RHS = U_R^dagger U_L  (sA),  T = V_R V_L  (sB)
    call la_zgemm('C', 'N', n, n, n, z1, udw_U(1,1,2,nfe), n, udw_U(1,1,1,nfe), n, z0, sA, n)
    call la_zgemm('N', 'N', n, n, n, z1, udw_V(1,1,2,nfe), n, udw_V(1,1,1,nfe), n, z0, sB, n)
#if defined(STABLOG)
    call k_tpup_log(n, sB, sA, udw_L(1,2,nfe), udw_L(1,1,nfe))
#else
    call k_tpup(n, sB, sA, udw_D(1,2,nfe), udw_D(1,1,nfe))
#endif
    ! NOTE: upstream's one-time "Large number encountered" warning of CGR is omitted here (the scale test of decompose warns instead).
    ! det(U_R^dagger U_L): from the tracked dets of U_L and U_R (cs(7), on the device, no transfer) and / or from the LU of the RHS.
    use_tr = (cgr_mode /= 1)
    use_lu = (cgr_mode /= 0)
    ph0 = z1
    if (use_tr) then
       if (.not. (w_det_ok(1,nfe) .and. w_det_ok(2,nfe))) call untracked_fail(nfe)
       call k_det_phase_u(udw_det(1,nfe), st_c, merge(1, 0, fault == 4))
    endif
    if (use_lu) then
       ! LU of the RHS (sA is not read after k_tpup, cgr_solve overwrites it).  The kernel returns the normalised phase and the raw
       ! determinant, the host selects (LU_PHASE_SAFE rule) and forms Phase = conj(det)/|det|.
#if defined(ALF_LA_HIP)
       ! GPU=amd: this LU (debug modes 1, 2 only) runs on host LAPACK: the device la_zgetrf here was followed by layout-dependent
       ! wrong results (KNOWN_ISSUES R2).  D->H of the RHS, host LU on an unmapped copy, H->D of the factors and pivots.
       call host_lu(n, info)
#else
       call la_zgetrf(n, n, sA, n, st_ipiv, info)
#endif
       if (info < 0) call lapack_fail('CGR: zgetrf', info)
       call k_det_phase(n, sA, n, st_ipiv, st_c)
       ! D->H: normalised LU phase and raw determinant (32 B)
       !$omp target update from(st_c(4:5))
       ph0 = conjg(gpu_lu_phase_select(st_c(5), st_c(4)))
       ph0 = ph0/abs(ph0)
    endif
    if (cgr_mode == 2) then
       ! check mode (debug): D->H of the tracked phase (16 B), compared with the LU phase on the host
       !$omp target update from(st_c(7:7))
       call check_phase(ph0, st_c(7))
    endif
    if (nvar == 1) then
       call cgr_solve(nfe, nvar, sB, ph0, use_tr)
    else
       call k_ct(n, sB, n, sC, n)                          ! T = T^dagger
       call cgr_solve(nfe, nvar, sC, ph0, use_tr)
    endif
  end subroutine dev_cgr

  !> Check mode (ALF_GPU_CGR_LU=2): |tracked phase - LU phase| counted and maximised; above PH_ABORT (or not finite) the run stops.
  Subroutine check_phase(p_lu, p_tr)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: p_lu, p_tr
    Real (Kind=Kind(0.d0)) :: dif
    Character (len=200) :: msg
    dif = abs(p_lu - p_tr)
    ph_ncmp = ph_ncmp + 1
    if (Control_Is_Finite(dif)) then
       if (dif > ph_maxdiff) ph_maxdiff = dif
       if (dif <= PH_ABORT) return
    endif
    write(msg,'(A,ES10.3,A,I0)') 'CGR phase check (ALF_GPU_CGR_LU=2): tracked det phase differs from the LU phase by ', dif, &
         & ' (limit 1e-8) at comparison ', ph_ncmp
    if (in_selftest) call selftest_fail('stab: '//trim(msg))
    write(error_unit,*)
    write(error_unit,*) 'Error in the GPU stabilisation: ', trim(msg)
    call info_line('GPU error: '//trim(msg))
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine check_phase

  !> Fail-fast: CGR needs det(U) of a slot whose det is not tracked (never a silent LU fallback)
  Subroutine untracked_fail(nfe)
    Implicit none
    Integer, intent(in) :: nfe
    write(error_unit,*)
    write(error_unit,*) 'Error in the GPU stabilisation: untracked U in CGR (flavour ', nfe, '): det(U) of udvl or udvr is not known.'
    write(error_unit,*) '       Every writer of a device U must set its det (decompose, reset, copy, gpu_stab_put_udv).'
    call info_line('GPU error: untracked U in CGR')
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine untracked_fail

!> @brief Test hook: marks slot (1 udvl, 2 udvr) of flavour nfe as untracked (the next tracked CGR must stop the run).
  Subroutine gpu_stab_untrack(slot, nfe)
    Implicit none
    Integer, intent(in) :: slot, nfe
    w_det_ok(slot,nfe) = .false.
  end subroutine gpu_stab_untrack

#if defined(ALF_LA_HIP)
  !> GPU=amd, debug modes ALF_GPU_CGR_LU = 1, 2: LU of the RHS sA on host LAPACK (KNOWN_ISSUES R2); sA and st_ipiv end as after
  !> la_zgetrf on the device.  The copies are not mapped, so linalg_mod's presence dispatch takes the host path.
  Subroutine host_lu(n, info)
    Implicit none
    Integer, intent(in) :: n
    Integer, intent(out) :: info
    Complex (Kind=Kind(0.d0)), Allocatable :: hA(:,:)
    Integer, Allocatable :: hp(:)
    allocate(hA(n,n), hp(n))
    !$omp target update from(sA)
    hA = sA
    call la_zgetrf(n, n, hA, n, hp, info)
    sA = hA
    st_ipiv(1:n) = hp
    !$omp target update to(sA, st_ipiv)
    deallocate(hA, hp)
  end subroutine host_lu
#endif

  !> The part of CGR after the matrix T = TPUP is built: QR, phase, solve.  T is destroyed.  Result in st_Gn, phase st_c(3).
  !> ph0: the LU phase (host), used if devph is .false.; devph: the tracked phase st_c(7) on the device is used instead.
  Subroutine cgr_solve(nfe, nvar, T, ph0, devph)
    Implicit none
    Integer, intent(in) :: nfe, nvar
    Complex (Kind=Kind(0.d0)), intent(inout) :: T(s_N,s_N)
    Complex (Kind=Kind(0.d0)), intent(in) :: ph0
    Logical, intent(in) :: devph
    Integer :: n, sgn, info
    Complex (Kind=Kind(0.d0)) :: ph1

    n = s_N
#if defined(STABLOG)
    ! STABLOG: host ZGEQP3 (upstream numerics): D->H of T, H->D of the factors, tau, ipvt (unavoidable: no device QRCP, 16 N^2 B each way)
    !$omp target update from(T)
    st_ipvt = 0
    call la_zgeqp3(n, n, T, n, st_ipvt, st_tau, st_work, s_lw, st_rw, info)
    if (info /= 0) call lapack_fail('CGR: zgeqp3', info)
    !$omp target update to(T, st_tau, st_ipvt)
    call fault_perm()
#else
    call la_zqrpp(n, n, T, n, st_ipvt, st_tau, st_work, s_lw, info)
    if (info /= 0) call lapack_fail('CGR: la_zqrpp (-101: NaN or Inf in the matrix)', info)
    call fault_perm()
#endif
    call k_diagabs(n, T, n, st_Dl)
    call k_divrows(n, T, n, st_Dl)
    ph1 = z1
    call Pivot_phase(ph1, st_ipvt, n)
    sgn = nint(real(ph1))
    call k_cgr_phase(n, T, n, st_tau, st_c, nvar, sgn, merge(1, 0, fault == 2), ph0, merge(1, 0, devph))
    if (nvar == 1) then
       call k_ct(n, udw_U(1,1,2,nfe), n, sA, n)             ! RHS = U_R^dagger
#if defined(STABLOG)
       call k_rowscale_log(n, sA, n, udw_L(1,2,nfe))
#endif
       call la_zunmqr('L', 'C', n, n, n, T, n, st_tau, sA, n, st_work, s_lw, info)
       call k_rowdiv(n, sA, n, st_Dl)
       call la_ztrsm('L', 'U', 'N', 'N', n, n, z1, T, n, sA, n)
       call k_iperm_rows(n, sA, n, st_ipvt, sC, n)           ! ZLAPMR backward  (sC is free: T is sB)
#if defined(STABLOG)
       call k_rowscale_log(n, sC, n, udw_L(1,1,nfe))
#endif
       call la_zgemm('N', 'N', n, n, n, z1, udw_U(1,1,1,nfe), n, sC, n, z0, st_Gn, n)
    else
       call k_copy_c(n*n, udw_U(1,1,1,nfe), sA)              ! RHS = U_L
#if defined(STABLOG)
       call k_colscale_log(n, sA, n, udw_L(1,1,nfe))
#endif
       call la_zunmqr('R', 'N', n, n, n, T, n, st_tau, sA, n, st_work, s_lw, info)
       call k_colinv(n, sA, n, st_Dl)
       call la_ztrsm('R', 'U', 'C', 'N', n, n, z1, T, n, sA, n)
       call k_iperm_cols(n, sA, n, st_ipvt, sB, n)           ! ZLAPMT backward  (sB is free: T is sC)
#if defined(STABLOG)
       call k_colscale_log(n, sB, n, udw_L(1,2,nfe))
#endif
       call la_zgemm('N', 'C', n, n, n, z1, sB, n, udw_U(1,1,2,nfe), n, z0, st_Gn, n)
    endif
    if (info /= 0) call lapack_fail('CGR: zunmqr', info)
  end subroutine cgr_solve

!====================================================================
!  LU phase (the rule of LU_PHASE_SAFE, mymats_mod)
!====================================================================
!--------------------------------------------------------------------
!> @brief
!> Kernel: from the LU factors A, ipiv of ZGETRF the raw determinant cs(5) = prod_i A(i,i) (times the pivot sign, as DET_C)
!> and the normalised phase cs(4) = prod_i A(i,i)/|A(i,i)| (times the pivot sign, modulus 1).
!> @details A zero pivot sets cs(4) = 0.  The raw product underflows for large N while |det| = 1 (fix/det-phase-underflow);
!> the normalised product does not.  Serial loop on the device, no transfer.
!> @param[in] n order; A(lda,n) LU factors; ipiv(n) pivots (LAPACK convention)
!> @param[inout] cs scalar array, entries 4 and 5 are written
!--------------------------------------------------------------------
  Subroutine k_det_phase(n, A, lda, ipiv, cs)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(in) :: A(lda,n)
    Integer, intent(in) :: ipiv(n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(8)
    Complex (Kind=Kind(0.d0)) :: raw, p
    Real (Kind=Kind(0.d0)) :: amod
    Integer :: i, sg
    Logical :: zero
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_det_phase(n, lda, A, ipiv, cs)
       return
    endif
#endif
    !$omp target
    raw = cmplx(1.d0, 0.d0, kind(0.d0))
    p = cmplx(1.d0, 0.d0, kind(0.d0))
    sg = 1
    zero = .false.
    do i = 1, n
       raw = raw*A(i,i)
       if (.not. zero) then
          amod = abs(A(i,i))
          if (amod == 0.d0) then
             zero = .true.
          else
             p = p*(A(i,i)/amod)
             if (ipiv(i) /= i) p = -p
          endif
       endif
       if (ipiv(i) /= i) sg = -sg
    enddo
    if (sg == -1) raw = -raw
    if (zero) p = cmplx(0.d0, 0.d0, kind(0.d0))
    cs(4) = p
    cs(5) = raw
    !$omp end target
  end subroutine k_det_phase

!--------------------------------------------------------------------
!> @brief
!> Host selection of the determinant phase: z_raw (raw product) if it is finite, of modulus in [1e-250, 1e250] and within
!> 1e-8 of the normalised product p; else p.  A zero pivot (p = 0) returns z_raw.
!> @details Same rule and constants as LU_PHASE_SAFE (mymats_mod); the test 52-gpu-stab compares the two bitwise.
!> @param[in] z_raw raw determinant (k_det_phase cs(5))
!> @param[in] p normalised phase (k_det_phase cs(4))
!--------------------------------------------------------------------
  Complex (Kind=Kind(0.d0)) Function gpu_lu_phase_select(z_raw, p)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: z_raw, p
    Real (Kind=Kind(0.d0)) :: x
    Real (Kind=Kind(0.d0)), parameter :: X_LOW = 1.d-250, X_HIGH = 1.d250, PHASE_TOL = 1.d-8
    gpu_lu_phase_select = z_raw
    if (abs(p) == 0.d0) return
    if (Control_Is_Finite(dble(z_raw)) .and. Control_Is_Finite(aimag(z_raw))) then
       x = abs(z_raw)
       if (x > X_LOW .and. x < X_HIGH) then
          if (abs(z_raw/x - p) < PHASE_TOL) return
       endif
    endif
    gpu_lu_phase_select = p
  end function gpu_lu_phase_select

!====================================================================
!  Determinant tracking (speedup S1)
!====================================================================
  !> det of the elementary reflector H = I - tau v v^H of ZGEQRF (v(1) = 1, |v|^2 = 2 Re(tau)/|tau|^2): -tau/conj(tau), modulus 1;
  !> 1 for tau = 0.  |tau| is 0 or in [1,2], so the division by |tau| cannot underflow.
  Complex (Kind=Kind(0.d0)) function refl_det(z)
    Implicit none
    !$omp declare target
    Complex (Kind=Kind(0.d0)), intent(in) :: z
    Complex (Kind=Kind(0.d0)) :: u
    if (z == cmplx(0.d0, 0.d0, kind(0.d0))) then
       refl_det = cmplx(1.d0, 0.d0, kind(0.d0))
    else
       u = z/abs(z)
       refl_det = -(u*u)
    endif
  end function refl_det

  !> d(1:nel) = 1 (det of a unit matrix)
  Subroutine k_set_one(nel, d)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)), intent(inout) :: d(nel)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nel
       d(i) = cmplx(1.d0, 0.d0, kind(0.d0))
    enddo
  end subroutine k_set_one

  !> d(1) = det(A)/|det A| for A = Q R from ZGEQRF (A holds R and the reflectors, tau their scalars): prod det(H_i) * prod R_ii/|R_ii|.
  !> The raw product of the R_ii is never formed (it under/overflows for large N).  A zero R_ii gives d = 0.  Serial loop on the device.
  Subroutine k_det_qr(n, A, lda, tau, d)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(in) :: A(lda,n), tau(n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: d(1:1)
    Complex (Kind=Kind(0.d0)) :: p
    Real (Kind=Kind(0.d0)) :: amod
    Integer :: i
    !$omp target
    p = cmplx(1.d0, 0.d0, kind(0.d0))
    do i = 1, n
       amod = abs(A(i,i))
       if (amod == 0.d0) then
          p = cmplx(0.d0, 0.d0, kind(0.d0))
       else
          p = p*(A(i,i)/amod)*refl_det(tau(i))
       endif
    enddo
    d(1) = p
    !$omp end target
  end subroutine k_det_qr

  !> cs(7) = conj(det)/|det| with det = det(U_R^dagger U_L) = conj(dR) dL, the tracked dets dw(1) = U_L, dw(2) = U_R (one tracked CGR
  !> phase, no LU).  bad = 1: fault injection (ALF_GPU_FAULT_STAB=det, the sign of the phase is flipped).
  Subroutine k_det_phase_u(dw, cs, bad)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: dw(2)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(8)
    Integer, intent(in) :: bad
    Complex (Kind=Kind(0.d0)) :: det
    !$omp target
    det = conjg(conjg(dw(2))*dw(1))
    det = det/abs(det)
    if (bad == 1) det = -det
    cs(7) = det
    !$omp end target
  end subroutine k_det_phase_u

!--------------------------------------------------------------------
!> @brief
!> H->D upload of the UDV stacks of the walker together with det(U) of every uploaded U (the only way the stacks may be re-uploaded).
!> @details gpu_put_udv, then a QR of each uploaded U on the device (udvst, udvl, udvr) gives its tracked det (k_det_qr).  A zero
!> or non-finite det leaves the slot untracked: the first CGR that needs it stops the run ("untracked U in CGR").
!> One D->H of the (2 + NSTM) N_FL_eff dets, to set the flags.
!> @param[in] walker walker holding udvst, udvl, udvr
!--------------------------------------------------------------------
  Subroutine gpu_stab_put_udv(walker)
    Implicit none
    Type (qmc_walker_t), intent(in) :: walker
    Integer :: nfe, slot, nst
    call gpu_put_udv(walker)
    if (s_proj) return                                    ! projector: CGRP takes its phase from the LU of S, no det tracking
    do nfe = 1, s_nfe
       do slot = 1, 2
          call upload_det(udw_U(1,1,slot,nfe), udw_det(slot,nfe))
       enddo
       do nst = 1, s_nstm
          call upload_det(udv_U(1,1,nst,nfe), udv_det(nst,nfe))
       enddo
    enddo
    ! D->H of the dets of the uploaded slots, to set the flags (once per run)
    !$omp target update from(udw_det, udv_det)
    do nfe = 1, s_nfe
       do slot = 1, 2
          w_det_ok(slot,nfe) = det_valid(udw_det(slot,nfe))
       enddo
       do nst = 1, s_nstm
          st_det_ok(nst,nfe) = det_valid(udv_det(nst,nfe))
       enddo
    enddo
  end subroutine gpu_stab_put_udv

  !> det of the uploaded U (n x n, device) into d(1): QR of a copy in sA
  Subroutine upload_det(U, d)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: U(s_N,s_N)
    Complex (Kind=Kind(0.d0)), intent(inout) :: d(1:1)
    integer :: info
    call k_copy_c(s_N*s_N, U, sA)
    call la_zgeqrf(s_N, s_N, sA, s_N, st_tau, st_work, s_lw, info)
    if (info /= 0) call lapack_fail('upload det: zgeqrf', info)
    call k_det_qr(s_N, sA, s_N, st_tau, d)
  end subroutine upload_det

  Logical function det_valid(d)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: d
    det_valid = Control_Is_Finite(dble(d)) .and. Control_Is_Finite(aimag(d)) .and. abs(d) > 0.d0
  end function det_valid

!====================================================================
!  Device kernels (plain arrays; synchronous target regions)
!====================================================================
  ! U(:,j) = U(:,j) * D(j)
  Subroutine k_scale_cols(n, U, ld, D)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_scale_cols(n, ld, U, D)
       return
    endif
#endif
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
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(n)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_diagabs(n, ld, U, D)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do i = 1, n
       D(i) = cmplx(abs(U(i,i)), 0.d0, kind(0.d0))
    enddo
  end subroutine k_diagabs

  ! U(i,j) = U(i,j) / D(i), j >= i   (D real valued)
  Subroutine k_divrows(n, U, ld, D)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: x
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_divrows(n, ld, U, D)
       return
    endif
#endif
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
  ! cs(6) = det(U_new) = prod_i det(H_i) * cs(1): U_new = Q diag(Phase,1,..,1) with Q = H_1 ... H_n (tau = the QR reflector scalars)
  Subroutine k_phase(n, U, ld, sgn, cj, cs, tau)
    Implicit none
    Integer, intent(in) :: n, ld, sgn
    Logical, intent(in) :: cj
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,n), tau(n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(8)
    Complex (Kind=Kind(0.d0)) :: ph, dq
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_phase(n, ld, U, sgn, cj, tau, cs)
       return
    endif
#endif
    !$omp target
    ph = cmplx(1.d0, 0.d0, kind(0.d0))
    dq = cmplx(1.d0, 0.d0, kind(0.d0))
    do i = 1, n
       ph = ph*U(i,i)
       dq = dq*refl_det(tau(i))
    enddo
    if (sgn < 0) ph = -ph
    if (cj) ph = conjg(ph)
    cs(1) = ph
    cs(2) = cmplx(1.d0, 0.d0, kind(0.d0)) / ph
    cs(6) = dq*ph
    !$omp end target
  end subroutine k_phase

  ! U(1,j) = cs(2) * U(1,j)   (ZSCAL of the first row of R with 1/Phase)
  Subroutine k_scale_row1(n, U, ld, cs)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: cs(8)
    Integer :: j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_scale_row1(n, ld, U, cs)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do j = 1, n
       U(1,j) = cs(2)*U(1,j)
    enddo
  end subroutine k_scale_row1

  ! U(i,1) = cs(1) * U(i,1)   (ZSCAL of the first column of Q with Phase)
  Subroutine k_scale_col1(n, U, ld, cs)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: cs(8)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_scale_col1(n, ld, U, cs)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do i = 1, n
       U(i,1) = cs(1)*U(i,1)
    enddo
  end subroutine k_scale_col1

  ! R(i,j) = U(i,j) for i <= j, else 0
  Subroutine k_triu(n, U, ld, R, ldr)
    Implicit none
    Integer, intent(in) :: n, ld, ldr
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: R(ldr,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_triu(n, ld, U, ldr, R)
       return
    endif
#endif
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
    Integer, intent(in) :: perm(n)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_perm(0, n, X, ldx, perm, W, ldw)
       return
    endif
#endif
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
    Integer, intent(in) :: perm(n)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_perm(1, n, X, ldx, perm, W, ldw)
       return
    endif
#endif
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
    Integer, intent(in) :: perm(n)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_perm(2, n, X, ldx, perm, W, ldw)
       return
    endif
#endif
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
    Integer, intent(in) :: perm(n)
    Complex (Kind=Kind(0.d0)), intent(in) :: X(ldx,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_perm(3, n, X, ldx, perm, W, ldw)
       return
    endif
#endif
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          W(i,perm(j)) = X(i,j)
       enddo
    enddo
  end subroutine k_iperm_cols

#if defined(STABLOG)
  ! W(i) = L(perm(i))
  Subroutine k_perm_vec(n, L, perm, W)
    Implicit none
    Integer, intent(in) :: n
    Integer, intent(in) :: perm(n)
    Real (Kind=Kind(0.d0)), intent(in) :: L(n)
    Real (Kind=Kind(0.d0)), intent(inout) :: W(n)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_perm_vec(n, L, perm, W)
       return
    endif
#endif
    !$omp target teams distribute parallel do
    do i = 1, n
       W(i) = L(perm(i))
    enddo
  end subroutine k_perm_vec
#endif

  ! B(i,j) = conjg(A(j,i))
  Subroutine k_ct(n, A, lda, B, ldb)
    Implicit none
    Integer, intent(in) :: n, lda, ldb
    Complex (Kind=Kind(0.d0)), intent(in) :: A(lda,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: B(ldb,n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_ct(n, A, lda, B, ldb)
       return
    endif
#endif
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_rowdiv(n, lda, A, D)
       return
    endif
#endif
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: sv
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_colinv(n, lda, A, D)
       return
    endif
#endif
    !$omp target teams distribute parallel do collapse(2) private(sv)
    do j = 1, n
       do i = 1, n
          sv = 1.d0 / real(D(j))
          A(i,j) = A(i,j) * sv
       enddo
    enddo
  end subroutine k_colinv

#if !defined(STABLOG)
  ! default CGR: T(i,j) = Dr(i) * T(i,j) * Dl(j) + RHS(i,j)
  Subroutine k_tpup(n, T, R, Dr, Dl)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: T(n,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: R(n,n), Dr(n), Dl(n)
    Integer :: i, j
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_tpup(n, T, R, Dr, Dl)
       return
    endif
#endif
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          T(i,j) = Dr(i) * T(i,j) * Dl(j)
          T(i,j) = T(i,j) + R(i,j)
       enddo
    enddo
  end subroutine k_tpup
#else
  ! STABLOG CGR: the four cases of cgr1_mod (exponent signs chosen so that no exp overflows)
  Subroutine k_tpup_log(n, T, R, Lr, Ll)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: T(n,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: R(n,n)
    Real (Kind=Kind(0.d0)), intent(in) :: Lr(n), Ll(n)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,n)
    Real (Kind=Kind(0.d0)), intent(in) :: L(n)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,n)
    Real (Kind=Kind(0.d0)), intent(in) :: L(n)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (L(j) > 0.d0) A(i,j) = cmplx(exp(-L(j)), 0.d0, kind(0.d0)) * A(i,j)
       enddo
    enddo
  end subroutine k_colscale_log
#endif

  ! cs(3) = ph0 * permutation sign * prod_i (R_ii phase) * reflector determinants   (cgr1_mod, after the QR); serial loop.
  ! ph0 = conj(det(U_R^dagger U_L))/|.| (host, from the LU), passed by value, or (dev = 1) the tracked phase cs(7) on the device.
  ! bad = 1: fault injection (phase sign flipped)
  Subroutine k_cgr_phase(n, T, ld, tau, cs, nvar, sgn, bad, ph0, dev)
    Implicit none
    Integer, intent(in) :: n, ld, nvar, sgn, bad, dev
    Complex (Kind=Kind(0.d0)), intent(in) :: T(ld,n), tau(n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: cs(8)
    Complex (Kind=Kind(0.d0)), intent(in) :: ph0
    Complex (Kind=Kind(0.d0)) :: ph, z
    Real (Kind=Kind(0.d0)) :: x
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_cgr_phase(n, ld, T, tau, cs, nvar, sgn, bad, ph0, dev)
       return
    endif
#endif
    ! map(to:): nvfortran 25.7 drops the imaginary part of an implicitly firstprivate complex scalar in a serial target region
    ! (KNOWN_ISSUES G5)
    !$omp target map(to: ph0)
    if (dev == 1) then
       ph = cs(7)
    else
       ph = ph0
    endif
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

  ! scale test: r(1) = number of NaN or Inf scales, r(2) = Re D(1), r(3) = Re D(n)
  Subroutine k_testscale(n, D, r)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Real (Kind=Kind(0.d0)), intent(inout) :: r(16)
    Integer :: i, nn
    Real (Kind=Kind(0.d0)) :: a, b
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_testscale(n, D, r, merge(1, 0, fault == 3))
       return
    endif
#endif
    nn = 0
    !$omp target teams distribute parallel do reduction(+:nn) private(a, b)
    do i = 1, n
       a = real(D(i)); b = aimag(D(i))
       if (is_nonfinite(a) .or. is_nonfinite(b)) nn = nn + 1
    enddo
    if (fault == 3) nn = nn + 1                             ! negative control (ALF_GPU_FAULT_STAB=nan)
    !$omp target
    r(1) = real(nn, kind(0.d0))
    r(2) = real(D(1))
    r(3) = real(D(n))
    !$omp end target
  end subroutine k_testscale

  ! precision check: r(4) = max |A-B|, r(5) = sum |A-B|, r(6)/r(7) = number of NaN or Inf entries in A / B (not in max/sum)
  Subroutine k_prec(n, A, B, r)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: A(n,n), B(n,n)
    Real (Kind=Kind(0.d0)), intent(inout) :: r(16)
    Integer :: i, j, na, nb
    Real (Kind=Kind(0.d0)) :: xm, xs, d, cs
    Logical :: fa, fb
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_prec(n, A, B, r)
       return
    endif
#endif
    xm = 0.d0; xs = 0.d0; na = 0; nb = 0
    ! the sum has a fixed order: one thread per column sums serially into st_pw(j), then one thread adds the n column
    ! sums in order (max and the integer counts are order independent and stay reductions)
    !$omp target teams distribute parallel do reduction(max:xm) reduction(+:na, nb) private(i, d, fa, fb, cs)
    do j = 1, n
       cs = 0.d0
       do i = 1, n
          fa = is_nonfinite(real(A(i,j))) .or. is_nonfinite(aimag(A(i,j)))
          fb = is_nonfinite(real(B(i,j))) .or. is_nonfinite(aimag(B(i,j)))
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

  ! true for NaN and +-Inf (bit test on the exponent; comparisons are folded away under -ffast-math)
  !> .true. if the scales of a UDV state (D, or L with STABLOG) are all finite (bit test, no FP exception)
  Logical function state_finite(st)
    Implicit none
    Class (UDV_State), intent(in) :: st
    Integer :: i
    state_finite = .true.
#if defined(STABLOG)
    do i = 1, size(st%L)
       if (is_nonfinite(st%L(i))) state_finite = .false.
    enddo
#else
    do i = 1, size(st%D)
       if (is_nonfinite(real(st%D(i))) .or. is_nonfinite(aimag(st%D(i)))) state_finite = .false.
    enddo
#endif
  end function state_finite

  Logical function is_nonfinite(x)
    Implicit none
    !$omp declare target
    Real (Kind=Kind(0.d0)), intent(in) :: x
    Integer(int64), parameter :: expmask = int(z'7FF0000000000000', int64)
    Integer(int64) :: b
    b = transfer(x, 0_int64)
    is_nonfinite = (iand(b, expmask) == expmask)
  end function is_nonfinite

#if defined(STABLOG)
  ! out(j) = log(|U(:,j)|_2) + L(j)   (the CPU's tmpnorm); scaled norm
  Subroutine k_lognorm(n, U, ld, L, out)
    Implicit none
    Integer, intent(in) :: n, ld
    Complex (Kind=Kind(0.d0)), intent(in) :: U(ld,n)
    Real (Kind=Kind(0.d0)), intent(in) :: L(n)
    Real (Kind=Kind(0.d0)), intent(inout) :: out(n)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: mx, s, v
    ! one thread per column, serial sums (fixed order; a device reduction combines in an unspecified order)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,n)
    Real (Kind=Kind(0.d0)), intent(inout) :: L(n)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
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
!  Hook and self-test
!====================================================================
!--------------------------------------------------------------------
!> @brief
!> Debug hook of Setup_gpu: arms the fault of ALF_GPU_FAULT_STAB (perm|phase|nan|det) and runs the self-test when
!> ALF_GPU_SELFTEST=3.  The self-test leaves the device slots dirty: afterwards they are reset to the identity states.
!--------------------------------------------------------------------
  Subroutine gpu_stab_hook()
    Implicit none
    Character (len=64) :: f
    Integer :: st, ib
    if (.not. gpu_stab_on) return
    if (gpu_env_get('ALF_GPU_FAULT_STAB', f)) then
       select case (trim(f))
       case ('perm');  fault = 1
       case ('phase'); fault = 2
       case ('nan');   fault = 3
       case ('det');   fault = 4
       end select
    endif
    if (gpu_env_is('ALF_GPU_SELFTEST', '3')) then
       if (s_proj) then
          call gpu_stab_selftest_proj()
       else
          call gpu_stab_selftest()
       endif
       do st = 1, s_nfe
          call gpu_stab_reset_l(st); call gpu_stab_reset_r(st)
          do ib = 1, s_nstm
             call gpu_stab_reset_st(ib, st)
          enddo
       enddo
    endif
  end subroutine gpu_stab_hook

  Subroutine say(line)
    Implicit none
    Character (len=*), intent(in) :: line
    call info_line(line)
    write(output_unit,'(1X,A)') trim(line)
  end subroutine say

  ! download / upload the state of a device slot to / from a host UDV_State object (self-test only)
  Subroutine dl_state(slot, nfe, ust)
    Implicit none
    Integer, intent(in) :: slot, nfe
    Class (UDV_State), intent(inout) :: ust
    ! self-test: D->H of the device state, the host oracle starts from it
    !$omp target update from(udw_U(:,:,slot,nfe), udw_V(:,:,slot,nfe))
    ust%U = udw_U(:,:,slot,nfe); ust%V = udw_V(:,:,slot,nfe)
#if defined(STABLOG)
    !$omp target update from(udw_L(:,slot,nfe))
    ust%L = udw_L(:,slot,nfe)
#else
    !$omp target update from(udw_D(:,slot,nfe))
    ust%D = udw_D(:,slot,nfe)
#endif
  end subroutine dl_state

  ! relative difference of two states: U, V (max |a-b| / max |b|), scales (componentwise relative, D; L absolute), and the
  ! matrix the state represents (U D V for side R, U D V^dagger for side L; default stabilisation only): invariant under the
  ! choice of the pivots
  Subroutine state_diff(a, b, dU, dV, dS, dP)
    Implicit none
    Class (UDV_State), intent(in) :: a, b
    Real (Kind=Kind(0.d0)), intent(out) :: dU, dV, dS, dP
    Integer :: i
    Real (Kind=Kind(0.d0)) :: d
    Complex (Kind=Kind(0.d0)), Allocatable :: Ma(:,:), Mb(:,:)
    dU = rel_diff(a%U, b%U)
    dV = rel_diff(a%V, b%V)
    dS = 0.d0
    dP = 0.d0
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
    allocate(Ma(size(a%U,1),size(a%V,1)), Mb(size(a%U,1),size(a%V,1)))
    Ma = represented(a); Mb = represented(b)
    dP = rel_diff(Ma, Mb)
    deallocate(Ma, Mb)
#endif
  end subroutine state_diff

#if !defined(STABLOG)
  function represented(a) result(M)
    Implicit none
    Class (UDV_State), intent(in) :: a
    Complex (Kind=Kind(0.d0)) :: M(size(a%U,1),size(a%V,1))
    Complex (Kind=Kind(0.d0)) :: UD(size(a%U,1),size(a%U,2))
    Integer :: j
    do j = 1, size(a%U,2)
       UD(:,j) = a%U(:,j)*a%D(j)
    enddo
    if (a%side == 'R' .or. a%side == 'r') then
       M = matmul(UD, a%V)
    else
       M = matmul(UD, conjg(transpose(a%V)))
    endif
  end function represented
#endif

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

  ! hp = a with the random matrix B = A + 1e-3 * noise
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
    Complex (Kind=Kind(0.d0)), intent(in) :: A(s_N,s_N), B(s_N,s_N)
    Real (Kind=Kind(0.d0)), intent(out) :: xmax_h, xmean_h, xmax_d, xsum_d
    Integer, intent(out) :: nna, nnb
    call COMPARE(A, B, xmax_h, xmean_h)
    st_Gn = A; st_Gt = B
    ! self-test: H->D of the two matrices of the reduction test
    !$omp target update to(st_Gn, st_Gt)
    call dev_prec(s_N, st_Gn, st_Gt, xmax_d, xsum_d, nna, nnb)
  end subroutine dev_check_prec

!--------------------------------------------------------------------
!> @brief
!> Self-test (ALF_GPU_SELFTEST=3): the fixed initial fields; one pass of the sweep.
!> @details
!> Left chain WRAPUL + decompose for NST = NSTM-1 .. 1, right chain WRAPUR + decompose for NST = 1 .. NSTM, and at every
!> stabilisation point CGR (NVAR = 1 and 2) from the states exchanged as in the sweep, each against the upstream host routines
!> (WRAPUL/WRAPUR/decompose/CGR on host UDV_State objects that hold the same input) and the precision reduction against
!> COMPARE.  Every block starts from the same input on both sides.  Compared: G and the CGR phase (tolerance 1e-10), the
!> matrix represented by the decomposed state U D V (default stabilisation; invariant under the choice of the pivots), and,
!> in a build whose host decomposition makes the same pivots (STABLOG, or -DQR_PREPIVOT), U, V and the scales one by one.
!> A failure ends the run with "GPU self-test FAILED".
!--------------------------------------------------------------------
  Subroutine gpu_stab_selftest()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-10, tolr = 1.d-12
    Class (UDV_State), Allocatable :: hl(:), hr(:), dl_(:), dr_(:)
    Integer :: nfe, nf, nst, nvar, N
    Integer :: NT1, NT
    Real (Kind=Kind(0.d0)) :: dU, dV, dS, dP, mu(2), mv(2), ms(2), mp(2), mg(2), mph(2), d, xmax_h, xmean_h, xmax_d, xsum_d
    Real (Kind=Kind(0.d0)) :: dpr(4)
    Complex (Kind=Kind(0.d0)), Allocatable :: Gh(:,:), Gd(:,:,:), Gp(:,:)
    Complex (Kind=Kind(0.d0)) :: ph_h, ph_d
    Integer(int64) :: seed
    Integer :: nna, nnb, nblk(2)
    Logical :: same_pivots
    Character (len=512) :: line
    Integer :: mode_save, ncmp_save
    Real (Kind=Kind(0.d0)) :: maxdiff_save

    N = s_N
    seed = 20251101_int64
    in_selftest = .true.
    ! the self-test runs in check mode whatever the module default is: the tracked phase is used and compared with the LU phase
    mode_save = cgr_mode; ncmp_save = ph_ncmp; maxdiff_save = ph_maxdiff
    cgr_mode = 2; ph_ncmp = 0; ph_maxdiff = 0.d0
#if defined(STABLOG) || defined(QR_PREPIVOT)
    same_pivots = .true.
#else
    same_pivots = .false.
#endif
    call say('GPU stabilization self-test (step 1d) start')
    allocate(hl(s_nfe), hr(s_nfe), dl_(s_nfe), dr_(s_nfe))
    allocate(Gh(N,N), Gp(N,N), Gd(N,N,N_FL))
    mu = 0.d0; mv = 0.d0; ms = 0.d0; mp = 0.d0; mg = 0.d0; mph = 0.d0; dpr = 0.d0; nblk = 0
    do nfe = 1, s_nfe
       call hl(nfe)%init(N, 'l'); call hr(nfe)%init(N, 'r')
       call dl_(nfe)%init(N, 'l'); call dr_(nfe)%init(N, 'r')
       call gpu_stab_reset_l(nfe); call gpu_stab_reset_r(nfe)
       call gpu_stab_reset_st(s_nstm, nfe)
    enddo
    ! ---- left chain
    do nst = s_nstm-1, 1, -1
       NT1 = s_nt(nst+1); NT = s_nt(nst)
       do nfe = 1, s_nfe
          call dl_state(1, nfe, hl(nfe))               ! the host oracle starts from the device input
          ! a non-finite device state fails here: the host decomposition would trap on it (-ffpe-trap, CI debug job)
          if (.not. state_finite(hl(nfe))) call selftest_fail('stab: non-finite scales in the device state (WRAPUL + decompose)')
       enddo
       call WRAPUL(NT1, NT, hl)
       call gpu_stab_wrapul(NT1, NT)
       do nfe = 1, s_nfe
          call dl_state(1, nfe, dl_(nfe))              ! device result for the comparison
          call state_diff(dl_(nfe), hl(nfe), dU, dV, dS, dP)
          call mx(mu(1), dU); call mx(mv(1), dV); call mx(ms(1), dS); call mx(mp(1), dP)
          call gpu_stab_store_l(nst, nfe)
       enddo
       nblk(1) = nblk(1) + 1
    enddo
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test WRAPUL+decompose over ', nblk(1), &
         & ' blocks: max rel diff U ', mu(1), ', V ', mv(1), ', scales ', ms(1), ', represented matrix ', mp(1), ' (tol ', tol
    call say(line)
    ! fail before the right chain uses the stored left states (a NaN in them would stop CGR with another message)
    if (.not. (mp(1) <= tol)) call selftest_fail('stab: WRAPUL + decompose (represented matrix)')
    if (same_pivots) then
       if (.not. (mu(1) <= tol .and. mv(1) <= tol .and. ms(1) <= tol)) call selftest_fail('stab: WRAPUL + decompose')
    endif
    ! ---- right chain with CGR at every stabilisation point
    do nst = 1, s_nstm
       NT1 = s_nt(nst-1); NT = s_nt(nst)
       do nfe = 1, s_nfe
          call dl_state(2, nfe, hr(nfe))
       enddo
       call WRAPUR(NT1, NT, hr)
       call gpu_stab_wrapur(NT1, NT)
       do nfe = 1, s_nfe
          call dl_state(2, nfe, dr_(nfe))
          call state_diff(dr_(nfe), hr(nfe), dU, dV, dS, dP)
          call mx(mu(2), dU); call mx(mv(2), dV); call mx(ms(2), dS); call mx(mp(2), dP)
       enddo
       nblk(2) = nblk(2) + 1
       ! CGR as in the sweep: udvl = udvst(NST), udvst(NST) = udvr (device), then compare with the host CGR on the same states
       do nfe = 1, s_nfe
          nf = Calc_Fl_map(nfe)
          call gpu_stab_exchange(nst, nfe, .true., hr(nfe))   ! (the host object is only written in the backward direction)
          call dl_state(1, nfe, hl(nfe))
          call dl_state(2, nfe, hr(nfe))
          do nvar = 1, 2
             call CGR(ph_h, nvar, Gh, hr(nfe), hl(nfe))
             call gpu_stab_cgr(nfe, nf, nvar, Gd, ph_d, .false.)
             d = rel_diff(Gd(:,:,nf), Gh)
             call mx(mg(nvar), d)
             call mx(mph(nvar), abs(ph_d - ph_h))
             ! precision reduction vs COMPARE: (new G, reference G) and (new G, G + noise of order 1e-3)
             call rand_perturb(Gh, Gp, seed)
             call dev_check_prec(Gd(:,:,nf), Gh, xmax_h, xmean_h, xmax_d, xsum_d, nna, nnb)
             call mx(dpr(1), abs(xmax_d - xmax_h) / max(xmax_h, tiny(1.d0)))
             call mx(dpr(2), abs(xsum_d/dble(N*N) - xmean_h) / max(xmean_h, tiny(1.d0)))
             call dev_check_prec(Gd(:,:,nf), Gp, xmax_h, xmean_h, xmax_d, xsum_d, nna, nnb)
             call mx(dpr(3), abs(xmax_d - xmax_h) / max(xmax_h, tiny(1.d0)))
             call mx(dpr(4), abs(xsum_d/dble(N*N) - xmean_h) / max(xmean_h, tiny(1.d0)))
             if (nna /= 0 .or. nnb /= 0) call selftest_fail('stab: non-finite count of the precision reduction')
          enddo
       enddo
    enddo
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test WRAPUR+decompose over ', nblk(2), &
         & ' blocks: max rel diff U ', mu(2), ', V ', mv(2), ', scales ', ms(2), ', represented matrix ', mp(2), ' (tol ', tol
    call say(line)
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test CGR at ', nblk(2), &
         & ' stabilization points x flavours: max rel diff G NVAR=1 ', mg(1), ', NVAR=2 ', mg(2), '; phase NVAR=1 ', mph(1), &
         & ', NVAR=2 ', mph(2), ' (tol ', tol
    call say(line)
    write(line,'(A,I0,A,ES10.3,A,ES8.1)') 'GPU stab self-test CGR phase, tracked det(U) against the LU phase over ', ph_ncmp, &
         & ' calls: max difference ', ph_maxdiff, ' (limit ', PH_ABORT
    call say(trim(line)//')')
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU stab self-test precision reduction vs COMPARE: rel diff max/mean ', &
         & dpr(1), ' / ', dpr(2), ' (G vs G_ref), ', max(dpr(3),dpr(4)), ' (G vs perturbed) (tol ', tolr
    call say(line)
    if (.not. (mp(1) <= tol .and. mp(2) <= tol)) call selftest_fail('stab: WRAPUL/WRAPUR + decompose (represented matrix)')
    if (same_pivots) then
       if (.not. (mu(1) <= tol .and. mv(1) <= tol .and. ms(1) <= tol .and. mu(2) <= tol .and. mv(2) <= tol .and. ms(2) <= tol)) &
            & call selftest_fail('stab: WRAPUL/WRAPUR + decompose')
    endif
    if (.not. (mg(1) <= tol .and. mg(2) <= tol .and. mph(1) <= tol .and. mph(2) <= tol)) call selftest_fail('stab: CGR')
    if (.not. (dpr(1) <= tolr .and. dpr(2) <= tolr .and. dpr(3) <= tolr .and. dpr(4) <= tolr)) &
         & call selftest_fail('stab: precision reduction')
    do nfe = 1, s_nfe
       call hl(nfe)%dealloc; call hr(nfe)%dealloc; call dl_(nfe)%dealloc; call dr_(nfe)%dealloc
    enddo
    deallocate(hl, hr, dl_, dr_, Gh, Gp, Gd)
    in_selftest = .false.
    cgr_mode = mode_save; ph_ncmp = ncmp_save; ph_maxdiff = maxdiff_save
    call say('GPU stab self-test PASSED')
  end subroutine gpu_stab_selftest

!====================================================================
!  Projector code (Projector = .T.): every U is Ndim x N_part (the first s_nc(nfe) columns of its slot), no V, no tracked det;
!  CGRP in place of CGR.  Finite-temperature runs never reach these routines (s_proj = .false.).
!====================================================================
!> @brief Projector setup: N_part per flavour (s_nc), WF_L%P and WF_R%P uploaded once (p_wf), CGRP buffers.
  Subroutine proj_setup()
    Implicit none
    Integer :: nfe, nf, np
    s_npm = 0
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe)
       s_nc(nfe) = size(WF_R(nf)%P, 2)
       if (size(WF_L(nf)%P, 2) /= s_nc(nfe) .or. size(WF_L(nf)%P, 1) /= s_N .or. size(WF_R(nf)%P, 1) /= s_N) then
          write(error_unit,*) 'GPU stabilisation (projector): WF_L and WF_R of flavour ', nf, ' differ in shape'
          call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
       endif
       s_npm = max(s_npm, s_nc(nfe))
    enddo
    allocate(p_wf(s_N, s_npm, 2, s_nfe), p_S(s_npm, s_npm), p_R(s_npm, s_N))
    p_wf = z0
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe); np = s_nc(nfe)
       p_wf(:,1:np,1,nfe) = WF_L(nf)%P
       p_wf(:,1:np,2,nfe) = WF_R(nf)%P
    enddo
    ! H->D once: the trial wave functions (resets of the sweeps)
    !$omp target enter data map(to: p_wf) map(alloc: p_S, p_R)
  end subroutine proj_setup

!> @brief Projector decompose of an Ndim x np state (decompose_UDV_state without V): QR, U = Q, column 1 times the phase.
!> @details Pre-pivoted QR (la_zqrpp; the CPU uses ZGEQP3, or a log-norm sort in STABLOG builds): another pivot order only
!> changes U by a matrix of positive determinant (Q R P^-1 times the phase of det R P^-1), so CGRP's G and phase are the same.
  Subroutine dev_decompose_p(side, np, U)
    Implicit none
    Character (len=1), intent(in) :: side
    Integer, intent(in) :: np
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(s_N,s_N)
    Integer :: n, sgn, info
    Complex (Kind=Kind(0.d0)) :: ph1

    n = s_N
    call la_zqrpp(n, np, U, n, st_ipvt, st_tau, st_work, s_lw, info)
    if (info /= 0) call lapack_fail('decompose (projector): la_zqrpp (-101: NaN or Inf in the matrix to decompose)', info)
    call fault_perm()
    call k_diagabs(np, U, n, st_Dl)
    call k_divrows(np, U, n, st_Dl)                       ! R(i,j) /= |R(i,i)| (QDRP_decompose)
    ph1 = z1
    call Pivot_phase(ph1, st_ipvt, np)
    sgn = nint(real(ph1))
    call k_phase(np, U, n, sgn, .not. (side == 'R' .or. side == 'r'), st_c, st_tau)
    call k_ident(n, sC, n)
    call la_zunmqr('L', 'N', n, np, np, U, n, st_tau, sC, n, st_work, s_lw, info)   ! Q (first np columns)
    if (info /= 0) call lapack_fail('decompose (projector): zunmqr', info)
    call k_copy_c(n*np, sC, U)
    call k_scale_col1(n, U, n, st_c)                      ! first column of U times the phase
  end subroutine dev_decompose_p

!> @brief CGRP (cgr1_mod) on the device states of flavour nfe: S = U_L^dagger U_R, G = 1 - U_R S^-1 U_L^dagger -> st_Gn,
!> phase of the LU of S -> st_c(3).
!> @details The LU of S (np x np) runs on host LAPACK (16 np^2 B each way): the phase is CGRP's loop bit for bit, the pivots
!> become one row permutation, applied while R = U_L^dagger is formed; the two triangular solves (= ZGETRS) and the GEMM run
!> on the device.
  Subroutine dev_cgrp(nfe)
    Implicit none
    Integer, intent(in) :: nfe
    Integer :: n, np, ld, i, j, info
    Complex (Kind=Kind(0.d0)), Allocatable :: hS(:,:)
    Integer, Allocatable :: hp(:)
    Complex (Kind=Kind(0.d0)) :: ph

    n = s_N; np = s_nc(nfe); ld = s_npm
    call la_zgemm('C', 'N', np, np, n, z1, udw_U(1,1,1,nfe), n, udw_U(1,1,2,nfe), n, z0, p_S, ld)
#if defined(ALF_LA_CUDA)
    call la_cuda_sync()                                   ! S complete before the D->H
#endif
    !$omp target update from(p_S)
    allocate(hS(np,np), hp(np))
    hS = p_S(1:np,1:np)                                   ! unmapped copy: linalg_mod takes the host path
    call la_zgetrf(np, np, hS, np, hp, info)
    ph = z1
    do i = 1, np
       if (hp(i) /= i) then
          ph = -ph * hS(i,i)/abs(hS(i,i))
       else
          ph =  ph * hS(i,i)/abs(hS(i,i))
       endif
    enddo
    if (fault == 2) ph = -ph                              ! negative control: flipped phase
    do i = 1, np
       st_perm(i) = i
    enddo
    do i = 1, np                                          ! ZLASWP: row i <-> row hp(i), in order
       j = st_perm(i); st_perm(i) = st_perm(hp(i)); st_perm(hp(i)) = j
    enddo
    if (fault == 1 .and. np > 1) then                     ! negative control: wrong row permutation (a wrong pivot sign of the
       j = st_perm(1); st_perm(1) = st_perm(2); st_perm(2) = j   ! decompose alone is a constant sign of every CGRP: invisible)
    endif
    p_S(1:np,1:np) = hS
    st_c(3) = ph
    !$omp target update to(p_S, st_c(3:3), st_perm(1:np))
    deallocate(hS, hp)
    ! R = P U_L^dagger (np x n)
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, np
          p_R(i,j) = conjg(udw_U(j, st_perm(i), 1, nfe))
       enddo
    enddo
    call la_ztrsm('L', 'L', 'N', 'U', np, n, z1, p_S, ld, p_R, ld)
    call la_ztrsm('L', 'U', 'N', 'N', np, n, z1, p_S, ld, p_R, ld)
    call k_ident(n, st_Gn, n)
    call la_zgemm('N', 'N', n, n, np, -z1, udw_U(1,1,2,nfe), n, p_R, ld, z1, st_Gn, n)
  end subroutine dev_cgrp

!> @brief Projector, Ltau = 1: D->H of the device udvr (slot 2) into the host walker%udvr, which the host Tau_p reads.
!> @param[inout] udvr host objects walker%udvr (U, N_part columns, side)
  Subroutine gpu_stab_export_r(udvr)
    Implicit none
    Class (UDV_State), Dimension(:), intent(inout) :: udvr
    Integer :: nfe
    do nfe = 1, s_nfe
       !$omp target update from(udw_U(:,1:s_nc(nfe),2,nfe))
       udvr(nfe)%U = udw_U(:,1:s_nc(nfe),2,nfe)
       udvr(nfe)%side = 'r'
    enddo
  end subroutine gpu_stab_export_r

!> @brief Self-test of the projector path (ALF_GPU_SELFTEST=3, Projector = .T.): left chain from WF_L, right chain from WF_R
!> with WRAPUL / WRAPUR + decompose on host and device, CGRP at every stabilisation point: G and phase against the host CGRP
!> (both are independent of the pivot choice of the decompositions).  ALF_GPU_FAULT_STAB=perm (CGRP row permutation) and phase
!> must make it fail.
  Subroutine gpu_stab_selftest_proj()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-10
    Class (UDV_State), Allocatable :: hl(:), hr(:), hst(:,:)
    Complex (Kind=Kind(0.d0)), Allocatable :: Gh(:,:), Gd(:,:,:)
    Complex (Kind=Kind(0.d0)) :: ph_h, ph_d
    Integer :: N, nfe, nf, nst, NT1, NT, nblk
    Real (Kind=Kind(0.d0)) :: mg, mph
    Character (len=256) :: line

    N = s_N
    in_selftest = .true.
    call say('GPU stabilization self-test (projector) start')
    allocate(hl(s_nfe), hr(s_nfe), hst(s_nstm, s_nfe), Gh(N,N), Gd(N,N,N_FL))
    Gd = z0
    do nfe = 1, s_nfe
       nf = Calc_Fl_map(nfe)
       call hl(nfe)%init(N, 'l', WF_L(nf)%P); call hr(nfe)%init(N, 'r', WF_R(nf)%P)
       do nst = 1, s_nstm
          call hst(nst,nfe)%alloc(N, s_nc(nfe))
       enddo
       hst(s_nstm,nfe) = hl(nfe)
       call gpu_stab_reset_l(nfe); call gpu_stab_reset_r(nfe)
       call gpu_stab_reset_st(s_nstm, nfe)
    enddo
    do nst = s_nstm-1, 1, -1
       NT1 = s_nt(nst+1); NT = s_nt(nst)
       call WRAPUL(NT1, NT, hl)
       call gpu_stab_wrapul(NT1, NT)
       do nfe = 1, s_nfe
          hst(nst,nfe) = hl(nfe)
          call gpu_stab_store_l(nst, nfe)
       enddo
    enddo
    mg = 0.d0; mph = 0.d0; nblk = 0
    do nst = 1, s_nstm
       NT1 = s_nt(nst-1); NT = s_nt(nst)
       call WRAPUR(NT1, NT, hr)
       call gpu_stab_wrapur(NT1, NT)
       do nfe = 1, s_nfe
          nf = Calc_Fl_map(nfe)
          call gpu_stab_exchange(nst, nfe, .true., hr(nfe))   ! device: udvl = udvst(nst) (host object not written)
          call CGRP(ph_h, Gh, hr(nfe), hst(nst,nfe))
          call gpu_stab_cgr(nfe, nf, 1, Gd, ph_d, .false.)
          call mx(mg, rel_diff(Gd(:,:,nf), Gh))
          call mx(mph, abs(ph_d - ph_h))
       enddo
       nblk = nblk + 1
    enddo
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES8.1,A)') 'GPU stab self-test projector CGRP at ', nblk, &
         & ' stabilization points x flavours: max rel diff G ', mg, ', phase ', mph, ' (tol ', tol, ')'
    call say(line)
    if (.not. (mg <= tol .and. mph <= tol)) call selftest_fail('stab: projector CGRP')
    do nfe = 1, s_nfe
       call hl(nfe)%dealloc; call hr(nfe)%dealloc
       do nst = 1, s_nstm
          call hst(nst,nfe)%dealloc
       enddo
    enddo
    deallocate(hl, hr, hst, Gh, Gd)
    in_selftest = .false.
    call say('GPU stab self-test PASSED')
  end subroutine gpu_stab_selftest_proj

End Module gpu_stab_mod
#endif

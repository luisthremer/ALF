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
!> GPU port, P2.4 step 1e: the time-displaced pass TAU_M (tau_m_mod) of the finite-temperature code on the device.
!> Compiled only with -DALF_GPU.
!>
!> @details
!> Called at the end of Sweep_down when gpu_tau_on (device stabilisation, Ltau = 1, device tau observables of gpu_obs_mod;
!> no Projector, no Langevin/HMC, no weight reconstruction: refused by gpu_supported).
!>  * GT0, G0T, GTT (gpu_state_mod GT0_d, G0T_d, GTT_d) and G00 (G00_d, here) stay on the device during the pass.
!>  * per slice, the propagation of TAU_M: GT0 = PROPR, G0T = PROPRM1, GTT = PROPRM1 then PROPR (device kernels of gpu_wrap_mod).
!>  * at every stabilisation point Stab_nt(NST): the right propagation WRAPUR + decompose of a fresh udvr (device slot 2 of
!>    gpu_stab_mod, reset at the start of the pass) and CGR2_2 (below) with the left state udvst(NST) read on the device (the D->H
!>    export of udvst of gpu_stab_exchange is switched off by gpu_stab_tau_dev), then the four Control_Precision_tau checks as
!>    device reductions (Control_Precision_tau_vals).
!>  * observables of every slice: gpu_obs_tau_measure on the device matrices.
!>
!> CGR2_2 on the device (dev_cgr22), the algorithm of cgr2_2_mod (default, STAB3 and STABLOG; STAB1/2 are refused), 2N x 2N, N = Ndim:
!>   MYU2 = U2^dagger, V1INV = V1^-1 (ZGETRF + ZGETRI as INV), D1m, D2m = D, or the STAB3 / STABLOG splitting, the branch D1(1) > D2(1)
!>   (L1(1) > L2(1) in STABLOG: two scalars D->H), HLPB2 from the four N x N blocks, HLPB1 = HLPB2^dagger, QR: default = column-norm
!>   pre-pivot + unpivoted ZGEQRF (la_zqrpp, as the device decomposition of gpu_stab_mod), STAB3, STABLOG = host ZGEQP3 with D->H / H->D of
!>   the matrix; D separated as |R_ii|; solve_extended_System = rows of [[UCT,0],[0,VINV]] permuted by IPVT, ZTRSM('L','U','C'), row
!>   scaling with conj(1/D), ZUNMQR('L','N'); get_blocks into G00 / G0T / GT0 / GTT.
!>
!> Transfers per stabilisation point and flavour: D1(1), D2(1) (32 B, or L1(1), L2(1) 16 B) D->H; precision numbers 4 x 32 B D->H;
!> STAB3 and STABLOG additionally the 2N x 2N matrix both ways for ZGEQP3.  Nothing per slice.
!>
!> Negative controls ALF_GPU_FAULT_TAU = perm (wrong permutation), prop (wrong propagation of GT0 at slice 2), prec (wrong precision
!> reduction), blocks (wrong block of GT0).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_tau_mod

  Use iso_fortran_env, only: output_unit, error_unit
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Ltrot, Calc_Fl_map, Projector, reconstruction_needed
  Use QMC_runtime_var, only: get_Ltau
  Use Control, only: Control_Precision_tau_vals
  Use MyMats, only: COMPARE
  Use UDV_State_mod, only: UDV_State
  Use wrapur_mod, only: WRAPUR
  Use cgr2_2_mod, only: CGR2_2
  Use tau_m_mod, only: PropR, PropRm1
  Use linalg_mod, only: la_zgemm, la_ztrsm, la_zunmqr, la_zgetrf, la_zgetri, la_zqrpp
#if defined(ALF_LA_CUDA)
  Use linalg_cuda_mod, only: la_cuda_sync
#endif
#if defined(STABLOG) || defined(STAB3)
  Use linalg_mod, only: la_zgeqp3
#endif
  Use qmc_walker_mod, only: qmc_walker_t
  Use gpu_state_mod, only: gpu_active, G_d, g_dev_ok, GT0_d, G0T_d, GTT_d, udv_U, udv_V, udv_D, udv_L, udw_U, udw_V, udw_D, udw_L, &
       &                   info_line, gpu_env_get, gpu_env_is, rel_diff, selftest_fail
  Use gpu_wrap_mod, only: gpu_mmthr, gpu_mmthl_m1, gpu_Vslice_R, gpu_Vslice_L
  Use gpu_stab_mod, only: gpu_stab_on, gpu_stab_tau_dev, gpu_stab_wrapur, gpu_stab_reset_r, dev_prec, k_ident, k_copy_c, k_ct, &
       &                  k_diagabs, k_divrows, say
#if defined(STAB3)
  Use gpu_stab_mod, only: k_rowscale_d
#endif
  Use gpu_obs_mod, only: gpu_obs_tau_active, gpu_obs_tau_measure
  Implicit none

  Private
  Public :: gpu_tau_on, gpu_tau_setup, gpu_tau_finalize, gpu_tau_pass, gpu_tau_selftest

  Logical, Save :: gpu_tau_on = .false.
  Integer, Save :: t_N = 0, t_N2 = 0, t_nfe = 0, t_lw = 0
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

  ! ---- own device buffers (allocated once in gpu_tau_setup) ----
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: G00_d(:,:,:)                 ! N x N x N_FL
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_H1(:,:), t_H2(:,:)          ! 2N x 2N: QR matrix / solution
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_V1i(:,:), t_U2h(:,:)        ! N x N: V1^-1, U2^dagger
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_old(:,:,:)                 ! N x N x 3: GTT, GT0, G0T before CGR2_2
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_GR(:,:,:)                  ! N x N x N_FL: G at time 0 (precision reference)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_tau(:), t_D3(:), t_Dm1(:), t_Dm2(:), t_work(:)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: t_Rs(:,:)                    ! N x 2 row scale factors (STABLOG)
  Integer, Allocatable, Save :: t_ipvt(:), t_ipiv(:)
#if defined(STABLOG) || defined(STAB3)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: t_rw(:)                      ! RWORK of ZGEQP3
#endif
  Integer, Save :: fault = 0                                                   ! 1 perm, 2 prop, 3 prec, 4 blocks
  Integer, Save :: n_pass = 0, n_cgr = 0
  Logical, Save :: prec_report = .true.                                      ! .false. in the self-test: the Control counters are not touched
  Real (Kind=Kind(0.d0)), Save :: last_prec(2,4) = 0.d0                         ! last reduction: (max, sum) of the four pairs of tau_stab_cgr
  Integer, Save :: prec_k = 0

Contains

!--------------------------------------------------------------------
!> @brief
!> Decides once whether the device tau pass is used, allocates its buffers, writes one info line.
!> @details
!> On iff device stabilisation, Ltau = 1 and the device tau observables are active. Called at the end of Setup_gpu after
!> gpu_obs_report. If it stays off, ALF_GPU_FAULT_TAU and ALF_GPU_SELFTEST=4 stop the run (they would act on nothing).
!--------------------------------------------------------------------
  Subroutine gpu_tau_setup(walker)
    Implicit none
    Type (qmc_walker_t), intent(in) :: walker
    Integer :: n, n2, lw, info
    Complex (Kind=Kind(0.d0)) :: wq(1)
    Character (len=256) :: line, f
    Character (len=96) :: reason

    gpu_tau_on = .false.
    reason = ''
    if (.not. gpu_active .or. .not. gpu_stab_on) then
       reason = 'device stabilisation inactive'
    elseif (get_Ltau() /= 1) then
       reason = 'Ltau /= 1'
    elseif (Projector .or. reconstruction_needed) then
       reason = 'Projector or reconstruction'
    elseif (.not. gpu_obs_tau_active()) then
       reason = 'device tau observables inactive'
    endif
    if (len_trim(reason) > 0) then
       call tau_off_refuse(trim(reason))
       if (gpu_active) call info_line('GPU tau pass (step 1e): OFF, reason: '//trim(reason))
       return
    endif
    t_N = Ndim; t_N2 = 2*Ndim; t_nfe = N_FL_eff
    n = t_N; n2 = t_N2
    allocate(t_H1(n2,n2), t_H2(n2,n2), t_V1i(n,n), t_U2h(n,n), t_old(n,n,3), t_GR(n,n,N_FL))
    allocate(t_tau(n2), t_D3(n2), t_Dm1(n), t_Dm2(n), t_Rs(n,2), t_ipvt(n2), t_ipiv(n), G00_d(n,n,N_FL))
    t_ipvt = 1; t_ipiv = 1
    ! device allocation once, no transfer; the first use of each buffer is a device write
    !$omp target enter data map(alloc: t_H1, t_H2, t_V1i, t_U2h, t_old, t_GR, t_tau, t_D3, t_Dm1, t_Dm2, t_Rs, t_ipvt, t_ipiv, G00_d)
    ! LAPACK work space: the largest of the queries of the calls of dev_cgr22
    call la_zqrpp(n2, n2, t_H1, n2, t_ipvt, t_tau, wq, -1, info)
    lw = int(real(wq(1)))
    call la_zunmqr('L', 'N', n2, n2, n2, t_H1, n2, t_tau, t_H2, n2, wq, -1, info)
    lw = max(lw, int(real(wq(1))))
    call la_zgetri(n, t_V1i, n, t_ipiv, wq, -1, info)
    lw = max(lw, int(real(wq(1))))
#if defined(STABLOG) || defined(STAB3)
    allocate(t_rw(2*n2))
    call la_zgeqp3(n2, n2, t_H1, n2, t_ipvt, t_tau, wq, -1, t_rw, info)
    lw = max(lw, int(real(wq(1))))
#endif
    t_lw = max(lw, 1)
    allocate(t_work(t_lw))
    !$omp target enter data map(alloc: t_work)
    if (gpu_env_get('ALF_GPU_FAULT_TAU', f)) then
       select case (trim(f))
       case ('perm');   fault = 1
       case ('prop');   fault = 2
       case ('prec');   fault = 3
       case ('blocks'); fault = 4
       end select
    endif
    n_pass = 0; n_cgr = 0
    gpu_tau_on = .true.
    gpu_stab_tau_dev = .true.
    write(line,'(A,I0,A,I0,A,I0)') 'GPU tau pass (step 1e): ACTIVE (propagation, WRAPUR, CGR2_2, precision check and observables on the device), N = ', &
         & n, ', Nstm = ', walker%NSTM, ', flavours ', t_nfe
#if defined(STABLOG)
    line = trim(line)//', CGR2_2 QR: host ZGEQP3 (STABLOG)'
#elif defined(STAB3)
    line = trim(line)//', CGR2_2 QR: host ZGEQP3, D+/D- split (STAB3)'
#else
    line = trim(line)//', CGR2_2 QR: pre-pivoted (default)'
#endif
    call info_line(trim(line))
  end subroutine gpu_tau_setup

!> The device tau pass is off: stop with an error if a variable is set that acts only in it (ALF_GPU_FAULT_TAU, ALF_GPU_SELFTEST=4),
!> else an rc 0 could not tell right from wrong.
  Subroutine tau_off_refuse(reason)
    Implicit none
    Character (len=*), intent(in) :: reason
    Character (len=256) :: v
    if (gpu_env_is('ALF_GPU_SELFTEST', '4')) &
         & call tau_off_fail('GPU tau pass (step 1e): OFF ('//reason//'), but ALF_GPU_SELFTEST=4 is set: the tau self-test cannot run (refused)')
    if (gpu_env_get('ALF_GPU_FAULT_TAU', v)) &
         & call tau_off_fail('GPU tau pass (step 1e): OFF ('//reason//'), but ALF_GPU_FAULT_TAU is set: it would act on nothing (refused)')
  end subroutine tau_off_refuse

  Subroutine tau_off_fail(text)
    Implicit none
    Character (len=*), intent(in) :: text
    call say(text)
    flush(output_unit)
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine tau_off_fail

!> @brief Writes the counters and releases the buffers.
  Subroutine gpu_tau_finalize()
    Implicit none
    Character (len=256) :: line
    if (.not. gpu_tau_on) return
    write(line,'(A,I0,A,I0)') 'GPU tau pass counters: passes ', n_pass, ', CGR2_2 calls ', n_cgr
    call info_line(trim(line))
    !$omp target exit data map(delete: t_H1, t_H2, t_V1i, t_U2h, t_old, t_GR, t_tau, t_D3, t_Dm1, t_Dm2, t_Rs, t_ipvt, t_ipiv, G00_d, t_work)
    deallocate(t_H1, t_H2, t_V1i, t_U2h, t_old, t_GR, t_tau, t_D3, t_Dm1, t_Dm2, t_Rs, t_ipvt, t_ipiv, G00_d, t_work)
#if defined(STABLOG) || defined(STAB3)
    deallocate(t_rw)
#endif
    gpu_tau_on = .false.
    gpu_stab_tau_dev = .false.
  end subroutine gpu_tau_finalize

!====================================================================
!  The pass (what TAU_M does, without the Langevin / HMC / reconstruction parts, which gpu_supported refuses)
!====================================================================
!> @brief One pass of TAU_M on the device: the four Green functions on all slices, observables, CGR2_2 at the stabilisation points.
  Subroutine gpu_tau_pass(walker)
    Implicit none
    Type (qmc_walker_t), intent(in) :: walker
    Integer :: nfe, NT, NT1, NTST, NST
    Real (Kind=Kind(0.d0)) :: Mc_step_Weight

    Mc_step_Weight = 1.d0
    call tau_begin(walker%GR)
    call gpu_obs_tau_measure(0, GT0_d, G0T_d, G00_d, GTT_d, walker%PHASE, Mc_step_Weight)
    NST = 1
    do NT = 0, Ltrot - 1
       NT1 = NT + 1
       call tau_slice(NT1)
       call gpu_obs_tau_measure(NT1, GT0_d, G0T_d, G00_d, GTT_d, walker%PHASE, Mc_step_Weight)
       if (walker%Stab_nt(NST) == NT1) then
          NTST = walker%Stab_nt(NST-1)
          call gpu_stab_wrapur(NTST, NT1)          ! WRAPUR(NTST, NT1, udvr) of TAU_M: device slot udvr, decompose, scale test
          do nfe = 1, t_nfe
             call tau_stab_cgr(nfe, NST)           ! CGR2_2 + the four precision checks
          enddo
          NST = NST + 1
       endif
    enddo
    call tau_end()
#if defined(ALF_LA_CUDA)
    call la_cuda_sync()                            ! cuSOLVER info of the tau pass
#endif
    n_pass = n_pass + 1
  end subroutine gpu_tau_pass

!> Start of the pass: G00 = GT0 = GTT = G, G0T = -(1 - G) (TAU_M, tau = 0), from G_d on the device.
  Subroutine tau_begin(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL)
    Integer :: nfe, nf
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       if (g_dev_ok(nf)) then
          call k_copy_c(t_N*t_N, G_d(1,1,nf), G00_d(1,1,nf))
       else
          ! H->D: the device copy of G is stale (a host writer ran since the last device update); once per flavour and pass
          G00_d(:,:,nf) = GR(:,:,nf)
          !$omp target update to(G00_d(:,:,nf))
       endif
       call k_tau_init(t_N, G00_d(1,1,nf), GT0_d(1,1,nf), GTT_d(1,1,nf), G0T_d(1,1,nf))
       call k_copy_c(t_N*t_N, G00_d(1,1,nf), t_GR(1,1,nf))
       call gpu_stab_reset_r(nfe)         ! fresh udvr (slot 2): the right propagation of TAU_M starts from the unit state
    enddo
  end subroutine tau_begin

  Subroutine tau_end()
    Implicit none
    integer :: nfe
    do nfe = 1, t_nfe
       call gpu_stab_reset_r(nfe)         ! leave slot 2 clean (the sweeps reset it before the next use anyway)
    enddo
  end subroutine tau_end

!> One time slice of TAU_M for all flavours: PROPR(GT0), PROPRM1(G0T), PROPRM1+PROPR(GTT) at slice nt1.
  Subroutine tau_slice(nt1)
    Implicit none
    Integer, intent(in) :: nt1
    Integer :: nfe, nf, N
    N = t_N
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       ! PROPR: B_T from the left, then Op_mmultR('n') for n = 1..N_op
       call gpu_mmthr(GT0_d(1,1,nf), N, N, nf, nt1)
       call gpu_Vslice_R(GT0_d(1,1,nf), N, N, nf, nt1, 'n', .false.)
       ! PROPRM1: B_T^-1 from the right, then Op_mmultL('n', sign = -1) for n = 1..N_op
       call gpu_mmthl_m1(G0T_d(1,1,nf), N, N, nf, nt1)
       call gpu_Vslice_L(G0T_d(1,1,nf), N, N, nf, nt1, 'n', -1, .false.)
       call gpu_mmthl_m1(GTT_d(1,1,nf), N, N, nf, nt1)
       call gpu_Vslice_L(GTT_d(1,1,nf), N, N, nf, nt1, 'n', -1, .false.)
       call gpu_mmthr(GTT_d(1,1,nf), N, N, nf, nt1)
       call gpu_Vslice_R(GTT_d(1,1,nf), N, N, nf, nt1, 'n', .false.)
       if (fault == 2 .and. nt1 == 2) call k_fault_scale(N, GT0_d(1,1,nf))    ! negative control: wrong propagation
    enddo
  end subroutine tau_slice

!====================================================================
!  Stabilisation point: CGR2_2 and the precision checks (per flavour)
!====================================================================
  Subroutine tau_stab_cgr(nfe, nst)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Integer :: nf, n

    n = t_N
    nf = Calc_Fl_map(nfe)
    ! HLP4/5/6 of TAU_M: the propagated matrices before CGR2_2 (device copies)
    call k_copy_c(n*n, GTT_d(1,1,nf), t_old(1,1,1))
    call k_copy_c(n*n, GT0_d(1,1,nf), t_old(1,1,2))
    call k_copy_c(n*n, G0T_d(1,1,nf), t_old(1,1,3))
    call dev_cgr22(nfe, nst, G00_d(1,1,nf), GT0_d(1,1,nf), GTT_d(1,1,nf), G0T_d(1,1,nf))
    n_cgr = n_cgr + 1
    ! Control_Precision_tau(GR, G00), (HLP4, GTT), (HLP5, GT0), (HLP6, G0T)
    call prec_pair(1, t_GR(1,1,nf), G00_d(1,1,nf))
    call prec_pair(2, t_old(1,1,1), GTT_d(1,1,nf))
    call prec_pair(3, t_old(1,1,2), GT0_d(1,1,nf))
    call prec_pair(4, t_old(1,1,3), G0T_d(1,1,nf))
  end subroutine tau_stab_cgr

!> Control_Precision_tau(A, B) as a device reduction (dev_prec of gpu_stab_mod: one D->H of four numbers, 32 B)
  Subroutine prec_pair(k, A, B)
    Implicit none
    Integer, intent(in) :: k
    Complex (Kind=Kind(0.d0)), intent(in) :: A(t_N,t_N), B(t_N,t_N)
    Real (Kind=Kind(0.d0)) :: xmax, xsum
    Integer :: nna, nnb
    call dev_prec(t_N, A, B, xmax, xsum, nna, nnb)
    if (fault == 3) xmax = xmax*(1.d0 + 1.d-6)            ! negative control: wrong reduction
    last_prec(1,k) = xmax; last_prec(2,k) = xsum
    if (prec_report) call Control_Precision_tau_vals(xmax, xsum, nna, nnb, t_N)
  end subroutine prec_pair

!> Stops the run after a failed LAPACK call
  Subroutine tau_fail(what, info)
    Implicit none
    Character (len=*), intent(in) :: what
    Integer, intent(in) :: info
    write(error_unit,*)
    write(error_unit,*) 'Error in the GPU tau pass: ', what, ', info = ', info
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine tau_fail

!====================================================================
!  CGR2_2 on the device
!====================================================================
!> CGR2_2(GRT0, GR00, GRTT, GR0T, udv2 = udvr (slot 2), udv1 = udvst(nst), LQ) of cgr2_2_mod for flavour nfe.  Outputs: N x N
!> device arrays.
  Subroutine dev_cgr22(nfe, nst, G00, GT0, GTT, G0T)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Complex (Kind=Kind(0.d0)), intent(inout) :: G00(t_N,t_N), GT0(t_N,t_N), GTT(t_N,t_N), G0T(t_N,t_N)
    Integer :: n, n2, br, info
    Real (Kind=Kind(0.d0)) :: d1, d2
    Integer :: kk

    n = t_N; n2 = t_N2
    ! the branch decision needs D1(1) and D2(1) (L1(1), L2(1)) on the host
#if defined(STABLOG)
    ! D->H: L1(1), L2(1), 16 B (branch of CGR2_2)
    !$omp target update from(udv_L(1:1,nst,nfe), udw_L(1:1,2,nfe))
    d1 = udv_L(1,nst,nfe); d2 = udw_L(1,2,nfe)
#else
    ! D->H: D1(1), D2(1), 32 B (branch of CGR2_2)
    !$omp target update from(udv_D(1:1,nst,nfe), udw_D(1:1,2,nfe))
    d1 = real(udv_D(1,nst,nfe)); d2 = real(udw_D(1,2,nfe))
#endif
    br = 2
    if (d1 > d2) br = 1
    ! MYU2 = U2^dagger
    call k_ct(n, udw_U(1,1,2,nfe), n, t_U2h, n)
    ! V1INV = V1^-1: ZGETRF + ZGETRI as INV
    call k_copy_c(n*n, udv_V(1,1,nst,nfe), t_V1i)
    call la_zgetrf(n, n, t_V1i, n, t_ipiv, info)
    if (info < 0) call tau_fail('CGR2_2: zgetrf', info)
    call la_zgetri(n, t_V1i, n, t_ipiv, t_work, t_lw, info)
    if (info /= 0) call tau_fail('CGR2_2: zgetri', info)
#if defined(STABLOG)
    ! scales smaller than 1 stay in D1m / D2m, larger ones move to the rows of V1INV / MYU2 (STABLOG branch of CGR2_2)
    call k_dm_log(n, udv_L(1,nst,nfe), t_Dm1, t_Rs(1,1))
    call k_dm_log(n, udw_L(1,2,nfe), t_Dm2, t_Rs(1,2))
    call k_rowscale_r(n, t_V1i, n, t_Rs(1,1))
    call k_rowscale_r(n, t_U2h, n, t_Rs(1,2))
    call k_build_H2(n, br, t_V1i, t_U2h, t_Dm1, udv_U(1,1,nst,nfe), t_Dm2, udw_V(1,1,2,nfe), t_H2)
#elif defined(STAB3)
    ! scales <= 1 stay in D1m / D2m, larger ones move (as 1/D) to the rows of V1INV / MYU2 (STAB3 branch of CGR2_2)
    call k_dm_d(n, udv_D(1,nst,nfe), t_Dm1)
    call k_dm_d(n, udw_D(1,2,nfe), t_Dm2)
    call k_rowscale_d(n, t_V1i, n, udv_D(1,nst,nfe))
    call k_rowscale_d(n, t_U2h, n, udw_D(1,2,nfe))
    call k_build_H2(n, br, t_V1i, t_U2h, t_Dm1, udv_U(1,1,nst,nfe), t_Dm2, udw_V(1,1,2,nfe), t_H2)
#else
    call k_build_H2(n, br, t_V1i, t_U2h, udv_D(1,nst,nfe), udv_U(1,1,nst,nfe), udw_D(1,2,nfe), udw_V(1,1,2,nfe), t_H2)
#endif
    call k_ct(n2, t_H2, n2, t_H1, n2)                      ! HLPB1 = conjg(transpose(HLPB2))
#if defined(STABLOG) || defined(STAB3)
    ! host ZGEQP3: D->H of HLPB1 (16 (2N)^2 B), H->D of the factors, tau and the permutation (unavoidable: no device QRCP)
    !$omp target update from(t_H1)
    t_ipvt = 0
    call la_zgeqp3(n2, n2, t_H1, n2, t_ipvt, t_tau, t_work, t_lw, t_rw, info)
    if (info /= 0) call tau_fail('CGR2_2: zgeqp3', info)
    !$omp target update to(t_H1, t_tau, t_ipvt)
#else
    ! column-norm pre-pivot + unpivoted QR
    call la_zqrpp(n2, n2, t_H1, n2, t_ipvt, t_tau, t_work, t_lw, info)
    if (info /= 0) call tau_fail('CGR2_2: la_zqrpp (-101: NaN or Inf in the matrix)', info)
#endif
    if (fault == 1) then                                    ! negative control: wrong permutation (host and device copy)
       kk = t_ipvt(1); t_ipvt(1) = t_ipvt(2); t_ipvt(2) = kk
       !$omp target update to(t_ipvt(1:2))
    endif
    call k_diagabs(n2, t_H1, n2, t_D3)                      ! separate off D = |R_ii|
    call k_divrows(n2, t_H1, n2, t_D3)
    ! solve_extended_System: rows of [[UCT,0],[0,VINV]] in the order of ZLAPMR(forward, IPVT), ZTRSM, scaling, ZUNMQR
    call k_build_rhs(n, br, t_ipvt, t_V1i, t_U2h, t_H2)
    call la_ztrsm('L', 'U', 'C', 'N', n2, n2, z1, t_H1, n2, t_H2, n2)
    call k_rowscale_inv(n2, t_H2, n2, t_D3)                 ! conj(1/D3(i)) * row i
    call la_zunmqr('L', 'N', n2, n2, n2, t_H1, n2, t_tau, t_H2, n2, t_work, t_lw, info)
    if (info /= 0) call tau_fail('CGR2_2: zunmqr', info)
    call k_blocks(n, br, fault == 4, t_H2, G00, GT0, GTT, G0T)
  end subroutine dev_cgr22

!====================================================================
!  Kernels (plain contiguous arrays, explicit extents)
!====================================================================
  Subroutine k_tau_init(n, G00, GT0, GTT, G0T)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: G00(n,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: GT0(n,n), GTT(n,n), G0T(n,n)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          GT0(i,j) = G00(i,j)
          GTT(i,j) = G00(i,j)
          if (i == j) then
             G0T(i,j) = -(cmplx(1.d0, 0.d0, kind(0.d0)) - G00(i,j))
          else
             G0T(i,j) = -(cmplx(0.d0, 0.d0, kind(0.d0)) - G00(i,j))
          endif
       enddo
    enddo
  end subroutine k_tau_init

  ! fault injection: A(1,1) *= 1 + 1e-6
  Subroutine k_fault_scale(n, A)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(n,n)
    !$omp target
    A(1,1) = A(1,1) * (1.d0 + 1.d-6)
    !$omp end target
  end subroutine k_fault_scale

#if defined(STABLOG)
  ! STABLOG: L(i) <= 0: Dm = exp(L), Rs = 1;  else Dm = 1, Rs = exp(-L)
  Subroutine k_dm_log(n, L, Dm, Rs)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(in) :: L(n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dm(n)
    Real (Kind=Kind(0.d0)), intent(inout) :: Rs(n)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       if (L(i) <= 0.d0) then
          Dm(i) = cmplx(exp(L(i)), 0.d0, kind(0.d0))
          Rs(i) = 1.d0
       else
          Dm(i) = cmplx(1.d0, 0.d0, kind(0.d0))
          Rs(i) = exp(-L(i))
       endif
    enddo
  end subroutine k_dm_log

  ! A(i,:) *= Rs(i)  (ZSCAL with a real factor; factors equal to 1 leave A unchanged)
  Subroutine k_rowscale_r(n, A, lda, Rs)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,n)
    Real (Kind=Kind(0.d0)), intent(in) :: Rs(n)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          A(i,j) = cmplx(Rs(i), 0.d0, kind(0.d0)) * A(i,j)
       enddo
    enddo
  end subroutine k_rowscale_r
#endif

#if defined(STAB3)
  ! STAB3: Dm = D if Re D <= 1, else 1
  Subroutine k_dm_d(n, D, Dm)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dm(n)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       if (real(D(i)) <= 1.d0) then
          Dm(i) = D(i)
       else
          Dm(i) = cmplx(1.d0, 0.d0, kind(0.d0))
       endif
    enddo
  end subroutine k_dm_d
#endif

  ! HLPB2 (2N x 2N) from the blocks, branch 1 (D1(1) > D2(1)): [[V1INV, D1m*U1^dagger], [-D2m*V2, MYU2]]; branch 2:
  ! [[MYU2, -D2m*V2], [D1m*U1^dagger, V1INV]]   (D1m(I) conjg(U1(J,I)) at (I,J), D2m(I) V2(I,J))
  Subroutine k_build_H2(n, br, V1i, U2h, D1m, U1, D2m, V2, H)
    Implicit none
    Integer, intent(in) :: n, br
    Complex (Kind=Kind(0.d0)), intent(in) :: V1i(n,n), U2h(n,n), D1m(n), U1(n,n), D2m(n), V2(n,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: H(2*n,2*n)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (br == 1) then
             H(i,j)     = V1i(i,j)
             H(i+n,j+n) = U2h(i,j)
             H(i,j+n)   =  D1m(i)*conjg(U1(j,i))
             H(i+n,j)   = -D2m(i)*V2(i,j)
          else
             H(i,j)     = U2h(i,j)
             H(i+n,j+n) = V1i(i,j)
             H(i,j+n)   = -D2m(i)*V2(i,j)
             H(i+n,j)   =  D1m(i)*conjg(U1(j,i))
          endif
       enddo
    enddo
  end subroutine k_build_H2

  ! right-hand side of solve_extended_System after ZLAPMR(forward, IPVT): row i = row ipvt(i) of [[UCT,0],[0,VINV]],
  ! (UCT, VINV) = (V1INV, MYU2) in branch 1, (MYU2, V1INV) in branch 2
  Subroutine k_build_rhs(n, br, ipvt, V1i, U2h, H)
    Implicit none
    Integer, intent(in) :: n, br
    Integer, intent(in) :: ipvt(2*n)
    Complex (Kind=Kind(0.d0)), intent(in) :: V1i(n,n), U2h(n,n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: H(2*n,2*n)
    Integer :: i, j, r
    !$omp target teams distribute parallel do collapse(2) private(r)
    do j = 1, 2*n
       do i = 1, 2*n
          r = ipvt(i)
          H(i,j) = cmplx(0.d0, 0.d0, kind(0.d0))
          if (r <= n .and. j <= n) then
             if (br == 1) then
                H(i,j) = V1i(r,j)
             else
                H(i,j) = U2h(r,j)
             endif
          else if (r > n .and. j > n) then
             if (br == 1) then
                H(i,j) = U2h(r-n,j-n)
             else
                H(i,j) = V1i(r-n,j-n)
             endif
          endif
       enddo
    enddo
  end subroutine k_build_rhs

  ! H(i,:) *= conj(1/D(i))   (D real valued: the factor is the real 1/Re D(i)); H is 2N x 2N, n = 2N
  Subroutine k_rowscale_inv(n, H, ldh, D)
    Implicit none
    Integer, intent(in) :: n, ldh
    Complex (Kind=Kind(0.d0)), intent(inout) :: H(ldh,n)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(n)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: sv
    !$omp target teams distribute parallel do collapse(2) private(sv)
    do j = 1, n
       do i = 1, n
          sv = 1.d0 / real(D(i))
          H(i,j) = H(i,j) * sv
       enddo
    enddo
  end subroutine k_rowscale_inv

  ! get_blocks: branch 1: GR00 = (1,1), GR0T = (1,2), GRT0 = (2,1), GRTT = (2,2); branch 2: GRTT = (1,1), GRT0 = (1,2),
  ! GR0T = (2,1), GR00 = (2,2)   (blocks of the 2N x 2N matrix).  swp: fault injection, GT0(1,1) + 1e-6.
  Subroutine k_blocks(n, br, swp, H, G00, GT0, GTT, G0T)
    Implicit none
    Integer, intent(in) :: n, br
    Logical, intent(in) :: swp
    Complex (Kind=Kind(0.d0)), intent(in) :: H(2*n,2*n)
    Complex (Kind=Kind(0.d0)), intent(inout) :: G00(n,n), GT0(n,n), GTT(n,n), G0T(n,n)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          if (br == 1) then
             G00(i,j) = H(i,j)
             G0T(i,j) = H(i,j+n)
             GT0(i,j) = H(i+n,j)
             GTT(i,j) = H(i+n,j+n)
          else
             GTT(i,j) = H(i,j)
             GT0(i,j) = H(i,j+n)
             G0T(i,j) = H(i+n,j)
             G00(i,j) = H(i+n,j+n)
          endif
          if (swp) then
             if (i == 1 .and. j == 1) then
                GT0(i,j) = GT0(i,j) + cmplx(1.d-6, 0.d0, kind(0.d0))
             endif
          endif
       enddo
    enddo
  end subroutine k_blocks

!====================================================================
!  Self-test (ALF_GPU_SELFTEST=4, run by the test executable 54-gpu-tau)
!====================================================================
!--------------------------------------------------------------------
!> @brief
!> One pass of the device tau pass against the upstream host code on the same state.
!> @details
!> Needs the device tau pass on and walker%udvst / walker%GR as set by Walker_init_stack / Walker_init_GR and handed to the
!> device by Setup_gpu_stack.  Oracle: PropR / PropRm1 (tau_m_mod), WRAPUR and CGR2_2 on host objects.  Per slice GT0, G0T, GTT
!> of the device against the host (max |A-B| / max |B| <= 1e-10); at every stabilisation point the four outputs of CGR2_2 the same
!> way, and the four precision numbers (max, mean) of the device reduction against COMPARE of the downloaded device matrices
!> (relative 1e-12: the reduction itself is tested, not the roundoff noise of the numbers).  A failure ends the run
!> (selftest_fail).  The faults ALF_GPU_FAULT_TAU = perm, prop, prec, blocks must fail here.
!--------------------------------------------------------------------
  Subroutine gpu_tau_selftest(walker)
    Implicit none
    Type (qmc_walker_t), intent(inout) :: walker
    Real (Kind=Kind(0.d0)), parameter :: tol_g = 1.d-10, tol_p = 1.d-12
    Complex (Kind=Kind(0.d0)), Allocatable :: hGT0(:,:,:), hG0T(:,:,:), hGTT(:,:,:), hG00(:,:,:), hold(:,:,:), dA(:,:), dB(:,:)
    Class (UDV_State), Allocatable :: udvr(:)
    Integer :: n, nfe, nf, NT, NT1, NTST, NST, I, J, k
    Real (Kind=Kind(0.d0)) :: r, xmax, xmean, worst_g, worst_p
    Character (len=160) :: msg

    if (.not. gpu_tau_on) call selftest_fail('tau self-test: the device tau pass is off')
    n = t_N
    prec_report = .false.
    worst_g = 0.d0; worst_p = 0.d0
    allocate(hGT0(n,n,N_FL), hG0T(n,n,N_FL), hGTT(n,n,N_FL), hG00(n,n,N_FL), hold(n,n,3), dA(n,n), dB(n,n))
    ! host: TAU_M at tau = 0
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       do J = 1, n
          do I = 1, n
             hG00(I,J,nf) = walker%GR(I,J,nf)
             hGT0(I,J,nf) = walker%GR(I,J,nf)
             hGTT(I,J,nf) = walker%GR(I,J,nf)
             if (I == J) then
                hG0T(I,J,nf) = -(z1 - walker%GR(I,J,nf))
             else
                hG0T(I,J,nf) = walker%GR(I,J,nf)
             endif
          enddo
       enddo
       ! G0T off the diagonal is -(0 - G) = +G, the diagonal -(1 - G)
    enddo
    call tau_begin(walker%GR)
    allocate(UDV_State :: udvr(t_nfe))
    do nfe = 1, t_nfe
       call udvr(nfe)%init(Ndim, 'r')
    enddo
    NST = 1
    do NT = 0, Ltrot - 1
       NT1 = NT + 1
       call PropR(hGT0, NT1)
       call PropRm1(hG0T, NT1)
       call PropRm1(hGTT, NT1)
       call PropR(hGTT, NT1)
       call tau_slice(NT1)
       do nfe = 1, t_nfe
          nf = Calc_Fl_map(nfe)
          ! D->H: device matrices of the slice, for the comparison with the host (self-test only)
          !$omp target update from(GT0_d(:,:,nf), G0T_d(:,:,nf), GTT_d(:,:,nf))
          r = max(rel_diff(GT0_d(:,:,nf), hGT0(:,:,nf)), rel_diff(G0T_d(:,:,nf), hG0T(:,:,nf)), rel_diff(GTT_d(:,:,nf), hGTT(:,:,nf)))
          worst_g = max(worst_g, r)
          if (.not. (r <= tol_g)) then
             write(msg,'(A,I0,A,I0,A,ES10.3)') 'tau pass: propagation slice ', NT1, ', flavour ', nf, ': relative difference ', r
             call selftest_fail(trim(msg))
          endif
       enddo
       if (walker%Stab_nt(NST) == NT1) then
          NTST = walker%Stab_nt(NST-1)
          call WRAPUR(NTST, NT1, udvr)
          call gpu_stab_wrapur(NTST, NT1)
          do nfe = 1, t_nfe
             nf = Calc_Fl_map(nfe)
             hold(:,:,1) = hGTT(:,:,nf); hold(:,:,2) = hGT0(:,:,nf); hold(:,:,3) = hG0T(:,:,nf)
             call CGR2_2(hGT0(:,:,nf), hG00(:,:,nf), hGTT(:,:,nf), hG0T(:,:,nf), udvr(nfe), walker%udvst(NST,nfe), Ndim)
             call tau_stab_cgr(nfe, NST)
             ! D->H: device outputs and the inputs of the precision pairs (self-test only)
             !$omp target update from(G00_d(:,:,nf), GT0_d(:,:,nf), GTT_d(:,:,nf), G0T_d(:,:,nf), t_GR(:,:,nf), t_old(:,:,:))
             r = max(rel_diff(G00_d(:,:,nf), hG00(:,:,nf)), rel_diff(GT0_d(:,:,nf), hGT0(:,:,nf)), &
                  & rel_diff(GTT_d(:,:,nf), hGTT(:,:,nf)), rel_diff(G0T_d(:,:,nf), hG0T(:,:,nf)))
             worst_g = max(worst_g, r)
             if (.not. (r <= tol_g)) then
                write(msg,'(A,I0,A,I0,A,ES10.3)') 'tau pass: CGR2_2 at stabilisation point ', NST, ', flavour ', nf, ': relative difference ', r
                call selftest_fail(trim(msg))
             endif
             ! precision numbers: device reduction against COMPARE of the same (downloaded) device matrices
             do k = 1, 4
                select case (k)
                case (1); dA = t_GR(:,:,nf);   dB = G00_d(:,:,nf)
                case (2); dA = t_old(:,:,1);   dB = GTT_d(:,:,nf)
                case (3); dA = t_old(:,:,2);   dB = GT0_d(:,:,nf)
                case default; dA = t_old(:,:,3); dB = G0T_d(:,:,nf)
                end select
                if (rel_diff(dA, dA) > 0.d0 .or. rel_diff(dB, dB) > 0.d0) call selftest_fail('tau pass: non-finite matrix in a precision pair')
                call COMPARE(dA, dB, xmax, xmean)
                r = max(abs(last_prec(1,k) - xmax)/max(xmax, tiny(1.d0)), abs(last_prec(2,k)/dble(n*n) - xmean)/max(xmean, tiny(1.d0)))
                worst_p = max(worst_p, r)
                if (.not. (r <= tol_p)) then
                   write(msg,'(A,I0,A,I0,A,ES10.3)') 'tau pass: precision reduction, pair ', k, ', flavour ', nf, ': relative difference ', r
                   call selftest_fail(trim(msg))
                endif
             enddo
          enddo
          NST = NST + 1
       endif
    enddo
    call tau_end()
    do nfe = 1, t_nfe
       call udvr(nfe)%dealloc
    enddo
    prec_report = .true.
    write(msg,'(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)') 'GPU tau self-test PASSED: largest relative difference of the matrices ', worst_g, &
         & ' (tolerance ', tol_g, '), of the precision reductions ', worst_p
    call say(trim(msg))
  end subroutine gpu_tau_selftest
End Module gpu_tau_mod
#endif

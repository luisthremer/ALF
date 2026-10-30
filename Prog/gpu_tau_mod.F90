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
!> GPU port, ticket 207 (ADR-001 section 4, ADR-002): the time-displaced pass TAU_M of the finite-temperature code on the
!> device.  Compiled only with -DALF_GPU.  Called from the guarded early branch in TAU_M (tau_m_mod) when gpu_tau_on.
!>
!> What runs on the device (Ltau = 1, no projector, no Langevin/HMC, no weight reconstruction: all refused by gpu_supported):
!>   * the four Green functions GT0, G0T, GTT, G00 (gpu_state_mod GT0_d, G0T_d, GTT_d and this module's G00_d; the arrays are the host shadows
!>     of the device arrays) stay on the device during the pass;
!>   * the per-slice propagation PROPR / PROPRM1 (ticket 204 kernels: hopping mmthr / mmthl_m1 and the fused Op_V slice
!>     kernels), in TAU_M's order: GT0 = PROPR, G0T = PROPRM1, GTT = PROPRM1 then PROPR;
!>   * at every stabilization point Stab_nt(NST): the right propagation WRAPUR + decompose of a fresh udvr state (the device slot
!>     udw(.,.,2,.) of gpu_stab_mod, ticket 206 blocks, reset at the start of the pass) and CGR2_2 (below) with the left state
!>     udvst(NST) read on the device from the slot udv(.,.,NST,.) (no D->H export of udvst any more), followed by the four
!>     Control_Precision_tau checks as device reductions (control_mod: Control_Precision_tau_vals).
!>   Observables (ham%OBSERT) and Hop_mod_Symm run on the host on the host shadows of the four matrices: INTERIM, removed by 208.
!>
!> CGR2_2 on the device (dev_cgr22), the algorithm of cgr2_2_mod (STAB1/2/3 excluded), 2N x 2N, N = Ndim:
!>   MYU2 = U2^dagger, V1INV = V1^-1 (LU: la_getrf + la_getrs on the identity; upstream INV is ZGETRF + ZGETRI), D1m, D2m = D (default)
!>   or the STABLOG splitting (D = exp(L) for L <= 0, else 1 and the rows of V1INV / MYU2 scaled by exp(-L)), the branch
!>   D1(1) > D2(1) (L1(1) > L2(1) in STABLOG: two scalars D->H), HLPB2 from the four N x N blocks, HLPB1 = HLPB2^dagger,
!>   QR: default = column-norm pre-pivot + unpivoted geqrf (la_qrpp, ADR-001 section 4, owner-approved), STABLOG = host ZGEQP3
!>   (la_geqp3_host); D separated as |R_ii|; solve_extended_System = rows of [[UCT,0],[0,VINV]] permuted by IPVT (ZLAPMR forward,
!>   built directly in permuted order), ZTRSM('L','U','C'), row scaling with conj(1/D), then Q * X: Q is formed explicitly (la_ormqr on
!>   the unit matrix), orthogonalised by one Newton-Schulz step Q <- Q (3 - Q^H Q)/2 (as in the 206 decompose) and applied by one GEMM
!>   (tau_ns; ALF_GPU_TAU_NS=0 switches back to ZUNMQR applied directly to X, as upstream).  Measured in ticket 207 (12x12, beta 4,
!>   LOG, 100 sweeps, same seeds): Precision tau mean without the step 7.99e-12, with it 6.75e-12 (CPU LOG 7.36e-12); without it the
!>   max over 400 sweeps was up to 2x the CPU LOG max (3 of 4 seed sets above 1.3x).  Replacing, one at a time, the device V1INV (LU
!>   solve), ZTRSM and ZUNMQR by their host LAPACK twins attributed the whole excess to cuSOLVER's ZUNMQR (V1INV and ZTRSM: no effect).
!>   get_blocks into G00 / G0T / GT0 / GTT.
!>
!> Transfers (every one is marked in the code; N = Ndim, F = number of flavours, L = Ltrot, S = number of stabilization points):
!>   per pass, INTERIM 207 (removed by 208):  H->D  GR, 16 N^2 B per flavour (host-authoritative G at time 0);
!>                                            D->H  GT0, G0T, GTT, 3 * 16 N^2 B per flavour and slice (for ham%OBSERT / Symm),
!>                                            D->H  G00, 16 N^2 B per flavour and stabilization point (G00 changes only there);
!>   per stabilization point and flavour, permanent: CGR2_2 default: 8 * 2N B (column norms, la_qrpp) D->H, 4 * 2N B
!>                                            (permutation) H->D, 32 B D->H (D1(1), D2(1), complex); 4 * 32 B D->H (precision numbers);
!>                                            STABLOG: 16 (2N)^2 B D->H and 16 (2N)^2 + 16 * 2N + 4 * 2N B H->D (la_geqp3_host),
!>                                            16 B D->H (L1(1), L2(1)), 4 * 32 B D->H (precision numbers);
!>   plus the WRAPUR block of gpu_stab_mod (206: 8N + 4N + 24 B per flavour and point).
!>   host-update path only (ALF_GPU_NOUPDATE=1): nsigma H->D once per pass (N_op * L * 16 B).
!>
!> Environment (ticket 212): every ALF_GPU_* variable is read through gpu_state_mod (gpu_env_get / gpu_env_is), which refuses all of the
!> variables below in a production run (no ALF_GPU_SELFTEST / _UPDTEST / _DEBUG), validates ALF_GPU_FAULT_TAU in every run and prints its
!> "GPU fault armed" line.
!> Self-test ALF_GPU_SELFTEST=4 (oracle: upstream WRAPUR / CGR2_2 and PROPR / PROPRM1 copies on host objects; the default path
!> needs a build with -DQR_PREPIVOT), timing ALF_GPU_TAUTIME=1, ALF_GPU_EXIT=1 ends the run after the hook;
!> ALF_GPU_FAULT_TAU=perm|prop|prec|blocks injects faults (negative controls; they act wherever the faulted step runs), and the graded
!> fault ALF_GPU_FAULT_TAU=eps:<x> multiplies, in every device-vs-host comparison of GT0 in the SELFTEST=4 self-test (pass and CGR2_2
!> alone), the element of largest |GT0| of the host copy of the device GT0 by (1 + x) (the device data and the run are unchanged);
!> ALF_GPU_TAUTEST_V=1 prints the per-slice differences of the self-test.  DEBUG AIDS (not used in production, no effect when unset): ALF_GPU_TAUDIAG=1|2 prints one "TAUDIAG" line per
!> stabilization point and flavour (device precision numbers next to a full host tau pass run in lockstep on the same trajectory, host
!> CGR2_2 on the device inputs, GTT chains with 1-ulp perturbed right states; 1 adds host twins of single CGR2_2 steps; slow, see
!> diag_post and docs/gpu/scripts/nividia/t207_diagtab.py / t207_diagstat.py), ALF_GPU_TAUDIAG_FROM=<pass> starts the perturbation
!> chains late, ALF_GPU_TAUSWAP=<mask> replaces CGR2_2 steps by host twins in the production calls (bits as swp below).
!> Finding (ticket 207, STABLOG precision tail): the tail of Precision tau is the GTT pair at the last stabilization points; its size
!> per event changes by up to ~10x under a 1-ulp perturbation of U of the right state, on the host as on the device; device / host
!> per event is unbiased with the same spread as two 1-ulp perturbed host runs (conditioning, not a device defect).
!> Host copies of device matrices go through upd_from_c / upd_to_c (dummy-argument views, as alf_la): a direct
!> "target update from(GT0_d(:,:,nf))" of a module allocatable in a small routine returned stale host data with nvfortran 26.9 -fast
!> (self-test: device GT0 unchanged after a propagation slice, no error), the same update through a dummy argument works.
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_tau_mod

  Use omp_lib
  Use iso_c_binding, only: c_int
  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Operator_mod
  Use Fields_mod
  Use Hamiltonian_main
  Use QMC_runtime_var
  Use Hop_mod
  Use MyMats, only: COMPARE, INV
  Use UDV_State_mod
  Use wrapur_mod, only: WRAPUR
  Use cgr2_2_mod, only: CGR2_2
  Use Control, only: Control_Precision_tau_vals
  Use alf_la_mod
  Use gpu_state_mod
  Use gpu_wrap_mod, only: gpu_mmthr, gpu_mmthl_m1, gpu_Vslice_R, gpu_Vslice_L
  Use gpu_stab_mod

  Implicit none

  Private
  Public :: gpu_tau_on, gpu_tau_setup, gpu_tau_finalize, gpu_tau_hook, gpu_tau_pass

  Logical, Save :: gpu_tau_on = .false.

  Integer, Save :: t_N = 0, t_N2 = 0, t_nfe = 0, t_nstm = 0
  Integer, Allocatable, Save :: t_nt(:)                       ! Stab_nt(0:NSTM)
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

  ! ---- own device buffers (allocated once in gpu_tau_setup) ----
  ! G00 at tau = 0 / after each CGR2_2 (N x N x N_FL, device + host shadow).  Here and not in gpu_state_mod next to GT0_d, G0T_d,
  ! GTT_d: with a fourth such array in gpu_state_mod, nvfortran 26.9 (-fast, STABLOG) crashed compiling gpu_stab_mod (fort1 heap
  ! corruption, signal 6; reproducible; -O0 compiles)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: G00_d(:,:,:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_H1(:,:), t_H2(:,:)          ! 2N x 2N: QR matrix / solution
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_V1i(:,:), t_U2h(:,:), t_W(:,:)   ! N x N: V1^-1, U2^dagger, LU work
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_Q(:,:), t_A(:,:), t_B(:,:)   ! 2N x 2N: explicit Q, Q^H Q, orthogonalised Q (tau_ns)
  Logical, Save :: tau_ns = .true.                                            ! Newton-Schulz corrected explicit Q in CGR2_2 (ALF_GPU_TAU_NS=0 switches it off)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_old(:,:,:)                 ! N x N x 3: GTT, GT0, G0T before CGR2_2
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_GR(:,:,:)                  ! N x N x N_FL: GR (precision reference, device copy)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_S(:,:,:)                   ! N x N x 4: scratch outputs of CGR2_2 (tests only)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: t_tau(:), t_D3(:), t_Dm1(:), t_Dm2(:)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: t_Rs(:,:)                    ! N x 2 row scale factors (STABLOG)
  Integer (c_int), Allocatable, Save :: t_ipvt(:), t_ipiv(:)
  ! ---- host side ----
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: G00_T(:,:,:), G0T_T(:,:,:), GT0_T(:,:,:), GTT_T(:,:,:)   ! Symm only

  ! ---- counters / test hooks ----
  Integer, Save :: fault = 0                                  ! 1 perm, 2 prop, 3 prec, 4 blocks
  Logical, Save :: prec_report = .true.                       ! .false. in the tests: the Control counters are not touched
  Real (Kind=Kind(0.d0)), Save :: last_prec(2,4) = 0.d0       ! last reduction: (max, sum) of the four pairs of tau_stab_cgr
  Integer(int64), Save :: n_pass = 0, n_branch(2) = 0, n_cgr = 0
  Integer(int64), Save :: b_h2d = 0, b_d2h = 0                ! interim bytes (this module's per-pass copies)
  ! ---- DEBUG AID (ALF_GPU_TAUDIAG=1, not used in production): per stabilization point and flavour, the device CGR2_2 against host
  ! chains on the same inputs; swp: host twins of single steps of dev_cgr22 (bit 1 V1INV, 2 ZTRSM, 4 Q by host ZUNMQR, 8 Q by device
  ! ZUNMQR without the Newton-Schulz step); ALF_GPU_TAUSWAP=<mask> sets the mask of the production calls (isolation experiments)
  Logical, Save :: diag_on = .false., diag_full = .true.        ! ALF_GPU_TAUDIAG=1: all variants; =2: no step twins (E) and no hD
  Integer, Parameter :: npc = 5                                     ! GTT chains (see diag_post)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: dg_p(:,:,:,:)      ! (N,N,N_FL,0:npc-1), propagated on the host
  Integer(int64), Save :: dg_seed = 20261029_int64
  Integer, Save :: dg_from = 0                                      ! ALF_GPU_TAUDIAG_FROM: perturbation chains from this pass on
  Class (UDV_State), Allocatable, Save :: dg_hp(:)
  Integer, Save :: swp = 0, swp_prod = 0
  Class (UDV_State), Allocatable, Save :: dg_hr(:), dg_hrd(:), dg_hw(:), dg_hl(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: dh00(:,:,:), dhT0(:,:,:), dhTT(:,:,:), dh0T(:,:,:), dg_c(:,:,:), dg_o(:,:,:), dg_ho(:,:,:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: dg_hwk(:)
  Integer(int64), Save :: dg_npt = 0
  Integer, Save :: last_br = 0

Contains

!====================================================================
!  Setup
!====================================================================
  Subroutine gpu_tau_setup(Stab_nt, NSTM)
    Implicit none
    Integer, intent(in) :: NSTM
    Integer, intent(in) :: Stab_nt(0:NSTM)
    Integer :: n, n2, lw
    Character (len=256) :: line

    gpu_tau_on = .false.
    if (.not. gpu_stab_on) return
    if (get_Ltau() /= 1) return
    if (Projector .or. reconstruction_needed) return
    t_N = Ndim; t_N2 = 2*Ndim; t_nfe = N_FL_eff; t_nstm = NSTM
    n = t_N; n2 = t_N2
    allocate(t_nt(0:NSTM)); t_nt = Stab_nt
    allocate(t_H1(n2,n2), t_H2(n2,n2), t_V1i(n,n), t_U2h(n,n), t_W(n,n), t_old(n,n,3), t_GR(n,n,N_FL), t_S(n,n,4))
    allocate(t_tau(n2), t_D3(n2), t_Dm1(n), t_Dm2(n), t_Rs(n,2), t_ipvt(n2), t_ipiv(n), G00_d(n,n,N_FL))
    if (gpu_env_is('ALF_GPU_TAU_NS', '0')) tau_ns = .false.        ! debug runs only (ticket 212 refuses it in production)
    if (tau_ns) then
       allocate(t_Q(n2,n2), t_A(n2,n2), t_B(n2,n2))
       !$omp target enter data map(alloc: t_Q, t_A, t_B)
    endif
    t_ipvt = 1; t_ipiv = 1
    if (Symm) allocate(G00_T(n,n,N_FL), G0T_T(n,n,N_FL), GT0_T(n,n,N_FL), GTT_T(n,n,N_FL))
    call diag_setup()
    ! device allocation once, no transfer; the first use of each buffer is a device write
    !$omp target enter data map(alloc: t_H1, t_H2, t_V1i, t_U2h, t_W, t_old, t_GR, t_S, t_tau, t_D3, t_Dm1, t_Dm2, t_Rs, t_ipvt, t_ipiv, G00_d)
    ! pre-grow the library work spaces for the 2N x 2N factorization and run the 2N QR once, so that nothing is allocated in the sweep
    call k_ident(n2, t_H1, n2)
    lw = la_geqrf_lwork(n2, n2, t_H1, n2)
    lw = max(lw, la_ormqr_lwork('L', LA_N, n2, n2, n2, t_H1, n2, t_tau, t_H2, n2))
    lw = max(lw, la_getrf_lwork(n, n, t_W, n))
    call la_reserve_work(lw)
    call la_qrpp(n2, n2, t_H1, n2, t_ipvt, t_tau)
    call la_sync()
    gpu_tau_on = .true.
    gpu_stab_tau_dev = .true.
    write(line,'(A,I0,A,I0,A,I0)') 'GPU tau pass (207): ACTIVE (propagation, WRAPUR, CGR2_2 and precision check on the device; observables on the host, interim), N = ', &
         & n, ', Nstm = ', NSTM, ', flavours ', t_nfe
#if defined(STABLOG)
    line = trim(line)//', CGR2_2 QR: host ZGEQP3 (STABLOG)'
#else
    line = trim(line)//', CGR2_2 QR: pre-pivoted (default)'
#endif
    call info_line(line)
  end subroutine gpu_tau_setup

  Subroutine gpu_tau_finalize()
    Implicit none
    Character (len=256) :: line
    if (.not. gpu_tau_on) return
    write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'GPU tau pass counters: passes ', n_pass, ', CGR2_2 calls ', n_cgr, &
         & ' (branch D1(1) > D2(1): ', n_branch(1), ', else: ', n_branch(2), '), interim bytes H->D ', b_h2d, ' and D->H ', b_d2h, ' (INTERIM 207)'
    call info_line(line)
    if (tau_ns) then
       !$omp target exit data map(delete: t_Q, t_A, t_B)
       deallocate(t_Q, t_A, t_B)
    endif
    !$omp target exit data map(delete: t_H1, t_H2, t_V1i, t_U2h, t_W, t_old, t_GR, t_S, t_tau, t_D3, t_Dm1, t_Dm2, t_Rs, t_ipvt, t_ipiv, G00_d)
    deallocate(t_H1, t_H2, t_V1i, t_U2h, t_W, t_old, t_GR, t_S, t_tau, t_D3, t_Dm1, t_Dm2, t_Rs, t_ipvt, t_ipiv, t_nt, G00_d)
    if (allocated(G00_T)) deallocate(G00_T, G0T_T, GT0_T, GTT_T)
    gpu_tau_on = .false.
    gpu_stab_tau_dev = .false.
  end subroutine gpu_tau_finalize

!====================================================================
!  The pass (what TAU_M does, without the Langevin / HMC / reconstruction parts, which gpu_supported refuses)
!====================================================================
  Subroutine gpu_tau_pass(GR, PHASE, NSTM, STAB_NT)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL), Phase
    Integer, intent(in) :: NSTM
    Integer, intent(in) :: STAB_NT(0:NSTM)
    Integer :: nfe, NT, NT1, NTST, NST
    Real (Kind=Kind(0.d0)) :: Mc_step_Weight

    Mc_step_Weight = 1.d0
    call tau_begin(GR)
    if (diag_on) call diag_begin(GR)
    call tau_obs(0, PHASE, Mc_step_Weight)
    NST = 1
    do NT = 0, LTROT - 1
       NT1 = NT + 1
       call tau_slice(NT1)
       if (diag_on) call diag_slice(NT1)
       call tau_download(.false.)
       call tau_obs(NT1, PHASE, Mc_step_Weight)
       if (STAB_NT(NST) == NT1) then
          NTST = STAB_NT(NST-1)
          if (diag_on) call diag_prewrap(NTST, NT1)
          call gpu_stab_wrapur(NTST, NT1)             ! WRAPUR(NTST, NT1, udvr) of TAU_M: device slot udvr, decompose, scale test
          do nfe = 1, t_nfe
             if (diag_on) call diag_pre(nfe, NST)
             call tau_stab_cgr(nfe, NST)              ! CGR2_2 + the four precision checks
             if (diag_on) call diag_post(nfe, NST, GR)
          enddo
          NST = NST + 1
       endif
    enddo
    call tau_end()
    n_pass = n_pass + 1
  end subroutine gpu_tau_pass

  !> Observables of one time slice on the host shadows of the device matrices (INTERIM 207: removed by 208).
  Subroutine tau_obs(nt, PHASE, Mc_step_Weight)
    Implicit none
    Integer, intent(in) :: nt
    Complex (Kind=Kind(0.d0)), intent(in) :: Phase
    Real (Kind=Kind(0.d0)), intent(in) :: Mc_step_Weight
    if (nt == 0) then
       if (Symm) then
          Call Hop_mod_Symm(G00_T,G00_d,nt)
          Call Hop_mod_Symm(GTT_T,GTT_d,nt)
          Call Hop_mod_Symm(G0T_T,G0T_d,nt)
          Call Hop_mod_Symm(GT0_T,GT0_d,nt)
          CALL ham%OBSERT(nt, GT0_T,G0T_T,G00_T,GTT_T, PHASE, Mc_step_Weight)
       else
          CALL ham%OBSERT(nt, GT0_d,G0T_d,G00_d,GTT_d, PHASE, Mc_step_Weight)
       endif
    else
       if (Symm) then
          Call Hop_mod_Symm(G00_T,G00_d,0)
          Call Hop_mod_Symm(GTT_T,GTT_d,nt)
          Call Hop_mod_Symm(G0T_T,G0T_d,Ltrot,nt)
          Call Hop_mod_Symm(GT0_T,GT0_d,nt,Ltrot)
          CALL ham%OBSERT(nt, GT0_T,G0T_T,G00_T,GTT_T,PHASE, Mc_step_Weight)
       else
          CALL ham%OBSERT(nt, GT0_d,G0T_d,G00_d,GTT_d,PHASE, Mc_step_Weight)
       endif
    endif
  end subroutine tau_obs

  !> Start of the pass: G00 = GT0 = GTT = GR, G0T = -(1 - GR) (TAU_M, tau = 0).  The host shadows get the CPU's loop (they feed
  !> OBSERT at tau = 0); the device gets GR by one H->D copy per flavour and derives the other three on the device.
  Subroutine tau_begin(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL)
    Integer :: nfe, nf, I, J
    Complex (Kind=Kind(0.d0)) :: Z
    call sync_fields(1, Ltrot)           ! host-update path only (ALF_GPU_NOUPDATE=1): nsigma%f -> nsigma_d, INTERIM 206; no-op with device updates
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       do J = 1, t_N
          do I = 1, t_N
             Z = cmplx(0.d0, 0.d0, kind(0.D0))
             if (I == J) Z = cmplx(1.d0, 0.d0, kind(0.d0))
             G00_d(I,J,nf) = GR(I,J,nf)
             GT0_d(I,J,nf) = GR(I,J,nf)
             GTT_d(I,J,nf) = GR(I,J,nf)
             G0T_d(I,J,nf) = -(Z - GR(I,J,nf))
          enddo
       enddo
       ! INTERIM 207: removed by 208: H->D of GR (16 N^2 B per flavour): GR is host-authoritative until the measurements move to the device
       call upd_to_c(G00_d(1,1,nf), t_N*t_N)
       b_h2d = b_h2d + 16_int64*t_N*t_N
       call k_tau_init(t_N, G00_d(1,1,nf), GT0_d(1,1,nf), GTT_d(1,1,nf), G0T_d(1,1,nf))
       call k_copy_c(t_N*t_N, G00_d(1,1,nf), t_GR(1,1,nf))
       call gpu_stab_reset_r(nfe)         ! fresh udvr (slot 2): the right propagation of TAU_M starts from the unit state
    enddo
  end subroutine tau_begin

  Subroutine tau_end()
    Implicit none
    integer :: nfe
    do nfe = 1, t_nfe
       call gpu_stab_reset_r(nfe)         ! leave slot 2 clean (main resets it before the next use anyway)
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
       call gpu_mmthr(GT0_d(1,1,nf), N, N, nf)
       call gpu_Vslice_R(GT0_d(1,1,nf), N, N, nf, nt1, 'n', .false.)
       ! PROPRM1: B_T^-1 from the right, then Op_mmultL('n', sign = -1) for n = 1..N_op
       call gpu_mmthl_m1(G0T_d(1,1,nf), N, N, nf)
       call gpu_Vslice_L(G0T_d(1,1,nf), N, N, nf, nt1, 'n', -1, .false.)
       call gpu_mmthl_m1(GTT_d(1,1,nf), N, N, nf)
       call gpu_Vslice_L(GTT_d(1,1,nf), N, N, nf, nt1, 'n', -1, .false.)
       call gpu_mmthr(GTT_d(1,1,nf), N, N, nf)
       call gpu_Vslice_R(GTT_d(1,1,nf), N, N, nf, nt1, 'n', .false.)
       if (fault == 2 .and. nt1 == 2) call k_fault_scale(GT0_d(1,1,nf))    ! negative control: wrong propagation
    enddo
  end subroutine tau_slice

  !> D->H of the matrices the host observables read: GT0, G0T, GTT (incl00 = .false.: every slice) and G00 (incl00 = .true.).
  Subroutine tau_download(incl00)
    Implicit none
    Logical, intent(in) :: incl00
    Integer :: nfe, nf
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       ! INTERIM 207: removed by 208: D->H of GT0, G0T, GTT (3 * 16 N^2 B per flavour and slice) for ham%OBSERT / Hop_mod_Symm on the host
       call upd_from_c(GT0_d(1,1,nf), t_N*t_N)
       call upd_from_c(G0T_d(1,1,nf), t_N*t_N)
       call upd_from_c(GTT_d(1,1,nf), t_N*t_N)
       b_d2h = b_d2h + 3_int64*16_int64*t_N*t_N
       if (incl00) then
          ! INTERIM 207: removed by 208: D->H of G00 (16 N^2 B per flavour), only needed after CGR2_2 (tests)
          call upd_from_c(G00_d(1,1,nf), t_N*t_N)
          b_d2h = b_d2h + 16_int64*t_N*t_N
       endif
    enddo
  end subroutine tau_download

  !> H->D of nel contiguous complex elements starting at A (same reason as upd_from_c)
  Subroutine upd_to_c(A, nel)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)) :: A(*)
    !$omp target update to(A(1:nel))
  end subroutine upd_to_c

  !> H->D of nel contiguous integer elements starting at A (same reason as upd_from_c)
  Subroutine upd_to_i(A, nel)
    Implicit none
    Integer, intent(in) :: nel
    Integer (c_int) :: A(*)
    !$omp target update to(A(1:nel))
  end subroutine upd_to_i

  !> D->H of nel contiguous real elements starting at A (same reason as upd_from_c)
  Subroutine upd_from_r(A, nel)
    Implicit none
    Integer, intent(in) :: nel
    Real (Kind=Kind(0.d0)) :: A(*)
    !$omp target update from(A(1:nel))
  end subroutine upd_from_r

  !> D->H of nel contiguous complex elements starting at A (a flat view, as in alf_la: a section of a mapped array).  A dummy
  !> argument, not a directive on the module array: nvfortran 26.9 -fast gave stale host data after a direct
  !> "target update from(GT0_d(:,:,nf))" in a small routine (see DISCREPANCIES.md).
  Subroutine upd_from_c(A, nel)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)) :: A(*)
    !$omp target update from(A(1:nel))
  end subroutine upd_from_c

!====================================================================
!  Stabilization point: CGR2_2 and the precision checks (per flavour)
!====================================================================
  Subroutine tau_stab_cgr(nfe, nst)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Integer :: nf, n, br

    n = t_N
    nf = Calc_Fl_map(nfe)
    ! HLP4/5/6 of TAU_M: the propagated matrices before CGR2_2 (device copies)
    call k_copy_c(n*n, GTT_d(1,1,nf), t_old(1,1,1))
    call k_copy_c(n*n, GT0_d(1,1,nf), t_old(1,1,2))
    call k_copy_c(n*n, G0T_d(1,1,nf), t_old(1,1,3))
    call dev_cgr22(nfe, nst, 1.d0, 0.d0, br, G00_d(1,1,nf), GT0_d(1,1,nf), GTT_d(1,1,nf), G0T_d(1,1,nf))
    n_cgr = n_cgr + 1
    n_branch(br) = n_branch(br) + 1
    ! Control_Precision_tau(GR, G00), (HLP4, GTT), (HLP5, GT0), (HLP6, G0T)
    call prec_pair(1, t_GR(1,1,nf), G00_d(1,1,nf))
    call prec_pair(2, t_old(1,1,1), GTT_d(1,1,nf))
    call prec_pair(3, t_old(1,1,2), GT0_d(1,1,nf))
    call prec_pair(4, t_old(1,1,3), G0T_d(1,1,nf))
    ! INTERIM 207: removed by 208: D->H of the new G00 (16 N^2 B), the host observables read it
    call upd_from_c(G00_d(1,1,nf), n*n)
    b_d2h = b_d2h + 16_int64*n*n
  end subroutine tau_stab_cgr

  !> Control_Precision_tau(A, B) as a device reduction (dev_prec of gpu_stab_mod: one D->H of four numbers, 32 B)
  Subroutine prec_pair(k, A, B)
    Implicit none
    Integer, intent(in) :: k
    Complex (Kind=Kind(0.d0)), intent(in) :: A(t_N,*), B(t_N,*)
    Real (Kind=Kind(0.d0)) :: xmax, xsum
    Integer :: nna, nnb
    call dev_prec(t_N, A, B, xmax, xsum, nna, nnb)
    if (fault == 3) xmax = xmax*(1.d0 + 1.d-6)            ! negative control: wrong reduction
    last_prec(1,k) = xmax; last_prec(2,k) = xsum
    if (prec_report) call Control_Precision_tau_vals(xmax, xsum, nna, nnb, t_N)
  end subroutine prec_pair

!====================================================================
!  CGR2_2 on the device
!====================================================================
!> CGR2_2(GRT0, GR00, GRTT, GR0T, udv2 = udvr (slot 2), udv1 = udvst(nst), LQ) of cgr2_2_mod for flavour nfe.  Outputs: N x N
!> device arrays.  s1 (default) / dl1 (STABLOG): test-only synthetic scaling of the left state, D1 -> s1 * D1 resp. L1 -> L1 + dl1
!> (1 and 0 in production), used to reach the other branch.  br = 1: D1(1) > D2(1) (L1(1) > L2(1)), 2: else.
  Subroutine dev_cgr22(nfe, nst, s1, dl1, br, G00, GT0, GTT, G0T)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Real (Kind=Kind(0.d0)), intent(in) :: s1, dl1
    Integer, intent(out) :: br
    Complex (Kind=Kind(0.d0)), intent(inout) :: G00(t_N,*), GT0(t_N,*), GTT(t_N,*), G0T(t_N,*)
    Integer :: n, n2
    Real (Kind=Kind(0.d0)) :: d1, d2
    Integer(c_int) :: kk

    n = t_N; n2 = t_N2
    ! the branch decision needs D1(1) and D2(1) (L1(1), L2(1)) on the host
#if defined(STABLOG)
    ! D->H: L1(1), L2(1), 16 B (branch of CGR2_2)
    call upd_from_r(udv_L(1,nst,nfe), 1)
    call upd_from_r(udw_L(1,2,nfe), 1)
    d1 = udv_L(1,nst,nfe) + dl1; d2 = udw_L(1,2,nfe)
#else
    ! D->H: D1(1), D2(1), 32 B (branch of CGR2_2)
    call upd_from_c(udv_D(1,nst,nfe), 1)
    call upd_from_c(udw_D(1,2,nfe), 1)
    d1 = real(udv_D(1,nst,nfe))*s1; d2 = real(udw_D(1,2,nfe))
#endif
    br = 2
    if (d1 > d2) br = 1
    last_br = br
    ! MYU2 = U2^dagger
    call k_ct(n, udw_U(1,1,2,nfe), n, t_U2h, n)
    ! V1INV = V1^-1: LU and a solve with the identity (upstream: ZGETRF + ZGETRI)
    call k_copy_c(n*n, udv_V(1,1,nst,nfe), t_W)
    call la_getrf(n, n, t_W, n, t_ipiv)
    call k_ident(n, t_V1i, n)
    call la_getrs(LA_N, n, n, t_W, n, t_ipiv, t_V1i, n)
    if (iand(swp, 1) /= 0) call sw_v1inv(nfe, nst)          ! DEBUG AID (diagnostic only): host ZGETRF + ZGETRI
#if defined(STABLOG)
    ! scales smaller than 1 stay in D1m / D2m, larger ones move to the rows of V1INV / MYU2 (STAB3 / STABLOG branch of CGR2_2)
    call k_dm_log(n, udv_L(1,nst,nfe), dl1, t_Dm1, t_Rs(1,1))
    call k_dm_log(n, udw_L(1,2,nfe), 0.d0, t_Dm2, t_Rs(1,2))
    call k_rowscale_r(n, t_V1i, n, t_Rs(1,1))
    call k_rowscale_r(n, t_U2h, n, t_Rs(1,2))
#else
    call k_dm(n, udv_D(1,nst,nfe), s1, t_Dm1)
    call k_dm(n, udw_D(1,2,nfe), 1.d0, t_Dm2)
#endif
    call k_build_H2(n, br, t_V1i, t_U2h, t_Dm1, udv_U(1,1,nst,nfe), t_Dm2, udw_V(1,1,2,nfe), t_H2)
    call k_ct(n2, t_H2, n2, t_H1, n2)                      ! HLPB1 = conjg(transpose(HLPB2))
#if defined(STABLOG)
    ! host ZGEQP3 (ADR-001 section 4): D->H of HLPB1 (16 (2N)^2 B), H->D of the factors, tau and the permutation (inside la_geqp3_host)
    call la_geqp3_host(n2, n2, t_H1, n2, t_ipvt, t_tau)
#else
    ! column-norm pre-pivot + unpivoted QR: 8 * 2N B D->H (norms), 4 * 2N B H->D (permutation) inside la_qrpp
    call la_qrpp(n2, n2, t_H1, n2, t_ipvt, t_tau)
#endif
    if (fault == 1) then                                    ! negative control: wrong permutation (host and device copy)
       kk = t_ipvt(1); t_ipvt(1) = t_ipvt(2); t_ipvt(2) = kk
       call upd_to_i(t_ipvt, 2)
    endif
    call k_diagabs(n2, t_H1, n2, t_D3)                      ! separate off D = |R_ii|
    call k_divrows(n2, t_H1, n2, t_D3)
    ! solve_extended_System: rows of [[UCT,0],[0,VINV]] in the order of ZLAPMR(forward, IPVT), ZTRSM, scaling, ZUNMQR
    call k_build_rhs(n, br, t_ipvt, t_V1i, t_U2h, t_H2)
    if (iand(swp, 2) /= 0) then
       call sw_trsm()                                       ! DEBUG AID (diagnostic only): host ZTRSM
    else
       call la_trsm('L', 'U', LA_C, 'N', n2, n2, z1, t_H1, n2, t_H2, n2)
    endif
    call k_rowscale_inv(n2, t_H2, n2, t_D3)                 ! conj(1/D3(i)) * row i
    if (iand(swp, 4) /= 0) then
       call sw_unmqr()                                      ! DEBUG AID (diagnostic only): host ZUNMQR applied to X
       call k_blocks(n, br, fault == 4, t_H2, G00, GT0, GTT, G0T)
    else if (iand(swp, 8) /= 0) then
       call la_ormqr('L', LA_N, n2, n2, n2, t_H1, n2, t_tau, t_H2, n2)   ! DEBUG AID: device ZUNMQR without the Newton-Schulz step
       call k_blocks(n, br, fault == 4, t_H2, G00, GT0, GTT, G0T)
    else if (tau_ns) then
       ! Q formed explicitly (ormqr on the unit matrix), one Newton-Schulz step Q <- Q (3 - Q^H Q)/2, then Q * X by one GEMM: cuSOLVER's
       ! ZUNMQR applied directly to X gave a Precision tau about 17 % (mean) and up to 2x (max) above the CPU's (ticket 207; same cause as
       ! the Q orthogonality defect found in 206, docs/gpu/DISCREPANCIES.md)
       call k_ident(n2, t_Q, n2)
       call la_ormqr('L', LA_N, n2, n2, n2, t_H1, n2, t_tau, t_Q, n2)
       call la_gemm(LA_C, LA_N, n2, n2, n2, z1, t_Q, n2, t_Q, n2, z0, t_A, n2)
       call k_ns2(n2, t_A)
       call la_gemm(LA_N, LA_N, n2, n2, n2, z1, t_Q, n2, t_A, n2, z0, t_B, n2)
       call la_gemm(LA_N, LA_N, n2, n2, n2, z1, t_B, n2, t_H2, n2, z0, t_A, n2)
       call k_blocks(n, br, fault == 4, t_A, G00, GT0, GTT, G0T)
    else
       call la_ormqr('L', LA_N, n2, n2, n2, t_H1, n2, t_tau, t_H2, n2)
       call k_blocks(n, br, fault == 4, t_H2, G00, GT0, GTT, G0T)
    endif
  end subroutine dev_cgr22

!====================================================================
!  DEBUG AID, not used in production: ALF_GPU_TAUDIAG=1 (ticket 207, STABLOG precision tail).  Per stabilization point and
!  flavour one line "TAUDIAG": the device precision numbers of the production CGR2_2, the same numbers of a full host tau pass
!  run in lockstep (host PROPR / PROPRM1 / WRAPUR / CGR2_2 from the same GR and the same device left states: what the host TAU_M
!  computes), of the host CGR2_2 on the device inputs (device udvr, device propagated matrices), of the host CGR2_2 on the host
!  one-block WRAPUR of the device udvr, and of dev_cgr22 with single steps replaced by their host twins (swp masks).  All the
!  transfers below are diagnostic D->H / H->D copies (host oracle on device data); none of them feeds back into the run.
!====================================================================
  Subroutine diag_setup()
    Implicit none
    Integer :: nfe, n, n2, info
    Complex (Kind=Kind(0.d0)) :: wq(1)
    Character (len=256) :: v
    Integer :: ios
    ! debug variables, read through the ticket 212 reader (refused in production runs); an unparsable value stops the run
    if (gpu_env_get('ALF_GPU_TAUSWAP', v)) then
       read(v,*,iostat=ios) swp_prod
       if (ios /= 0 .or. swp_prod < 0 .or. swp_prod > 15) call env_bad('ALF_GPU_TAUSWAP', v)
    endif
    swp = swp_prod
    if (gpu_env_get('ALF_GPU_TAUDIAG', v)) then
       if (trim(v) /= '1' .and. trim(v) /= '2' .and. trim(v) /= '0') call env_bad('ALF_GPU_TAUDIAG', v)
    endif
    diag_on = gpu_env_is('ALF_GPU_TAUDIAG', '1') .or. gpu_env_is('ALF_GPU_TAUDIAG', '2')
    diag_full = gpu_env_is('ALF_GPU_TAUDIAG', '1')
    if (gpu_env_get('ALF_GPU_TAUDIAG_FROM', v)) then
       read(v,*,iostat=ios) dg_from
       if (ios /= 0 .or. dg_from < 0) call env_bad('ALF_GPU_TAUDIAG_FROM', v)
    endif
    if (.not. (diag_on .or. swp_prod /= 0)) return
    n = t_N; n2 = t_N2
    call ZUNMQR('L', 'N', n2, n2, n2, t_H1, n2, t_tau, t_H2, n2, wq, -1, info)
    allocate(dg_hwk(max(1, int(dble(wq(1))))))
    if (swp_prod /= 0) call info_line('GPU tau pass: DEBUG ALF_GPU_TAUSWAP active (host twins of CGR2_2 steps)')
    if (.not. diag_on) return
    allocate(dg_hr(t_nfe), dg_hrd(t_nfe), dg_hw(t_nfe), dg_hl(t_nfe), dg_hp(1))
    call dg_hp(1)%init(n, 'r')
    do nfe = 1, t_nfe
       call dg_hr(nfe)%init(n, 'r'); call dg_hrd(nfe)%init(n, 'r'); call dg_hw(nfe)%init(n, 'r')
       call dg_hl(nfe)%alloc(n); dg_hl(nfe)%side = 'l'
    enddo
    allocate(dh00(n,n,N_FL), dhT0(n,n,N_FL), dhTT(n,n,N_FL), dh0T(n,n,N_FL), dg_c(n,n,4), dg_o(n,n,3), dg_ho(n,n,3))
    dh00 = z0; dhT0 = z0; dhTT = z0; dh0T = z0
    allocate(dg_p(n,n,N_FL,0:npc-1)); dg_p = z0
    call info_line('GPU tau pass: DEBUG ALF_GPU_TAUDIAG=1 active (host oracle chains at every stabilization point, slow)')
  end subroutine diag_setup

  !> unparsable value of a debug variable of this module: stop (stdout and info)
  Subroutine env_bad(name, v)
    Implicit none
    Character (len=*), intent(in) :: name, v
    call say('GPU debug value invalid: '//trim(name)//'='//trim(v))
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine env_bad

  Subroutine diag_begin(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL)
    Integer :: nfe
    call h_begin(GR, dh00, dhT0, dhTT, dh0T)
    do nfe = 1, t_nfe
       call dg_hr(nfe)%reset('r')
       call diag_pset(Calc_Fl_map(nfe), GR(:,:,Calc_Fl_map(nfe)))
    enddo
  end subroutine diag_begin

  Subroutine diag_slice(nt1)
    Implicit none
    Integer, intent(in) :: nt1
    Integer :: k
    call h_propr(dhT0, nt1)
    call h_proprm1(dh0T, nt1)
    call h_proprm1(dhTT, nt1)
    call h_propr(dhTT, nt1)
    do k = 0, npc-1
       if (k > 0 .and. n_pass+1 < dg_from) exit
       call h_proprm1(dg_p(:,:,:,k), nt1)
       call h_propr(dg_p(:,:,:,k), nt1)
    enddo
  end subroutine diag_slice

  ! GTT chains restart from G (tau = 0)
  Subroutine diag_pset(nf, G)
    Implicit none
    Integer, intent(in) :: nf
    Complex (Kind=Kind(0.d0)), intent(in) :: G(:,:)
    Integer :: k
    do k = 0, npc-1
       dg_p(:,:,nf,k) = G
    enddo
  end subroutine diag_pset

  ! before the device WRAPUR: host WRAPUR of the full host chain, and a host one-block WRAPUR of the current device udvr
  Subroutine diag_prewrap(ntst, nt1)
    Implicit none
    Integer, intent(in) :: ntst, nt1
    Integer :: nfe
    call la_sync()
    do nfe = 1, t_nfe
       call dl_state(1, 2, nfe, dg_hw(nfe))        ! diagnostic D->H: device udvr before this block
    enddo
    call WRAPUR(ntst, nt1, dg_hw)
    call WRAPUR(ntst, nt1, dg_hr)
  end subroutine diag_prewrap

  Subroutine diag_pre(nfe, nst)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Integer :: nf, n
    n = t_N; nf = Calc_Fl_map(nfe)
    call la_sync()
    call dl_state(1, 2, nfe, dg_hrd(nfe))           ! diagnostic D->H: device udvr after the device WRAPUR
    call dl_state(2, nst, nfe, dg_hl(nfe))          ! diagnostic D->H: left state udvst(nst)
    call upd_from_c(GTT_d(1,1,nf), n*n); call upd_from_c(GT0_d(1,1,nf), n*n); call upd_from_c(G0T_d(1,1,nf), n*n)
    dg_o(:,:,1) = GTT_d(:,:,nf); dg_o(:,:,2) = GT0_d(:,:,nf); dg_o(:,:,3) = G0T_d(:,:,nf)
    dg_ho(:,:,1) = dhTT(:,:,nf); dg_ho(:,:,2) = dhT0(:,:,nf); dg_ho(:,:,3) = dh0T(:,:,nf)
  end subroutine diag_pre

  ! precision numbers (max |A-B|) of the four pairs: (GR, G00), (oTT, GTT), (oT0, GT0), (o0T, G0T); c = (G00, GT0, GTT, G0T)
  Subroutine diag_p4(GR, o, c00, cT0, cTT, c0T, p)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(:,:), o(:,:,:), c00(:,:), cT0(:,:), cTT(:,:), c0T(:,:)
    Real (Kind=Kind(0.d0)), intent(out) :: p(4)
    Real (Kind=Kind(0.d0)) :: xm
    call COMPARE(GR, c00, p(1), xm)
    call COMPARE(o(:,:,1), cTT, p(2), xm)
    call COMPARE(o(:,:,2), cT0, p(3), xm)
    call COMPARE(o(:,:,3), c0T, p(4), xm)
  end subroutine diag_p4

  Subroutine diag_post(nfe, nst, GR)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL)
    Integer, parameter :: nv = 5
    Integer, parameter :: masks(nv) = (/ 1, 2, 4, 8, 7 /)
    Integer :: nf, n, n2, v, br, brd, i
    Real (Kind=Kind(0.d0)) :: pdev(4), pB(4), pC(4), pD(4), pE(4,nv), dCd, dU, dV, dS, d3mx, d3mn
    Real (Kind=Kind(0.d0)) :: l1(3), l2(3), pP(0:npc-1), xm
    Integer :: k
    Character (len=1024) :: line
    n = t_N; n2 = t_N2; nf = Calc_Fl_map(nfe)
    dg_npt = dg_npt + 1
    brd = last_br
    pdev = last_prec(1,:)
    call la_sync()
    call upd_from_c(t_D3, n2)                       ! diagnostic D->H: the 2N scales of the production CGR2_2
    d3mx = maxval(real(t_D3)); d3mn = minval(real(t_D3))
    ! device outputs of the production call
    call upd_from_c(GT0_d(1,1,nf), n*n); call upd_from_c(G0T_d(1,1,nf), n*n); call upd_from_c(GTT_d(1,1,nf), n*n)
    ! (B) full host chain (host TAU_M)
    call CGR2_2(dhT0(:,:,nf), dh00(:,:,nf), dhTT(:,:,nf), dh0T(:,:,nf), dg_hr(nfe), dg_hl(nfe), n)
    call diag_p4(GR(:,:,nf), dg_ho, dh00(:,:,nf), dhT0(:,:,nf), dhTT(:,:,nf), dh0T(:,:,nf), pB)
    ! (C) host CGR2_2 on the device inputs
    pC = 0.d0; dCd = 0.d0
    if (diag_full) then
    call CGR2_2(dg_c(:,:,2), dg_c(:,:,1), dg_c(:,:,3), dg_c(:,:,4), dg_hrd(nfe), dg_hl(nfe), n)
    call diag_p4(GR(:,:,nf), dg_o, dg_c(:,:,1), dg_c(:,:,2), dg_c(:,:,3), dg_c(:,:,4), pC)
    call mx(dCd, rel_diff(G00_d(:,:,nf), dg_c(:,:,1))); call mx(dCd, rel_diff(GT0_d(:,:,nf), dg_c(:,:,2)))
    call mx(dCd, rel_diff(GTT_d(:,:,nf), dg_c(:,:,3))); call mx(dCd, rel_diff(G0T_d(:,:,nf), dg_c(:,:,4)))
    endif
    ! GTT chains, each propagated on the host from its own CGR2_2 output of the previous point and compared with its own new GTT:
    !   0: the device output (device udvr), against the device's new GTT;
    !   1, 2: host CGR2_2 on the host udvr chain with U perturbed by ~1 ulp (fresh perturbation at each point);
    !   3, 4: host CGR2_2 on the device udvr with U perturbed by ~1 ulp.
    call COMPARE(dg_p(:,:,nf,0), GTT_d(:,:,nf), pP(0), xm)
    dg_p(:,:,nf,0) = GTT_d(:,:,nf)
    pP(1:) = 0.d0
    if (n_pass+1 >= dg_from) then
       do k = 1, npc-1
          if (k <= 2) then
             call perturb_state(dg_hr(nfe), dg_hp(1), dg_seed)
          else
             call perturb_state(dg_hrd(nfe), dg_hp(1), dg_seed)
          endif
          call CGR2_2(dg_c(:,:,2), dg_c(:,:,1), dg_c(:,:,3), dg_c(:,:,4), dg_hp(1), dg_hl(nfe), n)
          call COMPARE(dg_p(:,:,nf,k), dg_c(:,:,3), pP(k), xm)
          dg_p(:,:,nf,k) = dg_c(:,:,3)
       enddo
    endif
    pD = 0.d0; pE = 0.d0; dU = 0.d0; dV = 0.d0; dS = 0.d0
    if (diag_full) then
    ! (D) host CGR2_2 on the host one-block WRAPUR of the device udvr
    call CGR2_2(dg_c(:,:,2), dg_c(:,:,1), dg_c(:,:,3), dg_c(:,:,4), dg_hw(nfe), dg_hl(nfe), n)
    call diag_p4(GR(:,:,nf), dg_o, dg_c(:,:,1), dg_c(:,:,2), dg_c(:,:,3), dg_c(:,:,4), pD)
    call state_diff(dg_hrd(nfe), dg_hw(nfe), dU, dV, dS)
    ! (E) dev_cgr22 with host twins of single steps
    do v = 1, nv
       swp = masks(v)
       call dev_cgr22(nfe, nst, 1.d0, 0.d0, br, t_S(1,1,1), t_S(1,1,2), t_S(1,1,3), t_S(1,1,4))
       call la_sync()
       call upd_from_c(t_S(1,1,1), 4*n*n)            ! diagnostic D->H of the variant's outputs
       call diag_p4(GR(:,:,nf), dg_o, t_S(:,:,1), t_S(:,:,2), t_S(:,:,3), t_S(:,:,4), pE(:,v))
    enddo
    swp = swp_prod
    endif
    last_br = brd
#if defined(STABLOG)
    l1 = (/ dg_hl(nfe)%L(1), minval(dg_hl(nfe)%L), maxval(dg_hl(nfe)%L) /)
    l2 = (/ dg_hrd(nfe)%L(1), minval(dg_hrd(nfe)%L), maxval(dg_hrd(nfe)%L) /)
#else
    l1 = (/ log(real(dg_hl(nfe)%D(1))), log(minval(real(dg_hl(nfe)%D))), log(maxval(real(dg_hl(nfe)%D))) /)
    l2 = (/ log(real(dg_hrd(nfe)%D(1))), log(minval(real(dg_hrd(nfe)%D))), log(maxval(real(dg_hrd(nfe)%D))) /)
#endif
    write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,4ES9.2,A,4ES9.2,A,ES9.2,A,ES9.2,A,5ES9.2,A,5ES9.2,A,ES9.2,A,3ES9.2,A,3F9.2,A,3F9.2,A,2ES9.2)') &
         & 'TAUDIAG pt ', dg_npt, ' pass ', n_pass+1, ' nst ', nst, ' nfe ', nfe, ' br ', brd, &
         & ' dev', pdev, ' | hostfull', pB, ' | hC', maxval(pC), ' hD', maxval(pD), ' | E1,2,4,8,7', (maxval(pE(:,i)), i = 1, nv), &
         & ' | P', pP, ' | dCd', dCd, ' dUVS', dU, dV, dS, ' | L1', l1, ' L2', l2, ' D3', d3mx, d3mn
    call say(line)
  end subroutine diag_post

  ! host twins of single steps of dev_cgr22 (diagnostic only; every copy below is a diagnostic transfer)
  Subroutine sw_v1inv(nfe, nst)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Complex (Kind=Kind(0.d0)) :: z
    call la_sync()
    call upd_from_c(udv_V(1,1,nst,nfe), t_N*t_N)
    call INV(udv_V(:,:,nst,nfe), t_V1i, z)
    call upd_to_c(t_V1i, t_N*t_N)
  end subroutine sw_v1inv

  Subroutine sw_trsm()
    Implicit none
    integer :: n2
    n2 = t_N2
    call la_sync()
    call upd_from_c(t_H1, n2*n2); call upd_from_c(t_H2, n2*n2)
    call ZTRSM('L', 'U', 'C', 'N', n2, n2, z1, t_H1, n2, t_H2, n2)
    call upd_to_c(t_H2, n2*n2)
  end subroutine sw_trsm

  Subroutine sw_unmqr()
    Implicit none
    integer :: n2, info
    n2 = t_N2
    call la_sync()
    call upd_from_c(t_H1, n2*n2); call upd_from_c(t_H2, n2*n2); call upd_from_c(t_tau, n2)
    call ZUNMQR('L', 'N', n2, n2, n2, t_H1, n2, t_tau, t_H2, n2, dg_hwk, size(dg_hwk), info)
    call upd_to_c(t_H2, n2*n2)
  end subroutine sw_unmqr

!====================================================================
!  Device kernels (plain arrays; synchronous target regions)
!====================================================================
  ! G00 given: GT0 = GTT = G00, G0T = -(1 - G00)   (TAU_M at tau = 0)
  Subroutine k_tau_init(n, G00, GT0, GTT, G0T)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: G00(n,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: GT0(n,*), GTT(n,*), G0T(n,*)
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
  Subroutine k_fault_scale(A)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(*)
    !$omp target
    A(1) = A(1) * (1.d0 + 1.d-6)
    !$omp end target
  end subroutine k_fault_scale

  ! Dm(i) = D(i) * s     (default stabilization: D1m = D1, D2m = D2; s = 1 in production)
  Subroutine k_dm(n, D, s, Dm)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
    Real (Kind=Kind(0.d0)), intent(in) :: s
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dm(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
       Dm(i) = D(i) * s
    enddo
  end subroutine k_dm

  ! STABLOG: Le = L(i) + dl;  Le <= 0: Dm = exp(Le), Rs = 1;  else Dm = 1, Rs = exp(-Le)
  Subroutine k_dm_log(n, L, dl, Dm, Rs)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(in) :: L(*), dl
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dm(*)
    Real (Kind=Kind(0.d0)), intent(inout) :: Rs(*)
    Integer :: i
    Real (Kind=Kind(0.d0)) :: le
    !$omp target teams distribute parallel do private(le)
    do i = 1, n
       le = L(i) + dl
       if (le <= 0.d0) then
          Dm(i) = cmplx(exp(le), 0.d0, kind(0.d0))
          Rs(i) = 1.d0
       else
          Dm(i) = cmplx(1.d0, 0.d0, kind(0.d0))
          Rs(i) = exp(-le)
       endif
    enddo
  end subroutine k_dm_log

  ! A(i,:) *= Rs(i)  (ZSCAL with a real factor; factors equal to 1 leave A unchanged)
  Subroutine k_rowscale_r(n, A, lda, Rs)
    Implicit none
    Integer, intent(in) :: n, lda
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(lda,*)
    Real (Kind=Kind(0.d0)), intent(in) :: Rs(*)
    Integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
       do i = 1, n
          A(i,j) = cmplx(Rs(i), 0.d0, kind(0.d0)) * A(i,j)
       enddo
    enddo
  end subroutine k_rowscale_r

  ! HLPB2 (2N x 2N) from the blocks, branch 1 (D1(1) > D2(1)): [[V1INV, D1m*U1^dagger], [-D2m*V2, MYU2]]; branch 2:
  ! [[MYU2, -D2m*V2], [D1m*U1^dagger, V1INV]]   (D1m(I) conjg(U1(J,I)) at (I,J), D2m(I) V2(I,J))
  Subroutine k_build_H2(n, br, V1i, U2h, D1m, U1, D2m, V2, H)
    Implicit none
    Integer, intent(in) :: n, br
    Complex (Kind=Kind(0.d0)), intent(in) :: V1i(n,*), U2h(n,*), D1m(*), U1(n,*), D2m(*), V2(n,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: H(2*n,*)
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
    Integer(c_int), intent(in) :: ipvt(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: V1i(n,*), U2h(n,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: H(2*n,*)
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

  ! A <- 3/2 - A/2  (A = Q^H Q on input)
  Subroutine k_ns2(n, A)
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
  end subroutine k_ns2

  ! H(i,:) *= conj(1/D(i))   (D real valued: the factor is the real 1/Re D(i))
  Subroutine k_rowscale_inv(n, H, ldh, D)
    Implicit none
    Integer, intent(in) :: n, ldh
    Complex (Kind=Kind(0.d0)), intent(inout) :: H(ldh,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: D(*)
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
  ! GR0T = (2,1), GR00 = (2,2)   (blocks of the 2N x 2N matrix).  swp: fault injection, GT0 and G0T exchanged.
  Subroutine k_blocks(n, br, swp, H, G00, GT0, GTT, G0T)
    Implicit none
    Integer, intent(in) :: n, br
    Logical, intent(in) :: swp
    Complex (Kind=Kind(0.d0)), intent(in) :: H(2*n,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: G00(n,*), GT0(n,*), GTT(n,*), G0T(n,*)
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
!  Host oracle (copies of PROPR / PROPRM1 of tau_m_mod, which cannot be used here: tau_m_mod uses this module)
!====================================================================
  Subroutine h_propr(A, nt)
    Implicit none
    Complex (Kind=Kind(0.D0)), intent(INOUT) :: A(Ndim,Ndim,N_FL)
    Integer, INTENT(IN) :: nt
    Integer :: nf, nfe, n
    Do nfe = 1, N_FL_eff
       nf = Calc_Fl_map(nfe)
       Call Hop_mod_mmthr(A(:,:,nf), nf, nt)
       Do n = 1, Size(Op_V,1)
          Call Op_mmultR(A(:,:,nf), Op_V(n,nf), nsigma%f(n,nt), 'n', nt)
       ENDDO
    Enddo
  end subroutine h_propr

  Subroutine h_proprm1(A, nt)
    Implicit none
    Complex (Kind=Kind(0.D0)), intent(Inout) :: A(Ndim,Ndim,N_FL)
    Integer, intent(in) :: nt
    Integer :: nf, nfe, n, sgn
    sgn = -1
    do nfe = 1, N_FL_eff
       nf = Calc_Fl_map(nfe)
       Call Hop_mod_mmthl_m1(A(:,:,nf), nf, nt)
       Do n = 1, Size(Op_V,1)
          Call Op_mmultL(A(:,:,nf), Op_V(n,nf), nsigma%f(n,nt), 'n', nt, sgn)
       Enddo
    enddo
  end subroutine h_proprm1

!====================================================================
!  Hook: self-test, timing
!====================================================================
  Subroutine gpu_tau_hook()
    Implicit none
    Character (len=256) :: f
    Integer :: nfe, ib
    Logical :: did
    did = .false.
    if (.not. gpu_tau_on) return
    ! ticket 212: ALF_GPU_FAULT_TAU is refused in production runs and its value was validated (and the armed line printed) by
    ! gpu_env_check; eps:<x> acts in the SELFTEST=4 comparison only (gpu_env_fault_tau_eps)
    if (gpu_env_get('ALF_GPU_FAULT_TAU', f)) then
       select case (trim(f))
       case ('perm');   fault = 1
       case ('prop');   fault = 2
       case ('prec');   fault = 3
       case ('blocks'); fault = 4
       case default                  ! eps:<x>
       end select
    endif
    if (gpu_env_is('ALF_GPU_SELFTEST', '4')) then
       call gpu_tau_selftest()
       did = .true.
    endif
    if (gpu_env_is('ALF_GPU_TAUTIME', '1')) then
       call gpu_tau_timing()
       did = .true.
    endif
    if (did) then
       ! the tests leave the device slots dirty and the counters changed: restore the state before the hook
       do nfe = 1, t_nfe
          call gpu_stab_reset_l(nfe); call gpu_stab_reset_r(nfe)
          do ib = 1, t_nstm
             call gpu_stab_reset_st(ib, nfe)
          enddo
       enddo
       prec_report = .true.
       n_pass = 0; n_branch = 0; n_cgr = 0; b_h2d = 0; b_d2h = 0
    endif
    if (did .and. gpu_env_is('ALF_GPU_EXIT', '1')) then
       call gpu_tau_finalize()
       call gpu_stab_finalize()
       call gpu_state_finalize()
       write(output_unit,'(A)') ' GPU tau pass hook: ALF_GPU_EXIT=1, stopping after the test'
       stop
    endif
  end subroutine gpu_tau_hook

  !> The states of the sweep as left after the backward sweep: udvst(NST) = left propagation of the fixed fields (device),
  !> GR = G at time 0 (NVAR = 1 CGR of udvl and the unit udvr), also on the host: hst(NST,nfe) host copies of the udvst slots.
  Subroutine build_states(GR, hst)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(Ndim,Ndim,N_FL)
    Class (UDV_State), intent(inout) :: hst(:,:)
    Integer :: nfe, nf, nst
    Complex (Kind=Kind(0.d0)) :: ph
    do nfe = 1, t_nfe
       call gpu_stab_reset_l(nfe); call gpu_stab_reset_r(nfe)
       call gpu_stab_reset_st(t_nstm, nfe)
    enddo
    do nst = t_nstm-1, 1, -1
       call gpu_stab_wrapul(t_nt(nst+1), t_nt(nst))
       do nfe = 1, t_nfe
          call gpu_stab_store_l(nst, nfe)
       enddo
    enddo
    call gpu_stab_wrapul(t_nt(1), 0)
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       call gpu_stab_reset_r(nfe)
       call gpu_stab_cgr(nfe, nf, 1, GR, ph, .false.)
    enddo
    do nfe = 1, t_nfe
       do nst = 1, t_nstm
          call dl_state(2, nst, nfe, hst(nst,nfe))        ! self-test / timing: host copy of the device state for the oracle
          hst(nst,nfe)%side = 'l'
       enddo
    enddo
  end subroutine build_states

  Subroutine alloc_hst(hst)
    Implicit none
    Class (UDV_State), Allocatable, intent(inout) :: hst(:,:)
    Integer :: i, j
    allocate(hst(t_nstm, t_nfe))
    do j = 1, t_nfe
       do i = 1, t_nstm
          call hst(i,j)%alloc(t_N)
       enddo
    enddo
  end subroutine alloc_hst

  ! tau = 0 of TAU_M on host arrays (the CPU's loops)
  Subroutine h_begin(GR, h00, hT0, hTT, h0T)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(Ndim,Ndim,N_FL)
    Complex (Kind=Kind(0.d0)), intent(inout) :: h00(Ndim,Ndim,N_FL), hT0(Ndim,Ndim,N_FL), hTT(Ndim,Ndim,N_FL), h0T(Ndim,Ndim,N_FL)
    Integer :: nfe, nf, I, J
    Complex (Kind=Kind(0.d0)) :: Z
    do nfe = 1, N_FL_eff
       nf = Calc_Fl_map(nfe)
       do J = 1, Ndim
          do I = 1, Ndim
             Z = cmplx(0.d0, 0.d0, kind(0.D0))
             if (I == J) Z = cmplx(1.d0, 0.d0, kind(0.d0))
             h00(I,J,nf) = GR(I,J,nf)
             hT0(I,J,nf) = GR(I,J,nf)
             hTT(I,J,nf) = GR(I,J,nf)
             h0T(I,J,nf) = -(Z - GR(I,J,nf))
          enddo
       enddo
    enddo
  end subroutine h_begin

!--------------------------------------------------------------------
!> Self-test (ALF_GPU_SELFTEST=4): the fixed initial fields; one full tau pass, device vs the upstream host routines.
!>   * pass: lockstep device / host (independent trajectories from the same GR and the same udvst): GT0, G0T, GTT, G00 at every
!>     time slice, max relative difference (|a-b|_max / |b|_max) per matrix, overall and right after CGR2_2;
!>   * CGR2_2 alone: at every stabilization point and flavour, from the same UDV states (the device udvr copied to a host
!>     object), device vs host CGR2_2 (outputs); natural branch counts; two synthetic variants, D1 (L1) scaled so that the
!>     branches are exercised from either side;
!>   * precision reduction vs COMPARE on the same matrices; 1-ulp sensitivity (reported).
!--------------------------------------------------------------------
  Subroutine gpu_tau_selftest()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-10, tolr = 1.d-12
    Integer, parameter :: nvar = 3                         ! variants: 1 natural, 2 D1 * 1e3, 3 D1 * 1e-3
    Real (Kind=Kind(0.d0)), parameter :: ssc(3) = (/ 1.d0, 1.d3, 1.d-3 /)
    Class (UDV_State), Allocatable :: hst(:,:), hr(:), hrd(:), hsy(:), hp(:)
    Complex (Kind=Kind(0.d0)), Allocatable :: GR(:,:,:), h00(:,:,:), hT0(:,:,:), hTT(:,:,:), h0T(:,:,:)
    Complex (Kind=Kind(0.d0)), Allocatable :: o00(:,:,:), oT0(:,:,:), oTT(:,:,:), o0T(:,:,:)       ! device matrices before CGR2_2 (host copies)
    Complex (Kind=Kind(0.d0)), Allocatable :: c(:,:,:), cp(:,:,:), dS(:,:,:)       ! host CGR2_2 outputs: (N,N,4) = G00, GT0, GTT, G0T
    Integer :: nfe, nf, nst, nt, nt1, ntst, N, k, v, br
    Real (Kind=Kind(0.d0)) :: dU, dV, dS_, xmax_h, xmean_h
    Real (Kind=Kind(0.d0)) :: m_all(4), m_cgr(4), m_alone(4,nvar), m_state(3), m_prec, sens_h, sens_d, wslice(4)
    Integer :: nslices, nstab, nalone(2), nsyn(2,2), worst(4)
    Integer(int64) :: seed
    Character (len=512) :: line
    Character (len=4), parameter :: nm(4) = (/ 'G00 ', 'GT0 ', 'GTT ', 'G0T ' /)
    Real (Kind=Kind(0.d0)) :: dd(4)
    Real (Kind=Kind(0.d0)) :: dvec(4)
    Logical :: verbose

#if !defined(STABLOG) && !defined(QR_PREPIVOT)
    write(error_unit,*) ' ALF_GPU_SELFTEST=4 with the default stabilization needs a build with -DQR_PREPIVOT (ALF_FLAGS_EXT=-DQR_PREPIVOT before configure.sh): the host oracle'
    write(error_unit,*) ' must be the same pre-pivoted algorithm (ZGEQP3 gives other, equally valid, pivots).'
    call info_line('GPU tau self-test: needs a build with ALF_FLAGS_EXT=-DQR_PREPIVOT for the default stabilization')
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
    N = t_N
    seed = 20251107_int64
    call say('GPU tau self-test (ticket 207) start')
    allocate(GR(N,N,N_FL), h00(N,N,N_FL), hT0(N,N,N_FL), hTT(N,N,N_FL), h0T(N,N,N_FL))
    allocate(o00(N,N,N_FL), oT0(N,N,N_FL), oTT(N,N,N_FL), o0T(N,N,N_FL), c(N,N,4), cp(N,N,4), dS(N,N,4))
    GR = z0; h00 = z0; hT0 = z0; hTT = z0; h0T = z0; o00 = z0; oT0 = z0; oTT = z0; o0T = z0
    call alloc_hst(hst)
    allocate(hr(t_nfe), hrd(t_nfe), hsy(t_nfe), hp(t_nfe))
    do nfe = 1, t_nfe
       call hr(nfe)%init(N, 'r'); call hrd(nfe)%init(N, 'r'); call hsy(nfe)%alloc(N); call hp(nfe)%init(N, 'r')
       hsy(nfe)%side = 'l'
    enddo
    call build_states(GR, hst)
    prec_report = .false.
    verbose = gpu_env_is('ALF_GPU_TAUTEST_V', '1')
    m_all = 0.d0; m_cgr = 0.d0; m_alone = 0.d0; m_state = 0.d0; m_prec = 0.d0; sens_h = 0.d0; sens_d = 0.d0
    nslices = 0; nstab = 0; nalone = 0; nsyn = 0; worst = 0; wslice = 0.d0

    ! ---- the pass: device (tau_begin, tau_slice, gpu_stab_wrapur, dev_cgr22) against the host routines, slice by slice
    call tau_begin(GR)
    call h_begin(GR, h00, hT0, hTT, h0T)
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       call hr(nfe)%reset('r')
    enddo
    ! tau = 0: the device derives GT0, GTT, G0T from G00 (self-test: D->H of all four and comparison with the CPU's loops)
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       call upd_from_c(GT0_d(1,1,nf), N*N); call upd_from_c(G0T_d(1,1,nf), N*N)
       call upd_from_c(GTT_d(1,1,nf), N*N); call upd_from_c(G00_d(1,1,nf), N*N)
    enddo
    call cmp_all(h00, hT0, hTT, h0T, dd)
    call say('GPU tau self-test slice 0 (device init vs host loops): G00, GT0, GTT, G0T '//fmt4(dd))
    do k = 1, 4
       call mx(m_all(k), dd(k))
    enddo
    nst = 1
    do nt = 0, Ltrot - 1
       nt1 = nt + 1
       ! device and host propagation of one slice
       call tau_slice(nt1)
       call tau_download(.false.)
       call h_propr(hT0, nt1)
       call h_proprm1(h0T, nt1)
       call h_proprm1(hTT, nt1)
       call h_propr(hTT, nt1)
       call cmp_all(h00, hT0, hTT, h0T, dd)
       do k = 1, 4
          call mx(m_all(k), dd(k))
          if (dd(k) > wslice(k)) then
             wslice(k) = dd(k); worst(k) = nt1
          endif
       enddo
       nslices = nslices + 1
       if (verbose) then
          write(line,'(A,I4,A,4ES11.3)') 'GPU tau self-test slice ', nt1, ' (before a CGR2_2 at this slice): G00, GT0, GTT, G0T ', dd
          call say(line)
       endif
       if (t_nt(nst) == nt1) then
          ntst = t_nt(nst-1)
          nstab = nstab + 1
          ! right propagation: host chain (own objects) and device slot; states compared for information
          call WRAPUR(ntst, nt1, hr)
          call gpu_stab_wrapur(ntst, nt1)
          do nfe = 1, t_nfe
             call dl_state(1, 2, nfe, hrd(nfe))                   ! self-test: the device udvr as a host object (input of the CGR2_2 comparison)
             call state_diff(hrd(nfe), hr(nfe), dU, dV, dS_)
             call mx(m_state(1), dU); call mx(m_state(2), dV); call mx(m_state(3), dS_)
          enddo
          ! ---- CGR2_2 alone, per flavour: the same states on both sides
          do nfe = 1, t_nfe
             nf = Calc_Fl_map(nfe)
             do v = 1, nvar
                ! left state (host): udvst(nst), scaled in the synthetic variants
                hsy(nfe)%U = hst(nst,nfe)%U; hsy(nfe)%V = hst(nst,nfe)%V
#if defined(STABLOG)
                if (v == 1) then
                   hsy(nfe)%L = hst(nst,nfe)%L
                else
                   hsy(nfe)%L = hst(nst,nfe)%L + log(ssc(v))
                endif
#else
                hsy(nfe)%D = hst(nst,nfe)%D * ssc(v)
#endif
                call CGR2_2(c(:,:,2), c(:,:,1), c(:,:,3), c(:,:,4), hrd(nfe), hsy(nfe), N)
                call dev_cgr22(nfe, nst, ssc(v), log(ssc(v)), br, t_S(1,1,1), t_S(1,1,2), t_S(1,1,3), t_S(1,1,4))
                if (v == 1) then
                   nalone(br) = nalone(br) + 1
                else
                   nsyn(br,v-1) = nsyn(br,v-1) + 1
                endif
                ! self-test: D->H of the four scratch outputs for the comparison
                call upd_from_c(t_S(1,1,1), 4*N*N)
                call cmp_set(t_S, c, dvec, .true.)
                do k = 1, 4
                   call mx(m_alone(k,v), dvec(k))
                enddo
                if (v == 1) then
                   ! intrinsic sensitivity (reported, not a criterion): U_R perturbed by 2e-16 relative per entry
                   call perturb_state(hrd(nfe), hp(nfe), seed)
                   call CGR2_2(cp(:,:,2), cp(:,:,1), cp(:,:,3), cp(:,:,4), hp(nfe), hsy(nfe), N)
                   call cmp_set(cp, c, dvec, .false.)
                   do k = 1, 4
                      call mx(sens_h, dvec(k))
                   enddo
                   dS = t_S                                        ! the unperturbed device result (host copy)
                   call ul_state(2, nfe, hp(nfe))                  ! self-test: upload the perturbed udvr
                   call dev_cgr22(nfe, nst, 1.d0, 0.d0, br, t_S(1,1,1), t_S(1,1,2), t_S(1,1,3), t_S(1,1,4))
                   call upd_from_c(t_S(1,1,1), 4*N*N)
                   call cmp_set(t_S, dS, dvec, .false.)
                   do k = 1, 4
                      call mx(sens_d, dvec(k))
                   enddo
                   call ul_state(2, nfe, hrd(nfe))                 ! self-test: restore udvr
                endif
             enddo
          enddo
          ! ---- the pass itself: host chain CGR2_2 with the precision numbers vs the device reduction on the same matrices
          do nfe = 1, t_nfe
             nf = Calc_Fl_map(nfe)
             ! device matrices before CGR2_2 (host shadows are current for GT0, G0T, GTT, G00)
             o00(:,:,nf) = G00_d(:,:,nf); oT0(:,:,nf) = GT0_d(:,:,nf); oTT(:,:,nf) = GTT_d(:,:,nf); o0T(:,:,nf) = G0T_d(:,:,nf)
             Call CGR2_2(hT0(:,:,nf), h00(:,:,nf), hTT(:,:,nf), h0T(:,:,nf), hr(nfe), hst(nst,nfe), N)
             call tau_stab_cgr(nfe, nst)
             ! self-test: D->H of all four device matrices (production copies G00 only)
             call upd_from_c(GT0_d(1,1,nf), N*N); call upd_from_c(G0T_d(1,1,nf), N*N)
             call upd_from_c(GTT_d(1,1,nf), N*N); call upd_from_c(G00_d(1,1,nf), N*N)
             ! precision reduction vs COMPARE on the device's own matrices: (GR, G00new), (TTold, TTnew), (T0old, T0new), (0Told, 0Tnew)
             call COMPARE(GR(:,:,nf), G00_d(:,:,nf), xmax_h, xmean_h)
             call cmp_prec(1, xmax_h, xmean_h, m_prec)
             call COMPARE(oTT(:,:,nf), GTT_d(:,:,nf), xmax_h, xmean_h)
             call cmp_prec(2, xmax_h, xmean_h, m_prec)
             call COMPARE(oT0(:,:,nf), GT0_d(:,:,nf), xmax_h, xmean_h)
             call cmp_prec(3, xmax_h, xmean_h, m_prec)
             call COMPARE(o0T(:,:,nf), G0T_d(:,:,nf), xmax_h, xmean_h)
             call cmp_prec(4, xmax_h, xmean_h, m_prec)
          enddo
          call cmp_all(h00, hT0, hTT, h0T, dd)
          do k = 1, 4
             call mx(m_cgr(k), dd(k))
             call mx(m_all(k), dd(k))
          enddo
          nst = nst + 1
       endif
    enddo
    call tau_end()
    prec_report = .true.

    write(line,'(A,I0,A,I0,A)') 'GPU tau self-test pass over ', nslices, ' slices (', nstab, ' stabilization points): max rel diff device vs host (all slices incl. after CGR2_2):'
    call say(line)
    write(line,'(A,ES10.3,A,I0,A,ES10.3,A,I0,A,ES10.3,A,I0,A,ES10.3,A,I0,A)') '   G00 ', m_all(1), ' (worst slice ', worst(1), '), GT0 ', m_all(2), ' (', worst(2), '), GTT ', m_all(3), ' (', worst(3), '), G0T ', m_all(4), ' (', worst(4), ')'
    call say(line)
    write(line,'(A,4ES11.3,A,ES8.1)') 'GPU tau self-test right after CGR2_2 (chains independent): G00, GT0, GTT, G0T ', m_cgr, ' (tol ', tol
    call say(line)
    write(line,'(A,I0,A,I0,A,I0)') 'GPU tau self-test CGR2_2 alone: natural calls with D1(1) > D2(1): ', nalone(1), ', else: ', nalone(2), ', total ', nalone(1)+nalone(2)
    call say(line)
    write(line,'(A,4ES11.3)') 'GPU tau self-test CGR2_2 alone, natural states, max rel diff G00, GT0, GTT, G0T:   ', m_alone(:,1)
    call say(line)
    write(line,'(A,I0,A,I0,A,4ES11.3)') 'GPU tau self-test CGR2_2 alone, D1 scaled x1e3 (branches 1/2: ', nsyn(1,1), '/', nsyn(2,1), '):  ', m_alone(:,2)
    call say(line)
    write(line,'(A,I0,A,I0,A,4ES11.3)') 'GPU tau self-test CGR2_2 alone, D1 scaled x1e-3 (branches 1/2: ', nsyn(1,2), '/', nsyn(2,2), '): ', m_alone(:,3)
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU tau self-test (reported) udvr state device vs host chain: U ', m_state(1), ', V ', m_state(2), ', scales ', m_state(3)
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3)') 'GPU tau self-test (reported) intrinsic sensitivity of CGR2_2 to a 2e-16 relative (1 ulp) perturbation of U_R: host ', sens_h, ', device ', sens_d
    call say(line)
    write(line,'(A,ES10.3,A,ES8.1)') 'GPU tau self-test precision reduction vs COMPARE: max rel diff (max and mean) ', m_prec, ' (tol ', tolr
    call say(line)
    if (.not. (all(m_all <= tol))) call selftest_fail('tau: pass (propagation / WRAPUR / CGR2_2) vs host')
    if (.not. (all(m_alone <= tol))) call selftest_fail('tau: CGR2_2 alone vs host')
    if (.not. (m_prec <= tolr)) call selftest_fail('tau: precision reduction')
    if (nalone(1) + nalone(2) /= nstab*t_nfe) call selftest_fail('tau: CGR2_2 call count')
    do nfe = 1, t_nfe
       call hr(nfe)%dealloc; call hrd(nfe)%dealloc; call hsy(nfe)%dealloc; call hp(nfe)%dealloc
    enddo
    deallocate(GR, h00, hT0, hTT, h0T, o00, oT0, oTT, o0T, c, cp, dS, hr, hrd, hsy, hp, hst)
    call say('GPU tau self-test PASSED')
  end subroutine gpu_tau_selftest

  ! m = max(m, d); a non-finite d makes m huge (same as mx of gpu_stab_mod)
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

  Character (len=64) function fmt4(d)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(in) :: d(4)
    write(fmt4,'(4ES11.3)') d
  end function fmt4

  ! relative differences (device shadows vs host arrays) of the four matrices, flavours in the map only
  Subroutine cmp_all(h00, hT0, hTT, h0T, dd)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: h00(Ndim,Ndim,N_FL), hT0(Ndim,Ndim,N_FL), hTT(Ndim,Ndim,N_FL), h0T(Ndim,Ndim,N_FL)
    Real (Kind=Kind(0.d0)), intent(out) :: dd(4)
    Integer :: nfe, nf
    dd = 0.d0
    do nfe = 1, t_nfe
       nf = Calc_Fl_map(nfe)
       call mx(dd(1), rel_diff(G00_d(:,:,nf), h00(:,:,nf)))
       call mx(dd(2), rel_diff_f(GT0_d(:,:,nf), hT0(:,:,nf)))     ! graded fault (eps:<x>) on a copy, see rel_diff_f
       call mx(dd(3), rel_diff(GTT_d(:,:,nf), hTT(:,:,nf)))
       call mx(dd(4), rel_diff(G0T_d(:,:,nf), h0T(:,:,nf)))
    enddo
  end subroutine cmp_all

  ! relative differences of two (N,N,4) sets (G00, GT0, GTT, G0T); dev: a is the device result compared with the host (the graded
  ! fault acts on its GT0)
  Subroutine cmp_set(a, b, dvec, dev)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: a(:,:,:), b(:,:,:)
    Real (Kind=Kind(0.d0)), intent(out) :: dvec(4)
    Logical, intent(in) :: dev
    Integer :: k
    do k = 1, 4
       if (dev .and. k == 2) then
          dvec(k) = rel_diff_f(a(:,:,k), b(:,:,k))
       else
          dvec(k) = rel_diff(a(:,:,k), b(:,:,k))
       endif
    enddo
  end subroutine cmp_set

  ! rel_diff(A, B) of the host copy A of a device GT0 against the host B; with the graded fault ALF_GPU_FAULT_TAU=eps:<x> armed, the
  ! element of largest |A| of a copy of A is multiplied by (1 + x) first (self-test comparison only; A itself is not changed)
  Real (Kind=Kind(0.d0)) function rel_diff_f(A, B)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:), B(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: Af(:,:)
    Integer :: i, j, im, jm
    Real (Kind=Kind(0.d0)) :: x, a_, amax
    x = gpu_env_fault_tau_eps()
    if (x <= 0.d0) then
       rel_diff_f = rel_diff(A, B)
       return
    endif
    Af = A
    amax = -1.d0; im = 1; jm = 1
    do j = 1, size(Af,2)
       do i = 1, size(Af,1)
          a_ = abs(Af(i,j))
          if (a_ > amax) then
             amax = a_; im = i; jm = j
          endif
       enddo
    enddo
    Af(im,jm) = Af(im,jm) * (1.d0 + x)
    rel_diff_f = rel_diff(Af, B)
    deallocate(Af)
  end function rel_diff_f

  ! device reduction numbers of the pair k (last tau_stab_cgr) against COMPARE: relative difference of max and mean
  Subroutine cmp_prec(k, xmax_h, xmean_h, m)
    Implicit none
    Integer, intent(in) :: k
    Real (Kind=Kind(0.d0)), intent(in) :: xmax_h, xmean_h
    Real (Kind=Kind(0.d0)), intent(inout) :: m
    Real (Kind=Kind(0.d0)) :: dmax, dmean
    dmax  = abs(last_prec(1,k) - xmax_h) / max(xmax_h, tiny(1.d0))
    dmean = abs(last_prec(2,k)/dble(t_N*t_N) - xmean_h) / max(xmean_h, tiny(1.d0))
    call mx(m, dmax); call mx(m, dmean)
  end subroutine cmp_prec

  ! hp = a with every entry of U multiplied by (1 + 4e-16 * r), r in (-1/2, 1/2): +-2e-16, about one ulp
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

!====================================================================
!  Timing (ALF_GPU_TAUTIME=1): one tau pass without the observables, device vs host, same states
!====================================================================
  Subroutine gpu_tau_timing()
    Implicit none
    Class (UDV_State), Allocatable :: hst(:,:), hr(:)
    Complex (Kind=Kind(0.d0)), Allocatable :: GR(:,:,:), h00(:,:,:), hT0(:,:,:), hTT(:,:,:), h0T(:,:,:), c(:,:,:)
    Integer :: nfe, nf, nst, nt, nt1, ntst, k, N, br, nrep, nrep_h
    Integer (int64) :: c0, c1, crate, h2d0, d2h0
    Real (Kind=Kind(0.d0)) :: tg(5), tgi(5), tgc(5), tch(5), tcc(5), by_h2d, by_d2h
    Character (len=512) :: line
    Character (len=8) :: mode

    N = t_N
    nrep = 5
    allocate(GR(N,N,N_FL), h00(N,N,N_FL), hT0(N,N,N_FL), hTT(N,N,N_FL), h0T(N,N,N_FL), c(N,N,4))
    GR = z0; h00 = z0; hT0 = z0; hTT = z0; h0T = z0
    call alloc_hst(hst)
    allocate(hr(t_nfe))
    do nfe = 1, t_nfe
       call hr(nfe)%init(N, 'r')
    enddo
    call build_states(GR, hst)
    prec_report = .false.
    by_h2d = 0.d0; by_d2h = 0.d0
#if defined(STABLOG)
    mode = 'LOG'
#else
    mode = 'default'
#endif
    tg = 0.d0; tgi = 0.d0; tgc = 0.d0
    do k = 1, nrep
       ! (a) device only: no D->H of the matrices, no precision reductions
       call tau_begin(GR)
       call la_sync()
       call system_clock(c0, crate)
       nst = 1
       do nt = 0, Ltrot - 1
          nt1 = nt + 1
          call tau_slice(nt1)
          if (t_nt(nst) == nt1) then
             ntst = t_nt(nst-1)
             call gpu_stab_wrapur(ntst, nt1)
             do nfe = 1, t_nfe
                call tau_stab_cgr_nodl(nfe, nst)
             enddo
             nst = nst + 1
          endif
       enddo
       call la_sync()
       call system_clock(c1)
       tg(k) = dble(c1-c0)/dble(crate)
       call tau_end()
       ! (b) as in the run, without the observables: D->H of GT0, G0T, GTT per slice, G00 per stabilization point, precision reductions
       h2d0 = b_h2d; d2h0 = b_d2h
       call system_clock(c0, crate)
       call tau_begin(GR)
       nst = 1
       do nt = 0, Ltrot - 1
          nt1 = nt + 1
          call tau_slice(nt1)
          call tau_download(.false.)
          if (t_nt(nst) == nt1) then
             ntst = t_nt(nst-1)
             call gpu_stab_wrapur(ntst, nt1)
             do nfe = 1, t_nfe
                call tau_stab_cgr(nfe, nst)
             enddo
             nst = nst + 1
          endif
       enddo
       call la_sync()
       call system_clock(c1)
       tgi(k) = dble(c1-c0)/dble(crate)
       call tau_end()
       by_h2d = dble(b_h2d - h2d0); by_d2h = dble(b_d2h - d2h0)
    enddo
    ! (c) CGR2_2 alone at the last stabilization point from a realistic right state (the whole chain 0..Ltrot, block by block)
    do nfe = 1, t_nfe
       call gpu_stab_reset_r(nfe)
    enddo
    do nst = 1, t_nstm
       call gpu_stab_wrapur(t_nt(nst-1), t_nt(nst))
    enddo
    nst = t_nstm
    do k = 1, nrep
       call la_sync()
       call system_clock(c0, crate)
       do nfe = 1, t_nfe
          call dev_cgr22(nfe, nst, 1.d0, 0.d0, br, t_S(1,1,1), t_S(1,1,2), t_S(1,1,3), t_S(1,1,4))
       enddo
       call la_sync()
       call system_clock(c1)
       tgc(k) = dble(c1-c0)/dble(crate)
    enddo
    do nfe = 1, t_nfe
       call dl_state(1, 2, nfe, hr(nfe))                    ! timing: host copy of the right state for the host CGR2_2
    enddo
    ! host pass without observables (the CPU's PROPR / PROPRM1 / WRAPUR / CGR2_2), fewer repetitions
    nrep_h = 3
    if (N >= 800) nrep_h = 2
    tch = 0.d0; tcc = 0.d0
    do k = 1, nrep_h
       nst = t_nstm
       call system_clock(c0, crate)
       do nfe = 1, t_nfe
          Call CGR2_2(c(:,:,2), c(:,:,1), c(:,:,3), c(:,:,4), hr(nfe), hst(nst,nfe), N)
       enddo
       call system_clock(c1)
       tcc(k) = dble(c1-c0)/dble(crate)
    enddo
    do k = 1, nrep_h
       call h_begin(GR, h00, hT0, hTT, h0T)
       do nfe = 1, t_nfe
          call hr(nfe)%reset('r')
       enddo
       call system_clock(c0, crate)
       nst = 1
       do nt = 0, Ltrot - 1
          nt1 = nt + 1
          call h_propr(hT0, nt1)
          call h_proprm1(h0T, nt1)
          call h_proprm1(hTT, nt1)
          call h_propr(hTT, nt1)
          if (t_nt(nst) == nt1) then
             ntst = t_nt(nst-1)
             call WRAPUR(ntst, nt1, hr)
             do nfe = 1, t_nfe
                nf = Calc_Fl_map(nfe)
                Call CGR2_2(hT0(:,:,nf), h00(:,:,nf), hTT(:,:,nf), h0T(:,:,nf), hr(nfe), hst(nst,nfe), N)
             enddo
             nst = nst + 1
          endif
       enddo
       call system_clock(c1)
       tch(k) = dble(c1-c0)/dble(crate)
    enddo
    call sort5(tg); call sort5(tgi); call sort5(tgc)
    call sortn(tch, nrep_h); call sortn(tcc, nrep_h)
    write(line,'(A,I0,A,I0,A,A,A,I0,A,I0,A,I0,A)') 'GPU tau timing N = ', N, ', Ltrot = ', Ltrot, ', ', trim(mode), ', flavours ', t_nfe, &
         & ', stabilization points ', t_nstm, ', CGR2_2 QR 2N = ', t_N2, ' (s per tau pass without observables):'
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,F8.1)') '   GPU device only ', tg(3), ' | GPU with the interim D->H copies ', tgi(3), ' | CPU ', tch((nrep_h+1)/2), &
         & ' | speedup incl. copies ', tch((nrep_h+1)/2)/tgi(3)
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3,A,F8.1)') '   CGR2_2 all flavours, one point: GPU ', tgc(3), ' CPU ', tcc((nrep_h+1)/2), ' speedup ', tcc((nrep_h+1)/2)/tgc(3)
    call say(line)
    write(line,'(A,F14.0,A,F14.0,A)') '   interim transfers per pass: H->D ', by_h2d, ' B, D->H ', by_d2h, ' B'
    call say(line)
    do nfe = 1, t_nfe
       call hr(nfe)%dealloc
    enddo
    deallocate(GR, h00, hT0, hTT, h0T, c, hr, hst)
    prec_report = .true.
  end subroutine gpu_tau_timing

  ! as tau_stab_cgr without the interim D->H of G00 and without the precision reductions (device only timing)
  Subroutine tau_stab_cgr_nodl(nfe, nst)
    Implicit none
    Integer, intent(in) :: nfe, nst
    Integer :: nf, n, br
    n = t_N
    nf = Calc_Fl_map(nfe)
    call k_copy_c(n*n, GTT_d(1,1,nf), t_old(1,1,1))
    call k_copy_c(n*n, GT0_d(1,1,nf), t_old(1,1,2))
    call k_copy_c(n*n, G0T_d(1,1,nf), t_old(1,1,3))
    call dev_cgr22(nfe, nst, 1.d0, 0.d0, br, G00_d(1,1,nf), GT0_d(1,1,nf), GTT_d(1,1,nf), G0T_d(1,1,nf))
  end subroutine tau_stab_cgr_nodl

  Subroutine sort5(t)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: t(5)
    call sortn(t, 5)
  end subroutine sort5

  Subroutine sortn(t, m)
    Implicit none
    Integer, intent(in) :: m
    Real (Kind=Kind(0.d0)), intent(inout) :: t(*)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: x
    do i = 2, m
       x = t(i); j = i-1
       do while (j >= 1)
          if (t(j) <= x) exit
          t(j+1) = t(j); j = j-1
       enddo
       t(j+1) = x
    enddo
  end subroutine sortn

End Module gpu_tau_mod
#endif

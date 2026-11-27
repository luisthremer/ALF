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
!> GPU port, P2.4 step 1e: equal-time and time-displaced observables of the Hubbard model on the device.
!> Compiled only with -DALF_GPU.
!>
!> @details
!> Replaces ham%Obser (and, through gpu_obs_tau_measure, ham%OBSERT) of Hamiltonian_Hubbard_smod for the supported case:
!> N_FL = 1, or N_FL = 2 with the Obs_eq layout of the model (three or five observables), no time-dependent Symm, device
!> stabilisation active (G_d is authoritative).  Everything else keeps the host Obser; gpu_obs_report writes one info line
!> that says which (DEVICE or HOST with the reason).
!>
!> Scalars (Kin, Pot, Part, Ener): device reductions, three complex numbers back per call (48 B); the host accumulates them in
!> Obs_scal as Obser does.  Correlations (Green, SpinZ/SpinXY/SpinT or SpinSUN, Den): for every output element (imj, no_I, no_J)
!> of Obs_Latt one gather over the lattice sites I (the partner J at displacement imj is unique, table o_jc), accumulated on the
!> device in acc / acc0 (shape of Obs_Latt / Obs_Latt0); the sum over I is a team reduction, so the summation order differs
!> from the host at the 1e-13 level.  Local observable (ham_h0 /= 0, N_FL = 2): gather kernel into acc_loc(N,Norb).
!> gpu_obs_flush adds the accumulators into the host arrays once per bin (before Pr_obs) and zeroes them.
!>
!> Weights of the interaction energy: U on the orbitals no <= Norb/2 and U2 above, for the bilayer lattices only (the rule of
!> Obser, KNOWN_ISSUES H7); U everywhere else.
!> Symm: the measurement is made on a symmetrised copy of G (A G B, built once from the hopping elements by gpu_symm_build);
!> G_d itself is not changed.  Static hopping only: a time-dependent Symm half step declines the device observables.
!>
!> Time-displaced observables (Ltau = 1): gpu_obs_tau_measure works on the device matrices GT0, G0T, G00, GTT of the tau pass
!> (gpu_tau_mod), accumulated in tacc / tacc0, also with Symm.
!>
!> Negative control ALF_GPU_FAULT=obs: the device weight of U is multiplied by 1 + 1e-3.
!> Debug switch ALF_GPU_NOOBS=1 keeps the host Obser (A/B tests).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_obs_mod
  Use iso_fortran_env, only: output_unit
  Use Natural_Constants, only: Eps_small
  Use Lattices_v3, only: Lattice, Unit_cell
  Use Predefined_Hoppings, only: Hopping_Matrix_type, inquire_hop, Generic_hopping
  Use Hamiltonian_main, only: Obs_scal, Obs_eq, Obs_tau, Obs_local
  Use linalg_mod, only: la_zgemm
#if defined(ALF_LA_HIP)
  Use alf_hipk_mod, only: hk_enabled, hk_corr, hk_scalars
#endif
  Use gpu_state_mod, only: gpu_active, G_d, g_dev_ok, td_symm_active, info_line, gpu_env_is
#if defined(ALF_LA_ONEMKL)
  Use iso_c_binding, only: c_ptr, c_int, c_double, c_loc
  Use iso_fortran_env, only: error_unit
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use linalg_onemkl_mod, only: om_get_queue
#endif
  Use gpu_stab_mod, only: gpu_stab_on
  Use gpu_wrap_mod, only: gpu_symm_build
  Implicit none
  Private
  Public :: gpu_obs_setup, gpu_obs_report, gpu_obs_active, gpu_obs_measure, gpu_obs_flush, gpu_obs_finalize
  Public :: gpu_obs_tau_active, gpu_obs_tau_measure

  Logical, Save :: o_set = .false.         ! model data stored (gpu_obs_setup called)
  Logical, Save :: o_ok = .false.          ! the model is supported
  Logical, Save :: o_active = .false.      ! decided once by gpu_obs_report
  Logical, Save :: o_dev = .false.         ! device arrays allocated
  Logical, Save :: o_off = .false.         ! ALF_GPU_NOOBS=1
  Logical, Save :: o_fault = .false.       ! ALF_GPU_FAULT=obs
  Logical, Save :: o_symm = .false., o_h0 = .false.
  Logical, Save :: o_tau = .false., o_tdev = .false.     ! tau-resolved observables present / device arrays allocated
  Logical, Save :: symm_built = .false.
  Character (len=200), Save :: o_reason = 'gpu_obs_setup not called'
#if defined(ALF_LA_ONEMKL)
  Logical, Save :: o_corr_sycl = .false.   ! ALFSH_OBSCORR=1: corr as a SYCL kernel (libalfsb alfsb_corr) on the OpenMP interop queue
  Logical, Parameter :: OBSCORR_DEFAULT = .true.    ! when ALFSH_OBSCORR is not set (an explicit 0|1 always wins); on since gpu-intel-speed (Q4b)
  Interface
    Integer (c_int) Function alfsb_corr(q, kind, k, n, norb, nfl, ns, sq3, wr, wi, ld, g, jc, inv, acc, acc0) Bind(C, name='alfsb_corr')
      Import :: c_ptr, c_int, c_double
      Type (c_ptr), Value :: q, g, jc, inv, acc, acc0
      Integer (c_int), Value :: kind, k, n, norb, nfl, ld
      Real (c_double), Value :: ns, sq3, wr, wi
    End Function
  End Interface
#endif
  Integer, Save :: o_N = 0, o_Norb = 0, o_NFL = 0, o_NSUN = 0, o_Ndim = 0, o_nk = 0, o_nobs = 0, o_dec = 1, o_Nt = 0
  Integer, Save :: o_kind(5) = 0           ! per Obs_eq: 1 Green, 2 SpinSUN, 3 Den, 4 SpinMz (fills Obs_eq(k:k+2))
  Integer, Allocatable, Save :: o_inv(:,:), o_jc(:,:), hk_i(:), hk_j(:), hk_f(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: hk_c(:)
  Real (Kind=Kind(0.d0)), Allocatable, Save :: o_uo(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: acc(:,:,:,:,:)    ! (N,1,Norb,Norb,nobs)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: acc0(:,:)        ! (Norb,nobs)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: acc_loc(:,:)     ! (N,Norb) local observable (ham_h0 /= 0)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: tacc(:,:,:,:,:)  ! (N,Nt,Norb,Norb,nobs)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: tacc0(:,:)       ! (Norb,nobs)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: Asy(:,:,:), Bsy(:,:,:), Tmp(:,:), Gw(:,:,:), Ts(:,:,:,:)

Contains

!--------------------------------------------------------------------
!> @brief
!> Stores the model data of the Hubbard model (host only; the device arrays are made at the first measurement).
!> @details
!> Called from Alloc_obs (needs Obs_eq, Obs_tau, Obs_local allocated before).  Sets o_ok and, if the model is declined,
!> o_reason.  May be called again (the previous data are dropped).
!> @param[in] bilayer .true. for BILAYER_SQUARE / BILAYER_HONEYCOMB (U2 on the upper orbitals)
!> @param[in] h0 the pinning field ham_h0 (local observable if /= 0)
!> @param[in] Symm_ the model symmetrises G before measuring
!--------------------------------------------------------------------
  Subroutine gpu_obs_setup(Latt, Latt_unit, List, Invlist, hop, N_SUN, N_FL, ham_U, ham_U2, bilayer, h0, Symm_)
    Implicit none
    Type (Lattice), intent(in) :: Latt
    Type (Unit_cell), intent(in) :: Latt_unit
    Integer, intent(in) :: List(:,:), Invlist(:,:), N_SUN, N_FL
    Type (Hopping_Matrix_type), allocatable, intent(in) :: hop(:)
    Real (Kind=Kind(0.d0)), intent(in) :: ham_U, ham_U2, h0
    Logical, intent(in) :: bilayer, Symm_
    Integer :: I, J, no_I, nf, Nb, n_1, n_2, nz, k, mode, N_Phi, pass
    Real (Kind=Kind(0.d0)) :: Phi_X, Phi_Y
    Logical :: Bulk
    Complex (Kind=Kind(0.d0)) :: Z
    Integer, Allocatable :: cnt(:,:)
    Character (len=16) :: s
#if defined(ALF_LA_ONEMKL)
    Integer :: st_env
#endif

    call gpu_obs_finalize()
    o_set = .true.; o_ok = .false.; o_active = .false.; o_tau = .false.
    call get_environment_variable('ALF_GPU_NOOBS', s)
    o_off = (trim(s) == '1')
    o_fault = gpu_env_is('ALF_GPU_FAULT', 'obs')
#if defined(ALF_LA_ONEMKL)
    ! ALFSH_OBSCORR=1: the correlation sums as one SYCL work-group per output element (alfsb_corr, the SYCL port of the HIP
    ! k_corr); the summation order over I differs from the OpenMP kernels (rounding level)
    s = ''
    call get_environment_variable('ALFSH_OBSCORR', s, status=st_env)
    o_corr_sycl = OBSCORR_DEFAULT
    if (st_env == 0) o_corr_sycl = (trim(s) == '1')
#endif
    o_symm = Symm_
    o_h0 = (abs(h0) >= Eps_small)
    o_reason = 'N_FL is not 1 or 2'
    if (N_FL /= 1 .and. N_FL /= 2) return
    o_reason = 'Obs_eq is not allocated'
    if (.not. allocated(Obs_eq)) return
    o_reason = 'Obs_eq layout is not the one of the Hubbard model'
    if (size(Obs_eq,1) /= merge(5, 3, N_FL == 2)) return
    o_reason = 'ham_h0 /= 0 needs N_FL = 2 and Obs_local'
    if (o_h0) then
       if (N_FL /= 2) return
       if (.not. allocated(Obs_local)) return
       if (size(Obs_local,1) /= 1) return
    endif
    o_N = Latt%N; o_Norb = Latt_unit%Norb; o_NFL = N_FL; o_NSUN = N_SUN; o_Ndim = size(List,1)
    o_reason = 'Ndim /= N*Norb'
    if (o_Ndim /= o_N*o_Norb) return
    allocate(o_inv(o_N,o_Norb), o_jc(o_N,o_N), o_uo(o_Norb))
    o_inv = Invlist(1:o_N, 1:o_Norb)
    ! partner table: o_jc(I, imj) = the J with imj(I,J) = imj; it must be a bijection per I
    allocate(cnt(o_N,o_N)); cnt = 0; o_jc = 0
    o_reason = 'imj is not a bijection per site'
    do I = 1, o_N
       do J = 1, o_N
          k = Latt%imj(I,J)
          if (k < 1 .or. k > o_N) then; deallocate(o_inv, o_jc, o_uo, cnt); return; endif
          o_jc(I,k) = J; cnt(I,k) = cnt(I,k) + 1
       enddo
    enddo
    if (any(cnt /= 1)) then; deallocate(o_inv, o_jc, o_uo, cnt); return; endif
    deallocate(cnt)
    ! weights of the interaction energy, as in Obser
    o_uo = ham_U
    if (bilayer) then
       do no_I = o_Norb/2 + 1, o_Norb
          o_uo(no_I) = ham_U2
       enddo
    endif
    o_dec = 1; if (N_FL == 2) o_dec = 2
    ! kinetic energy as a sparse linear functional of GRC: Z_kin = sum_k hk_c(k) GRC(hk_i(k), hk_j(k), hk_f(k)); the code of
    ! Predefined_Hoppings_Compute_Kin (no pinning): mode 1 = T_Loc only, else bonds + T_Loc
    mode = inquire_hop(hop)
    nz = 0
    do pass = 1, 2
       k = 0
       if (mode == 1) then
          do nf = 1, N_FL
             do I = 1, o_N
                do no_I = 1, o_Norb
                   k = k + 1
                   if (pass == 2) then
                      hk_i(k) = o_inv(I,no_I); hk_j(k) = o_inv(I,no_I); hk_f(k) = nf; hk_c(k) = hop(nf)%T_Loc(no_I)
                   endif
                enddo
             enddo
          enddo
       elseif (mode /= 0) then
          do nf = 1, N_FL
             N_Phi = hop(nf)%N_Phi; Phi_X = hop(nf)%Phi_X; Phi_Y = hop(nf)%Phi_Y; Bulk = hop(nf)%Bulk
             do I = 1, o_N
                do Nb = 1, hop(nf)%N_bonds
                   no_I = hop(nf)%list(Nb,1)
                   n_1 = hop(nf)%list(Nb,3); n_2 = hop(nf)%list(Nb,4)
                   J = Latt%nnlist(I,n_1,n_2)
                   Z = hop(nf)%T(Nb) * Generic_hopping(I, no_I, n_1, n_2, hop(nf)%list(Nb,2), N_Phi, Phi_X, Phi_Y, Bulk, Latt, Latt_unit)
                   k = k + 1
                   if (pass == 2) then
                      hk_i(k) = o_inv(I,no_I); hk_j(k) = o_inv(J,hop(nf)%list(Nb,2)); hk_f(k) = nf; hk_c(k) = Z
                   endif
                   k = k + 1
                   if (pass == 2) then
                      hk_i(k) = o_inv(J,hop(nf)%list(Nb,2)); hk_j(k) = o_inv(I,no_I); hk_f(k) = nf; hk_c(k) = conjg(Z)
                   endif
                enddo
                do no_I = 1, o_Norb
                   k = k + 1
                   if (pass == 2) then
                      hk_i(k) = o_inv(I,no_I); hk_j(k) = o_inv(I,no_I); hk_f(k) = nf; hk_c(k) = hop(nf)%T_Loc(no_I)
                   endif
                enddo
             enddo
          enddo
       endif
       if (pass == 1) then
          nz = k
          allocate(hk_i(max(nz,1)), hk_j(max(nz,1)), hk_f(max(nz,1)), hk_c(max(nz,1)))
       endif
    enddo
    o_nk = nz
    if (N_FL == 2) then
       o_nobs = 5; o_kind = (/ 1, 4, 4, 4, 3 /)
    else
       o_nobs = 3; o_kind(1:3) = (/ 1, 2, 3 /)
    endif
    o_reason = ''
    o_ok = .true.
    ! tau-resolved observables (allocated by Alloc_obs only for Ltau = 1): the same observables, one per Obs_eq
    if (allocated(Obs_tau)) then
       if (size(Obs_tau,1) == o_nobs) then
          o_Nt = size(Obs_tau(1)%Obs_Latt,2)
          o_tau = .true.
       endif
    endif
  End Subroutine gpu_obs_setup

!--------------------------------------------------------------------
!> @brief
!> Decides once whether the device observables are used and writes the one info line about it.
!> @details
!> Called at the end of Setup_gpu (device stabilisation set up, so gpu_stab_on is final).  Device observables need the
!> device stabilisation (G_d is then authoritative) and no time-dependent Symm half step.
!--------------------------------------------------------------------
  Subroutine gpu_obs_report()
    Implicit none
    Character (len=400) :: line
    o_active = .false.
    if (.not. o_set) then
       line = 'GPU observables (step 1e): HOST, reason: model has no device observables (Hubbard only)'
    elseif (.not. o_ok) then
       line = 'GPU observables (step 1e): HOST, reason: '//trim(o_reason)
    elseif (o_off) then
       line = 'GPU observables (step 1e): HOST, reason: ALF_GPU_NOOBS=1 (debug)'
    elseif (.not. gpu_active .or. .not. gpu_stab_on) then
       line = 'GPU observables (step 1e): HOST, reason: device stabilisation inactive'
    elseif (td_symm_active .and. o_symm) then
       line = 'GPU observables (step 1e): HOST, reason: time-dependent Symm half step'
    else
       o_active = .true.
       write(line,'(A,I0,A,I0,A,L1,A,L1,A,L1,A)') 'GPU observables (step 1e): DEVICE (N_FL ', o_NFL, ', N_SUN ', o_NSUN, &
            & ', Symm ', o_symm, ', local h0 ', o_h0, ', tau-resolved ', o_tau, ')'
#if defined(ALF_LA_ONEMKL)
       if (o_corr_sycl) then
          line = trim(line)//', ALFSH_OBSCORR=1 (SYCL correlation kernel)'
       else
          line = trim(line)//', ALFSH_OBSCORR=0'
       endif
#endif
    endif
    call info_line(trim(line))
  End Subroutine gpu_obs_report

!> .true. when the device measurement replaces ham%Obser
  Logical Function gpu_obs_active()
    Implicit none
    gpu_obs_active = o_set .and. o_ok .and. o_active .and. gpu_stab_on
  End Function gpu_obs_active

!> .true. when the device tau-resolved measurement replaces ham%OBSERT (checked by gpu_tau_mod)
  Logical Function gpu_obs_tau_active()
    Implicit none
    gpu_obs_tau_active = gpu_obs_active() .and. o_tau
  End Function gpu_obs_tau_active

!> out = Symm(inn) = A inn B (A, B built once from the hopping elements); inn and out must be different arrays
  Subroutine symm_apply(nd, inn, out, nf)
    Implicit none
    Integer, intent(in) :: nd, nf
    Complex (Kind=Kind(0.d0)), intent(in) :: inn(nd,nd)
    Complex (Kind=Kind(0.d0)), intent(inout) :: out(nd,nd)
    Complex (Kind=Kind(0.d0)), parameter :: cone = (1.d0,0.d0), czero = (0.d0,0.d0)
    Integer :: kb
    if (.not. symm_built) then
       allocate(Asy(nd,nd,o_NFL), Bsy(nd,nd,o_NFL), Tmp(nd,nd))
       !$omp target enter data map(alloc: Asy, Bsy, Tmp)
       do kb = 1, o_NFL
          call gpu_symm_build(Asy(1,1,kb), Bsy(1,1,kb), nd, nd, kb)
       enddo
       symm_built = .true.
    endif
    call la_zgemm('N', 'N', nd, nd, nd, cone, Asy(1,1,nf), nd, inn, nd, czero, Tmp, nd)
    call la_zgemm('N', 'N', nd, nd, nd, cone, Tmp, nd, Bsy(1,1,nf), nd, czero, out, nd)
  End Subroutine symm_apply

  Subroutine alloc_dev()
    Implicit none
    allocate(acc(o_N,1,o_Norb,o_Norb,o_nobs), acc0(o_Norb,o_nobs), acc_loc(o_N,o_Norb))
    acc = (0.d0, 0.d0); acc0 = (0.d0, 0.d0); acc_loc = (0.d0, 0.d0)
    if (o_fault) o_uo(1) = o_uo(1)*(1.d0 + 1.d-3)       ! negative control ALF_GPU_FAULT=obs
    ! one-time upload of the constant tables and the zeroed accumulators
    !$omp target enter data map(to: o_inv, o_jc, o_uo, hk_i, hk_j, hk_f, hk_c, acc, acc0, acc_loc)
    o_dev = .true.
  End Subroutine alloc_dev

!--------------------------------------------------------------------
!> @brief
!> ham%Obser on the device.
!> @details
!> Reads G_d (device-authoritative with the device stabilisation); the host GR is copied to G_d only if the device copy of a
!> flavour is stale.  Ntau is not used (static Symm only).
!--------------------------------------------------------------------
  Subroutine gpu_obs_measure(GR, Phase, Ntau, Mc_step_weight)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(o_Ndim,o_Ndim,o_NFL)
    Complex (Kind=Kind(0.d0)), intent(in) :: Phase
    Integer, intent(in) :: Ntau
    Real (Kind=Kind(0.d0)), intent(in) :: Mc_step_weight
    Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZW
    Integer :: I, nf

    if (.not. o_dev) call alloc_dev()
    do nf = 1, o_NFL
       if (.not. g_dev_ok(nf)) then
          ! H->D: the device copy of G is stale (a host writer ran since the last device update); once, not per slice
          G_d(:,:,nf) = GR(:,:,nf)
          !$omp target update to(G_d(:,:,nf))
          g_dev_ok(nf) = .true.
       endif
    enddo
    ZP = Phase/Real(Phase, kind(0.D0))
    ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
    ZS = ZS*Mc_step_weight
    ZW = ZP*ZS
    do I = 1, size(Obs_scal,1)
       Obs_scal(I)%N = Obs_scal(I)%N + 1
       Obs_scal(I)%Ave_sign = Obs_scal(I)%Ave_sign + Real(ZS, kind(0.d0))
    enddo
    do I = 1, size(Obs_eq,1)
       Obs_eq(I)%N = Obs_eq(I)%N + 1
       Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + Real(ZS, kind(0.d0))
    enddo
    if (o_h0) then
       Obs_local(1)%N = Obs_local(1)%N + 1
       Obs_local(1)%Ave_sign = Obs_local(1)%Ave_sign + Real(ZS, kind(0.d0))
    endif
    if (o_symm) then                  ! measure on Symm(G), G_d stays as it is
       if (.not. allocated(Gw)) then
          allocate(Gw(o_Ndim,o_Ndim,o_NFL))
          !$omp target enter data map(alloc: Gw)
       endif
       do nf = 1, o_NFL
          call symm_apply(o_Ndim, G_d(1,1,nf), Gw(1,1,nf), nf)
       enddo
       call measure_on(o_Ndim, o_NFL, Gw, ZP, ZS, ZW)
    else
       call measure_on(o_Ndim, o_NFL, G_d, ZP, ZS, ZW)
    endif
  End Subroutine gpu_obs_measure

  !> scalars, local observable and correlations of the equal-time matrix Gm
  Subroutine measure_on(nd, nfl, Gm, ZP, ZS, ZW)
    Implicit none
    Integer, intent(in) :: nd, nfl
    Complex (Kind=Kind(0.d0)), intent(in) :: Gm(nd,nd,nfl)
    Complex (Kind=Kind(0.d0)), intent(in) :: ZP, ZS, ZW
    Complex (Kind=Kind(0.d0)) :: Zkin, Zpot, Zrho, zc, zt
    Real (Kind=Kind(0.d0)) :: wr, wi, sq3
    Integer :: I, nf, k, no, i1, J1, N, Norb, NFL_, NSUN, dec, kind_k
#if defined(ALF_LA_HIP)
    Real (Kind=Kind(0.d0)) :: hsc(6)
#endif

    N = o_N; Norb = o_Norb; NFL_ = o_NFL; NSUN = o_NSUN; dec = o_dec; sq3 = sqrt(3.d0)
    wr = real(ZW); wi = aimag(ZW)
    ! ---- scalars: GRC(i,j,nf) = delta_ij - G(j,i,nf)
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_scalars(o_nk, hk_i, hk_j, hk_f, hk_c, N, Norb, NFL_, dec, nd, Gm, o_inv, o_uo, hsc)
       Zkin = cmplx(hsc(1), hsc(2), kind(0.d0))*dble(NSUN)
       Zpot = cmplx(hsc(3), hsc(4), kind(0.d0))
       Zrho = cmplx(hsc(5), hsc(6), kind(0.d0))
    else
#endif
    Zkin = (0.d0, 0.d0)
    !$omp target teams distribute parallel do reduction(+:Zkin) private(i1,J1,nf,zc)
    do k = 1, o_nk
       i1 = hk_i(k); J1 = hk_j(k); nf = hk_f(k)
       zc = -Gm(J1,i1,nf)
       if (i1 == J1) zc = zc + (1.d0, 0.d0)
       Zkin = Zkin + hk_c(k)*zc
    enddo
    Zkin = Zkin*dble(NSUN)
    Zpot = (0.d0, 0.d0)
    !$omp target teams distribute parallel do collapse(2) reduction(+:Zpot) private(i1,zc,zt)
    do I = 1, N
       do no = 1, Norb
          i1 = o_inv(I,no)
          zc = (1.d0, 0.d0) - Gm(i1,i1,1)
          zt = (1.d0, 0.d0) - Gm(i1,i1,dec)
          Zpot = Zpot + zc*zt*o_uo(no)
       enddo
    enddo
    Zrho = (0.d0, 0.d0)
    !$omp target teams distribute parallel do collapse(3) reduction(+:Zrho) private(i1)
    do nf = 1, NFL_
       do I = 1, N
          do no = 1, Norb
             i1 = o_inv(I,no)
             Zrho = Zrho + ((1.d0, 0.d0) - Gm(i1,i1,nf))
          enddo
       enddo
    enddo
    ! D->H: three complex scalars per call (48 B) enter the host sums of Obs_scal
#if defined(ALF_LA_HIP)
    endif
#endif
    Zrho = Zrho*dble(NSUN)
    Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*ZP*ZS
    Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Zpot*ZP*ZS
    Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zrho*ZP*ZS
    Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

    if (o_h0) then
       call k_local(N, Norb, nd, nfl, dec, wr, wi, Gm)
    endif
    ! ---- correlations: one gather per output element (imj, no, no2)
    do k = 1, o_nobs
       kind_k = o_kind(k)
       if (kind_k == 4 .and. k /= 2) cycle        ! SpinMz fills Obs_eq(2:4) from k = 2
       if (kind_k == 4 .and. NSUN /= 1) cycle     ! host Predefined_Obs_eq_SpinMz_measure measures only N_SUN = 1 (it just counts otherwise)
       call k_corr(k, kind_k, N, Norb, NFL_, dble(NSUN), sq3, wr, wi, nd, Gm)
    enddo
  End Subroutine measure_on

  !> local observable Obs_local(1) = (GRC(i,i,1) - GRC(i,i,dec)) ZP ZS, gather into acc_loc
  Subroutine k_local(N, Norb, nd, nfl, dec, wr, wi, Gm)
    Implicit none
    Integer, intent(in) :: N, Norb, nd, nfl, dec
    Real (Kind=Kind(0.d0)), intent(in) :: wr, wi
    Complex (Kind=Kind(0.d0)), intent(in) :: Gm(nd,nd,nfl)
    Integer :: I, no, i1
    Complex (Kind=Kind(0.d0)) :: zc, zt
    !$omp target teams distribute parallel do collapse(2) private(i1,zc,zt)
    do no = 1, Norb
       do I = 1, N
          i1 = o_inv(I,no)
          zc = (1.d0, 0.d0) - Gm(i1,i1,1)
          zt = (1.d0, 0.d0) - Gm(i1,i1,dec)
          acc_loc(I,no) = acc_loc(I,no) + (zc - zt)*cmplx(wr,wi,kind(0.d0))
       enddo
    enddo
  End Subroutine k_local

#if defined(ALF_LA_ONEMKL)
  !> corr on the SYCL queue of the OpenMP offload runtime (libalfsb alfsb_corr; ported from fork/intel-batched ee682f09):
  !> device addresses by use_device_addr, the same arguments as hk_corr of alf_hipk_mod
  Subroutine sb_corr(kind_k, k, n, norb, nfl, ns, sq3, wr, wi, ld, g, jc, inv, acc_, acc0_)
    Implicit none
    Integer, intent(in) :: kind_k, k, n, norb, nfl, ld
    Real (Kind=Kind(0.d0)), intent(in) :: ns, sq3, wr, wi
    Complex (Kind=Kind(0.d0)), Target :: g(*), acc_(*), acc0_(*)
    Integer, Target :: jc(*), inv(*)
    Integer (c_int) :: rc
    Type (c_ptr) :: qh
    qh = om_get_queue()
    !$omp target data use_device_addr(g, jc, inv, acc_, acc0_)
    rc = alfsb_corr(qh, int(kind_k,c_int), int(k,c_int), int(n,c_int), int(norb,c_int), int(nfl,c_int), ns, sq3, wr, wi, &
         int(ld,c_int), c_loc(g(1)), c_loc(jc(1)), c_loc(inv(1)), c_loc(acc_(1)), c_loc(acc0_(1)))
    !$omp end target data
    if (rc /= 0) then
       write(error_unit,*) 'gpu_obs: alfsb_corr (ALFSH_OBSCORR=1) failed, rc = ', rc
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    endif
  End Subroutine sb_corr
#endif

  !> equal-time correlators of kind kind_k into acc(:,:,:,:,k..)
  Subroutine k_corr(k, kind_k, N, Norb, NFL, ns, sq3, wr, wi, nd, Gm)
    Implicit none
    Integer, intent(in) :: k, kind_k, N, Norb, NFL, nd
    Real (Kind=Kind(0.d0)), intent(in) :: ns, sq3, wr, wi
    Complex (Kind=Kind(0.d0)), intent(in) :: Gm(nd,nd,NFL)
    Integer :: m, imj, no, no2, I, J, i1, j1, nf
    Complex (Kind=Kind(0.d0)) :: s1, s2, s3, z, zi, zj, zxy, zz, a, b, zs0

#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_corr(kind_k, k, N, Norb, NFL, ns, sq3, wr, wi, nd, Gm, o_jc, o_inv, acc, acc0)
       return
    endif
#endif
#if defined(ALF_LA_ONEMKL)
    if (o_corr_sycl) then
       call sb_corr(kind_k, k, N, Norb, NFL, ns, sq3, wr, wi, nd, Gm, o_jc, o_inv, acc, acc0)
       return
    endif
#endif
    if (kind_k == 1) then            ! Green
       !$omp target teams distribute private(imj,no,no2,s1)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,nf,z,a) reduction(+:s1)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             z = (0.d0, 0.d0)
             do nf = 1, NFL
                a = -Gm(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                z = z + a
             enddo
             s1 = s1 + z*ns
          enddo
          !$omp end parallel do
          acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
    elseif (kind_k == 2) then        ! SpinSUN (N_FL = 1)
       !$omp target teams distribute private(imj,no,no2,s1)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,a) reduction(+:s1)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             a = -Gm(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
             s1 = s1 + a*Gm(i1,j1,1)*ns
          enddo
          !$omp end parallel do
          acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
    elseif (kind_k == 3) then        ! Den
       !$omp target teams distribute private(imj,no,no2,s1)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,nf,z,zi,zj,a) reduction(+:s1)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             zi = (0.d0, 0.d0); zj = (0.d0, 0.d0); z = (0.d0, 0.d0)
             do nf = 1, NFL
                zi = zi + ((1.d0, 0.d0) - Gm(i1,i1,nf))
                zj = zj + ((1.d0, 0.d0) - Gm(j1,j1,nf))
                a = -Gm(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                z = z + a*Gm(i1,j1,nf)
             enddo
             s1 = s1 + (zi*ns)*(zj*ns) + z*ns
          enddo
          !$omp end parallel do
          acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute private(I,i1,zi,nf,zs0)
       do no = 1, Norb
          zs0 = (0.d0, 0.d0)
          !$omp parallel do private(I,i1,zi,nf) reduction(+:zs0)
          do I = 1, N
             i1 = o_inv(I,no); zi = (0.d0, 0.d0)
             do nf = 1, NFL
                zi = zi + ((1.d0, 0.d0) - Gm(i1,i1,nf))
             enddo
             zs0 = zs0 + zi*ns
          enddo
          !$omp end parallel do
          acc0(no,k) = acc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
       enddo
    else                             ! SpinMz (N_FL = 2, N_SUN = 1): Obs_eq(k) = SpinZ, k+1 = SpinXY, k+2 = SpinT
       !$omp target teams distribute private(imj,no,no2,s1,s2,s3)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,zxy,zz,a,b,zi,zj) reduction(+:s1,s2,s3)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             a = -Gm(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)      ! GRC(I1,J1,1)
             b = -Gm(j1,i1,2); if (i1 == j1) b = b + (1.d0, 0.d0)      ! GRC(I1,J1,2)
             zxy = a*Gm(i1,j1,2) + b*Gm(i1,j1,1)
             zi = ((1.d0, 0.d0) - Gm(i1,i1,2)) - ((1.d0, 0.d0) - Gm(i1,i1,1))     ! GRC(I1,I1,2) - GRC(I1,I1,1)
             zj = ((1.d0, 0.d0) - Gm(j1,j1,2)) - ((1.d0, 0.d0) - Gm(j1,j1,1))
             zz = a*Gm(i1,j1,1) + b*Gm(i1,j1,2) + zi*zj
             s1 = s1 + zz
             s2 = s2 + zxy
             s3 = s3 + (2.d0*zxy + zz)/3.d0
          enddo
          !$omp end parallel do
          acc(imj,1,no,no2,k)   = acc(imj,1,no,no2,k)   + s1*cmplx(wr,wi,kind(0.d0))
          acc(imj,1,no,no2,k+1) = acc(imj,1,no,no2,k+1) + s2*cmplx(wr,wi,kind(0.d0))
          acc(imj,1,no,no2,k+2) = acc(imj,1,no,no2,k+2) + s3*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute private(I,i1,zi,zs0)
       do no = 1, Norb
          zs0 = (0.d0, 0.d0)
          !$omp parallel do private(I,i1,zi) reduction(+:zs0)
          do I = 1, N
             i1 = o_inv(I,no)
             zi = ((1.d0, 0.d0) - Gm(i1,i1,2)) - ((1.d0, 0.d0) - Gm(i1,i1,1))
             zs0 = zs0 + zi
          enddo
          !$omp end parallel do
          acc0(no,k)   = acc0(no,k)   + zs0*cmplx(wr,wi,kind(0.d0))
          acc0(no,k+2) = acc0(no,k+2) + zs0*cmplx(wr,wi,kind(0.d0))/sq3
       enddo
    endif
  End Subroutine k_corr

!--------------------------------------------------------------------
!> @brief
!> Adds the device accumulators to the host Obs_* (before Pr_obs normalises them) and zeroes them.
!--------------------------------------------------------------------
  Subroutine gpu_obs_flush()
    Implicit none
    Integer :: k
    if (o_tdev) then
       ! D->H once per bin: Pr_obs and the bin files read the host arrays
       !$omp target update from(tacc, tacc0)
       do k = 1, o_nobs
          Obs_tau(k)%Obs_Latt = Obs_tau(k)%Obs_Latt + tacc(:,:,:,:,k)
          Obs_tau(k)%Obs_Latt0 = Obs_tau(k)%Obs_Latt0 + tacc0(:,k)
       enddo
       tacc = (0.d0, 0.d0); tacc0 = (0.d0, 0.d0)
       ! H->D: zero the device accumulators again
       !$omp target update to(tacc, tacc0)
    endif
    if (.not. o_dev) return
    ! D->H once per bin: Pr_obs and the bin files read the host arrays
    !$omp target update from(acc, acc0, acc_loc)
    do k = 1, o_nobs
       Obs_eq(k)%Obs_Latt = Obs_eq(k)%Obs_Latt + acc(:,:,:,:,k)
       Obs_eq(k)%Obs_Latt0 = Obs_eq(k)%Obs_Latt0 + acc0(:,k)
    enddo
    if (o_h0) Obs_local(1)%Obs_Latt(:,1,:) = Obs_local(1)%Obs_Latt(:,1,:) + acc_loc
    acc = (0.d0, 0.d0); acc0 = (0.d0, 0.d0); acc_loc = (0.d0, 0.d0)
    ! H->D: zero the device accumulators again
    !$omp target update to(acc, acc0, acc_loc)
  End Subroutine gpu_obs_flush

!> @brief Releases the device arrays and the host tables.
  Subroutine gpu_obs_finalize()
    Implicit none
    if (allocated(Gw)) then
       !$omp target exit data map(delete: Gw)
       deallocate(Gw)
    endif
    if (allocated(Ts)) then
       !$omp target exit data map(delete: Ts)
       deallocate(Ts)
    endif
    if (symm_built) then
       !$omp target exit data map(delete: Asy, Bsy, Tmp)
       deallocate(Asy, Bsy, Tmp)
       symm_built = .false.
    endif
    if (o_tdev) then
       !$omp target exit data map(delete: tacc, tacc0)
       deallocate(tacc, tacc0)
       o_tdev = .false.
    endif
    if (o_dev) then
       !$omp target exit data map(delete: o_inv, o_jc, o_uo, hk_i, hk_j, hk_f, hk_c, acc, acc0, acc_loc)
       deallocate(acc, acc0, acc_loc)
       o_dev = .false.
    endif
    if (allocated(o_inv)) deallocate(o_inv)
    if (allocated(o_jc)) deallocate(o_jc)
    if (allocated(o_uo)) deallocate(o_uo)
    if (allocated(hk_i)) deallocate(hk_i, hk_j, hk_f, hk_c)
    o_set = .false.; o_ok = .false.; o_active = .false.; o_tau = .false.
  End Subroutine gpu_obs_finalize

!--------------------------------------------------------------------
!> @brief
!> ham%OBSERT on the device: GT0, G0T, G00, GTT are the device-resident matrices of the tau pass (no D->H per slice).
!> @details
!> nt = 0 counts the measurement once per pass, as the Predefined_Obs_tau_* routines.  With Symm the four matrices are
!> symmetrised into Ts first.
!--------------------------------------------------------------------
  Subroutine gpu_obs_tau_measure(nt, GT0, G0T, G00, GTT, Phase, Mc_step_weight)
    Implicit none
    Integer, intent(in) :: nt
    Complex (Kind=Kind(0.d0)), intent(in) :: GT0(o_Ndim,o_Ndim,o_NFL), G0T(o_Ndim,o_Ndim,o_NFL), G00(o_Ndim,o_Ndim,o_NFL), GTT(o_Ndim,o_Ndim,o_NFL)
    Complex (Kind=Kind(0.d0)), intent(in) :: Phase
    Real (Kind=Kind(0.d0)), intent(in) :: Mc_step_weight
    Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZW
    Integer :: I, k, nf
    Real (Kind=Kind(0.d0)) :: rs, wr, wi

    if (.not. o_tdev) then
       allocate(tacc(o_N,o_Nt,o_Norb,o_Norb,o_nobs), tacc0(o_Norb,o_nobs))
       tacc = (0.d0, 0.d0); tacc0 = (0.d0, 0.d0)
       if (.not. o_dev) call alloc_dev()
       ! one-time upload of the zeroed accumulators
       !$omp target enter data map(to: tacc, tacc0)
       o_tdev = .true.
    endif
    ZP = Phase/Real(Phase, kind(0.D0))
    ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
    ZS = ZS*Mc_step_weight
    ZW = ZP*ZS
    wr = real(ZW); wi = aimag(ZW)
    rs = Real(ZS, kind(0.d0))
    if (nt == 0) then
       do I = 1, o_nobs
          Obs_tau(I)%N = Obs_tau(I)%N + 1
          Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + rs
       enddo
    endif
    if (o_symm) then
       if (.not. allocated(Ts)) then
          allocate(Ts(o_Ndim,o_Ndim,o_NFL,4))
          !$omp target enter data map(alloc: Ts)
       endif
       do nf = 1, o_NFL
          call symm_apply(o_Ndim, GT0(1,1,nf), Ts(1,1,nf,1), nf)
          call symm_apply(o_Ndim, G0T(1,1,nf), Ts(1,1,nf,2), nf)
          call symm_apply(o_Ndim, G00(1,1,nf), Ts(1,1,nf,3), nf)
          call symm_apply(o_Ndim, GTT(1,1,nf), Ts(1,1,nf,4), nf)
       enddo
       do k = 1, o_nobs
          if (o_kind(k) == 4 .and. k /= 2) cycle
          call k_tcorr(k, o_kind(k), nt, wr, wi, o_N, o_Norb, o_NFL, o_NSUN, o_Ndim, &
               & Ts(1,1,1,1), Ts(1,1,1,2), Ts(1,1,1,3), Ts(1,1,1,4))
       enddo
    else
       do k = 1, o_nobs
          if (o_kind(k) == 4 .and. k /= 2) cycle
          call k_tcorr(k, o_kind(k), nt, wr, wi, o_N, o_Norb, o_NFL, o_NSUN, o_Ndim, GT0, G0T, G00, GTT)
       enddo
    endif
  End Subroutine gpu_obs_tau_measure

  !> time-displaced correlators of kind kind_k at slice nt into tacc(:,nt+1,:,:,k..)
  Subroutine k_tcorr(k, kind_k, nt, wr, wi, N, Norb, NFL, NSUN, nd, GT0, G0T, G00, GTT)
    Implicit none
    Integer, intent(in) :: k, kind_k, nt, N, Norb, NFL, NSUN, nd
    Real (Kind=Kind(0.d0)), intent(in) :: wr, wi
    Complex (Kind=Kind(0.d0)), intent(in) :: GT0(nd,nd,NFL), G0T(nd,nd,NFL), G00(nd,nd,NFL), GTT(nd,nd,NFL)
    Integer :: m, imj, no, no2, I, J, i1, j1, nf, it
    Complex (Kind=Kind(0.d0)) :: s1, s2, s3, z, zi, zj, zxy, zz, zs0
    Real (Kind=Kind(0.d0)) :: ns, rn, sq3

    it = nt + 1
    ns = dble(NSUN); rn = 1.d0/dble(NFL); sq3 = sqrt(3.d0)
    if (kind_k == 1) then            ! Green: sum_nf GT0 / N_FL
#if defined(ALF_LA_HIP)
       !$omp target teams distribute parallel do private(imj,no,no2,s1,I,J,i1,j1,nf,z)
#else
       !$omp target teams distribute private(imj,no,no2,s1)
#endif
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
#if !defined(ALF_LA_HIP)
          !$omp parallel do private(I,J,i1,j1,nf,z) reduction(+:s1)
#endif
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             z = (0.d0, 0.d0)
             do nf = 1, NFL
                z = z + GT0(i1,j1,nf)
             enddo
             s1 = s1 + z*rn
          enddo
#if !defined(ALF_LA_HIP)
          !$omp end parallel do
#endif
          tacc(imj,it,no,no2,k) = tacc(imj,it,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
    elseif (kind_k == 2) then        ! SpinSUN (N_FL = 1): - N_SUN G0T(J1,I1,1) GT0(I1,J1,1)
#if defined(ALF_LA_HIP)
       !$omp target teams distribute parallel do private(imj,no,no2,s1,I,J,i1,j1)
#else
       !$omp target teams distribute private(imj,no,no2,s1)
#endif
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
#if !defined(ALF_LA_HIP)
          !$omp parallel do private(I,J,i1,j1) reduction(+:s1)
#endif
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             s1 = s1 - ns*G0T(j1,i1,1)*GT0(i1,j1,1)
          enddo
#if !defined(ALF_LA_HIP)
          !$omp end parallel do
#endif
          tacc(imj,it,no,no2,k) = tacc(imj,it,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
    elseif (kind_k == 3) then        ! Den
#if defined(ALF_LA_HIP)
       !$omp target teams distribute parallel do private(imj,no,no2,s1,I,J,i1,j1,nf,z,zi,zj)
#else
       !$omp target teams distribute private(imj,no,no2,s1)
#endif
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
#if !defined(ALF_LA_HIP)
          !$omp parallel do private(I,J,i1,j1,nf,z,zi,zj) reduction(+:s1)
#endif
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             zi = (0.d0, 0.d0); zj = (0.d0, 0.d0); z = (0.d0, 0.d0)
             do nf = 1, NFL
                zi = zi + ((1.d0, 0.d0) - GTT(i1,i1,nf))
                zj = zj + ((1.d0, 0.d0) - G00(j1,j1,nf))
                z = z - G0T(j1,i1,nf)*GT0(i1,j1,nf)
             enddo
             s1 = s1 + (zi*ns)*(zj*ns) + z*ns
          enddo
#if !defined(ALF_LA_HIP)
          !$omp end parallel do
#endif
          tacc(imj,it,no,no2,k) = tacc(imj,it,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute parallel do private(I,i1,zi,nf,zs0)
       do no = 1, Norb
          zs0 = (0.d0, 0.d0)
          do I = 1, N
             i1 = o_inv(I,no); zi = (0.d0, 0.d0)
             do nf = 1, NFL
                zi = zi + ((1.d0, 0.d0) - GTT(i1,i1,nf))
             enddo
             zs0 = zs0 + zi*ns
          enddo
          tacc0(no,k) = tacc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
       enddo
    else                             ! SpinMz (N_FL = 2): k = SpinZ, k+1 = SpinXY, k+2 = SpinT
#if defined(ALF_LA_HIP)
       !$omp target teams distribute parallel do private(imj,no,no2,s1,s2,s3,I,J,i1,j1,zxy,zz)
#else
       !$omp target teams distribute private(imj,no,no2,s1,s2,s3)
#endif
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0)
#if !defined(ALF_LA_HIP)
          !$omp parallel do private(I,J,i1,j1,zxy,zz) reduction(+:s1,s2,s3)
#endif
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             zz  = (GTT(i1,i1,1) - GTT(i1,i1,2))*(G00(j1,j1,1) - G00(j1,j1,2)) &
                  &  - G0T(j1,i1,1)*GT0(i1,j1,1) - G0T(j1,i1,2)*GT0(i1,j1,2)
             zxy = - G0T(j1,i1,1)*GT0(i1,j1,2) - G0T(j1,i1,2)*GT0(i1,j1,1)
             s1 = s1 + zz
             s2 = s2 + zxy
             s3 = s3 + (2.d0*zxy + zz)/3.d0
          enddo
#if !defined(ALF_LA_HIP)
          !$omp end parallel do
#endif
          tacc(imj,it,no,no2,k)   = tacc(imj,it,no,no2,k)   + s1*cmplx(wr,wi,kind(0.d0))
          tacc(imj,it,no,no2,k+1) = tacc(imj,it,no,no2,k+1) + s2*cmplx(wr,wi,kind(0.d0))
          tacc(imj,it,no,no2,k+2) = tacc(imj,it,no,no2,k+2) + s3*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute parallel do private(I,i1,zs0)
       do no = 1, Norb
          zs0 = (0.d0, 0.d0)
          do I = 1, N
             i1 = o_inv(I,no)
             zs0 = zs0 + (GTT(i1,i1,2) - GTT(i1,i1,1))
          enddo
          tacc0(no,k)   = tacc0(no,k)   + zs0*cmplx(wr,wi,kind(0.d0))
          tacc0(no,k+2) = tacc0(no,k+2) + zs0*cmplx(wr,wi,kind(0.d0))/sq3
       enddo
    endif
  End Subroutine k_tcorr
End Module gpu_obs_mod
#endif

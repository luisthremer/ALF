!  Copyright (C) 2016 - 2022 The ALF project contributors
!
!  This file is part of the Algorithms for Lattice Fermions package (GPU port, ticket 208-2).
!
!  ALF is free software under the GNU General Public License, version 3 or later (see license.GPL).

!--------------------------------------------------------------------
!> @brief
!> Equal-time observables of the Hubbard model on the device (ticket 208-2).
!>
!> @details
!> Replaces ham%Obser of Hamiltonian_Hubbard_smod for the supported case (N_FL = 1, or N_FL = 2 with N_SUN = 1; |ham_h0| < Eps_small; no
!> Symm; no Projector; device updates active).  The Green function is read from G_d (device-authoritative since 208-1), so no D->H of G
!> is needed per time slice.
!>
!> Scalars (kin, pot, rho): device reductions, one D->H of three complex numbers per slice (48 B); the host accumulates them in Obs_scal
!> exactly as Obser does.  Correlations (Green, SpinZ/SpinXY/SpinT or SpinSUN, Den): for every output element (imj, no_I, no_J) of
!> Obs_Latt a gather over the lattice sites I (the partner site J of I at displacement imj is unique, table o_jc), accumulated on the
!> device in acc_k (shape of Obs_Latt) and acc0_k (shape of Obs_Latt0); gpu_obs_flush adds them into the host arrays once per bin
!> (before Pr_obs) and zeroes them.  No atomics: the sum over I is a sequential loop per output element.
!>
!> The formulas are those of Predefined_Obs_eq_{Green,SpinSUN,SpinMz,Den}_measure and of Obser (GRC(i,j,nf) = delta_ij - GR(j,i,nf)).
!>
!> Environment: ALF_GPU_NOOBS=1 keeps the host Obser (A/B tests).
!> Note (nvfortran): the sign/phase weight enters the device loops as two real scalars (wr, wi); a complex scalar there was read as garbage.
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_obs_mod
  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Natural_Constants, only: Eps_small
  Use Lattices_v3
  Use Predefined_Hoppings
  Use Hamiltonian_main, only: Obs_scal, Obs_eq, Obs_tau
  Use gpu_state_mod
  Use gpu_update_mod, only: gpu_upd_enabled
  Use gpu_update_mod, only: gpu_tacc
  Use gpu_wrap_mod, only: gpu_symm, gpu_symm_build
  Use alf_la_mod, only: la_gemm, LA_N
  Implicit none
  Private
  Logical, Save :: o_symm = .false.
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: Gsv(:,:,:), Ts(:,:,:,:), Asy(:,:,:), Bsy(:,:,:), Tmp(:,:)
  Logical, Save :: symm_built = .false.
  Integer, Save :: k_b = 0
  Public :: gpu_obs_setup, gpu_obs_active, gpu_obs_measure, gpu_obs_flush, gpu_obs_finalize, gpu_obs_symm
  Public :: gpu_obs_tau_active, gpu_obs_tau_measure

  Logical, Save :: o_set = .false.         ! model data stored (gpu_obs_setup called)
  Logical, Save :: o_ok = .false.          ! the model is supported
  Logical, Save :: o_dev = .false.         ! device arrays allocated
  Logical, Save :: o_off = .false.         ! ALF_GPU_NOOBS=1
  Integer, Save :: o_N, o_Norb, o_NFL, o_NSUN, o_Ndim, o_nk, o_nobs
  Integer, Save :: o_kind(8)               ! per Obs_eq: 1 Green, 2 SpinSUN, 3 Den, 4 SpinMz (fills Obs_eq(k:k+2))
  Integer, Allocatable, Save :: o_inv(:,:), o_jc(:,:), hk_i(:), hk_j(:), hk_f(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: hk_c(:)
  Real (Kind=Kind(0.d0)), Allocatable, Save :: o_uo(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: acc(:,:,:,:,:)    ! (N,1,Norb,Norb,nobs)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: acc0(:,:)        ! (Norb,nobs)
  Integer, Save :: o_dec
  Logical, Save :: o_tau = .false., o_tdev = .false.     ! tau-resolved observables (ticket 208-3)
  Integer, Save :: o_Nt
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: tacc(:,:,:,:,:)   ! (N,Nt,Norb,Norb,nobs)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: tacc0(:,:)        ! (Norb,nobs)

Contains

  !> Store the model data of the Hubbard model (host only; the device arrays are made at the first measurement).
  Subroutine gpu_obs_setup(Latt, Latt_unit, List, Invlist, hop, N_SUN, N_FL, ham_U, ham_U2, bilayer, h0, Symm_, Projector_)
    Implicit none
    Type (Lattice), intent(in) :: Latt
    Type (Unit_cell), intent(in) :: Latt_unit
    Integer, intent(in) :: List(:,:), Invlist(:,:), N_SUN, N_FL
    Type (Hopping_Matrix_type), allocatable, intent(in) :: hop(:)
    Real (Kind=Kind(0.d0)), intent(in) :: ham_U, ham_U2, h0
    Logical, intent(in) :: bilayer, Symm_, Projector_
    Integer :: I, J, I1, no_I, nf, Nb, n_1, n_2, nz, k, mode
    Integer :: N_Phi, pass
    Real (Kind=Kind(0.d0)) :: Phi_X, Phi_Y
    Logical :: Bulk
    Complex (Kind=Kind(0.d0)) :: Z
    Integer, Allocatable :: cnt(:,:)
    Character (len=16) :: s

    o_set = .true.; o_ok = .false.; o_tau = .false.
    call get_environment_variable('ALF_GPU_NOOBS', s)
    o_off = (trim(s) == '1')
    o_symm = Symm_
    if (Projector_) return          ! Symm is handled on the device (gpu_symm)
    if (abs(h0) >= Eps_small) return
    if (N_FL /= 1 .and. N_FL /= 2) return
    if (N_FL == 2 .and. N_SUN /= 1) return
    if (size(Obs_eq,1) /= merge(5, 3, N_FL == 2)) return
    o_N = Latt%N; o_Norb = Latt_unit%Norb; o_NFL = N_FL; o_NSUN = N_SUN; o_Ndim = size(List,1)
    if (o_Ndim /= o_N*o_Norb) return
    allocate(o_inv(o_N,o_Norb), o_jc(o_N,o_N), o_uo(o_Norb))
    o_inv = Invlist(1:o_N, 1:o_Norb)
    ! partner table: o_jc(I, imj) = the J with imj(I,J) = imj; it must be a bijection per I
    allocate(cnt(o_N,o_N)); cnt = 0; o_jc = 0
    do I = 1, o_N
       do J = 1, o_N
          k = Latt%imj(I,J)
          if (k < 1 .or. k > o_N) then; deallocate(o_inv, o_jc, o_uo, cnt); return; endif
          o_jc(I,k) = J; cnt(I,k) = cnt(I,k) + 1
       enddo
    enddo
    if (any(cnt /= 1)) then; deallocate(o_inv, o_jc, o_uo, cnt); return; endif
    deallocate(cnt)
    o_uo = ham_U
    if (bilayer .and. o_Norb >= 2) o_uo(2) = ham_U2
    o_dec = 1; if (N_FL == 2) o_dec = 2
    ! kinetic energy as a sparse linear functional of GRC: Z_kin = sum_k hk_c(k) GRC(hk_i(k), hk_j(k), hk_f(k)); the code of
    ! Predefined_Hoppings_Compute_Kin (no pinning), mode 1: T_Loc only, else bonds + T_Loc
    mode = inquire_hop(hop)
    nz = 0
    do pass = 1, 2
       k = 0
       if (mode == 1) then
          do nf = 1, N_FL
             do I = 1, o_N
                do no_I = 1, o_Norb
                   k = k + 1
                   if (pass == 2) then; hk_i(k) = o_inv(I,no_I); hk_j(k) = o_inv(I,no_I); hk_f(k) = nf; hk_c(k) = hop(nf)%T_Loc(no_I); endif
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
                   if (pass == 2) then; hk_i(k) = o_inv(I,no_I); hk_j(k) = o_inv(J,hop(nf)%list(Nb,2)); hk_f(k) = nf; hk_c(k) = Z; endif
                   k = k + 1
                   if (pass == 2) then; hk_i(k) = o_inv(J,hop(nf)%list(Nb,2)); hk_j(k) = o_inv(I,no_I); hk_f(k) = nf; hk_c(k) = conjg(Z); endif
                enddo
                do no_I = 1, o_Norb
                   k = k + 1
                   if (pass == 2) then; hk_i(k) = o_inv(I,no_I); hk_j(k) = o_inv(I,no_I); hk_f(k) = nf; hk_c(k) = hop(nf)%T_Loc(no_I); endif
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
       o_nobs = 5; o_kind(1:5) = (/ 1, 4, 4, 4, 3 /)
    else
       o_nobs = 3; o_kind(1:3) = (/ 1, 2, 3 /)
    endif
    o_ok = .true.
    ! tau-resolved observables (allocated by Alloc_obs only for Ltau = 1): the same observables, one per Obs_eq
    if (allocated(Obs_tau)) then
       if (size(Obs_tau,1) == o_nobs) then
          o_Nt = size(Obs_tau(1)%Obs_Latt,2)
          o_tau = .true.
       endif
    endif
  End Subroutine gpu_obs_setup

  !> .true. when the device measurement also applies the symmetrisation Hop_mod_Symm (the host must then skip it)
  Logical Function gpu_obs_symm()
    Implicit none
    gpu_obs_symm = o_symm .and. gpu_obs_active()
  End Function gpu_obs_symm

  !> out = Symm(in) = A in B  (A, B built once from the hop elements); in and out must be different arrays
  Subroutine symm_apply(inn, out, nf)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: inn(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: out(*)
    Integer, intent(in) :: nf
    Complex (Kind=Kind(0.d0)), parameter :: cone = (1.d0,0.d0), czero = (0.d0,0.d0)
    if (.not. symm_built) then
       allocate(Asy(o_Ndim,o_Ndim,o_NFL), Bsy(o_Ndim,o_Ndim,o_NFL), Tmp(o_Ndim,o_Ndim))
       !$omp target enter data map(alloc: Asy, Bsy, Tmp)
       do k_b = 1, o_NFL
          call gpu_symm_build(Asy(1,1,k_b), Bsy(1,1,k_b), o_Ndim, o_Ndim, k_b)
       enddo
       symm_built = .true.
    endif
    call la_gemm(LA_N, LA_N, o_Ndim, o_Ndim, o_Ndim, cone, Asy(1,1,nf), o_Ndim, inn, o_Ndim, czero, Tmp, o_Ndim)
    call la_gemm(LA_N, LA_N, o_Ndim, o_Ndim, o_Ndim, cone, Tmp, o_Ndim, Bsy(1,1,nf), o_Ndim, czero, out, o_Ndim)
  end subroutine symm_apply

  Subroutine k_cpy(nel, a, b)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)), intent(in) :: a(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: b(*)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nel
       b(i) = a(i)
    enddo
  end subroutine k_cpy

  !> .true. when the device measurement replaces ham%Obser
  Logical Function gpu_obs_active()
    Implicit none
    gpu_obs_active = o_set .and. o_ok .and. (.not. o_off) .and. gpu_active .and. gpu_upd_enabled()
  End Function gpu_obs_active

  Subroutine alloc_dev()
    Implicit none
    allocate(acc(o_N,1,o_Norb,o_Norb,o_nobs), acc0(o_Norb,o_nobs))
    acc = (0.d0, 0.d0); acc0 = (0.d0, 0.d0)
    !$omp target enter data map(to: o_inv, o_jc, o_uo, hk_i, hk_j, hk_f, hk_c, acc, acc0)
    o_dev = .true.
  End Subroutine alloc_dev

  !> ham%Obser on the device.  GR (host) is read only when the device copy of a flavour is not current.
  Subroutine gpu_obs_measure(GR, Phase, Ntau, Mc_step_weight)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(:,:,:)
    Complex (Kind=Kind(0.d0)), intent(in) :: Phase
    Integer, intent(in) :: Ntau
    Real (Kind=Kind(0.d0)), intent(in) :: Mc_step_weight
    Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZW, Zkin, Zpot, Zrho, zc, zt
    Integer :: I, nf, k, no, no2, imj, I1, J1, kd, N, Norb, NFL, NSUN, dec, kind_k
    Real (Kind=Kind(0.d0)) :: rs, sq3

    if (.not. o_dev) call alloc_dev()
    do nf = 1, o_NFL          ! the device copy must be current (it is, except after a host write of GR)
       if (.not. g_dev_ok(nf)) then
          G_d(:,:,nf) = GR(:,:,nf)
          !$omp target update to(G_d(:,:,nf))
          g_dev_ok(nf) = .true.
       endif
    enddo
    if (o_symm) then                  ! measure on the symmetrised G: save G_d, transform in place, restore at the end
       if (.not. allocated(Gsv)) then
          allocate(Gsv(o_Ndim,o_Ndim,o_NFL))
          !$omp target enter data map(alloc: Gsv)
       endif
       do nf = 1, o_NFL
          call k_cpy(o_Ndim*o_Ndim, G_d(1,1,nf), Gsv(1,1,nf))
          call symm_apply(Gsv(1,1,nf), G_d(1,1,nf), nf)
       enddo
    endif
    ZP = Phase/Real(Phase, kind(0.D0))
    ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
    ZS = ZS*Mc_step_weight
    ZW = ZP*ZS
    N = o_N; Norb = o_Norb; NFL = o_NFL; NSUN = o_NSUN; dec = o_dec; sq3 = sqrt(3.d0)
    rs = Real(ZS, kind(0.d0))
    do I = 1, size(Obs_scal,1)
       Obs_scal(I)%N = Obs_scal(I)%N + 1
       Obs_scal(I)%Ave_sign = Obs_scal(I)%Ave_sign + rs
    enddo
    do I = 1, size(Obs_eq,1)
       Obs_eq(I)%N = Obs_eq(I)%N + 1
       Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + rs
    enddo

    ! ---- scalars: GRC(i,j,nf) = delta_ij - G(j,i,nf)
    Zkin = (0.d0, 0.d0)
    !$omp target teams distribute parallel do reduction(+:Zkin) private(i1,J1,nf,zc)
    do k = 1, o_nk
       i1 = hk_i(k); J1 = hk_j(k); nf = hk_f(k)
       zc = -G_d(J1,i1,nf)
       if (i1 == J1) zc = zc + (1.d0, 0.d0)
       Zkin = Zkin + hk_c(k)*zc
    enddo
    Zkin = Zkin*dble(NSUN)
    Zpot = (0.d0, 0.d0)
    !$omp target teams distribute parallel do collapse(2) reduction(+:Zpot) private(i1,zc,zt)
    do I = 1, N
       do no = 1, Norb
          i1 = o_inv(I,no)
          zc = (1.d0, 0.d0) - G_d(i1,i1,1)
          zt = (1.d0, 0.d0) - G_d(i1,i1,dec)
          Zpot = Zpot + zc*zt*o_uo(no)
       enddo
    enddo
    Zrho = (0.d0, 0.d0)
    !$omp target teams distribute parallel do collapse(3) reduction(+:Zrho) private(i1)
    do nf = 1, NFL
       do I = 1, N
          do no = 1, Norb
             i1 = o_inv(I,no)
             Zrho = Zrho + ((1.d0, 0.d0) - G_d(i1,i1,nf))
          enddo
       enddo
    enddo
    Zrho = Zrho*dble(NSUN)
    Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zkin*ZP*ZS
    Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Zpot*ZP*ZS
    Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zrho*ZP*ZS
    Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

    ! ---- correlations: one gather per output element (imj, no, no2)
    do k = 1, o_nobs
       kind_k = o_kind(k)
       if (kind_k == 4 .and. k /= 2) cycle        ! SpinMz fills Obs_eq(2:4) from k = 2
       call corr(k, kind_k, ZW, N, Norb, NFL, NSUN, sq3)
    enddo
    if (o_symm) then
       do nf = 1, o_NFL
          call k_cpy(o_Ndim*o_Ndim, Gsv(1,1,nf), G_d(1,1,nf))
       enddo
    endif
  End Subroutine gpu_obs_measure

  Subroutine corr(k, kind_k, ZW, N, Norb, NFL, NSUN, sq3)
    Implicit none
    Integer, intent(in) :: k, kind_k, N, Norb, NFL, NSUN
    Complex (Kind=Kind(0.d0)), intent(in) :: ZW
    Real (Kind=Kind(0.d0)), intent(in) :: sq3
    Integer :: m, imj, no, no2, I, J, i1, j1, nf
    Complex (Kind=Kind(0.d0)) :: s1, s2, s3, z, zi, zj, zxy, zz, g11, g22, a, b, zs0
    Real (Kind=Kind(0.d0)) :: ns, wr, wi

    ns = dble(NSUN)
    wr = real(ZW); wi = aimag(ZW)
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
                      a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
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
                   a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
                   s1 = s1 + a*G_d(i1,j1,1)*ns
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
                      zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
                      zj = zj + ((1.d0, 0.d0) - G_d(j1,j1,nf))
                      a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                      z = z + a*G_d(i1,j1,nf)
                   enddo
                   s1 = s1 + (zi*ns)*(zj*ns) + z*ns
                enddo
                !$omp end parallel do
                acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute parallel do private(I,i1,zi,nf,zs0)
       do no = 1, Norb
          zs0 = (0.d0, 0.d0)
          do I = 1, N
             i1 = o_inv(I,no); zi = (0.d0, 0.d0)
             do nf = 1, NFL
                zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
             enddo
             zs0 = zs0 + zi*ns
          enddo
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
                   ! GRC1 = c11 = delta - G(j1,i1,1) etc.
                   a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)      ! GRC(I1,J1,1)
                   b = -G_d(j1,i1,2); if (i1 == j1) b = b + (1.d0, 0.d0)      ! GRC(I1,J1,2)
                   zxy = a*G_d(i1,j1,2) + b*G_d(i1,j1,1)
                   zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))     ! GRC(I1,I1,2) - GRC(I1,I1,1)
                   zj = ((1.d0, 0.d0) - G_d(j1,j1,2)) - ((1.d0, 0.d0) - G_d(j1,j1,1))
                   zz = a*G_d(i1,j1,1) + b*G_d(i1,j1,2) + zi*zj
                   s1 = s1 + zz
                   s2 = s2 + zxy
                   s3 = s3 + (2.d0*zxy + zz)/3.d0
                enddo
                !$omp end parallel do
                acc(imj,1,no,no2,k)   = acc(imj,1,no,no2,k)   + s1*cmplx(wr,wi,kind(0.d0))
                acc(imj,1,no,no2,k+1) = acc(imj,1,no,no2,k+1) + s2*cmplx(wr,wi,kind(0.d0))
                acc(imj,1,no,no2,k+2) = acc(imj,1,no,no2,k+2) + s3*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute parallel do private(I,i1,zi,zs0)
       do no = 1, Norb
          zs0 = (0.d0, 0.d0)
          do I = 1, N
             i1 = o_inv(I,no)
             zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
             zs0 = zs0 + zi
          enddo
          acc0(no,k)   = acc0(no,k)   + zs0*cmplx(wr,wi,kind(0.d0))
          acc0(no,k+2) = acc0(no,k+2) + zs0*cmplx(wr,wi,kind(0.d0))/sq3
       enddo
    endif
  End Subroutine corr

  !> Add the device accumulators to the host Obs_eq (before Pr_obs normalises them) and zero them
  Subroutine gpu_obs_flush()
    Implicit none
    Integer :: k
    if (o_tdev) then
       !$omp target update from(tacc, tacc0)
       do k = 1, o_nobs
          Obs_tau(k)%Obs_Latt = Obs_tau(k)%Obs_Latt + tacc(:,:,:,:,k)
          Obs_tau(k)%Obs_Latt0 = Obs_tau(k)%Obs_Latt0 + tacc0(:,k)
       enddo
       tacc = (0.d0, 0.d0); tacc0 = (0.d0, 0.d0)
       !$omp target update to(tacc, tacc0)
    endif
    if (.not. o_dev) return
    !$omp target update from(acc, acc0)
    do k = 1, o_nobs
       Obs_eq(k)%Obs_Latt = Obs_eq(k)%Obs_Latt + acc(:,:,:,:,k)
       Obs_eq(k)%Obs_Latt0 = Obs_eq(k)%Obs_Latt0 + acc0(:,k)
    enddo
    acc = (0.d0, 0.d0); acc0 = (0.d0, 0.d0)
    !$omp target update to(acc, acc0)
  End Subroutine gpu_obs_flush

  Subroutine gpu_obs_finalize()
    Implicit none
    if (allocated(Gsv)) then
       !$omp target exit data map(delete: Gsv)
       deallocate(Gsv)
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
       !$omp target exit data map(delete: o_inv, o_jc, o_uo, hk_i, hk_j, hk_f, hk_c, acc, acc0)
       deallocate(acc, acc0)
       o_dev = .false.
    endif
  End Subroutine gpu_obs_finalize

  !> .true. when the device tau-resolved measurement replaces ham%OBSERT (not with Symm / Projector: refused by gpu_obs_setup)
  Logical Function gpu_obs_tau_active()
    Implicit none
    gpu_obs_tau_active = o_set .and. o_ok .and. o_tau .and. (.not. o_off) .and. gpu_active
  End Function gpu_obs_tau_active

  !> ham%OBSERT on the device: GT0, G0T, G00, GTT are the device-resident matrices of the tau pass (no D->H per slice)
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
       !$omp target enter data map(to: tacc, tacc0)
       o_tdev = .true.
    endif
    ZP = Phase/Real(Phase, kind(0.D0))
    ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
    ZS = ZS*Mc_step_weight
    ZW = ZP*ZS
    wr = real(ZW); wi = aimag(ZW)
    rs = Real(ZS, kind(0.d0))
    if (nt == 0) then                       ! counters once per pass, as the Predefined_Obs_tau_* routines
       do I = 1, o_nobs
          Obs_tau(I)%N = Obs_tau(I)%N + 1
          Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + rs
       enddo
    endif
    call gpu_tacc(6, -1.d0)
    if (o_symm) then                  ! symmetrised copies of the four matrices (Hop_mod_Symm of tau_obs on the host before)
       if (.not. allocated(Ts)) then
          allocate(Ts(o_Ndim,o_Ndim,o_NFL,4))
          !$omp target enter data map(alloc: Ts)
       endif
       do nf = 1, o_NFL
          call symm_apply(GT0(1,1,nf), Ts(1,1,nf,1), nf)
          call symm_apply(G0T(1,1,nf), Ts(1,1,nf,2), nf)
          call symm_apply(G00(1,1,nf), Ts(1,1,nf,3), nf)
          call symm_apply(GTT(1,1,nf), Ts(1,1,nf,4), nf)
       enddo
       do k = 1, o_nobs
          if (o_kind(k) == 4 .and. k /= 2) cycle
          call tcorr(k, o_kind(k), nt, wr, wi, Ts(:,:,:,1), Ts(:,:,:,2), Ts(:,:,:,3), Ts(:,:,:,4))
       enddo
    else
       do k = 1, o_nobs
          if (o_kind(k) == 4 .and. k /= 2) cycle
          call gpu_tacc(6+k, -1.d0)
          call tcorr(k, o_kind(k), nt, wr, wi, GT0, G0T, G00, GTT)
          call gpu_tacc(6+k, 1.d0)
       enddo
    endif
    call gpu_tacc(6, 1.d0)
  End Subroutine gpu_obs_tau_measure

  Subroutine tcorr(k, kind_k, nt, wr, wi, GT0, G0T, G00, GTT)
    Implicit none
    Integer, intent(in) :: k, kind_k, nt
    Real (Kind=Kind(0.d0)), intent(in) :: wr, wi
    Complex (Kind=Kind(0.d0)), intent(in) :: GT0(o_Ndim,o_Ndim,o_NFL), G0T(o_Ndim,o_Ndim,o_NFL), G00(o_Ndim,o_Ndim,o_NFL), GTT(o_Ndim,o_Ndim,o_NFL)
    Integer :: m, imj, no, no2, I, J, i1, j1, nf, N, Norb, NFL, it
    Complex (Kind=Kind(0.d0)) :: s1, s2, s3, z, zi, zj, zxy, zz, zs0
    Real (Kind=Kind(0.d0)) :: ns, rn, sq3

    N = o_N; Norb = o_Norb; NFL = o_NFL; it = nt + 1
    ns = dble(o_NSUN); rn = 1.d0/dble(o_NFL); sq3 = sqrt(3.d0)
    if (kind_k == 1) then            ! Green: sum_nf GT0 / N_FL
       !$omp target teams distribute private(imj,no,no2,s1)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,nf,z) reduction(+:s1)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             z = (0.d0, 0.d0)
             do nf = 1, NFL
                z = z + GT0(i1,j1,nf)
             enddo
             s1 = s1 + z*rn
          enddo
          !$omp end parallel do
          tacc(imj,it,no,no2,k) = tacc(imj,it,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
    elseif (kind_k == 2) then        ! SpinSUN (N_FL = 1): - N_SUN G0T(J1,I1,1) GT0(I1,J1,1)
       !$omp target teams distribute private(imj,no,no2,s1)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1) reduction(+:s1)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             s1 = s1 - ns*G0T(j1,i1,1)*GT0(i1,j1,1)
          enddo
          !$omp end parallel do
          tacc(imj,it,no,no2,k) = tacc(imj,it,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
       enddo
    elseif (kind_k == 3) then        ! Den
       !$omp target teams distribute private(imj,no,no2,s1)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,nf,z,zi,zj) reduction(+:s1)
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
          !$omp end parallel do
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
    else                             ! SpinMz (N_FL = 2, N_SUN = 1): k = SpinZ, k+1 = SpinXY, k+2 = SpinT
       !$omp target teams distribute private(imj,no,no2,s1,s2,s3)
       do m = 0, N*Norb*Norb - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0)
          !$omp parallel do private(I,J,i1,j1,zxy,zz) reduction(+:s1,s2,s3)
          do I = 1, N
             J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
             zz  = (GTT(i1,i1,1) - GTT(i1,i1,2))*(G00(j1,j1,1) - G00(j1,j1,2)) &
                  &  - G0T(j1,i1,1)*GT0(i1,j1,1) - G0T(j1,i1,2)*GT0(i1,j1,2)
             zxy = - G0T(j1,i1,1)*GT0(i1,j1,2) - G0T(j1,i1,2)*GT0(i1,j1,1)
             s1 = s1 + zz
             s2 = s2 + zxy
             s3 = s3 + (2.d0*zxy + zz)/3.d0
          enddo
          !$omp end parallel do
          tacc(imj,it,no,no2,k)   = tacc(imj,it,no,no2,k)   + s1*cmplx(wr,wi,kind(0.d0))
          tacc(imj,it,no,no2,k+1) = tacc(imj,it,no,no2,k+1) + s2*cmplx(wr,wi,kind(0.d0))
          tacc(imj,it,no,no2,k+2) = tacc(imj,it,no,no2,k+2) + s3*cmplx(wr,wi,kind(0.d0))
       enddo
       !$omp target teams distribute parallel do private(I,i1,zi,zs0)
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
  End Subroutine tcorr
End Module gpu_obs_mod
#endif

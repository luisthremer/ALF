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
!> GPU port, ticket 203 (ADR-002 section 2): flat, device-resident mirror of everything the sweep needs that
!> does not change during a run, plus the device buffers later tickets fill.  Compiled only with -DALF_GPU.
!>
!> Plain contiguous arrays only (no derived types on the device).  Host objects (Op_T, Op_V, nsigma) stay the
!> source of truth; this module copies from them once at setup (ticket 203: nothing in the sweep uses the mirror).
!>
!> Flat layout (all offsets are 0-based "number of entries before", Fortran arrays are 1-based):
!>  hopping (per flavour nf, elements in the order stored in Hop_mod, active elements only, i.e. those the Hop_mod
!>  routines do not skip):
!>    hop_start(nf)       first element of flavour nf, elements hop_start(nf) : hop_start(nf+1)-1
!>    hop_meta(1:4,e)     kind (1 real, 2 complex), nd = Ndim_hop, poff, moff
!>    hop_P(poff+1:poff+nd)                  projector (site list) of element e
!>    kind 1: hop_Rval(moff + (k-1)*nd*nd + i + (j-1)*nd)  = matrix k (i,j), real
!>    kind 2: hop_Cval(moff + (k-1)*nd*nd + i + (j-1)*nd)  = matrix k (i,j), complex
!>    k = 1 mat (e^{-dtau T}), 2 invmat, 3 mat_1D2 (half step), 4 invmat_1D2.  Full nd x nd arrays.  OpTTypes_mod
!>    symmetrises only i <= j.  Which triangle ALF reads (Mat_subroutines_mod.F90):
!>      - complex, any nd: ZSLHEMM, nd <= 8 hand-written loops on the upper triangle completed by conjg, nd > 8 BLAS
!>        ZHEMM 'U' (upper triangle only);
!>      - real, nd <= 8: ZDSLSYMM hand-written loops, upper triangle completed by copy;
!>      - real, nd > 8 (dense hopping): ZDSLSYMM calls ZLARCM with the FULL stored matrix (the commented-out ZHEMM
!>        calls there are not active), so the raw lower triangle (not symmetrised, differs ~1e-16) is what ALF applies.
!>    build_hopping therefore fills the lower triangle from the upper one except for real elements with nd > 8, and
!>    the flat arrays encode exactly the operator ALF applies in every case.
!>    NOTE for ticket 204: for nd <= 8 ALF uses hand-written loops (sum order l = 1..nd); for nd > 8 it uses BLAS
!>    (ZHEMM / DGEMM via ZLARCM), so a device kernel is not bit-identical there (~1e-15; kernel bound 1e-10).
!>    Order of application = order of the Hop_mod routines (see gpu_selftest): mmthr nc=N..1 (mat, from the left),
!>    mmthr_m1 nc=1..N (invmat, left), mmthl nc=1..N (mat, right), mmthlc nc=1..N (mat, left),
!>    mmthl_m1 nc=N..1 (invmat, right), Symm nc=N..1 (mat_1D2 left then invmat_1D2 right).
!>  Op_V (diagonal operators of type 1 or 2 with constant coupling):
!>    V_meta(1:4,n,nf)    N (number of sites), poff (into V_P), eoff (into V_exp), skip flag
!>    V_P(poff+1:poff+N)  sites; V_type(n) field type (same for all flavours, checked in gpu_supported)
!>    V_exp(eoff + i + (sp+t)*N), i=1..N, sp=-t..t = E_exp(i,sp) of Operator_mod, always the real table (the factor
!>    Op_mmultL/R/Op_Wrapup/Op_Wrapdo apply for field value sp); sp=0 is unused by ALF and set to 1.
!>    skip flag = 1 iff |g| < Eps_machine: Op_mmultL/R then return early (kernels for L/R must honour the flag),
!>    Op_Wrapup/Op_Wrapdo do NOT check g and always apply the table (wrap kernels must ignore the flag).
!>  Fields: nsigma_d(N_op,Ltrot) = nsigma%f, nsigma_t(N_op) = nsigma%t.
!>
!>  Notes for ticket 205 (updates): the determinant ratio contains exp(Z1*Op_V%alpha) (upgrade_mod.F90:193) and
!>  alpha is not in the mirror; nsigma_d is a copy made at setup and goes stale as soon as updates run (205 must keep
!>  it in sync); N and P of an operator may differ between flavours (V_meta is per (n,nf)).  gpu_supported refuses
!>  runs with Propose_S0 or reconstruction_needed (ham%weight_reconstruction).
!>  Ticket 205 (done, gpu_update_mod): alpha is used on the host (Op_V%alpha); nsigma_d is updated per accepted flip (host
!>  copy + one scalar H->D); per-flavour N/P are read from the host Op_V; gpu_supported also refuses N_non_zero > 1 (the
!>  device update is rank 1).
!>
!>  Ticket 204 additions (kernels live in gpu_wrap_mod):
!>    hop_gstart(nf), hop_gptr(g)   checkerboard groups: maximal runs of consecutive active elements of one flavour
!>                        whose site lists P are pairwise disjoint and nd <= 8 (an element with nd > 8 forms a group of
!>                        its own).  Group g = elements hop_gptr(g) : hop_gptr(g+1)-1; flavour nf = groups
!>                        hop_gstart(nf) : hop_gstart(nf+1)-1.  Elements of one group touch disjoint rows (left) /
!>                        columns (right), so they can run in one kernel and give the sequential result; the same
!>                        partition serves both application directions (a group is order independent inside).
!>                        hop_big(nf) = flavour has an element with nd > 8.
!>    hop_aux(1:2,e)      (xoff, pident): pident = 1 iff P = (1,2,...,nd); xoff = offset (0-based) into hop_At/hop_Ac for
!>                        real elements with nd > 8, else -1.
!>    hop_At(xoff + (k-1)*nd*nd + j + (i-1)*nd) = A_k(i,j) (transposed copy, real), hop_Ac(xoff + (k-1)*nd*nd + i +
!>                        (j-1)*nd) = A_k(i,j) as complex; k = 1 (mat), 2 (invmat) only.  Used by the dense left
!>                        multiplication variants of gpu_wrap_mod.
!>    complex elements with nd > 8: the imaginary part of the diagonal is zeroed (ZHEMM, which ALF calls there, ignores
!>                        it; the nd <= 8 loops use it, but a Hermitian matrix has none).
!>    vr_ptr(r,nf), vr_ent(1:2,p)  per matrix index r (site) and flavour: the entries p = vr_ptr(r,nf) : vr_ptr(r+1,nf)-1
!>                        list (n, i) = (operator, position in P) with V_P(poff(n,nf)+i) = r, in ALF's operator order
!>                        (n ascending).  vr_ptr has Ndim+1 rows.
!>  Offsets are default (int32) integers: arrays must stay below 2**31-1 entries (dense hopping needs 4*N*N*N_FL
!>  reals, i.e. fine up to N ~ 20000).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_state_mod

  Use omp_lib
  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Natural_Constants, only: Eps_machine
  Use Operator_mod
  Use Fields_mod
  Use Hamiltonian_main
  Use QMC_runtime_var
  Use Hop_mod
  Use Mat_subroutines
#if !defined(ALF_LA_ONEMKL)
  Use alf_la_mod, only: la_init, la_finalize
#endif

  Implicit none

  Private
  Public :: gpu_state_setup, gpu_state_finalize, gpu_supported, gpu_active, gpu_reason, gpu_scratch_lwork
  Public :: hop_start, hop_meta, hop_P, hop_Rval, hop_Cval
  Public :: hop_gstart, hop_gptr, hop_big, hop_aux, hop_At, hop_Ac, vr_ptr, vr_ent, g_ngroups
  Public :: rnd, rand_mat, rel_diff, is_finite, finite_c, selftest_fail, info_line, same_bits_c
  Public :: V_meta, V_type, V_P, V_exp, nsigma_d, nsigma_t
  Public :: G_d, udv_U, udv_V, udv_D, udv_L, udw_U, udw_V, udw_D, udw_L
  Public :: GT0_d, G0T_d, GTT_d, G00_d, scr_c, scr_vec, scr_r, scr_work, scr_ipiv
  Public :: upd_cnt

  Logical, Save :: gpu_active = .false.
  Character (len=256), Save :: gpu_reason = 'gpu_state_setup not called'

  Integer, Save :: g_NFL = 0, g_Ndim = 0, g_Nop = 0, g_Ltrot = 0, g_NSTM = 0, g_NFLeff = 0
  Integer, Save :: g_nel = 0, g_ngroups = 0

  ! ---- constant data, built once from the host objects ----
  Integer,                 Allocatable, Save :: hop_start(:), hop_meta(:,:), hop_P(:)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: hop_Rval(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: hop_Cval(:)
  Integer,                 Allocatable, Save :: V_meta(:,:,:), V_type(:), V_P(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: V_exp(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: nsigma_d(:,:)
  Integer,                 Allocatable, Save :: nsigma_t(:)
  ! ticket 204: groups, dense auxiliaries, per-site operator lists (see header)
  Integer,                 Allocatable, Save :: hop_gstart(:), hop_gptr(:), hop_aux(:,:), vr_ptr(:,:), vr_ent(:,:)
  Logical,                 Allocatable, Save :: hop_big(:)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: hop_At(:)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: hop_Ac(:)

  ! ---- buffers allocated on the device only (not filled in ticket 203) ----
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: G_d(:,:,:)                    ! (Ndim,Ndim,N_FL)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: udv_U(:,:,:,:), udv_V(:,:,:,:)   ! (Ndim,Ndim,NSTM,N_FL_eff)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: udv_D(:,:,:)                  ! (Ndim,NSTM,N_FL_eff)  default
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: udv_L(:,:,:)                  ! (Ndim,NSTM,N_FL_eff)  STABLOG
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: udw_U(:,:,:,:), udw_V(:,:,:,:)   ! (Ndim,Ndim,2,N_FL_eff) udvl, udvr
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: udw_D(:,:,:)                  ! (Ndim,2,N_FL_eff)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: udw_L(:,:,:)                  ! (Ndim,2,N_FL_eff)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: GT0_d(:,:,:), G0T_d(:,:,:), GTT_d(:,:,:), G00_d(:,:,:)   ! only if Ltau = 1 (G00_d: ticket 207)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: scr_c(:,:,:), scr_vec(:,:), scr_work(:)
  Real    (Kind=Kind(0.d0)), Allocatable, Save :: scr_r(:,:)
  Integer,                 Allocatable, Save :: scr_ipiv(:)
  Logical, Save :: have_tau = .false.
  ! ticket 205 counters (gpu_update_mod), reported by gpu_state_finalize: 1 site proposals reaching the ratio, 2 accepted,
  ! 3 site-wrap launches, 4 update launches, 5 ratio-block D->H elements, 6 field H->D scalars, 7 interim G H->D copies,
  ! 8 interim G D->H copies, 9 interim bytes (both directions), 10 WRAPGRUP/DO calls
  Integer(int64), Save :: upd_cnt(10) = 0

  Character (len=256), Save :: info_file = 'info'
  Logical, Save :: info_root = .true.

Contains

!--------------------------------------------------------------------
!> QR/LU work space (complex entries) per scratch buffer: 64 blocks of Ndim, provisional size for tickets 204-206.
!--------------------------------------------------------------------
  Integer function gpu_scratch_lwork(Ndim)
    Implicit none
    Integer, intent(in) :: Ndim
    gpu_scratch_lwork = 64*Ndim + 64
  end function gpu_scratch_lwork

!--------------------------------------------------------------------
  Subroutine info_line(text)
    Implicit none
    Character (len=*), intent(in) :: text
    if (info_root) then
       Open (Unit=50, file=trim(info_file), status="unknown", position="append")
       Write(50,'(1X,A)') trim(text)
       Close(50)
    endif
  end subroutine info_line

!--------------------------------------------------------------------
!> Decides whether the run can use the GPU path.  Returns .true. if yes; else reason says why not.
!--------------------------------------------------------------------
  Logical function gpu_supported(reason)
    Implicit none
    Character (len=*), intent(out) :: reason

    Character (len=16) :: env
    Integer :: n, nf, nc, st, nh, ek
    Logical :: act
    Integer, allocatable :: P(:)
    Real    (Kind=Kind(0.d0)), allocatable :: rm(:,:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: cm(:,:,:)

    gpu_supported = .false.
    reason = ''
    env = ''
    call get_environment_variable('ALF_GPU_DISABLE', env, status=st)
    if (st == 0 .and. trim(env) == '1') then
       reason = 'disabled by ALF_GPU_DISABLE=1'; return
    endif
#if defined(MPI)
    reason = 'MPI build: GPU path not supported with MPI yet'; return
#endif
#if defined(TEMPERING)
    reason = 'tempering is not supported on the GPU path'; return
#endif
#if defined(STAB1) || defined(STAB2) || defined(STAB3)
    reason = 'STAB1/STAB2/STAB3 stabilization is not supported on the GPU path (default and LOG only, ticket 206)'; return
#endif
    if (Projector) then
       reason = 'projector code is not supported on the GPU path'; return
    endif
    if (get_Langevin() .or. get_HMC()) then
       reason = 'Langevin/HMC is not supported on the GPU path'; return
    endif
    if (get_Global_moves() .or. get_Global_tau_moves() .or. get_N_Global_tau() > 0) then
       reason = 'global moves are not supported on the GPU path'; return
    endif
    if (get_Propose_S0()) then
       reason = 'Propose_S0 is not supported on the GPU path'; return
    endif
    if (reconstruction_needed) then
       reason = 'Calc_Fl with weight reconstruction (ham%weight_reconstruction) is not supported on the GPU path'; return
    endif
    if (omp_get_num_devices() < 1) then
       reason = 'no OpenMP offload device found'; return
    endif
    do nf = 1, N_FL
       do n = 1, size(Op_V,1)
          if (.not. Op_V(n,nf)%diag) then
             write(reason,'(A,A,A,I0,A,I0,A)') 'Hamiltonian ', trim(get_ham_name()), &
                  & ': Op_V(', n, ',', nf, ') is not diagonal'
             return
          endif
          if (Op_V(n,nf)%type /= 1 .and. Op_V(n,nf)%type /= 2) then
             write(reason,'(A,A,A,I0,A,I0,A,I0,A)') 'Hamiltonian ', trim(get_ham_name()), &
                  & ': Op_V(', n, ',', nf, ') has field type ', Op_V(n,nf)%type, ' (only 1 and 2 supported)'
             return
          endif
          if (Op_V(n,nf)%type /= Op_V(n,1)%type) then
             write(reason,'(A,A,A,I0,A,I0,A)') 'Hamiltonian ', trim(get_ham_name()), &
                  & ': field type of Op_V(', n, ',', nf, ') differs from flavour 1'
             return
          endif
          if (Op_V(n,nf)%N_non_zero > 1) then
             ! ticket 205: the device update implements the rank-1 case (N_non_zero <= 1) only
             write(reason,'(A,A,A,I0,A,I0,A,I0,A)') 'Hamiltonian ', trim(get_ham_name()), &
                  & ': Op_V(', n, ',', nf, ') has N_non_zero = ', Op_V(n,nf)%N_non_zero, ' (device update supports <= 1)'
             return
          endif
          if (Op_V(n,nf)%get_g_t_alloc()) then
             write(reason,'(A,A,A,I0,A,I0,A)') 'Hamiltonian ', trim(get_ham_name()), &
                  & ': Op_V(', n, ',', nf, ') has a time-dependent coupling g_t'
             return
          endif
       enddo
    enddo
    nc = Hop_mod_gpu_ncheck()
    do nf = 1, N_FL
       do n = 1, nc
          call Hop_mod_gpu_element(nf, n, ek, nh, act, P, rm, cm)
          if (ek /= 1 .and. ek /= 2) then
             write(reason,'(A,A,A,I0,A,I0,A)') 'Hamiltonian ', trim(get_ham_name()), &
                  & ': hopping element (', n, ',', nf, ') is time dependent or of unknown type'
             return
          endif
       enddo
    enddo
    gpu_supported = .true.
    reason = 'all checks passed'
  end function gpu_supported

!--------------------------------------------------------------------
  Subroutine gpu_state_setup(nsigma_in, NSTM, file_info, root)
    Implicit none
    Type (Fields),     intent(in) :: nsigma_in
    Integer,           intent(in) :: NSTM
    Character (len=*), intent(in) :: file_info
    Logical,           intent(in) :: root

    Integer(int64) :: c0, c1, crate
    Real (Kind=Kind(0.d0)) :: tsetup
    Character (len=256) :: reason, line
    Character (len=16)  :: env
    Integer :: st
    Logical :: ok

    info_file = file_info
    info_root = root
    call system_clock(c0, crate)
    ok = gpu_supported(reason)
    gpu_reason = reason
    if (.not. ok) then
       gpu_active = .false.
       call info_line('GPU path (gpu_state): INACTIVE, reason: '//trim(reason))
       env = ''
       call get_environment_variable('ALF_GPU_SELFTEST', env, status=st)
       if (st == 0 .and. len_trim(env) > 0 .and. trim(env) /= '0') then
          write(error_unit,'(A)') ' ALF_GPU_SELFTEST is set but the GPU path is inactive: '//trim(reason)
          call info_line('GPU self-test requested but GPU path inactive: '//trim(reason))
          call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
       endif
       return
    endif

    g_NFL = N_FL; g_Ndim = Ndim; g_Nop = size(Op_V,1); g_Ltrot = Ltrot; g_NSTM = NSTM; g_NFLeff = N_FL_eff
    have_tau = ( get_Ltau() == 1 )

    call build_hopping()
    call build_groups()
    call build_V()
    call build_vrows()
    allocate( nsigma_d(g_Nop, g_Ltrot), nsigma_t(g_Nop) )
    nsigma_d = nsigma_in%f
    nsigma_t = nsigma_in%t
    call alloc_buffers()
    call map_to_device()
#if !defined(ALF_LA_ONEMKL)
    call la_init()      ! alf_la once, after the OpenMP runtime owns the device (ticket 205; was lazy in gpu_wrap_mod)
#endif
    gpu_active = .true.
    call system_clock(c1)
    tsetup = real(c1-c0,kind(0.d0))/real(crate,kind(0.d0))
    write(line,'(A,I0,A)') 'GPU path (gpu_state): ACTIVE (mirror built, used by the sweep), OpenMP device ', &
         & omp_get_default_device()
    call info_line(line)
    write(line,'(A,I0,A,I0,A,I0)') 'GPU mirror: hopping elements ', g_nel, ', Op_V operators ', g_Nop*g_NFL, &
         & ', omp devices ', omp_get_num_devices()
    call info_line(line)
    write(line,'(A,I0,A,I0,A,L1)') 'GPU mirror: checkerboard groups ', g_ngroups, ' (flavour 1: ', hop_gstart(2)-hop_gstart(1), &
         & ' groups), dense element(s) in flavour 1: ', hop_big(1)
    call info_line(line)
    call report_layout()
    ! timing of setup includes the host->device copies, not the self-test
    write(line,'(A,F10.4,A)') 'GPU mirror setup time (incl. host->device copies): ', tsetup, ' s'
    call info_line(line)
    write(output_unit,'(A,F10.4,A)') ' GPU mirror setup time (incl. host->device copies): ', tsetup, ' s'

    env = ''
    call get_environment_variable('ALF_GPU_SELFTEST', env, status=st)
    ! ALF_GPU_SELFTEST=1 (any value except empty/0): self-test, then the run continues.  Negative controls, which must
    ! fail with a nonzero exit: =fault (hop_Rval(1) *= 1+1e-9), =fault_c (hop_Cval(1) += 1e-9), =fault_v (V_exp(1) *=
    ! 1+1e-9), =nan (hop_Rval(1) = NaN), =nan_v (V_exp(1) = NaN); all applied to the host copy after the upload.
    ! ALF_GPU_SELFTEST=3 (ticket 206, stabilization self-test in gpu_stab_mod) and =4 (ticket 207, tau pass, gpu_tau_mod) do not run
    ! this mirror self-test: O(N^3) host loops
    if (st == 0 .and. len_trim(env) > 0 .and. trim(env) /= '0' .and. trim(env) /= '3' .and. trim(env) /= '4') then
       select case (trim(env))
       case ('fault');   hop_Rval(1) = hop_Rval(1)*(1.d0 + 1.d-9)
       case ('fault_c'); hop_Cval(1) = hop_Cval(1) + cmplx(1.d-9, 0.d0, kind(0.d0))
       case ('fault_v'); V_exp(1) = V_exp(1)*(1.d0 + 1.d-9)
       case ('nan');     hop_Rval(1) = transfer(9221120237041090560_int64, 0.d0)
       case ('nan_v');   V_exp(1) = cmplx(transfer(9221120237041090560_int64, 0.d0), 0.d0, kind(0.d0))
       end select
       call gpu_selftest()
    endif
  end subroutine gpu_state_setup

!--------------------------------------------------------------------
  Subroutine build_hopping()
    Implicit none
    Integer :: nf, nc, e, nel, npv, nrv, ncv, ek, nh, k, i, j
    Logical :: act
    Integer, allocatable :: P(:)
    Real    (Kind=Kind(0.d0)), allocatable :: rm(:,:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: cm(:,:,:)
    Integer :: ncheck

    ncheck = Hop_mod_gpu_ncheck()
    allocate(hop_start(N_FL+1))
    nel = 0; npv = 0; nrv = 0; ncv = 0
    do nf = 1, N_FL
       do nc = 1, ncheck
          call Hop_mod_gpu_element(nf, nc, ek, nh, act, P, rm, cm)
          if (act) then
             nel = nel + 1; npv = npv + nh
             if (ek == 1) nrv = nrv + 4*nh*nh
             if (ek == 2) ncv = ncv + 4*nh*nh
          endif
       enddo
    enddo
    g_nel = nel
    allocate(hop_meta(4,max(nel,1)), hop_P(max(npv,1)), hop_Rval(max(nrv,1)), hop_Cval(max(ncv,1)))
    hop_meta = 0; hop_P = 0; hop_Rval = 0.d0; hop_Cval = cmplx(0.d0,0.d0,kind(0.d0))
    e = 0; npv = 0; nrv = 0; ncv = 0
    do nf = 1, N_FL
       hop_start(nf) = e + 1
       do nc = 1, ncheck
          call Hop_mod_gpu_element(nf, nc, ek, nh, act, P, rm, cm)
          if (.not. act) cycle
          ! Lower triangle from the upper one = what ZSLHEMM/ZDSLSYMM ('U') see.  Exception (measured, see header):
          ! real elements with nh > 8 are applied by ZDSLSYMM through ZLARCM with the FULL stored matrix, so keep it as is.
          do k = 1, 4
             if (ek == 1 .and. nh > 8) exit
             do j = 1, nh
                do i = j+1, nh
                   if (ek == 1) then
                      rm(i,j,k) = rm(j,i,k)
                   else
                      cm(i,j,k) = conjg(cm(j,i,k))
                   endif
                enddo
             enddo
          enddo
          if (ek == 2 .and. nh > 8) then
             ! ZHEMM (ALF, nd > 8) ignores Im(diagonal); ticket 204 applies these elements with ZGEMM on the full matrix
             do k = 1, 4
                do i = 1, nh
                   cm(i,i,k) = cmplx(real(cm(i,i,k)), 0.d0, kind(0.d0))
                enddo
             enddo
          endif
          e = e + 1
          hop_meta(1,e) = ek
          hop_meta(2,e) = nh
          hop_meta(3,e) = npv
          hop_P(npv+1:npv+nh) = P
          npv = npv + nh
          if (ek == 1) then
             hop_meta(4,e) = nrv
             do k = 1, 4
                hop_Rval(nrv+1:nrv+nh*nh) = reshape(rm(:,:,k), (/nh*nh/))
                nrv = nrv + nh*nh
             enddo
          else
             hop_meta(4,e) = ncv
             do k = 1, 4
                hop_Cval(ncv+1:ncv+nh*nh) = reshape(cm(:,:,k), (/nh*nh/))
                ncv = ncv + nh*nh
             enddo
          endif
       enddo
    enddo
    hop_start(N_FL+1) = e + 1
  end subroutine build_hopping

!--------------------------------------------------------------------
!> Checkerboard groups and dense auxiliaries (ticket 204), from the flat element list.
!--------------------------------------------------------------------
  Subroutine build_groups()
    Implicit none
    Integer :: nf, e, e1, e2, ng, pass, nd, po, l, nx, xo, k, i, j, mo
    Integer, allocatable :: mark(:), gp(:)
    Logical :: ok

    allocate(hop_gstart(g_NFL+1), hop_big(g_NFL), hop_aux(2,max(g_nel,1)), mark(max(g_Ndim,1)))
    hop_big = .false.
    ! pass 1 counts the groups, pass 2 stores them
    ng = 0
    do pass = 1, 2
       if (pass == 2) then
          allocate(gp(ng+1)); gp = 0
          g_ngroups = ng
          ng = 0
       endif
       mark = 0
       do nf = 1, g_NFL
          e1 = hop_start(nf); e2 = hop_start(nf+1) - 1
          if (pass == 2) hop_gstart(nf) = ng + 1
          e = e1
          do while (e <= e2)
             ng = ng + 1
             if (pass == 2) gp(ng) = e
             nd = hop_meta(2,e)
             if (nd > 8) then
                hop_big(nf) = .true.
                e = e + 1
                cycle
             endif
             ! small group: extend while the next element is small and disjoint from all sites of this group
             do
                po = hop_meta(3,e)
                do l = 1, nd
                   mark(hop_P(po+l)) = ng
                enddo
                e = e + 1
                if (e > e2) exit
                if (hop_meta(2,e) > 8) exit
                ok = .true.
                po = hop_meta(3,e)
                do l = 1, hop_meta(2,e)
                   if (mark(hop_P(po+l)) == ng) ok = .false.
                enddo
                if (.not. ok) exit
                nd = hop_meta(2,e)
             enddo
          enddo
       enddo
       if (pass == 2) then
          hop_gstart(g_NFL+1) = ng + 1
          gp(ng+1) = g_nel + 1
       endif
    enddo
    allocate(hop_gptr(ng+1)); hop_gptr = gp; deallocate(gp)

    ! dense auxiliaries
    nx = 0
    do e = 1, g_nel
       if (hop_meta(1,e) == 1 .and. hop_meta(2,e) > 8) nx = nx + 2*hop_meta(2,e)**2
    enddo
    allocate(hop_At(max(nx,1)), hop_Ac(max(nx,1)))
    hop_At = 0.d0; hop_Ac = cmplx(0.d0,0.d0,kind(0.d0))
    xo = 0
    do e = 1, g_nel
       nd = hop_meta(2,e); po = hop_meta(3,e)
       hop_aux(1,e) = -1
       hop_aux(2,e) = 1
       do l = 1, nd
          if (hop_P(po+l) /= l) hop_aux(2,e) = 0
       enddo
       if (hop_meta(1,e) == 1 .and. nd > 8) then
          hop_aux(1,e) = xo
          do k = 1, 2
             mo = hop_meta(4,e) + (k-1)*nd*nd
             do j = 1, nd
                do i = 1, nd
                   hop_At(xo + (k-1)*nd*nd + j + (i-1)*nd) = hop_Rval(mo + i + (j-1)*nd)
                   hop_Ac(xo + (k-1)*nd*nd + i + (j-1)*nd) = cmplx(hop_Rval(mo + i + (j-1)*nd), 0.d0, kind(0.d0))
                enddo
             enddo
          enddo
          xo = xo + 2*nd*nd
       endif
    enddo
    deallocate(mark)
  end subroutine build_groups

!--------------------------------------------------------------------
!> Per matrix index and flavour: the (operator, position) entries in operator order (ticket 204).
!--------------------------------------------------------------------
  Subroutine build_vrows()
    Implicit none
    Integer :: n, nf, i, r, np, p, N1, po
    Integer, allocatable :: cnt(:)

    allocate(vr_ptr(g_Ndim+1, g_NFL))
    allocate(cnt(g_Ndim))
    np = 0
    do nf = 1, g_NFL
       do n = 1, g_Nop
          np = np + V_meta(1,n,nf)
       enddo
    enddo
    allocate(vr_ent(2, max(np,1))); vr_ent = 0
    p = 1
    do nf = 1, g_NFL
       cnt = 0
       do n = 1, g_Nop
          N1 = V_meta(1,n,nf); po = V_meta(2,n,nf)
          do i = 1, N1
             r = V_P(po+i)
             cnt(r) = cnt(r) + 1
          enddo
       enddo
       do r = 1, g_Ndim
          vr_ptr(r,nf) = p
          p = p + cnt(r)
       enddo
       vr_ptr(g_Ndim+1,nf) = p
       ! fill in operator order (n ascending), running position per index
       do r = 1, g_Ndim
          cnt(r) = vr_ptr(r,nf)
       enddo
       do n = 1, g_Nop
          N1 = V_meta(1,n,nf); po = V_meta(2,n,nf)
          do i = 1, N1
             r = V_P(po+i)
             vr_ent(1,cnt(r)) = n
             vr_ent(2,cnt(r)) = i
             cnt(r) = cnt(r) + 1
          enddo
       enddo
    enddo
    deallocate(cnt)
  end subroutine build_vrows

!--------------------------------------------------------------------
  Subroutine build_V()
    Implicit none
    Integer :: n, nf, nt, nsum, ne, np, ne0
    Complex (Kind=Kind(0.d0)), allocatable :: E(:,:)

    allocate(V_meta(4,g_Nop,N_FL), V_type(g_Nop))
    np = 0; ne = 0
    do nf = 1, N_FL
       do n = 1, g_Nop
          nt = Op_V(n,nf)%type
          V_meta(1,n,nf) = Op_V(n,nf)%N
          V_meta(2,n,nf) = np
          V_meta(3,n,nf) = ne
          V_meta(4,n,nf) = v_skip(Op_V(n,nf))
          np = np + Op_V(n,nf)%N
          ne = ne + Op_V(n,nf)%N*(2*nt+1)
       enddo
    enddo
    allocate(V_P(max(np,1)), V_exp(max(ne,1)))
    V_P = 0; V_exp = cmplx(1.d0,0.d0,kind(0.d0))
    do n = 1, g_Nop
       V_type(n) = Op_V(n,1)%type
    enddo
    do nf = 1, N_FL
       do n = 1, g_Nop
          nt = Op_V(n,nf)%type
          nsum = Op_V(n,nf)%N
          V_P(V_meta(2,n,nf)+1:V_meta(2,n,nf)+nsum) = Op_V(n,nf)%P(1:nsum)
          call Op_get_E_exp(Op_V(n,nf), E)          ! private table of Operator_mod, same numbers ALF applies
          ne0 = V_meta(3,n,nf)
          V_exp(ne0+1:ne0+nsum*(2*nt+1)) = reshape(E, (/nsum*(2*nt+1)/))
          V_exp(ne0+nt*nsum+1:ne0+(nt+1)*nsum) = cmplx(1.d0,0.d0,kind(0.d0))   ! sp = 0: unused
          deallocate(E)
       enddo
    enddo
  end subroutine build_V

!--------------------------------------------------------------------
  ! 1 iff Op_mmultL/R return early for this operator (constant coupling).
  Integer function v_skip(Op)
    Implicit none
    Type (Operator), intent(in) :: Op
    v_skip = 0
    if ( abs(Op%g) < Eps_machine ) v_skip = 1
  end function v_skip

!--------------------------------------------------------------------
  Subroutine alloc_buffers()
    Implicit none
    Integer :: lw
    allocate( G_d(g_Ndim,g_Ndim,g_NFL) )
    allocate( udv_U(g_Ndim,g_Ndim,g_NSTM,g_NFLeff), udv_V(g_Ndim,g_Ndim,g_NSTM,g_NFLeff) )
    allocate( udw_U(g_Ndim,g_Ndim,2,g_NFLeff), udw_V(g_Ndim,g_Ndim,2,g_NFLeff) )
#if defined(STABLOG)
    allocate( udv_L(g_Ndim,g_NSTM,g_NFLeff), udw_L(g_Ndim,2,g_NFLeff) )
    allocate( udv_D(1,1,1), udw_D(1,1,1) )
#else
    allocate( udv_D(g_Ndim,g_NSTM,g_NFLeff), udw_D(g_Ndim,2,g_NFLeff) )
    allocate( udv_L(1,1,1), udw_L(1,1,1) )
#endif
    if (have_tau) allocate( GT0_d(g_Ndim,g_Ndim,g_NFL), G0T_d(g_Ndim,g_Ndim,g_NFL), GTT_d(g_Ndim,g_Ndim,g_NFL), G00_d(g_Ndim,g_Ndim,g_NFL) )
    lw = gpu_scratch_lwork(g_Ndim)
    allocate( scr_c(g_Ndim,g_Ndim,2), scr_vec(g_Ndim,4), scr_r(g_Ndim,2), scr_ipiv(g_Ndim), scr_work(lw) )
  end subroutine alloc_buffers

!--------------------------------------------------------------------
  Subroutine map_to_device()
    Implicit none
    ! One-time host->device copies of the constant data (unavoidable: the device cannot read the host objects).
    !$omp target enter data map(to: hop_start, hop_meta, hop_P, hop_Rval, hop_Cval)
    !$omp target enter data map(to: V_meta, V_type, V_P, V_exp)
    !$omp target enter data map(to: nsigma_d, nsigma_t)
    ! ticket 204: group table, dense auxiliaries and per-site operator lists (constant data, uploaded once)
    !$omp target enter data map(to: hop_gptr, hop_aux, hop_At, hop_Ac, vr_ptr, vr_ent)
    ! Buffers that later tickets fill: device allocation only, no transfer.
    !$omp target enter data map(alloc: G_d, udv_U, udv_V, udw_U, udw_V)
    !$omp target enter data map(alloc: udv_D, udv_L, udw_D, udw_L)
    !$omp target enter data map(alloc: scr_c, scr_vec, scr_r, scr_ipiv, scr_work)
    if (have_tau) then
       !$omp target enter data map(alloc: GT0_d, G0T_d, GTT_d, G00_d)
    endif
  end subroutine map_to_device

!--------------------------------------------------------------------
  Subroutine gpu_state_finalize()
    Implicit none
    if (.not. gpu_active) return
    if (upd_cnt(10) > 0) call report_update_counts()
#if !defined(ALF_LA_ONEMKL)
    call la_finalize()
#endif
    if (have_tau) then
       !$omp target exit data map(delete: GT0_d, G0T_d, GTT_d, G00_d)
       deallocate(GT0_d, G0T_d, GTT_d, G00_d)
    endif
    !$omp target exit data map(delete: scr_c, scr_vec, scr_r, scr_ipiv, scr_work)
    !$omp target exit data map(delete: udv_D, udv_L, udw_D, udw_L)
    !$omp target exit data map(delete: G_d, udv_U, udv_V, udw_U, udw_V)
    !$omp target exit data map(delete: hop_gptr, hop_aux, hop_At, hop_Ac, vr_ptr, vr_ent)
    !$omp target exit data map(delete: nsigma_d, nsigma_t)
    !$omp target exit data map(delete: V_meta, V_type, V_P, V_exp)
    !$omp target exit data map(delete: hop_start, hop_meta, hop_P, hop_Rval, hop_Cval)
    deallocate(scr_c, scr_vec, scr_r, scr_ipiv, scr_work, udv_D, udv_L, udw_D, udw_L)
    deallocate(G_d, udv_U, udv_V, udw_U, udw_V, nsigma_d, nsigma_t)
    deallocate(V_meta, V_type, V_P, V_exp, hop_start, hop_meta, hop_P, hop_Rval, hop_Cval)
    deallocate(hop_gstart, hop_gptr, hop_big, hop_aux, hop_At, hop_Ac, vr_ptr, vr_ent)
    gpu_active = .false.
  end subroutine gpu_state_finalize

!--------------------------------------------------------------------
!> Ticket 205: transfer and launch counts of the device updates (see upd_cnt), per run.
  Subroutine report_update_counts()
    Implicit none
    Character (len=256) :: line
    write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'GPU update counts: WRAPGRUP/DO calls ', upd_cnt(10), ', site proposals tested ', &
         & upd_cnt(1), ', accepted ', upd_cnt(2), ', site-wrap launches ', upd_cnt(3), ', update launches ', upd_cnt(4)
    call info_line(line)
    write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'GPU update transfers: ratio-block D->H elements ', upd_cnt(5), &
         & ', field H->D scalars ', upd_cnt(6), ', interim G copies H->D ', upd_cnt(7), ' and D->H ', upd_cnt(8), &
         & ' (', upd_cnt(9), ' bytes, INTERIM 205)'
    call info_line(line)
  end subroutine report_update_counts

!--------------------------------------------------------------------
  Subroutine report_layout()
    Implicit none
    Integer(int64) :: tot
    tot = 0
    call lay('hop_start', 'int', shape(hop_start), 4_int64*size(hop_start), tot)
    call lay('hop_meta', 'int', shape(hop_meta), 4_int64*size(hop_meta), tot)
    call lay('hop_P', 'int', shape(hop_P), 4_int64*size(hop_P), tot)
    call lay('hop_Rval', 'real', shape(hop_Rval), 8_int64*size(hop_Rval), tot)
    call lay('hop_Cval', 'complex', shape(hop_Cval), 16_int64*size(hop_Cval), tot)
    call lay('V_meta', 'int', shape(V_meta), 4_int64*size(V_meta), tot)
    call lay('V_type', 'int', shape(V_type), 4_int64*size(V_type), tot)
    call lay('V_P', 'int', shape(V_P), 4_int64*size(V_P), tot)
    call lay('V_exp', 'complex', shape(V_exp), 16_int64*size(V_exp), tot)
    call lay('nsigma_d', 'complex', shape(nsigma_d), 16_int64*size(nsigma_d), tot)
    call lay('nsigma_t', 'int', shape(nsigma_t), 4_int64*size(nsigma_t), tot)
    call lay('hop_gptr', 'int', shape(hop_gptr), 4_int64*size(hop_gptr), tot)
    call lay('hop_aux', 'int', shape(hop_aux), 4_int64*size(hop_aux), tot)
    call lay('hop_At', 'real', shape(hop_At), 8_int64*size(hop_At), tot)
    call lay('hop_Ac', 'complex', shape(hop_Ac), 16_int64*size(hop_Ac), tot)
    call lay('vr_ptr', 'int', shape(vr_ptr), 4_int64*size(vr_ptr), tot)
    call lay('vr_ent', 'int', shape(vr_ent), 4_int64*size(vr_ent), tot)
    call lay('G_d', 'complex', shape(G_d), 16_int64*size(G_d), tot)
    call lay('udv_U', 'complex', shape(udv_U), 16_int64*size(udv_U), tot)
    call lay('udv_V', 'complex', shape(udv_V), 16_int64*size(udv_V), tot)
#if defined(STABLOG)
    call lay('udv_L', 'real', shape(udv_L), 8_int64*size(udv_L), tot)
    call lay('udw_L', 'real', shape(udw_L), 8_int64*size(udw_L), tot)
#else
    call lay('udv_D', 'complex', shape(udv_D), 16_int64*size(udv_D), tot)
    call lay('udw_D', 'complex', shape(udw_D), 16_int64*size(udw_D), tot)
#endif
    call lay('udw_U', 'complex', shape(udw_U), 16_int64*size(udw_U), tot)
    call lay('udw_V', 'complex', shape(udw_V), 16_int64*size(udw_V), tot)
    if (have_tau) then
       call lay('GT0_d', 'complex', shape(GT0_d), 16_int64*size(GT0_d), tot)
       call lay('G0T_d', 'complex', shape(G0T_d), 16_int64*size(G0T_d), tot)
       call lay('GTT_d', 'complex', shape(GTT_d), 16_int64*size(GTT_d), tot)
       call lay('G00_d', 'complex', shape(G00_d), 16_int64*size(G00_d), tot)
    endif
    call lay('scr_c', 'complex', shape(scr_c), 16_int64*size(scr_c), tot)
    call lay('scr_vec', 'complex', shape(scr_vec), 16_int64*size(scr_vec), tot)
    call lay('scr_r', 'real', shape(scr_r), 8_int64*size(scr_r), tot)
    call lay('scr_ipiv', 'int', shape(scr_ipiv), 4_int64*size(scr_ipiv), tot)
    call lay('scr_work', 'complex', shape(scr_work), 16_int64*size(scr_work), tot)
    block
      Character (len=256) :: line
      write(line,'(A,I0,A)') 'GPU mirror total device bytes: ', tot, ' (arrays listed above)'
      call info_line(line)
    end block
  end subroutine report_layout

  Subroutine lay(name, typ, shp, nbytes, tot)
    Implicit none
    Character (len=*), intent(in) :: name, typ
    Integer, intent(in) :: shp(:)
    Integer(int64), intent(in) :: nbytes
    Integer(int64), intent(inout) :: tot
    Character (len=256) :: line, s
    Integer :: i
    s = ''
    do i = 1, size(shp)
       if (i > 1) s = trim(s)//','
       write(s(len_trim(s)+1:),'(I0)') shp(i)
    enddo
    write(line,'(A,A,A,A,A,A,A,I0)') 'GPU mirror array ', trim(name), ' (', trim(s), ') ', trim(typ), ' bytes ', nbytes
    call info_line(line)
    tot = tot + nbytes
  end subroutine lay

!====================================================================
!  Self-test (ALF_GPU_SELFTEST=1)
!====================================================================
  Real (Kind=Kind(0.d0)) function rnd(s)
    Implicit none
    Integer(int64), intent(inout) :: s
    s = mod(1103515245_int64*s + 12345_int64, 2147483648_int64)
    rnd = real(s,kind(0.d0))/2147483648.d0 - 0.5d0
  end function rnd

  Subroutine rand_mat(M, s)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(out) :: M(:,:)
    Integer(int64), intent(inout) :: s
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: re, im
    do j = 1, size(M,2)
       do i = 1, size(M,1)
          re = rnd(s)
          im = rnd(s)
          M(i,j) = cmplx(re, im, kind(0.d0))
       enddo
    enddo
  end subroutine rand_mat

  ! NaN/Inf test on the bit pattern (exponent all ones); immune to fast-math options of the compilers
  Logical function is_finite(x)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(in) :: x
    is_finite = ( iand(ishft(transfer(x, 0_int64), -52), 2047_int64) /= 2047_int64 )
  end function is_finite

  Logical function finite_c(z)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: z
    finite_c = is_finite(real(z)) .and. is_finite(aimag(z))
  end function finite_c

  ! max |A-B| / max |B|; any non-finite entry in A or B gives huge(1.d0), so that no tolerance test can pass
  Real (Kind=Kind(0.d0)) function rel_diff(A, B)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:), B(:,:)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: d, md, mb
    md = 0.d0; mb = 0.d0
    do j = 1, size(A,2)
       do i = 1, size(A,1)
          if (.not. (finite_c(A(i,j)) .and. finite_c(B(i,j)))) then
             rel_diff = huge(1.d0); return
          endif
          d = abs(A(i,j)-B(i,j)); if (d > md) md = d
          d = abs(B(i,j));        if (d > mb) mb = d
       enddo
    enddo
    rel_diff = md / max(mb, tiny(1.d0))
  end function rel_diff

  ! number of non-finite entries of a real(8)/complex(8) vector
  Integer function n_nonfinite_r(n, x)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(in) :: x(n)
    Integer :: i
    n_nonfinite_r = 0
    do i = 1, n
       if (.not. is_finite(x(i))) n_nonfinite_r = n_nonfinite_r + 1
    enddo
  end function n_nonfinite_r

  Integer function n_nonfinite_c(n, x)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: x(n)
    Integer :: i
    n_nonfinite_c = 0
    do i = 1, n
       if (.not. finite_c(x(i))) n_nonfinite_c = n_nonfinite_c + 1
    enddo
  end function n_nonfinite_c

  Complex (Kind=Kind(0.d0)) function a_el(e, k, i, j)
    Implicit none
    Integer, intent(in) :: e, k, i, j
    Integer :: nd, mo
    nd = hop_meta(2,e); mo = hop_meta(4,e) + (k-1)*nd*nd + i + (j-1)*nd
    if (hop_meta(1,e) == 1) then
       a_el = cmplx(hop_Rval(mo), 0.d0, kind(0.d0))
    else
       a_el = hop_Cval(mo)
    endif
  end function a_el

  ! M(P,:) = A_k * M(P,:)  (own loops on the flat arrays, host)
  Subroutine el_left(e, k, M)
    Implicit none
    Integer, intent(in) :: e, k
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: W(:,:)
    Complex (Kind=Kind(0.d0)) :: z
    Integer :: nd, po, i, j, l
    nd = hop_meta(2,e); po = hop_meta(3,e)
    allocate(W(nd,size(M,2)))
    do j = 1, size(M,2)
       do i = 1, nd
          W(i,j) = M(hop_P(po+i),j)
       enddo
    enddo
    do j = 1, size(M,2)
       do i = 1, nd
          z = cmplx(0.d0,0.d0,kind(0.d0))
          do l = 1, nd
             z = z + a_el(e,k,i,l)*W(l,j)
          enddo
          M(hop_P(po+i),j) = z
       enddo
    enddo
    deallocate(W)
  end subroutine el_left

  ! M(:,P) = M(:,P) * A_k
  Subroutine el_right(e, k, M)
    Implicit none
    Integer, intent(in) :: e, k
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: W(:,:)
    Complex (Kind=Kind(0.d0)) :: z
    Integer :: nd, po, i, j, l
    nd = hop_meta(2,e); po = hop_meta(3,e)
    allocate(W(size(M,1),nd))
    do j = 1, nd
       do i = 1, size(M,1)
          W(i,j) = M(i,hop_P(po+j))
       enddo
    enddo
    do j = 1, nd
       do i = 1, size(M,1)
          z = cmplx(0.d0,0.d0,kind(0.d0))
          do l = 1, nd
             z = z + W(i,l)*a_el(e,k,l,j)
          enddo
          M(i,hop_P(po+j)) = z
       enddo
    enddo
    deallocate(W)
  end subroutine el_right

  ! Same Mat_subroutines routine Hop_mod calls, fed with slices of the flat arrays (tests offsets, order, data).
  Subroutine el_lib(e, k, side, M)
    Implicit none
    Integer, intent(in) :: e, k
    Character (len=1), intent(in) :: side
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer :: nd, po, mo
    nd = hop_meta(2,e); po = hop_meta(3,e); mo = hop_meta(4,e) + (k-1)*nd*nd
    if (hop_meta(1,e) == 1) then
       call ZDSLSYMM(side, 'U', nd, size(M,1), size(M,2), hop_Rval(mo+1:mo+nd*nd), hop_P(po+1:po+nd), M)
    else
       call ZSLHEMM (side, 'U', nd, size(M,1), size(M,2), hop_Cval(mo+1:mo+nd*nd), hop_P(po+1:po+nd), M)
    endif
  end subroutine el_lib

  ! Apply the elements of flavour nf in the given order and variant.  mode: 1 own loops, 2 library on slices
  Subroutine flat_apply(nf, op, M, mode)
    Implicit none
    Integer, intent(in) :: nf, op, mode
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer :: e, e1, e2
    e1 = hop_start(nf); e2 = hop_start(nf+1) - 1
    select case (op)
    case (1)   ! Hop_mod_mmthr:    left, mat, descending
       do e = e2, e1, -1
          call app(e, 1, 'L', M, mode)
       enddo
    case (2)   ! Hop_mod_mmthr_m1: left, invmat, ascending
       do e = e1, e2
          call app(e, 2, 'L', M, mode)
       enddo
    case (3)   ! Hop_mod_mmthl:    right, mat, ascending
       do e = e1, e2
          call app(e, 1, 'R', M, mode)
       enddo
    case (4)   ! Hop_mod_mmthlc:   left, mat, ascending
       do e = e1, e2
          call app(e, 1, 'L', M, mode)
       enddo
    case (5)   ! Hop_mod_mmthl_m1: right, invmat, descending
       do e = e2, e1, -1
          call app(e, 2, 'R', M, mode)
       enddo
    case (6)   ! Hop_mod_Symm:     per element left mat_1D2 then right invmat_1D2, descending
       do e = e2, e1, -1
          call app(e, 3, 'L', M, mode)
          call app(e, 4, 'R', M, mode)
       enddo
    end select
  contains
    Subroutine app(e, k, side, M, mode)
      Integer, intent(in) :: e, k, mode
      Character (len=1), intent(in) :: side
      Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
      if (mode == 1) then
         if (side == 'L') then
            call el_left(e, k, M)
         else
            call el_right(e, k, M)
         endif
      else
         call el_lib(e, k, side, M)
      endif
    end subroutine app
  end subroutine flat_apply

  Subroutine selftest_fail(msg)
    Implicit none
    Character (len=*), intent(in) :: msg
    write(error_unit,'(A,A)') ' GPU self-test FAILED: ', msg
    call info_line('GPU self-test FAILED: '//trim(msg))
    call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
  end subroutine selftest_fail

  Subroutine gpu_selftest()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-14
    Integer(int64) :: s
    Integer :: nf, nfe, op, n, sp, nt, i, nbad, nsp
    Type (Operator) :: Otmp
    Real (Kind=Kind(0.d0)) :: dw(2), dskip
    Complex (Kind=Kind(0.d0)), allocatable :: M0(:,:), M1(:,:), M2(:,:), In3(:,:,:), Out3(:,:,:)
    Real (Kind=Kind(0.d0)) :: dloop(6), dlib(6), dv(4), d, dnl, dnr
    Character (len=256) :: line
    Character (len=1) :: cop
    Character (len=10), parameter :: opname(6) = (/ 'mmthr     ', 'mmthr_m1  ', 'mmthl     ', 'mmthlc    ', &
         & 'mmthl_m1  ', 'Symm      ' /)

    s = 20251028_int64
    allocate(M0(Ndim,Ndim), M1(Ndim,Ndim), M2(Ndim,Ndim))
    dloop = 0.d0; dlib = 0.d0

    ! ---- inputs must be finite (max/maxval would drop NaN) ----
    nbad = n_nonfinite_r(size(hop_Rval), hop_Rval) + n_nonfinite_c(size(hop_Cval), hop_Cval) &
         & + n_nonfinite_c(size(V_exp), V_exp) + n_nonfinite_c(size(nsigma_d), nsigma_d)
    write(line,'(A,I0,A)') 'GPU self-test non-finite entries in mirror arrays: ', nbad, ' (expected 0)'
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    if (nbad /= 0) call selftest_fail('non-finite value in the mirror')

    ! ---- hopping: flat arrays vs Hop_mod routines ----
    do nf = 1, N_FL
       do op = 1, 5
          call rand_mat(M0, s)
          M1 = M0
          select case (op)
          case (1); call Hop_mod_mmthr   (M1, nf, 1)
          case (2); call Hop_mod_mmthr_m1(M1, nf, 1)
          case (3); call Hop_mod_mmthl   (M1, nf, 1)
          case (4); call Hop_mod_mmthlc  (M1, nf, 1)
          case (5); call Hop_mod_mmthl_m1(M1, nf, 1)
          end select
          M2 = M0; call flat_apply(nf, op, M2, 1); dloop(op) = max(dloop(op), rel_diff(M2, M1))
          M2 = M0; call flat_apply(nf, op, M2, 2); dlib(op)  = max(dlib(op),  rel_diff(M2, M1))
       enddo
    enddo
    ! Symm: Out = e^{-dtau T/2} In e^{dtau T/2} for the calculated flavours
    allocate(In3(Ndim,Ndim,N_FL), Out3(Ndim,Ndim,N_FL))
    do nf = 1, N_FL
       call rand_mat(M0, s); In3(:,:,nf) = M0
    enddo
    call Hop_mod_Symm(Out3, In3, 1)
    do nfe = 1, N_FL_eff
       nf = Calc_Fl_map(nfe)
       M2 = In3(:,:,nf); call flat_apply(nf, 6, M2, 1); dloop(6) = max(dloop(6), rel_diff(M2, Out3(:,:,nf)))
       M2 = In3(:,:,nf); call flat_apply(nf, 6, M2, 2); dlib(6)  = max(dlib(6),  rel_diff(M2, Out3(:,:,nf)))
    enddo
    deallocate(In3, Out3)
    do op = 1, 6
       write(line,'(A,A,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU self-test hopping ', opname(op), ' max rel diff: own loops ', &
            & dloop(op), ', same library routine ', dlib(op), ' (tol ', tol
       call info_line(line); write(output_unit,'(1X,A)') trim(line)
       if (.not. (dloop(op) <= tol .and. dlib(op) <= tol)) call selftest_fail('hopping '//trim(opname(op)))
    enddo

    ! ---- Op_V: flat table vs Op_mmultL/R ----
    dv = 0.d0; nsp = 0
    do nf = 1, N_FL
       do n = 1, g_Nop
          nt = V_type(n)
          do sp = -nt, nt
             if (sp == 0) cycle
             nsp = nsp + 1
             do i = 1, 2
                if (i == 1) cop = 'n'
                if (i == 2) cop = 'c'
                call rand_mat(M0, s)
                ! right side: rows scaled, Op_mmultR (Mat = Op * Mat)
                M1 = M0; call Op_mmultR(M1, Op_V(n,nf), cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), cop, 1)
                M2 = M0; call flat_V_R(n, nf, sp, cop, M2, V_meta(4,n,nf))
                dv(2*(i-1)+1) = max(dv(2*(i-1)+1), rel_diff(M2, M1))
                ! left side: columns scaled, Op_mmultL with sign +1 and -1
                M1 = M0; call Op_mmultL(M1, Op_V(n,nf), cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), cop, 1, 1)
                M2 = M0; call flat_V_L(n, nf, sp, cop, M2, V_meta(4,n,nf))
                d = rel_diff(M2, M1)
                M1 = M0; call Op_mmultL(M1, Op_V(n,nf), cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), cop, 1, -1)
                M2 = M0; call flat_V_L(n, nf, -sp, cop, M2, V_meta(4,n,nf))
                d = max(d, rel_diff(M2, M1))
                dv(2*(i-1)+2) = max(dv(2*(i-1)+2), d)
             enddo
          enddo
       enddo
    enddo
    ! wrap meaning (Op_Wrapup/Op_Wrapdo, N_type 1, no skip) and skip meaning (Op_mmultL/R early return for tiny g):
    ! a shallow copy of an operator with g set to 1e-20 (pointer components shared, only the scalar g changes)
    dw = 0.d0; dskip = 0.d0
    do nf = 1, N_FL
       do n = 1, g_Nop
          nt = V_type(n)
          do sp = -nt, nt
             if (sp == 0) cycle
             call rand_mat(M0, s)
             M1 = M0; call Op_Wrapup(M1, Op_V(n,nf), cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), Ndim, 1, 1)
             M2 = M0; call flat_V_wrap(n, nf, sp, 1, M2); dw(1) = max(dw(1), rel_diff(M2, M1))
             M1 = M0; call Op_Wrapdo(M1, Op_V(n,nf), cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), Ndim, 1, 1)
             M2 = M0; call flat_V_wrap(n, nf, sp, -1, M2); dw(2) = max(dw(2), rel_diff(M2, M1))
          enddo
       enddo
    enddo
    Otmp = Op_V(1,1)
    Otmp%g = cmplx(1.d-20, 0.d0, kind(0.d0))
    nt = V_type(1)
    if (v_skip(Otmp) /= 1 .or. v_skip(Op_V(1,1)) /= V_meta(4,1,1)) call selftest_fail('skip flag')
    do sp = -nt, nt
       if (sp == 0) cycle
       call rand_mat(M0, s)
       M1 = M0; call Op_mmultR(M1, Otmp, cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), 'n', 1)
       M2 = M0; call flat_V_R(1, 1, sp, 'n', M2, v_skip(Otmp)); dskip = max(dskip, rel_diff(M2, M1))
       if (rel_diff(M1, M0) /= 0.d0) call selftest_fail('Op_mmultR did not return early for tiny g')
       M1 = M0; call Op_mmultL(M1, Otmp, cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), 'n', 1, 1)
       M2 = M0; call flat_V_L(1, 1, sp, 'n', M2, v_skip(Otmp)); dskip = max(dskip, rel_diff(M2, M1))
       ! wrap ignores g: Op_Wrapup on the tiny-g copy applies the (unchanged) table
       M1 = M0; call Op_Wrapup(M1, Otmp, cmplx(real(sp,kind(0.d0)),0.d0,kind(0.d0)), Ndim, 1, 1)
       M2 = M0; call flat_V_wrap(1, 1, sp, 1, M2); dskip = max(dskip, rel_diff(M2, M1))
       if (abs(Op_V(1,1)%g) > 0.d0 .and. rel_diff(M1, M0) == 0.d0) call selftest_fail('Op_Wrapup skipped for tiny g')
    enddo
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU self-test Op_V wrap: Op_Wrapup ', dw(1), ', Op_Wrapdo ', dw(2), &
         & '; tiny-g copy (L/R skip honoured, wrap ignores it) ', dskip, ' (tol ', tol
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    if (.not. (dw(1) <= tol .and. dw(2) <= tol .and. dskip <= tol)) call selftest_fail('Op_V wrap/skip')

    dnr = max(dv(1), dv(3)); dnl = max(dv(2), dv(4))
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU self-test Op_V (', nsp, &
         & ' operator/field cases, N and C) max rel diff: right (Op_mmultR) ', dnr, &
         & ', left (Op_mmultL, both signs) ', dnl, ' (tol ', tol
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    if (.not. (dnr <= tol .and. dnl <= tol)) call selftest_fail('Op_V table')

    ! ---- device round trip: target update from == host-built array, bit for bit ----
    nbad = 0
    call rt_i('hop_start', size(hop_start), hop_start, .false., nbad)
    call rt_i('hop_meta', size(hop_meta), hop_meta, .false., nbad)
    call rt_i('hop_P', size(hop_P), hop_P, .false., nbad)
    call rt_r('hop_Rval', size(hop_Rval), hop_Rval, .false., nbad)
    call rt_c('hop_Cval', size(hop_Cval), hop_Cval, .false., nbad)
    call rt_i('V_meta', size(V_meta), V_meta, .false., nbad)
    call rt_i('V_type', size(V_type), V_type, .false., nbad)
    call rt_i('V_P', size(V_P), V_P, .false., nbad)
    call rt_c('V_exp', size(V_exp), V_exp, .false., nbad)
    call rt_c('nsigma_d', size(nsigma_d), nsigma_d, .false., nbad)
    call rt_i('nsigma_t', size(nsigma_t), nsigma_t, .false., nbad)
    ! buffers: fill a pattern, copy to the device, poison the host, copy back, compare, then zero both
    call rt_c('G_d', size(G_d), G_d, .true., nbad)
    call rt_c('udv_U', size(udv_U), udv_U, .true., nbad)
    call rt_c('udv_V', size(udv_V), udv_V, .true., nbad)
    call rt_c('udw_U', size(udw_U), udw_U, .true., nbad)
    call rt_c('udw_V', size(udw_V), udw_V, .true., nbad)
#if defined(STABLOG)
    call rt_r('udv_L', size(udv_L), udv_L, .true., nbad)
    call rt_r('udw_L', size(udw_L), udw_L, .true., nbad)
#else
    call rt_c('udv_D', size(udv_D), udv_D, .true., nbad)
    call rt_c('udw_D', size(udw_D), udw_D, .true., nbad)
#endif
    if (have_tau) then
       call rt_c('GT0_d', size(GT0_d), GT0_d, .true., nbad)
       call rt_c('G0T_d', size(G0T_d), G0T_d, .true., nbad)
       call rt_c('GTT_d', size(GTT_d), GTT_d, .true., nbad)
       call rt_c('G00_d', size(G00_d), G00_d, .true., nbad)
    endif
    call rt_c('scr_c', size(scr_c), scr_c, .true., nbad)
    call rt_c('scr_vec', size(scr_vec), scr_vec, .true., nbad)
    call rt_r('scr_r', size(scr_r), scr_r, .true., nbad)
    call rt_i('scr_ipiv', size(scr_ipiv), scr_ipiv, .true., nbad)
    call rt_c('scr_work', size(scr_work), scr_work, .true., nbad)
    write(line,'(A,I0,A)') 'GPU self-test device round trip: arrays with bit mismatch: ', nbad, ' (expected 0)'
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    if (nbad /= 0) call selftest_fail('device round trip')

    call device_kernel_check()
    deallocate(M0, M1, M2)
    call info_line('GPU self-test PASSED')
    write(output_unit,'(A)') ' GPU self-test PASSED'
  end subroutine gpu_selftest

  ! Rows of M scaled: Mat(P(i),:) *= E_exp(i,sp)   (Op_mmultR)
  Subroutine flat_V_R(n, nf, sp, cop, M, skip)
    Implicit none
    Integer, intent(in) :: n, nf, sp, skip
    Character (len=1), intent(in) :: cop
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer :: i, N1, po, eo, nt
    Complex (Kind=Kind(0.d0)) :: z
    if (skip == 1) return          ! Op_mmultR returns early for |g| < Eps_machine
    N1 = V_meta(1,n,nf); po = V_meta(2,n,nf); eo = V_meta(3,n,nf); nt = V_type(n)
    do i = 1, N1
       z = V_exp(eo + i + (sp+nt)*N1)
       if (cop == 'c' .or. cop == 'C') z = conjg(z)
       M(V_P(po+i),:) = z*M(V_P(po+i),:)
    enddo
  end subroutine flat_V_R

  ! Columns of M scaled: Mat(:,P(i)) *= E_exp(i,sp)   (Op_mmultL)
  Subroutine flat_V_L(n, nf, sp, cop, M, skip)
    Implicit none
    Integer, intent(in) :: n, nf, sp, skip
    Character (len=1), intent(in) :: cop
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer :: i, N1, po, eo, nt
    Complex (Kind=Kind(0.d0)) :: z
    if (skip == 1) return          ! Op_mmultL returns early for |g| < Eps_machine
    N1 = V_meta(1,n,nf); po = V_meta(2,n,nf); eo = V_meta(3,n,nf); nt = V_type(n)
    do i = 1, N1
       z = V_exp(eo + i + (sp+nt)*N1)
       if (cop == 'c' .or. cop == 'C') z = conjg(z)
       M(:,V_P(po+i)) = z*M(:,V_P(po+i))
    enddo
  end subroutine flat_V_L

  ! Op_Wrapup (dir = +1): rows by table(sp), then columns by table(-sp); Op_Wrapdo (dir = -1): the reverse.  No g check.
  Subroutine flat_V_wrap(n, nf, sp, dir, M)
    Implicit none
    Integer, intent(in) :: n, nf, sp, dir
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer :: i, N1, po, eo, nt
    N1 = V_meta(1,n,nf); po = V_meta(2,n,nf); eo = V_meta(3,n,nf); nt = V_type(n)
    do i = 1, N1
       M(V_P(po+i),:) = V_exp(eo + i + (dir*sp+nt)*N1) * M(V_P(po+i),:)
    enddo
    do i = 1, N1
       M(:,V_P(po+i)) = V_exp(eo + i + (-dir*sp+nt)*N1) * M(:,V_P(po+i))
    enddo
  end subroutine flat_V_wrap

  ! A trivial device kernel that reads the mirror without any map clause (arrays must be present).
  Subroutine device_kernel_check()
    Implicit none
    Integer :: i, isum_h, isum_d
    Real (Kind=Kind(0.d0)) :: rs_h, rs_d, errr
    Character (len=256) :: line
    isum_h = sum(hop_P) + sum(V_P) + sum(nsigma_t)
    rs_h = sum(hop_Rval) + sum(real(V_exp))
    if (.not. is_finite(rs_h)) call selftest_fail('non-finite host checksum')
    isum_d = 0; rs_d = 0.d0
    !$omp target teams distribute parallel do reduction(+:isum_d) map(tofrom: isum_d)
    do i = 1, size(hop_P)
       isum_d = isum_d + hop_P(i)
    enddo
    !$omp target teams distribute parallel do reduction(+:isum_d) map(tofrom: isum_d)
    do i = 1, size(V_P)
       isum_d = isum_d + V_P(i)
    enddo
    !$omp target teams distribute parallel do reduction(+:isum_d) map(tofrom: isum_d)
    do i = 1, size(nsigma_t)
       isum_d = isum_d + nsigma_t(i)
    enddo
    !$omp target teams distribute parallel do reduction(+:rs_d) map(tofrom: rs_d)
    do i = 1, size(hop_Rval)
       rs_d = rs_d + hop_Rval(i)
    enddo
    !$omp target teams distribute parallel do reduction(+:rs_d) map(tofrom: rs_d)
    do i = 1, size(V_exp)
       rs_d = rs_d + real(V_exp(i))
    enddo
    if (.not. is_finite(rs_d)) call selftest_fail('non-finite device checksum')
    errr = abs(rs_d - rs_h) / max(abs(rs_h), tiny(1.d0))
    write(line,'(A,I0,A,I0,A,ES10.3)') 'GPU self-test device kernel: integer checksum host ', isum_h, ' device ', isum_d, &
         & ', real checksum rel diff ', errr
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    if (isum_h /= isum_d .or. .not. (errr <= 1.d-12)) call selftest_fail('device kernel checksum')
  end subroutine device_kernel_check

!--------------------------------------------------------------------
! Round-trip helpers.  fill=.false.: x holds the host-built data; save it, poison x, "target update from", compare.
! fill=.true. (device-only buffers): fill pattern, "target update to", poison, "from", compare, then zero both.
!--------------------------------------------------------------------
  Subroutine rt_c(name, n, x, fill, nbad)
    Implicit none
    Character (len=*), intent(in) :: name
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: x(n)
    Logical, intent(in) :: fill
    Integer, intent(inout) :: nbad
    Complex (Kind=Kind(0.d0)), allocatable :: c(:)
    Integer :: i
    Logical :: same
    allocate(c(n))
    if (fill) then
       do i = 1, n
          x(i) = cmplx(1.d0 + mod(i,977), -1.d0 - mod(i,31)*0.5d0, kind(0.d0))
       enddo
       !$omp target update to(x)
    endif
    c = x
    x = cmplx(1.2345d-300, -7.6543d-300, kind(0.d0))
    !$omp target update from(x)
    same = same_bits_c(n, x, c)
    if (.not. same) then
       nbad = nbad + 1
       write(error_unit,*) 'round trip mismatch: ', name
    endif
    if (fill) then
       x = cmplx(0.d0, 0.d0, kind(0.d0))
       !$omp target update to(x)
    endif
    deallocate(c)
  end subroutine rt_c

  Subroutine rt_r(name, n, x, fill, nbad)
    Implicit none
    Character (len=*), intent(in) :: name
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(inout) :: x(n)
    Logical, intent(in) :: fill
    Integer, intent(inout) :: nbad
    Real (Kind=Kind(0.d0)), allocatable :: c(:)
    Integer :: i
    allocate(c(n))
    if (fill) then
       do i = 1, n
          x(i) = 1.d0 + mod(i,977)*0.25d0
       enddo
       !$omp target update to(x)
    endif
    c = x
    x = 1.2345d-300
    !$omp target update from(x)
    if (.not. same_bits(n, x, c)) then
       nbad = nbad + 1
       write(error_unit,*) 'round trip mismatch: ', name
    endif
    if (fill) then
       x = 0.d0
       !$omp target update to(x)
    endif
    deallocate(c)
  end subroutine rt_r

  Subroutine rt_i(name, n, x, fill, nbad)
    Implicit none
    Character (len=*), intent(in) :: name
    Integer, intent(in) :: n
    Integer, intent(inout) :: x(n)
    Logical, intent(in) :: fill
    Integer, intent(inout) :: nbad
    Integer, allocatable :: c(:)
    Integer :: i
    allocate(c(n))
    if (fill) then
       do i = 1, n
          x(i) = 1 + mod(i,977)
       enddo
       !$omp target update to(x)
    endif
    c = x
    x = -123456789
    !$omp target update from(x)
    if (any(x /= c)) then
       nbad = nbad + 1
       write(error_unit,*) 'round trip mismatch: ', name
    endif
    if (fill) then
       x = 0
       !$omp target update to(x)
    endif
    deallocate(c)
  end subroutine rt_i

  Logical function same_bits_c(n, a, b)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: a(n), b(n)
    Integer :: i
    same_bits_c = .true.
    do i = 1, n
       if (transfer(real(a(i)), 0_int64) /= transfer(real(b(i)), 0_int64) .or. &
           transfer(aimag(a(i)), 0_int64) /= transfer(aimag(b(i)), 0_int64)) then
          same_bits_c = .false.; return
       endif
    enddo
  end function same_bits_c

  ! bit-for-bit comparison of n real(8) values
  Logical function same_bits(n, a, b)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(in) :: a(n), b(n)
    Integer :: i
    same_bits = .true.
    do i = 1, n
       if (transfer(a(i), 0_int64) /= transfer(b(i), 0_int64)) then
          same_bits = .false.; return
       endif
    enddo
  end function same_bits

End Module gpu_state_mod
#endif

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
!>   Rank-k updates (N_non_zero > 1) are not implemented: gpu_supported refuses such models (the run stays on the CPU).
!>
!> INTERIM 205 (removed by 206/208): G is copied H->D on entry of the GPU WRAPGRUP/DO and D->H on exit, because the
!> stabilization, the tau pass and the measurements still run on the host.  The host GR is not read between entry and exit.
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

  Implicit none

  Private
  Public :: gpu_upd_host_only, gpu_upd_enabled, gpu_upd_override, gpu_upd_resync_fields, gpu_upd_check_fields, gpu_upd_selftest_due, gpu_upd_fault_set, gpu_upd_fault_clear
  Public :: gpu_wrapgrup, gpu_wrapgrdo, gpu_trace_on, gpu_trace_begin, gpu_trace_end

  Logical, Save :: checked = .false., enabled = .false., forced = .false., forced_val = .false.
  Logical, Save :: fault_upd = .false., fault_nsig = .false., test_host_only = .false.
  Logical, Save :: sel_done = .false., sel_checked = .false.
  Logical, Save :: tr_checked = .false., tr_on = .false.
  Integer, Save :: tr_unit = -1
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: tr_f(:)
  Type (Fields), Save :: nsigma_new
  Logical, Save :: nsigma_new_made = .false.

Contains

!--------------------------------------------------------------------
!> .true. iff the device update path replaces WRAPGRUP/WRAPGRDO: GPU path active and ALF_GPU_NOUPDATE /= 1.
  Logical function gpu_upd_enabled()
    Implicit none
    Character (len=16) :: env
    Integer :: st
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
    Character (len=16) :: env
    Integer :: st
    gpu_upd_selftest_due = .false.
    if (sel_done) return
    if (.not. sel_checked) then
       sel_checked = .true.
       sel_done = .not. gpu_env_is('ALF_GPU_UPDTEST', '1')
    endif
    if (.not. sel_done) then
       sel_done = .true.
       gpu_upd_selftest_due = .true.
    endif
  end function gpu_upd_selftest_due

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
  Subroutine g_in(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(:,:,:)
    Integer :: nf_eff, nf
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       G_d(:,:,nf) = GR(:,:,nf)
       ! INTERIM 205: removed by 206/208 (G lives on the device once stabilization and measurements do).  The sweep's
       ! stabilization / tau pass / measurements still read and write the host GR.
       !$omp target update to(G_d(:,:,nf))
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
       ! INTERIM 205: removed by 206/208 (see g_in)
       !$omp target update from(G_d(:,:,nf))
       GR(:,:,nf) = G_d(:,:,nf)
       upd_cnt(8) = upd_cnt(8) + 1
       upd_cnt(9) = upd_cnt(9) + 16_int64*Ndim*Ndim
    enddo
  end subroutine g_out

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

!====================================================================
!  Host-side orchestration
!====================================================================
  Subroutine site_wrap(nf, n, t, dir)
    Implicit none
    Integer, intent(in) :: nf, n, t, dir
    if (V_meta(1,n,nf) == 1) then
       call k_site_wrap(G_d(1,1,nf), Ndim, Ndim, n, t, dir, Op_V(n,nf)%P(1), V_meta(3,n,nf), V_type(n))
       upd_cnt(3) = upd_cnt(3) + 1
    else
       call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, t, dir, 1)
       upd_cnt(3) = upd_cnt(3) + 2
    endif
  end subroutine site_wrap

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
    Integer :: nf, nf_eff, op_dim_nf, p1

    toggle = .false.
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
       else
          p1 = Op_V(n_op,nf)%P(1)
          ! The ratio needs one element of G on the host, G(P(1),P(1),nf): 16 bytes D->H per flavour and site (unavoidable
          ! while the decision is on the host, ADR-001 section 5).
          !$omp target update from(G_d(p1,p1,nf))
          upd_cnt(5) = upd_cnt(5) + 1
          gpp(nf_eff) = G_d(p1,p1,nf)
          myexp = exp( Z1* Op_V(n_op,nf)%E(1) )
          Z = myexp - 1.d0
          Delta(nf_eff) = Z
          Mat11 = - Z * gpp(nf_eff)
          Mat11 = myexp + Mat11
          D_mat = Mat11
       endif
       Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
    Enddo

    Ratiotot = Product(Ratio)
    nf = 1
    Ratiotot = (Ratiotot**dble(N_SUN)) * nsigma_new%Gama(1,1)/nsigma%Gama(n_op,nt)
    Ratiotot = Ratiotot * Prev_Ratiotot
    Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )

    upd_cnt(1) = upd_cnt(1) + 1
    if ( Weight > ranf_wrap() )  Then
       toggle = .true.
       upd_cnt(2) = upd_cnt(2) + 1
       Phase = Phase * Ratiotot/sqrt(Ratiotot*conjg(Ratiotot))
       Do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          if (Op_V(n_op,nf)%N_non_zero > 0) then
             p1 = Op_V(n_op,nf)%P(1)
             ! u = Delta, v(P) = 1 - G(P,P): x_v = u/(1 + v u), ZGERU alpha = -x_v(P)
             vii = 1.d0 - gpp(nf_eff)
             xv  = Delta(nf_eff)/(1.d0 + vii*Delta(nf_eff))
             Zg  = - xv
             if (fault_upd) Zg = Zg * (1.d0 + 1.d-6)        ! negative control (ALF_GPU_FAULT=upd)
             call k_upd_prep(G_d(1,1,nf), Ndim, Ndim, p1, scr_vec(1,1), scr_vec(1,2))
             call k_upd_geru(G_d(1,1,nf), Ndim, Ndim, real(Zg), aimag(Zg), scr_vec(1,1), scr_vec(1,2))
             upd_cnt(4) = upd_cnt(4) + 2
          endif
       enddo
       ! Flip the spin: host field, and the device copy that the wrap kernels read
       nsigma%f(n_op,nt) = nsigma_new%f(1,1)
       nsigma_d(n_op,nt) = nsigma%f(n_op,nt)
       if (.not. fault_nsig) then
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
    Logical :: Acc, toggle1

    if (N_Global_tau > 0) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    upd_cnt(10) = upd_cnt(10) + 1
    NTAU1 = NTAU + 1
    call g_in(GR)
    ! Wrap up with B(tau1)
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       call gpu_mmthr   (G_d(1,1,nf), Ndim, Ndim, nf)
       call gpu_mmthl_m1(G_d(1,1,nf), Ndim, Ndim, nf)
    Enddo
    Do n = Nt_sequential_start,Nt_sequential_end
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          call site_wrap(nf, n, ntau1, 1)           ! Op_Wrapup, N_Type 1 (reads nsigma_d(n,ntau1), the old field)
       enddo
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
       If ( T0_proposal > ranf_wrap() ) Then
          call site_update(n, ntau1, PHASE, HS_new, S0_ratio, T0_Proposal_ratio, Acc)
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
       ! Op_Wrapup with N_Type 2: nothing for diagonal operators
    Enddo
    call g_out(GR)
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
    Logical :: Acc, toggle1

    if (N_Global_tau > 0) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    upd_cnt(10) = upd_cnt(10) + 1
    call g_in(GR)
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       ! Op_Wrapdo with N_Type 2: nothing for diagonal operators
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
       If ( T0_proposal > ranf_wrap() ) Then
          call site_update(n, ntau, PHASE, HS_new, S0_ratio, T0_Proposal_ratio, Acc)
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
       ! the field has changed after the update: Op_Wrapdo N_Type 1 reads the new nsigma_d(n,ntau)
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          call site_wrap(nf, n, ntau, -1)
       enddo
    Enddo
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       call gpu_mmthl   (G_d(1,1,nf), Ndim, Ndim, nf)
       call gpu_mmthr_m1(G_d(1,1,nf), Ndim, Ndim, nf)
    enddo
    call g_out(GR)
  end Subroutine gpu_wrapgrdo

end Module gpu_update_mod
#endif

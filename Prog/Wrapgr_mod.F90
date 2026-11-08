!  Copyright (C) 2016 - 2026 The ALF project
! 
!  This file is part of the ALF project.
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
!
!  Modified 2026 on the delayed/submatrix research branch (feature/delayed-submatrix): hooks of
!  the low-rank (delayed/submatrix) sequential update in the wraps and the local-update loop.
!  See docs/algorithm-notes.md.

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles global updates on a single time slice
!
!--------------------------------------------------------------------

Module Wrapgr_mod


  Use Hamiltonian_main
  Use MyMats 
  Use Operator_mod
  Use Control
  Use Random_Wrap
  Use Fields_mod
  Use Hamiltonian_main
  Use Hop_mod
  use upgrade_mod
#ifdef ALF_GPU
  use gpu_state_mod,  only: gpu_active, gpu_state_finalize, info_line, rel_diff, selftest_fail, gpu_env_get, gpu_env_is
  use gpu_update_mod
#endif
  use Lowrank_update_mod, only: Lowrank_initialised, Lowrank_init, Lowrank_active, Lowrank_eligible, &
       &                        Lowrank_flush, Lowrank_wrap_hook, Upgrade_lowrank, Lowrank_assert_flushed, &
       &                        Lowrank_cross_slice_active, Lowrank_hop_up, Lowrank_hop_down, &
       &                        Lowrank_autotune_slice_begin, Lowrank_autotune_slice_end
#ifdef ALF_LOWRANK_DEBUG
  use Lowrank_update_mod, only: Lowrank_debug_check
#endif
  use QMC_runtime_var, only: get_Sequential_scheme, get_N_delay, get_Lowrank_cross_slice, &
       &                     get_Nt_sequential_start, get_Nt_sequential_end
  use Langevin_HMC_mod

  Implicit none

  INTERFACE wrapgr_sort
    MODULE PROCEDURE  wrapgr_sort_value, wrapgr_sort_langevin
  END INTERFACE
  
  !> Privat 
  Complex (Kind=Kind(0.d0)),  private, allocatable ::  GR_ST(:,:,:)
  
Contains
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Allocate, Deallocate space
!--------------------------------------------------------------------
  Subroutine Wrapgr_alloc
    Implicit none
    Allocate (GR_ST(Ndim,Ndim,N_FL) )
  end Subroutine Wrapgr_alloc
  
  Subroutine  Wrapgr_dealloc
    Implicit none
    deallocate ( GR_ST )
  end Subroutine Wrapgr_dealloc

!--------------------------------------------------------------------
#ifdef ALF_GPU
  RECURSIVE &
#endif
  SUBROUTINE WRAPGRUP(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
             &        Sequential_MALA, Delta_t_MALA_sequential, Max_Force_MALA_sequential, &
             &        N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given the green function matrix GR at time NTAU   at m = 0, the routine   propagates 
!> it to time  NTAU + 1 and carries  out an update of the fields at time NTAU+1
!> NTAU: [0:LTROT-1]
!
!--------------------------------------------------------------------
    Implicit none
    
    ! Arguments
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
    INTEGER, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0, Sequential_MALA
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau, N_Global_tau_MALA
    Real    (kind=Kind(0.d0)), intent(in) :: Delta_t_MALA_sequential, Delta_t_MALA_global_tau
    real    (kind=kind(0.d0)), intent(in) :: Max_Force_MALA_sequential, Max_Force_MALA_global_tau

    !Local 
    Integer :: nf, nf_eff, N_Type, NTAU1,n, m
    Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
    Complex (Kind=Kind(0.d0)) :: force_old, force_new, phase_st, nsigma_st
    Real    (Kind=Kind(0.d0)) :: T0_proposal,  T0_Proposal_ratio,  S0_ratio
    real    (kind=kind(0.d0)) :: force_0_old, force_0_new, weight
    real    (kind=kind(0.d0)) :: delta_t_running_old, Xmax, Delta_t_running_new
    Character (Len=64)        :: Mode
    Logical                   :: Acc, toggle1, L_MALA
    Logical                   :: LR_n

#ifdef ALF_GPU
    ! GPU port, ticket 205: device-resident G (INTERIM copies per call), host decision, same RNG order
    if (gpu_upd_enabled()) then
       if (gpu_upd_selftest_due()) call Wrapgr_gpu_selftest(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
       if (gpu_trace_on()) call gpu_trace_begin(NTAU+1)
       call gpu_wrapgrup(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
       if (gpu_trace_on()) call gpu_trace_end('U',NTAU+1,PHASE)
       return
    endif
    if (gpu_trace_on()) call gpu_trace_begin(NTAU+1)     ! host path of a GPU build: same trace
#endif
    ! Low-rank scheme (Sequential_scheme / N_delay / Lowrank_cross_slice): initialised on the first call
    if (.not. Lowrank_initialised()) call Lowrank_init(get_Sequential_scheme(), get_N_delay(), &
         &                                             get_Lowrank_cross_slice(), N_Global_tau, &
         &                                             get_Nt_sequential_start(), get_Nt_sequential_end())
    ! Cross-slice delay (spec D10): the pending state of the previous slice may enter here; it is
    ! carried through the hopping wrap below (X <- B_T X, Y^T <- Y^T B_T^{-1}) instead of being flushed
    if (.not. Lowrank_cross_slice_active()) call Lowrank_assert_flushed("WRAPGRUP entry")

    ! Only used for MALA fields, where they are set before use.
    phase_st  = cmplx(1.d0,0.d0,kind(0.d0))
    nsigma_st = cmplx(0.d0,0.d0,kind(0.d0))

    ! Wrap up, upgrade ntau1.  with B^{1}(tau1) 
    NTAU1 = NTAU + 1
    if (Lowrank_cross_slice_active()) call Lowrank_hop_up(GR, ntau1)
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       CALL HOP_MOD_mmthr   (GR(:,:,nf), nf, ntau1 )
       CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf, ntau1 )
    Enddo
#ifdef ALF_LOWRANK_DEBUG
    if (Lowrank_cross_slice_active()) call Lowrank_debug_check(GR, 0, ntau1, "Lowrank_hop_up")
#endif
    ! Automatic N_delay with warm-up timing (N_delay = -1, sweeps 2-3 only; otherwise a no-op): candidate of
    ! this call, clock started after the hopping wrap; stopped after the slice-end flush (Global_tau excluded),
    ! so up- and down-calls time the same part (vertices + flushes)
    call Lowrank_autotune_slice_begin()
    Do n = Nt_sequential_start,Nt_sequential_end
       L_MALA = Sequential_MALA .and. Op_V(n,1)%type == 3
       LR_n = Lowrank_active() .and. Lowrank_eligible(n) .and. .not. L_MALA
       ! Ineligible vertex: apply the pending low-rank updates before G0 is wrapped
       if (Lowrank_active() .and. .not. LR_n) call Lowrank_flush(GR)
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          HS_Field =  nsigma%f(n,ntau1)
          N_type = 1
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
       enddo
       ! Eligible vertex: carry the pending state through the type-1 wrap (old field)
       if (LR_n) call Lowrank_wrap_hook(n, ntau1, HS_Field, 1, +1, GR)
#ifdef ALF_LOWRANK_DEBUG
       if (LR_n) call Lowrank_debug_check(GR, n, ntau1, "Lowrank_wrap_hook (up, 1)")
#endif
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       if (L_MALA) then
          nsigma_st = nsigma%f(n,ntau1)
          phase_st = Phase
          Gr_st = Gr
          N_type = 2
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Op_Wrapup(Gr_st(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
          enddo
          force_old = calculate_force(n,ntau1,Gr_st)
          Call ham%Ham_Langevin_HMC_S0_single( force_0_old,n, ntau1)
          call Control_MALA_sequential(force_old, force_0_old)
          Xmax = abs(dble(force_old))
          if (abs(force_0_old) > Xmax) Xmax = abs(force_0_old)
          delta_t_running_old = Delta_t_MALA_sequential
          if (Xmax > Max_Force_MALA_sequential) delta_t_running_old = Max_Force_MALA_sequential*Delta_t_MALA_sequential/Xmax
          HS_New =  nsigma%f(n,ntau1)  -  ( force_0_old +  &
              &  real( Phase*force_old,kind(0.d0)) / Real(Phase,kind(0.d0)) ) * delta_t_running_old + &
              &  sqrt( 2.d0 * delta_t_running_old) * rang_wrap()
       else
          Hs_New =   nsigma%flip(n,ntau1)
       Endif
       S0_ratio          = ham%S0(n,ntau1, Hs_New )
       if ( Propose_S0 .and.   Op_V(n,1)%type == 1)   then
          T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
          T0_Proposal_ratio = 1.d0/S0_ratio
       Endif
       If ( T0_proposal > ranf_wrap() ) Then
          !Write(6,*) 'Hi', n, Op_V(n,nf)%type, T0_Proposal_ratio, S0_ratio
          if (L_MALA) then
             mode = "Intermediate"
          else
             mode = "Final"
          endif
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          if (LR_n) then
             Call Upgrade_lowrank(GR,n,ntau1,PHASE,HS_new, S0_ratio,T0_Proposal_ratio, Acc)
          else
             Call Upgrade2(GR,n,ntau1,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode ) 
          endif
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
       do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          N_type =  2
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
       enddo
       ! Eligible vertex: carry the pending state through the type-2 wrap (no-op for diagonal flavours)
       if (LR_n) call Lowrank_wrap_hook(n, ntau1, HS_Field, 2, +1, GR)
#ifdef ALF_LOWRANK_DEBUG
       if (LR_n) call Lowrank_debug_check(GR, n, ntau1, "Lowrank_wrap_hook (up, 2)")
#endif

       if (L_MALA) then
          phase = Phase * Prev_Ratiotot/sqrt(Prev_Ratiotot*conjg(Prev_Ratiotot))
          Call ham%Ham_Langevin_HMC_S0_single( force_0_new,n,ntau1)
          force_new = calculate_force(n,ntau1,GR)
          call Control_MALA_sequential(force_new, force_0_new)
          Xmax = abs(dble(force_new))
          if (abs(force_0_new) > Xmax) Xmax = abs(force_0_new)
          Delta_t_running_new = Delta_t_MALA_sequential
          if (Xmax > Max_Force_MALA_sequential) Delta_t_running_new = Max_Force_MALA_sequential*Delta_t_MALA_sequential/Xmax

          t0_Proposal_ratio = sqrt(delta_t_running_old/Delta_t_running_new) * exp(-0.25d0/Delta_t_running_new * (Abs(nsigma_st - hs_new + &
            & Delta_t_running_new*(force_0_new +  real( Phase*force_new,kind(0.d0)) / Real(Phase,kind(0.d0))) )**2) + 0.25d0/delta_t_running_old * ( &
            & Abs(hs_new - nsigma_st + delta_t_running_old*(force_0_old + real( phase_st*force_old,kind(0.d0)) / Real(phase_st,kind(0.d0)))  )**2 ) )

          weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase_st * Prev_Ratiotot, kind=Kind(0.d0))/real(Phase_st,kind=Kind(0.d0)) )

          if (weight > ranf_wrap()) then
             acc = .true.
          else
             acc = .false.
             nsigma%f(n,ntau1) = nsigma_st
             Gr = Gr_st
             phase = phase_st
          endif

          Call Control_upgrade(acc)
          Call Control_upgrade_eff(acc)

       endif
    Enddo
    ! G must be the full Green function for the Global_tau moves and on return; with the cross-slice
    ! delay (refused with N_Global_tau > 0) the pending state is returned and main.F90 flushes it
    ! before every read of G
    if (Lowrank_active() .and. .not. Lowrank_cross_slice_active()) call Lowrank_flush(GR)
    call Lowrank_autotune_slice_end()

    if (Nt_sequential_end >= Nt_sequential_start) then
      m         = Nt_sequential_end
    else
      m = 0
    endif
    If ( N_Global_tau > 0 ) then
       !if ( Nt_sequential_start >  Nt_sequential_end ) m = Nt_sequential_start
       Call Wrapgr_Random_update(GR,m,ntau1, PHASE, N_Global_tau )
    Endif
    If ( N_Global_tau_MALA > 0 ) then
       !if ( Nt_sequential_start >  Nt_sequential_end ) m = Nt_sequential_start
       Call Wrapgr_Langevin_update(GR,m,ntau1, PHASE, N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau )
    Endif
    if ( N_Global_tau > 0 .or. N_Global_tau_MALA > 0 ) Call Wrapgr_PlaceGR(GR,m, Size(OP_V,1), ntau1)
#ifdef ALF_GPU
    if (gpu_trace_on()) call gpu_trace_end('U',NTAU+1,PHASE)
#endif

    if (.not. Lowrank_cross_slice_active()) call Lowrank_assert_flushed("WRAPGRUP")

  END SUBROUTINE WRAPGRUP


!--------------------------------------------------------------------    
#ifdef ALF_GPU
  RECURSIVE &
#endif
  SUBROUTINE WRAPGRDO(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
             &        Sequential_MALA, Delta_t_MALA_sequential, Max_Force_MALA_sequential, &
             &        N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given the green function matrix GR at time NTAU   at m = size(OP_V,1), the routine   carries out an 
!> update of the fields at time NTAU and  propagates  it to time NTAU-1
!> NTAU: [LTROT:1]
!
!--------------------------------------------------------------------    
    Implicit None
    
    ! Given GREEN at time NTAU => GREEN at time NTAU - 1,
    ! Upgrade NTAU  [LTROT:1]
    
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    Integer, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0, Sequential_MALA
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau, N_Global_tau_MALA
    real    (kind=kind(0.d0)), intent(in) :: Delta_t_MALA_sequential, Delta_t_MALA_global_tau
    real    (kind=kind(0.d0)), intent(in) :: Max_Force_MALA_sequential, Max_Force_MALA_global_tau
    
    ! Local
    Integer :: nf, nf_eff, N_Type, n, m, m1
    Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
    Complex (Kind=Kind(0.d0)) :: force_old, force_new, phase_st, nsigma_st
    Real    (Kind=Kind(0.d0)) :: T0_proposal,  T0_Proposal_ratio,  S0_ratio
    real    (kind=kind(0.d0)) :: force_0_old, force_0_new, weight
    real    (kind=kind(0.d0)) :: delta_t_running_old, Xmax, Delta_t_running_new
    Character (Len=64)        :: Mode
    Logical                   :: Acc, toggle1, L_MALA
    Logical                   :: LR_n

#ifdef ALF_GPU
    ! GPU port, ticket 205: see WRAPGRUP
    if (gpu_upd_enabled()) then
       if (gpu_trace_on()) call gpu_trace_begin(NTAU)
       call gpu_wrapgrdo(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
       if (gpu_trace_on()) call gpu_trace_end('D',NTAU,PHASE)
       return
    endif
    if (gpu_trace_on()) call gpu_trace_begin(NTAU)
#endif

    ! Low-rank scheme (Sequential_scheme / N_delay / Lowrank_cross_slice): initialised on the first call
    if (.not. Lowrank_initialised()) call Lowrank_init(get_Sequential_scheme(), get_N_delay(), &
         &                                             get_Lowrank_cross_slice(), N_Global_tau, &
         &                                             get_Nt_sequential_start(), get_Nt_sequential_end())
    ! Cross-slice delay: the pending state of the previous WRAPGRDO (already carried through its
    ! hopping wrap) may enter here
    if (.not. Lowrank_cross_slice_active()) call Lowrank_assert_flushed("WRAPGRDO entry")

    ! G is flushed here when N_Global_tau > 0: the previous WRAPGRUP/WRAPGRDO returned a flushed G (with the
    ! cross-slice key a pending state may enter WRAPGRDO, but the key is refused with N_Global_tau > 0)
    ! Only used for MALA fields, where they are set before use.
    phase_st  = cmplx(1.d0,0.d0,kind(0.d0))
    nsigma_st = cmplx(0.d0,0.d0,kind(0.d0))

    m         = Size(OP_V,1)
    If ( N_Global_tau > 0 ) then 
       Call Wrapgr_Random_update(GR,m,ntau, PHASE, N_Global_tau )
    Endif
    ! Automatic N_delay with warm-up timing (no-op otherwise): timed from here (after Global_tau) to the slice-end flush,
    ! before the hopping wrap
    call Lowrank_autotune_slice_begin()

    If ( N_Global_tau_MALA > 0 ) then
       !Write(6,*) 'Call Ran_up ', m,ntau
       Call Wrapgr_Langevin_update(GR,m,ntau, PHASE, N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau )
    Endif

    !HERE
    If (Nt_sequential_end >=  Nt_sequential_start) then 
       m1 = Nt_sequential_end
    else
       m1 = 0
    endif
    if ( N_Global_tau > 0 .or. N_Global_tau_MALA > 0 ) Call Wrapgr_PlaceGR(GR,m, m1, ntau)

    
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       L_MALA = Sequential_MALA .and. Op_V(n,1)%type == 3
       LR_n = Lowrank_active() .and. Lowrank_eligible(n) .and. .not. L_MALA
       ! Ineligible vertex: apply the pending low-rank updates before G0 is wrapped
       if (Lowrank_active() .and. .not. LR_n) call Lowrank_flush(GR)
       if (L_MALA) then
          nsigma_st = nsigma%f(n,ntau)
          phase_st = phase
          force_old = calculate_force(n,ntau,Gr)
       endif
       N_type = 2
       nf = 1
       HS_Field = nsigma%f(n,ntau) 
       do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type,ntau)
       enddo
       ! Eligible vertex: carry the pending state through the type-2 wrap (no-op for diagonal flavours)
       if (LR_n) call Lowrank_wrap_hook(n, ntau, HS_Field, 2, -1, GR)
#ifdef ALF_LOWRANK_DEBUG
       if (LR_n) call Lowrank_debug_check(GR, n, ntau, "Lowrank_wrap_hook (down, 2)")
#endif
       !Write(6,*) 'Upgrade : ', ntau,n 
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       if (L_MALA)  then
          Gr_st = Gr
          Call ham%Ham_Langevin_HMC_S0_single( force_0_old, n,ntau)
          call Control_MALA_sequential(force_old, force_0_old)
          Xmax = abs(dble(force_old))
          if (abs(force_0_old) > Xmax) Xmax = abs(force_0_old)
          delta_t_running_old = Delta_t_MALA_sequential
          if (Xmax > Max_Force_MALA_sequential) delta_t_running_old = Max_Force_MALA_sequential*Delta_t_MALA_sequential/Xmax
          HS_New =  nsigma%f(n,ntau)  -  ( force_0_old +  &
            &  real( Phase*force_old,kind(0.d0)) / Real(Phase,kind(0.d0)) ) * delta_t_running_old + &
            &  sqrt( 2.d0 * delta_t_running_old) * rang_wrap()
       else
          HS_new            = nsigma%flip(n,ntau)
       endif
       S0_ratio          = ham%S0(n,ntau,HS_new)
       if ( Propose_S0 .and.   Op_V(n,1)%type == 1)   then
           T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
           T0_Proposal_ratio = 1.d0/S0_ratio
       Endif
       If ( T0_proposal > ranf_wrap() ) Then
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          if (LR_n) then
             Call Upgrade_lowrank(GR,n,ntau,PHASE,HS_new, S0_ratio,T0_Proposal_ratio, Acc)
          else
             if (L_MALA) then
                mode = "Intermediate"
                Call Upgrade2(GR,n,ntau,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode )
                phase = Phase * Prev_Ratiotot/sqrt(Prev_Ratiotot*conjg(Prev_Ratiotot))

                N_type = 2
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
                enddo

                Call ham%Ham_Langevin_HMC_S0_single( force_0_new,n,ntau)
                force_new = calculate_force(n,ntau,GR)
                call Control_MALA_sequential(force_new, force_0_new)
                Xmax = abs(dble(force_new))
                if (abs(force_0_new) > Xmax) Xmax = abs(force_0_new)
                Delta_t_running_new = Delta_t_MALA_sequential
                if (Xmax > Max_Force_MALA_sequential) Delta_t_running_new = Max_Force_MALA_sequential*Delta_t_MALA_sequential/Xmax

                t0_Proposal_ratio = sqrt(delta_t_running_old/Delta_t_running_new) * exp(-0.25d0/Delta_t_running_new * (Abs(nsigma_st - hs_new + &
                    & Delta_t_running_new*(force_0_new +  real( Phase*force_new,kind(0.d0)) / Real(Phase,kind(0.d0))) )**2) + 0.25d0/delta_t_running_old *(  &
                    & Abs(hs_new - nsigma_st + delta_t_running_old*(force_0_old + real( phase_st*force_old,kind(0.d0)) / Real(phase_st,kind(0.d0)))  )**2 ) )

                weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase_st * Prev_Ratiotot, kind=Kind(0.d0))/real(Phase_st,kind=Kind(0.d0)) )

                if (weight > ranf_wrap()) then
                   acc = .true.
                   N_type = 2
                   do nf_eff = 1,N_FL_eff
                      nf=Calc_Fl_map(nf_eff)
                      Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type,ntau)
                   enddo
                else
                   acc = .false.
                   nsigma%f(n,ntau) = nsigma_st
                   Gr = Gr_st
                   phase = phase_st
                endif

                Call Control_upgrade(acc)
                Call Control_upgrade_eff(acc)

             else
                mode = "Final"
                Call Upgrade2(GR,n,ntau,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode )
             endif
          endif
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
     
       !Call Upgrade(GR,n,ntau,PHASE,Op_V(n,1)%N_non_zero) 
       ! The spin has changed after the upgrade!
       nf = 1
       HS_Field = nsigma%f(n,ntau)  
       N_type = 1
       do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type, ntau )
       enddo
       ! Eligible vertex: carry the pending state through the type-1 un-wrap (new field)
       if (LR_n) call Lowrank_wrap_hook(n, ntau, HS_Field, 1, -1, GR)
#ifdef ALF_LOWRANK_DEBUG
       if (LR_n) call Lowrank_debug_check(GR, n, ntau, "Lowrank_wrap_hook (down, 1)")
#endif
    enddo
    ! G must be the full Green function before the hopping propagation and on return. Cross-slice
    ! delay: instead, carry the pending state through the hopping wrap (X <- B_T^{-1} X, Y^T <- Y^T B_T);
    ! main.F90 flushes it before every read of G
    if (Lowrank_active()) then
       if (Lowrank_cross_slice_active()) then
          call Lowrank_hop_down(GR, ntau)
       else
          call Lowrank_flush(GR)
       endif
    endif
    call Lowrank_autotune_slice_end()
    DO nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       Call Hop_mod_mmthl   (GR(:,:,nf), nf,ntau)
       Call Hop_mod_mmthr_m1(GR(:,:,nf), nf,ntau)
    enddo
#ifdef ALF_LOWRANK_DEBUG
    if (Lowrank_cross_slice_active()) call Lowrank_debug_check(GR, 0, ntau, "Lowrank_hop_down")
#endif

    if (.not. Lowrank_cross_slice_active()) call Lowrank_assert_flushed("WRAPGRDO")
#ifdef ALF_GPU
    if (gpu_trace_on()) call gpu_trace_end('D',NTAU,PHASE)
#endif
    
  end SUBROUTINE WRAPGRDO
  

#ifdef ALF_GPU
!--------------------------------------------------------------------
!> GPU port, ticket 205: one-slice kernel test (ALF_GPU_UPDTEST=1, at the first WRAPGRUP call, i.e. with the true initial
!> G of the run).  From the same G, fields, Phase and RNG state: host WRAPGRUP(ntau) then WRAPGRDO(ntau+1) vs the device
!> versions.  Three trials with different RNG states.  Required: identical fields (accept/reject pattern) after the up and
!> after the down step, G relative difference (max-norm) <= 1e-10, |Phase difference| <= 1e-12, device fields in sync.  The
!> state (G, fields, Phase, RNG) is restored afterwards; the Control acceptance counters are not (use ALF_GPU_EXIT=1).
!--------------------------------------------------------------------
  Subroutine Wrapgr_gpu_selftest(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
    Use iso_fortran_env, only: output_unit
    Implicit none
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
    INTEGER, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

    Real (Kind=Kind(0.d0)), parameter :: tolG = 1.d-10, tolP = 1.d-12
    Integer :: K, itr, i, ntau1, nf_eff, nf, nbadp, nbadf, nup, ndo, nup_g, ndo_g, nsync, nacc_up, nacc_do
    Integer, allocatable :: seed0(:), seed1(:)
    Complex (Kind=Kind(0.d0)), allocatable :: f0(:,:), cfu(:,:), cfd(:,:), gfu(:,:), gfd(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: GR0(:,:,:), cGu(:,:,:), cGd(:,:,:)
    Complex (Kind=Kind(0.d0)) :: PH0, cPu, cPd
    Real (Kind=Kind(0.d0)) :: dGu, dGd, dPu, dPd, x, mGu, mGd, mPu, mPd
    Character (len=16) :: env
    Character (len=512) :: line
    ! (ticket 212: no environment locals left)
    ! (ticket 212: the ALF_GPU_* variables are read by gpu_state_mod; these comment lines keep the line count of this
    ! #ifdef ALF_GPU block equal to the base, so that the line numbers of the code after it, and with them the
    ! __LINE__ values of the GPU=none build, stay unchanged)
    Logical :: ok

    env = ''
    if (gpu_env_get('ALF_GPU_FAULT', env)) call gpu_upd_fault_set(trim(env))
    ntau1 = NTAU + 1
    call random_seed(size=K)
    allocate(seed0(K), seed1(K))
    call random_seed(get=seed0)
    allocate(f0(size(nsigma%f,1),size(nsigma%f,2)), cfu(size(nsigma%f,1),size(nsigma%f,2)), cfd(size(nsigma%f,1),size(nsigma%f,2)))
    allocate(gfu(size(nsigma%f,1),size(nsigma%f,2)), gfd(size(nsigma%f,1),size(nsigma%f,2)))
    allocate(GR0(size(GR,1),size(GR,2),size(GR,3)), cGu(size(GR,1),size(GR,2),size(GR,3)), cGd(size(GR,1),size(GR,2),size(GR,3)))
    f0 = nsigma%f; GR0 = GR; PH0 = PHASE
    mGu = 0.d0; mGd = 0.d0; mPu = 0.d0; mPd = 0.d0
    ok = .true.; nacc_up = 0; nacc_do = 0
    write(line,'(A,I0,A,I0)') 'GPU update self-test (ticket 205): one slice, WRAPGRUP ntau = ', NTAU, ' then WRAPGRDO ntau = ', ntau1
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    do itr = 1, 3
       call random_seed(put=seed0)
       do i = 1, (itr-1)*(1000*itr+3)       ! different RNG states: 0, 2003, 6006 draws skipped
          x = ranf_wrap()
       enddo
       call random_seed(get=seed1)
       ! ---- host
       call gpu_upd_override(.true., .false.)
       nsigma%f = f0; GR = GR0; PHASE = PH0
       call gpu_g_host_wrote()
       CALL WRAPGRUP(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, .false., 0.d0, 0.d0, 0, 0.d0, 0.d0)
       cfu = nsigma%f; cGu = GR; cPu = PHASE
       CALL WRAPGRDO(GR,ntau1,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, .false., 0.d0, 0.d0, 0, 0.d0, 0.d0)
       cfd = nsigma%f; cGd = GR; cPd = PHASE
       ! ---- device
       call random_seed(put=seed1)
       nsigma%f = f0; GR = GR0; PHASE = PH0
       call gpu_g_host_wrote()
       call gpu_upd_resync_fields()
       call gpu_upd_override(.true., .not. gpu_upd_host_only())
       CALL WRAPGRUP(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, .false., 0.d0, 0.d0, 0, 0.d0, 0.d0)
       call gpu_g_host(GR)
       gfu = nsigma%f
       dGu = 0.d0
       do nf_eff = 1, N_FL_eff
          nf = Calc_Fl_map(nf_eff)
          dGu = max(dGu, rel_diff(GR(:,:,nf), cGu(:,:,nf)))      ! rel_diff is huge for non-finite entries
       enddo
       dPu = abs(PHASE - cPu)
       CALL WRAPGRDO(GR,ntau1,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, .false., 0.d0, 0.d0, 0, 0.d0, 0.d0)
       call gpu_g_host(GR)
       gfd = nsigma%f
       dGd = 0.d0
       do nf_eff = 1, N_FL_eff
          nf = Calc_Fl_map(nf_eff)
          dGd = max(dGd, rel_diff(GR(:,:,nf), cGd(:,:,nf)))
       enddo
       dPd = abs(PHASE - cPd)
       nsync = gpu_upd_check_fields()
       nbadp = count(gfu /= cfu) + count(gfd /= cfd)
       nup = count(cfu(:,ntau1) /= f0(:,ntau1)); ndo = count(cfd(:,ntau1) /= cfu(:,ntau1))
       nup_g = count(gfu(:,ntau1) /= f0(:,ntau1)); ndo_g = count(gfd(:,ntau1) /= gfu(:,ntau1))
       write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU update self-test trial ', itr, &
            & ': accepted host up/down ', nup, '/', ndo, ', device ', nup_g, '/', ndo_g, ', differing field entries ', nbadp, &
            & ', stale device fields ', nsync, ', G rel diff up ', dGu, ' down ', dGd, ', Phase diff up ', dPu, ' down ', dPd
       call info_line(line); write(output_unit,'(1X,A)') trim(line)
       if (.not. (nbadp == 0 .and. nsync == 0 .and. dGu <= tolG .and. dGd <= tolG .and. dPu <= tolP .and. dPd <= tolP)) ok = .false.
       if (nup /= nup_g .or. ndo /= ndo_g) ok = .false.
       nacc_up = nacc_up + nup; nacc_do = nacc_do + ndo
       if (.not. (dGu <= mGu)) mGu = dGu
       if (.not. (dGd <= mGd)) mGd = dGd
       if (.not. (dPu <= mPu)) mPu = dPu
       if (.not. (dPd <= mPd)) mPd = dPd
    enddo
    ! the test is vacuous (and the fault controls blind) without accepted flips in both directions
    if (nacc_up == 0 .or. nacc_do == 0) then
       call info_line('GPU update self-test: no accepted flip in the up or the down step over all trials')
       write(output_unit,'(A)') ' GPU update self-test: no accepted flip in the up or the down step over all trials'
       ok = .false.
    endif
    ! restore the state
    call gpu_upd_override(.false., .false.)
    call gpu_upd_fault_clear()
    call random_seed(put=seed0)
    nsigma%f = f0; GR = GR0; PHASE = PH0
    call gpu_g_host_wrote()
    call gpu_upd_resync_fields()
    write(line,'(A,ES10.3,A,ES10.3,A,ES8.1,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU update self-test max over trials: G rel diff up ', mGu, &
         & ' down ', mGd, ' (tol ', tolG, '), Phase diff up ', mPu, ' down ', mPd, ' (tol ', tolP
    call info_line(line); write(output_unit,'(1X,A)') trim(line)
    if (.not. ok) call selftest_fail('one-slice update test (accept pattern, G, Phase or device fields differ)')
    call info_line('GPU update self-test PASSED'); write(output_unit,'(A)') ' GPU update self-test PASSED'
    if (gpu_env_is('ALF_GPU_EXIT', '1')) then
       call gpu_state_finalize()
       write(output_unit,'(A)') ' GPU update self-test: ALF_GPU_EXIT=1, stopping after the test'
       stop
    endif
    deallocate(seed0, seed1, f0, cfu, cfd, gfu, gfd, GR0, cGu, cGd)
  end Subroutine Wrapgr_gpu_selftest
#endif

!--------------------------------------------------------------------
  Subroutine  Wrapgr_PlaceGR(GR,m,m1,ntau)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> The Green function on a given time slice reads
!> G(tau) = (  1 + B(tau) B(tau-1,0)  B(beta,tau))^(-1) with B(tau) =   U_n e^(d_n) U_n^(dag) .... U_1 e^(V_1) U_1^(dag)  e^(-dtau H_t)
!> On input you have 
!> G(tau,m)  = [ 1 + U_m e^(d_m) U_m^(dag) U_m^(dag) ... U_1 e^(V_1) U_1^(dag) e^(-dtau H_t) B(tau-1,0) 
!>                   B(Beta,tau)  U_n e^(d_n) U_n^(dag) ...U_(m+1) e^(d_(m+1)) U_(m+1)^(dag) U_(m+1) ] 
!> On output you have 
!> G(tau,  m1 )  
!> 
!> Note that m,m1 in [0,n]
!--------------------------------------------------------------------
    
    Implicit none

    !Arguments
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer, INTENT(IN) :: m,m1, ntau

    !Local 
    Integer :: n, nf, nf_eff, N_Type 
    Complex (Kind=Kind(0.d0)) :: HS_Field

    If (m == m1)  then 
       return
    elseif  ( m1 > m  ) then
       !Write(6,*) "Wrapup from ",  m + 1, "to",  m1, " on tau=",  ntau
       Do n = m+1,m1
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             HS_Field = nsigma%f(n,ntau)
             N_type = 1
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
          enddo
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             N_type =  2
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
          enddo
       Enddo
    elseif  (m1 < m ) then
       !Write(6,*) "Wrapdo from ",  m, "to",  m1 + 1 
       Do n =  m, m1+1 ,-1 
          N_type = 2
          nf = 1
          HS_Field = nsigma%f(n,ntau)
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type,ntau)
          enddo
          !Write(6,*) 'Upgrade : ', ntau,n 
          nf = 1
          HS_Field= nsigma%f(n,ntau)
          N_type = 1
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type, ntau )
          enddo
       enddo
    endif
    
  end Subroutine Wrapgr_PlaceGR



!--------------------------------------------------------------------
  Subroutine  Wrapgr_Random_update( GR,m,ntau, PHASE, N_Global_tau )
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> On input: 
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot   to tau = 1
!> 
!> The routine calls  
!> Global_move_tau(T0_Proposal_ratio, Flip_list, Flip_length,ntau,m,direction)
!> in the Hamiltonian module and then carries out the update  
!> 
!> On output
!> 
!> Flip_length==1  
!>        Green function is on  GR(tau,Flip_list(1) +1 )  if direction = u 
!>        Green function is on  GR(tau,Flip_list(1) -1 )  if direction = d
!>        This is valid if the move has or has not been accepted. 
!>
!> Flip_length > 1 
!>        Let m_min = min(Flip_list), m_max = max(Flip_list)
!>        direction = u -->  On output Green on m_max is accepted. Green is on m_min if not accepted. 
!>        direction = d -->  On output Green on m_min if accepted. Green is on m_max if not accepted.
!--------------------------------------------------------------------
        
    Implicit none

    ! Arguments 
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer,           INTENT(INOUT) :: m
    Integer,           INTENT(IN)    :: ntau, N_Global_tau 
    Complex  (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    


    ! Space for local variables
    Integer                   :: n, Flip_length, nf, nf_eff, N_Type, ng_c, Flip_count
    Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, S0_ratio
    COMPLEX (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New 
    Logical                   :: Acc
    Character (Len=64)        :: Mode
    Integer,      allocatable :: Flip_list(:)
    Complex (Kind=Kind(0.d0)), allocatable :: Flip_value(:), Flip_value_st(:)
    Real    (Kind=Kind(0.d0)) :: Zero = 10D-8

    Allocate ( Flip_list(Size(Op_V,1)), Flip_value(Size(Op_V,1)), Flip_value_st(Size(Op_V,1)) )

    Do ng_c = 1,N_Global_tau
       ! New configuration
       Call ham%Global_move_tau(T0_Proposal_ratio, S0_ratio,  Flip_list, Flip_length,Flip_value,ntau )
       !Write(6,*)  "Calling global move",  m, Flip_list(1), nsigma(Flip_list(1),ntau),Flip_value(1)
       If ( T0_Proposal_ratio  >  Zero )  Then
          ! Order the list
          Call wrapgr_sort(Flip_length,Flip_list,Flip_value)
          If ( Flip_length > 1 ) then
             Do Flip_count = 1, Flip_length-1 
                Flip_value_st(Flip_count)  = nsigma%f( Flip_list(Flip_count), ntau  )
             Enddo
          endif
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          !Write(6,*) "-----", Flip_length
          do Flip_count = 1,Flip_length
             n = Flip_list(Flip_count)
             !Write(6,*)  "PlaceGR",  m, n-1,ntau
             Call Wrapgr_PlaceGR(GR,m, n-1, ntau)
             !Write(6,*)  "Back from PlaceGR",  m, n-1,ntau
             If ( Flip_count == 1 .and. Flip_length > 1 ) GR_st = Gr
             Do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                HS_Field = nsigma%f(n,ntau)
                N_type = 1
                Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
             enddo
             nf = 1
             If (Flip_count <  Flip_length)  then 
                mode = "Intermediate"
                HS_new = Flip_value(Flip_count)
                Call Upgrade2(GR,n,ntau,PHASE, HS_new , &
                     &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode ) 
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   N_type =  2
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
                enddo
             else
                !Write(6,*)  "Call Up mode final", n,ntau
                mode = "Final"
                HS_new = Flip_value(Flip_count)
                Call Upgrade2(GR,n,ntau,PHASE,HS_new, &
                     &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode ) 
                !Write(6,*)  "Back from up mode final", n,ntau
                !Write(6,*)  "Acceptance", Acc
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   N_type =  2
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
                enddo
             endif
             m = n
          enddo
          If ( .not. Acc .and. Flip_length > 1 ) then
             Gr = Gr_st
             m = Flip_list(1) - 1
             Do Flip_count = 1, Flip_length-1 
                nsigma%f( Flip_list(Flip_count), ntau  ) = Flip_value_st(Flip_count)
             Enddo
          Endif
          !If (Acc) Call Hamiltonian_Print(Ntau)
       endif
    Enddo
    
    Deallocate ( Flip_list, Flip_value, Flip_value_st )
    

  end Subroutine Wrapgr_Random_update


!--------------------------------------------------------------------
  Subroutine  Wrapgr_Langevin_update( GR,m,ntau, PHASE, N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force )
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> On input:
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot   to tau = 1
!>
!> The routine calls
!> Global_MALA_move_tau(Flip_list, Flip_length,ntau)
!> in the Hamiltonian module and then carries out the update
!>
!> On output
!>
!> Flip_length==1
!>        Green function is on  GR(tau,Flip_list(1) +1 )  if direction = u
!>        Green function is on  GR(tau,Flip_list(1) -1 )  if direction = d
!>        This is valid if the move has or has not been accepted.
!>
!> Flip_length > 1
!>        Let m_min = min(Flip_list), m_max = max(Flip_list)
!>        direction = u -->  On output Green on m_max is accepted. Green is on m_min if not accepted.
!>        direction = d -->  On output Green on m_min if accepted. Green is on m_max if not accepted.
!--------------------------------------------------------------------
        
    Implicit none

    ! Arguments
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer,           INTENT(INOUT) :: m
    Integer,           INTENT(IN)    :: ntau, N_Global_tau_MALA
    Complex  (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    Real    (kind=Kind(0.d0)), intent(in) :: Delta_t_MALA_global_tau
    real    (kind=kind(0.d0)), intent(in) :: Max_Force
    


    ! Space for local variables
    Integer                   :: n, Flip_length, nf, nf_eff, N_Type, ng_c, Flip_count
    Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, S0_ratio, force_0
    COMPLEX (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New, Phase_st
    Logical                   :: Acc
    Character (Len=64)        :: Mode
    Integer,      allocatable :: Flip_list(:)
    Complex (Kind=Kind(0.d0)), allocatable :: Flip_value(:), Flip_value_st(:), forces_old(:), forces_new(:)
    Real    (Kind=Kind(0.d0)) :: weight
    real    (kind=kind(0.d0)), allocatable :: Forces_0_old(:), Forces_0_new(:)
    Real    (Kind=Kind(0.d0)) :: X, Xmax, delta_t_running_old, Delta_t_running_new

    Allocate ( Flip_list(Size(Op_V,1)), Flip_value(Size(Op_V,1)), Flip_value_st(Size(Op_V,1)) )
    Allocate ( Forces_old(Size(Op_V,1)), Forces_new(Size(Op_V,1)), Forces_0_old(Size(Op_V,1)), Forces_0_new(Size(Op_V,1)))

    Do ng_c = 1,N_Global_tau_MALA
       Phase_st = Phase
       ! New configuration
       Call ham%Global_MALA_move_tau(Flip_list, Flip_length,ntau )
       !Write(6,*)  "Calling global move",  m, Flip_list(1), nsigma(Flip_list(1),ntau)
       ! Order the list
       Call wrapgr_sort(Flip_length,Flip_list)
       Do Flip_count = 1, Flip_length
          Flip_value_st(Flip_count)  = nsigma%f( Flip_list(Flip_count), ntau  )
       Enddo
       !Write(6,*) "-----", Flip_length
       !Calculate forces with current nsigma
       Xmax = 0.d0
       do Flip_count = 1,Flip_length
          n = Flip_list(Flip_count)
          !Write(6,*)  "PlaceGR",  m, n-1,ntau
          Call Wrapgr_PlaceGR(GR,m, n, ntau)
          !Write(6,*)  "Back from PlaceGR",  m, n-1,ntau
          If ( Flip_count == 1 ) GR_st = Gr
          forces_old(Flip_count)   = calculate_force(n,ntau,GR)
          call ham%Ham_Langevin_HMC_S0_single (force_0,n,ntau)
          forces_0_old(Flip_count) = force_0
          X = abs(dble(forces_old(Flip_count)))
          if (X > Xmax) Xmax = X
          X = abs(force_0)
          if (X > Xmax) Xmax = X
          m = n
       enddo
       delta_t_running_old = Delta_t_MALA_global_tau
       If (Xmax > Max_Force) delta_t_running_old = Max_Force*Delta_t_MALA_global_tau/Xmax

       call Control_MALA_Global_tau(Forces_old, forces_0_old, flip_length)

       do Flip_count = 1,Flip_length
          n = Flip_list(Flip_count)
          flip_value(Flip_count) = nsigma%f(n,ntau) - ( forces_0_old(flip_count) +  &
            &  real( Phase*forces_old(flip_count),kind(0.d0)) / Real(Phase,kind(0.d0)) ) * delta_t_running_old + &
            &  sqrt( 2.d0 * delta_t_running_old) * rang_wrap()
       enddo

       !Update Greens function
       S0_ratio = 1.d0
       Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
       do Flip_count = 1,Flip_length
          n = Flip_list(Flip_count)
          !Write(6,*)  "PlaceGR",  m, n-1,ntau
          Call Wrapgr_PlaceGR(GR,m, n-1, ntau)
          !Write(6,*)  "Back from PlaceGR",  m, n-1,ntau
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             HS_Field = nsigma%f(n,ntau)
             N_type = 1
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
          enddo
          nf = 1
          mode = "Intermediate"
          HS_new = Flip_value(Flip_count)
          S0_ratio = S0_ratio * ham%S0(n,ntau,hs_new)
          Call Upgrade2(GR,n,ntau,PHASE, HS_new , &
                &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode )
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             N_type =  2
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
           enddo
           m = n
        enddo
        phase = Phase * Prev_Ratiotot/sqrt(Prev_Ratiotot*conjg(Prev_Ratiotot))

       !Calculate forces with new nsigma
       Xmax = 0.d0
       do Flip_count = 1, Flip_length
          n = Flip_list(Flip_count)
          !Write(6,*)  "PlaceGR",  m, n-1,ntau
          Call Wrapgr_PlaceGR(GR,m, n, ntau)
          !Write(6,*)  "Back from PlaceGR",  m, n-1,ntau
          forces_new(Flip_count)   = calculate_force(n,ntau,GR)
          call ham%Ham_Langevin_HMC_S0_single (force_0,n,ntau)
          forces_0_new(Flip_count) = force_0
          X = abs(dble(forces_new(Flip_count)))
          if (X > Xmax) Xmax = X
          X = abs(force_0)
          if (X > Xmax) Xmax = X
          m = n
       enddo
       Delta_t_running_new = Delta_t_MALA_global_tau
       If (Xmax > Max_Force) Delta_t_running_new = Max_Force*Delta_t_MALA_global_tau/Xmax

       call Control_MALA_Global_tau(Forces_new, forces_0_new, flip_length)

       t0_proposal_ratio = 1.d0
       do Flip_count = 1, Flip_length
          n = Flip_list(Flip_count)
          t0_Proposal_ratio = t0_proposal_ratio * sqrt(delta_t_running_old/Delta_t_running_new) * exp(-0.25d0/Delta_t_running_new *  &
              & (Abs(Flip_value_st(Flip_count) - Flip_value(Flip_count) + Delta_t_running_new*(forces_0_new(Flip_count) + &
              & real( Phase*forces_new(Flip_count),kind(0.d0))    / Real(Phase,kind(0.d0)))     )**2 ) + 0.25d0/delta_t_running_old * &
              & (Abs(Flip_value(Flip_count) - Flip_value_st(Flip_count) + delta_t_running_old*(forces_0_old(Flip_count) + &
              & real( phase_st*forces_old(Flip_count),kind(0.d0)) / Real(phase_st,kind(0.d0)))  )**2 ) )
       enddo
      
       weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase_st * Prev_Ratiotot, kind=Kind(0.d0))/real(Phase_st,kind=Kind(0.d0)) )

       if (weight > ranf_wrap()) then
          acc = .true.
       else
          acc = .false.
          Gr = Gr_st
          phase = phase_st
          m = Flip_list(1)
          Do Flip_count = 1, Flip_length
             nsigma%f( Flip_list(Flip_count), ntau  ) = Flip_value_st(Flip_count)
          Enddo
       endif

       Call Control_upgrade(acc)
       Call Control_upgrade_eff(acc)

    Enddo
    
    Deallocate ( Flip_list, Flip_value, Flip_value_st )
    

  end Subroutine Wrapgr_Langevin_update

!----------------------------------------------------------------------------
  subroutine Wrapgr_sort_value(Flip_length,Flip_list,Flip_value)

    ! Arguments
    Integer, INTENT(IN) :: Flip_length
    Integer, INTENT(INOUT), allocatable :: Flip_list(:)
    Complex   (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: Flip_value(:)
    
    ! Local
    integer :: swaps            ! number of swaps made in one pass
    integer :: nc               ! loop variable
    integer :: temp             ! temporary holder for making swap
    Complex (Kind=Kind(0.d0))      :: X

    
    if ( Flip_length == 1 ) return 
    
    !Write(6,*) 'Before sort'
    !DO nc = 1,Flip_length
    !   Write(6,*) Flip_list(nc),  Flip_value(nc)
    !Enddo

    do 
       swaps      = 0           ! Initially, we've made no swaps
       do nc = 1, (Flip_length - 1)
          if ( Flip_list(nc)  >  Flip_list(nc+1) ) then
             temp              = Flip_list(nc  ) 
             Flip_list(nc)     = Flip_list(nc+1) 
             Flip_list(nc+1)   = temp
             X                   = Flip_value(nc   ) 
             Flip_value(nc)    = Flip_value(nc+1 ) 
             Flip_value(nc+1)  = X
             swaps             = swaps + 1
          end if
       end do
       if ( swaps == 0 ) exit ! do count swaps
    end do
    
    !Write(6,*) 'After sort'
    !DO nc = 1,Flip_length
    !   Write(6,*) Flip_list(nc),  Flip_value(nc)
    !Enddo
    
  end subroutine Wrapgr_sort_value

!----------------------------------------------------------------------------
  subroutine Wrapgr_sort_langevin(Flip_length,Flip_list)

    ! Arguments
    Integer, INTENT(INOUT) :: Flip_length
    Integer, INTENT(INOUT), allocatable :: Flip_list(:)
    
    ! Local
    integer :: swaps            ! number of swaps made in one pass
    integer :: nc               ! loop variable
    integer :: temp             ! temporary holder for making swap
    integer :: length

    !Check if all fields are continuous and remove non-continuous fields from flip_list
    length = 0
    do nc = 1, Flip_length
       if (Op_V(flip_list(nc),1)%type == 3) then
          length = length + 1
          flip_list(length) = flip_list(nc)
       endif
    enddo
    flip_length = length
    
    if ( Flip_length == 1 ) return
    
    !Write(6,*) 'Before sort'
    !DO nc = 1,Flip_length
    !   Write(6,*) Flip_list(nc)
    !Enddo

    do
       swaps      = 0           ! Initially, we've made no swaps
       do nc = 1, (Flip_length - 1)
          if ( Flip_list(nc)  >  Flip_list(nc+1) ) then
             temp              = Flip_list(nc  )
             Flip_list(nc)     = Flip_list(nc+1)
             Flip_list(nc+1)   = temp
             swaps             = swaps + 1
          end if
       end do
       if ( swaps == 0 ) exit ! do count swaps
    end do
    
    !Write(6,*) 'After sort'
    !DO nc = 1,Flip_length
    !   Write(6,*) Flip_list(nc)
    !Enddo
    
  end subroutine Wrapgr_sort_langevin
  
!----------------------------------------------------------------------------
  
  Subroutine Wrapgr_Test(Gr,ntau)
    
    ! Arguments 
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer :: ntau
    
    ! Local 
    Integer :: m, m1, N_size, nf, nf_eff, nth
    Real (Kind=kind(0.d0)) :: Xmean, Xmax

    !Input is the Green function on time slice tau 
    N_size = size(OP_V,1)
    GR_ST = GR
    m = N_Size
    m1 = 0
    DO nth = 1,10
       call Wrapgr_PlaceGR(GR,m,m1,ntau)        
       Write(6,*) m, m1
       m = m1
       m1 =  nranf(N_Size)-1
    enddo
    call Wrapgr_PlaceGR(GR,m,N_Size,ntau)        
    Write(6,*) m, N_size
    
    Xmax = 0.d0
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       Call COMPARE(GR_st(:,:,nf),GR(:,:,nf),XMAX,XMEAN)
    Enddo
    Write(6,*)  'Compare Global_tau_mod_Test ', Xmax
    Deallocate ( GR_ST )
    
    
  End Subroutine Wrapgr_Test

end Module Wrapgr_mod

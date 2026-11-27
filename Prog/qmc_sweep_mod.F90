!  Copyright (C) 2026 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

Module qmc_sweep_mod

  !--------------------------------------------------------------------
  !> @author
  !> ALF-project
  !
  !> @brief
  !> Sequential sweeps of one walker: local updates slice by slice
  !> (WRAPGRUP/WRAPGRDO), stabilisation at the slices Stab_nt, and the
  !> equal-time (and, for Ltau = 1, time-displaced) measurements.
  !>
  !> @details
  !> The measurement calls (ham%Obser, Tau_m, Tau_p) are the points where a
  !> GPU build needs the Green function on the host. With the device
  !> stabilisation (gpu_stab_on) the UDV states and G stay on the device: the
  !> sweeps call gpu_stab_* in place of WRAPUR/WRAPUL, the exchange with the
  !> storage and Recompute_GR_phase, and fetch G (gpu_g_host) before each host
  !> reader of walker%GR.
  !--------------------------------------------------------------------

  Use qmc_walker_mod, only: qmc_walker_t, Recompute_GR_phase, Reset_udv
  Use Hamiltonian_main, only: ham, Ltrot, Thtrot, Symm, Projector, reconstruction_needed, N_FL_eff
  Use Wrapgr_mod, only: WRAPGRUP, WRAPGRDO
  Use wrapur_mod, only: WRAPUR
  Use wrapul_mod, only: WRAPUL
  Use Tau_m_mod, only: Tau_m
  Use Tau_p_mod, only: Tau_p
  Use Hop_mod, only: Hop_mod_Symm
#ifdef ALF_GPU
  Use gpu_stab_mod, only: gpu_stab_on, gpu_stab_wrapur, gpu_stab_wrapul, gpu_stab_exchange, gpu_stab_recompute, &
       &                  gpu_stab_reset_l, gpu_stab_reset_r, gpu_stab_reset_st, gpu_stab_export_r
  Use gpu_update_mod, only: gpu_g_host
  Use gpu_obs_mod, only: gpu_obs_active
  Use gpu_tau_mod, only: gpu_tau_on, gpu_tau_pass
#endif
  Use QMC_runtime_var, only: get_Propose_S0, get_Nt_sequential_start, get_Nt_sequential_end, &
       &                     get_N_Global_tau, get_Sequential_MALA, get_Delta_t_MALA_sequential, &
       &                     get_Max_Force_MALA_sequential, get_N_Global_tau_MALA, &
       &                     get_delta_t_MALA_global_tau, get_Max_Force_MALA_global_tau, &
       &                     get_LOBS_ST, get_LOBS_EN, get_Ltau, get_Nwrap
#ifdef ALF_TIMERS
  Use Timers_mod, only: timer_push, timer_pop, T_STAB, T_MEAS
#endif

  Implicit none
  private

  public :: Sweep_up, Sweep_down

  contains

    !--------------------------------------------------------------------
    !> @brief
    !> Sweep from time slice 0 up to Ltrot: local updates and equal-time
    !> measurements; at each stabilisation slice the right propagation is
    !> stored and G is recomputed.
    !--------------------------------------------------------------------
    Subroutine Sweep_up(walker)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer :: NST, NTAU, NTAU1, NT1, NVAR, nf_eff
      Real (Kind=Kind(0.d0)) :: Mc_step_weight

      ! Propagation from 1 to Ltrot
      ! Set the right storage to 1
      Call Reset_udv(walker%udvr, 'r')
#ifdef ALF_GPU
      If (gpu_stab_on) then
         Do nf_eff = 1, N_FL_eff
            Call gpu_stab_reset_r(nf_eff)
         ENDDO
      Endif
#endif

      NST = 1
      DO NTAU = 0, LTROT-1
         NTAU1 = NTAU + 1
         CALL WRAPGRUP(walker%GR,NTAU,walker%PHASE,get_Propose_S0(), get_Nt_sequential_start(), get_Nt_sequential_end(), &    
              &        get_N_Global_tau(), &
              &        get_Sequential_MALA(), get_Delta_t_MALA_sequential(), get_Max_Force_MALA_sequential(), &
              &        get_N_Global_tau_MALA(), get_delta_t_MALA_global_tau(), get_Max_Force_MALA_global_tau())

         If (NTAU1 == walker%Stab_nt(NST) ) then
#ifdef ALF_TIMERS
            call timer_push(T_STAB)
#endif
            NT1 = walker%Stab_nt(NST-1)
#ifdef ALF_GPU
            If (gpu_stab_on) then
               CALL gpu_stab_wrapur(NT1, NTAU1)
               Do nf_eff = 1, N_FL_eff
                  Call gpu_stab_exchange(NST, nf_eff, .true., walker%udvst(NST, nf_eff))
               ENDDO
               NVAR = 1
               IF (NTAU1 .GT. LTROT/2) NVAR = 2
               Call gpu_stab_recompute(walker, NVAR)
            else
#endif
            CALL WRAPUR(NT1, NTAU1, walker%udvr)
            Do nf_eff = 1, N_FL_eff
               ! Read from storage left propagation from LTROT to  NTAU1
               walker%udvl(nf_eff) = walker%udvst(NST, nf_eff)
               ! Write in storage right prop from 1 to NTAU1
               walker%udvst(NST, nf_eff) = walker%udvr(nf_eff)
            ENDDO
            NVAR = 1
            IF (NTAU1 .GT. LTROT/2) NVAR = 2
            Call Recompute_GR_phase(walker, NVAR)
#ifdef ALF_GPU
            Endif
#endif
#ifdef ALF_TIMERS
            call timer_pop(T_STAB)
#endif
            NST = NST + 1
         ENDIF

         IF (NTAU1.GE. get_LOBS_ST() .AND. NTAU1.LE. get_LOBS_EN() ) THEN
#ifdef ALF_TIMERS
            call timer_push(T_MEAS)
#endif
            !Call  Global_tau_mod_Test(Gr,ntau1)
            !Stop
            !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
            !write(*,*) "Phase before obser : ",phase
            Mc_step_weight = 1.d0
#ifdef ALF_GPU
            If (gpu_obs_active()) then
               ! device observables: G_d is read in place, Symm is applied on the device (the host GR is read only if G_d is stale)
               CALL ham%Obser( walker%GR, walker%PHASE, Ntau1, Mc_step_weight )
            else
            ! G is device-authoritative: D->H of the stale flavours before the host readers below
            If (gpu_stab_on) Call gpu_g_host(walker%GR)
#endif
            If (Symm) then
               Call Hop_mod_Symm(walker%GR_Tilde,walker%GR,ntau1)
               !reconstruction of NOT calculated block!!!
               If (reconstruction_needed) Call ham%GR_reconstruction( walker%GR_Tilde )
               CALL ham%Obser( walker%GR_Tilde, walker%PHASE, Ntau1, Mc_step_weight )
            else
               !reconstruction of NOT calculated block!!!
               If (reconstruction_needed) Call ham%GR_reconstruction( walker%GR )
               CALL ham%Obser( walker%GR, walker%PHASE, Ntau1, Mc_step_weight  )
            endif
#ifdef ALF_GPU
            Endif
#endif
#ifdef ALF_TIMERS
            call timer_pop(T_MEAS)
#endif
         ENDIF
      ENDDO
    end Subroutine Sweep_up

    !--------------------------------------------------------------------
    !> @brief
    !> Sweep from Ltrot down to time slice 0: local updates and equal-time
    !> measurements; at each stabilisation slice the left propagation is
    !> stored and G is recomputed (projector: time-displaced measurement
    !> around Thtrot). Ends with G on slice 0 and, for Ltau = 1, the
    !> time-displaced measurements (Tau_m, or Tau_p in the projector code).
    !--------------------------------------------------------------------
    Subroutine Sweep_down(walker)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer :: NST, NTAU, NTAU1, NT, NT1, NVAR, nf_eff
      Real (Kind=Kind(0.d0)) :: Mc_step_weight

      Call Reset_udv(walker%udvl, 'l')
#ifdef ALF_GPU
      If (gpu_stab_on) then
         Do nf_eff = 1, N_FL_eff
            Call gpu_stab_reset_l(nf_eff)
         ENDDO
      Endif
#endif

      NST = walker%NSTM-1
      DO NTAU = LTROT,1,-1
         NTAU1 = NTAU - 1
         CALL WRAPGRDO(walker%GR,NTAU, walker%PHASE,get_Propose_S0(),get_Nt_sequential_start(), get_Nt_sequential_end(), & 
              &        get_N_Global_tau(), &
              &        get_Sequential_MALA(), get_Delta_t_MALA_sequential(), get_Max_Force_MALA_sequential(), &
              &        get_N_Global_tau_MALA(), get_delta_t_MALA_global_tau(), get_Max_Force_MALA_global_tau())
         IF (NTAU1.GE. get_LOBS_ST() .AND. NTAU1.LE. get_LOBS_EN() ) THEN
#ifdef ALF_TIMERS
            call timer_push(T_MEAS)
#endif
            !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
            !write(*,*) "Phase before obser : ",phase
            Mc_step_weight = 1.d0
#ifdef ALF_GPU
            If (gpu_obs_active()) then
               ! device observables: G_d is read in place, Symm is applied on the device (the host GR is read only if G_d is stale)
               CALL ham%Obser( walker%GR, walker%PHASE, Ntau1, Mc_step_weight )
            else
            ! G is device-authoritative: D->H of the stale flavours before the host readers below
            If (gpu_stab_on) Call gpu_g_host(walker%GR)
#endif
            If (Symm) then
               Call Hop_mod_Symm(walker%GR_Tilde,walker%GR,ntau1)
               !reconstruction of NOT calculated block!!!
               If (reconstruction_needed) Call ham%GR_reconstruction( walker%GR_Tilde )
               CALL ham%Obser( walker%GR_Tilde, walker%PHASE, Ntau1, Mc_step_weight )
            else
               !reconstruction of NOT calculated block!!!
               If (reconstruction_needed) Call ham%GR_reconstruction( walker%GR )
               CALL ham%Obser( walker%GR, walker%PHASE, Ntau1,Mc_step_weight )
            endif
#ifdef ALF_GPU
            Endif
#endif
#ifdef ALF_TIMERS
            call timer_pop(T_MEAS)
#endif
         ENDIF
         IF ( walker%Stab_nt(NST) == NTAU1 .AND. NTAU1.NE.0 ) THEN
#ifdef ALF_TIMERS
            call timer_push(T_STAB)
#endif
            NT1 = walker%Stab_nt(NST+1)
#ifdef ALF_GPU
            If (gpu_stab_on) then
               CALL gpu_stab_wrapul(NT1, NTAU1)
               do nf_eff = 1,N_FL_eff
                  ! with Ltau = 1 the new udvst(NST) is also copied to the host object (read by TAU_M)
                  Call gpu_stab_exchange(NST, nf_eff, .false., walker%udvst(NST, nf_eff))
               ENDDO
               NVAR = 1
               IF (NTAU1 .GT. LTROT/2) NVAR = 2
               Call gpu_stab_recompute(walker, NVAR)
               ! projector: the host Tau_p below reads udvr and GR (udvst came over at the exchange)
               IF( get_Ltau() == 1 .and. Projector .and. walker%Stab_nt(NST)<=THTROT+1 .and. THTROT+1<walker%Stab_nt(NST+1) ) then
                  Call gpu_stab_export_r(walker%udvr)
                  Call gpu_g_host(walker%GR)
               endif
            else
#endif
            !Write(6,*) 'Wrapul : ', NT1, NTAU1
            CALL WRAPUL(NT1, NTAU1, walker%udvl)
            !Write(6,*)  'Write UL, read UR ', NTAU1, NST
            do nf_eff = 1,N_FL_eff
               ! Read from store the right prop. from 1 to LTROT/NWRAP-1
               walker%udvr(nf_eff) = walker%udvst(NST, nf_eff)
               ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
               walker%udvst(NST, nf_eff) = walker%udvl(nf_eff)
            ENDDO
            NVAR = 1
            IF (NTAU1 .GT. LTROT/2) NVAR = 2
            Call Recompute_GR_phase(walker, NVAR)
#ifdef ALF_GPU
            Endif
#endif
#ifdef ALF_TIMERS
            call timer_pop(T_STAB)
#endif
            IF( get_Ltau() == 1 .and. Projector .and. walker%Stab_nt(NST)<=THTROT+1 .and. THTROT+1<walker%Stab_nt(NST+1) ) then
#ifdef ALF_TIMERS
               call timer_push(T_MEAS)
#endif
               Call tau_p ( walker%udvl, walker%udvr, walker%udvst, walker%GR, walker%PHASE, walker%NSTM, walker%STAB_NT, NST,  get_LOBS_ST(), get_LOBS_EN() )
#ifdef ALF_TIMERS
               call timer_pop(T_MEAS)
#endif
            endif
            NST = NST -1
         ENDIF
      ENDDO

      !Calculate and compare green functions on time slice 0.
#ifdef ALF_TIMERS
      call timer_push(T_STAB)
#endif
      NT1 = walker%Stab_nt(0)
      NT  = walker%Stab_nt(1)
#ifdef ALF_GPU
      If (gpu_stab_on) then
         CALL gpu_stab_wrapul(NT, NT1)
         Do nf_eff = 1, N_FL_eff
            Call gpu_stab_reset_r(nf_eff)
         ENDDO
         Call Reset_udv(walker%udvr, 'r')
         NVAR = 1
         Call gpu_stab_recompute(walker, NVAR)
         NST =  walker%NSTM
         Do nf_eff = 1, N_FL_eff
            Call gpu_stab_reset_st(NST, nf_eff)
         ENDDO
         Call Reset_udv(walker%udvst(NST,:), 'l')
         ! projector: the host Tau_p after the sweep reads GR (udvr was reset on the host above)
         IF ( get_LTAU() == 1 .and. Projector .and. walker%Stab_nt(1) > THTROT+1 ) Call gpu_g_host(walker%GR)
      else
#endif
      CALL WRAPUL(NT, NT1, walker%udvl)

      Call Reset_udv(walker%udvr, 'r')
      NVAR = 1
      Call Recompute_GR_phase(walker, NVAR)
      NST =  walker%NSTM
      Call Reset_udv(walker%udvst(NST,:), 'l')
#ifdef ALF_GPU
      Endif
#endif
#ifdef ALF_TIMERS
      call timer_pop(T_STAB)
      call timer_push(T_MEAS)
#endif

      IF ( get_Ltau() == 1 .and. .not. Projector ) then
#ifdef ALF_GPU
         If (gpu_tau_on) then
            ! device tau pass: G_d and udvst stay on the device, observables are measured there
            Call gpu_tau_pass(walker)
         else
         ! TAU_M reads GR (and udvst, exported at the exchange) on the host
         If (gpu_stab_on) Call gpu_g_host(walker%GR)
#endif
         Call TAU_M( walker%udvst, walker%GR, walker%PHASE, walker%NSTM, get_Nwrap(), walker%STAB_NT, get_LOBS_ST(), get_LOBS_EN() )
#ifdef ALF_GPU
         Endif
#endif
      endif
      ! When Nwrap > Thtrot+1, the smallest Stab_nt(1)=Nwrap already exceeds
      ! THTROT+1, so the backward-loop condition Stab_nt(NST)<=THTROT+1 was
      ! never met and tau_p was not called above.  At this point GR has just
      ! been recomputed at time 0 = Stab_nt(0) via CGR (line ~991), and
      ! udvl/udvr are set up such that CGRP reproduces GR.  The storage
      ! udvst(1..Nstm) contains fresh left propagations from the backward
      ! sweep, so calling tau_p with NST_IN=0 is correct here.
      IF ( get_LTAU() == 1 .and. Projector .and. walker%Stab_nt(1) > THTROT+1 ) then
         Call tau_p ( walker%udvl, walker%udvr, walker%udvst, walker%GR, walker%PHASE, walker%NSTM, walker%STAB_NT, 0, get_LOBS_ST(), get_LOBS_EN() )
      endif
#ifdef ALF_TIMERS
      call timer_pop(T_MEAS)
#endif
    end Subroutine Sweep_down

end Module qmc_sweep_mod

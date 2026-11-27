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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.

!--------------------------------------------------------------------
!> @brief
!> P2.6 (p26-langevin): the Langevin force pass and the storage reset of Langevin_HMC_mod on the device path.
!> @details
!> gpu_lang_forces replaces Langevin_HMC_Forces when the device stabilisation is on (gpu_stab_on): G is wrapped upward over all
!> slices on the device (hopping, then per operator segment the Op_Wrapup N_Type 1 and 2 of gpu_wrap_mod), the force of every
!> field of type 3 is formed on the device right after its segment is wrapped (Wrapgrup_Forces: Z = sum O(I,J) (delta_IJ -
!> G(P(J),P(I))) + alpha per flavour, F = - sum_nf g Z N_SUN), and the stabilisation points use the stored left propagations
!> read only.  The forces leave the device once per pass (N_op x Ltrot); G only for the host observables.
!> gpu_lang_reset replaces Langevin_HMC_Reset_storage: H->D of the fields (changed on the host by the update), storage and
!> G on time slice 0 rebuilt on the device; with Ltau = 1 the storage and G are copied to the host for Tau_m.
!> The field update, S0, clipping and every random number stay on the host in the CPU order (Langevin_HMC_update).
!> Operators of one segment (v_seg kind 1) have pairwise disjoint P, so wrapping a later operator of the segment does not
!> change G(P_n,P_n) of an earlier one: the sequential result of the CPU loop, bit for bit in the force inputs.
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_langevin_mod

  Use iso_fortran_env, only: int64
  Use Hamiltonian_main, only: ham, Ndim, N_FL, N_FL_eff, N_SUN, Calc_Fl_map, Ltrot, Op_V, nsigma, Symm, &
       &                      reconstruction_needed
  Use Operator_mod, only: Op_phase
  Use Control, only: Control_PrecisionP
  Use Hop_mod, only: Hop_mod_Symm
  Use QMC_runtime_var, only: get_Ltau
  Use UDV_State_mod, only: UDV_State
  Use gpu_state_mod, only: G_d, g_dev_ok, g_host_ok, nsigma_d, V_meta, V_type, V_g, V_P, v_sstart, v_seg, v_hasnd, info_line
  Use gpu_wrap_mod, only: gpu_mmthr, gpu_mmthl_m1, gpu_Vop_wrap, gpu_nd_wrap_seg
  Use gpu_stab_mod, only: gpu_stab_wrapur, gpu_stab_wrapul, gpu_stab_store_l, gpu_stab_reset_l, gpu_stab_reset_r, &
       &                  gpu_stab_reset_st, gpu_stab_cgr, gpu_stab_load_l, gpu_stab_get_st
  Use gpu_update_mod, only: gpu_g_host
  Use gpu_obs_mod, only: gpu_obs_active
  Use gpu_lang_hook_mod, only: gpu_lang_forces_p, gpu_lang_reset_p, gpu_lang_finalize_p

  Implicit none
  Private
  Public :: gpu_lang_forces, gpu_lang_reset, gpu_lang_finalize, gpu_lang_register

  Logical, Save :: lg_ready = .false.
  Integer, Save :: lg_nop = 0
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: lg_O(:), lg_alpha(:,:), lg_Z(:,:), lg_F(:,:)
  Integer, Allocatable, Save :: lg_ooff(:,:)
  Integer(int64), Save :: lg_npass = 0, lg_nreset = 0

Contains

!--------------------------------------------------------------------
!> Langevin_HMC_mod's device hooks (gpu_lang_hook_mod); called once the device stabilisation is set up (gpu_stab_on)
  Subroutine gpu_lang_register()
    Implicit none
    gpu_lang_forces_p   => gpu_lang_forces
    gpu_lang_reset_p    => gpu_lang_reset
    gpu_lang_finalize_p => gpu_lang_finalize
  end subroutine gpu_lang_register

!--------------------------------------------------------------------
!> device copies of Op_V%O (N x N per operator and flavour), alpha, and the force buffers
  Subroutine lg_setup()
    Implicit none
    Integer :: n, nf, no, N1
    if (lg_ready) return
    lg_nop = size(Op_V,1)
    allocate(lg_ooff(lg_nop,N_FL), lg_alpha(lg_nop,N_FL), lg_Z(lg_nop,N_FL), lg_F(lg_nop,Ltrot))
    no = 0
    do nf = 1, N_FL
       do n = 1, lg_nop
          lg_ooff(n,nf) = no
          no = no + Op_V(n,nf)%N**2
          lg_alpha(n,nf) = Op_V(n,nf)%alpha
       enddo
    enddo
    allocate(lg_O(max(no,1)))
    lg_O = cmplx(0.d0,0.d0,kind(0.d0))
    do nf = 1, N_FL
       do n = 1, lg_nop
          N1 = Op_V(n,nf)%N
          if (Op_V(n,1)%type == 3) lg_O(lg_ooff(n,nf)+1:lg_ooff(n,nf)+N1*N1) = reshape(Op_V(n,nf)%O(1:N1,1:N1), (/N1*N1/))
       enddo
    enddo
    lg_Z = cmplx(0.d0,0.d0,kind(0.d0)); lg_F = cmplx(0.d0,0.d0,kind(0.d0))
    !$omp target enter data map(to: lg_O, lg_ooff, lg_alpha, lg_Z, lg_F)
    lg_ready = .true.
    call info_line(' GPU Langevin force pass: ACTIVE (forces on the device, one D->H of the forces per pass)')
  end subroutine lg_setup

  Subroutine gpu_lang_finalize()
    Implicit none
    character(len=160) :: line
    if (.not. lg_ready) return
    !$omp target exit data map(delete: lg_O, lg_ooff, lg_alpha, lg_Z, lg_F)
    deallocate(lg_O, lg_ooff, lg_alpha, lg_Z, lg_F)
    lg_ready = .false.
    write(line,'(A,I0,A,I0)') ' GPU Langevin: device force passes ', lg_npass, ', device storage resets ', lg_nreset
    call info_line(trim(line))
  end subroutine gpu_lang_finalize

!--------------------------------------------------------------------
!> H->D of the host G for the computed flavours whose device copy is stale (start of the run, after a host write)
  Subroutine lg_g_in(GR)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: GR(:,:,:)
    Integer :: nf_eff, nf
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       if (g_dev_ok(nf)) cycle
       G_d(:,:,nf) = GR(:,:,nf)
       !$omp target update to(G_d(:,:,nf))
       g_dev_ok(nf) = .true.
    enddo
  end subroutine lg_g_in

!--------------------------------------------------------------------
!> Z(n,nf) = sum_{I,J} O(I,J) (delta_IJ - G(P(J),P(I))) + alpha for the type-3 operators n1..n2 of flavour nf (Wrapgrup_Forces,
!> same order of the sums).  Thread = operator.
  Subroutine k_force_z(G, ld, nf, n1, n2)
    Implicit none
    Integer, intent(in) :: ld, nf, n1, n2
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,Ndim)
    Integer :: n, i, j, nsz, po, oo
    Complex (Kind=Kind(0.d0)) :: z, z1
    !$omp target teams distribute parallel do private(i, j, nsz, po, oo, z, z1)
    do n = n1, n2
       if (V_type(n) == 3) then
          nsz = V_meta(1,n,nf); po = V_meta(2,n,nf); oo = lg_ooff(n,nf)
          z = cmplx(0.d0,0.d0,kind(0.d0))
          do i = 1, nsz
             do j = 1, nsz
                z1 = cmplx(0.d0,0.d0,kind(0.d0))
                if (i == j) z1 = cmplx(1.d0,0.d0,kind(0.d0))
                z = z + lg_O(oo + i + (j-1)*nsz) * (z1 - G(V_P(po+j), V_P(po+i)))
             enddo
          enddo
          lg_Z(n,nf) = z + lg_alpha(n,nf)
       endif
    enddo
  end subroutine k_force_z

!> F(n,nt) = - sum_nf g(n,nf) Z(n,nf) N_SUN (flavours in order), 0 for the fields of other types
  Subroutine k_force_f(nt, nsun)
    Implicit none
    Integer, intent(in) :: nt, nsun
    Integer :: n, nf
    Complex (Kind=Kind(0.d0)) :: f
    !$omp target teams distribute parallel do private(nf, f)
    do n = 1, lg_nop
       f = cmplx(0.d0,0.d0,kind(0.d0))
       if (V_type(n) == 3) then
          do nf = 1, N_FL
             f = f - V_g(n,nf) * lg_Z(n,nf) * cmplx(real(nsun,kind(0.d0)), 0.d0, kind(0.d0))
          enddo
       endif
       lg_F(n,nt) = f
    enddo
  end subroutine k_force_f

!--------------------------------------------------------------------
!> Wrapgrup_Forces of one slice on G_d: hopping of every computed flavour, then the operators segment by segment with the
!> force inputs Z after each segment, then the forces of the slice
  Subroutine lg_wrap_forces(nt1)
    Implicit none
    Integer, intent(in) :: nt1
    Integer :: nf_eff, nf, k, n
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       call gpu_mmthr   (G_d(1,1,nf), Ndim, Ndim, nf, nt1)
       call gpu_mmthl_m1(G_d(1,1,nf), Ndim, Ndim, nf, nt1)
    enddo
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       if (.not. v_hasnd(nf)) then
          do n = 1, lg_nop                                   ! table operators only (no field of type 3 in this flavour)
             call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, nt1, 1, 1)
          enddo
          cycle
       endif
       do k = v_sstart(nf), v_sstart(nf+1) - 1
          if (v_seg(1,k) == 0) then
             do n = v_seg(2,k), v_seg(3,k)
                call gpu_Vop_wrap(G_d(1,1,nf), Ndim, n, nf, nt1, 1, 1)
             enddo
          else
             call gpu_nd_wrap_seg(G_d(1,1,nf), Ndim, v_seg(2,k), v_seg(3,k), nf, nt1, 1, 1)
             call gpu_nd_wrap_seg(G_d(1,1,nf), Ndim, v_seg(2,k), v_seg(3,k), nf, nt1, 1, 2)
             call k_force_z(G_d(1,1,nf), Ndim, nf, v_seg(2,k), v_seg(3,k))
          endif
       enddo
    enddo
    call k_force_f(nt1, N_SUN)
  end subroutine lg_wrap_forces

!--------------------------------------------------------------------
!> Langevin_HMC_Forces on the device.  On input G (host GR or G_d) on slice 0 and the device storage full of left propagations.
!> On output the forces in Forces(:,:), Phase updated at the stabilisation points, the equal-time observables measured
!> (Calc_Obser_eq) with weight Delta_t_running; G_d on slice Ltrot (as the CPU routine leaves GR).
  Subroutine gpu_lang_forces(Forces, Phase, GR, GR_Tilde, Stab_nt, LOBS_ST, LOBS_EN, Calc_Obser_eq, Delta_t_running)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: Forces(:,:), Phase
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:), GR_Tilde(:,:,:)
    Integer, intent(in) :: Stab_nt(0:), LOBS_ST, LOBS_EN
    Logical, intent(in) :: Calc_Obser_eq
    Real (Kind=Kind(0.d0)), intent(in) :: Delta_t_running
    Integer :: nf_eff, nf, NST, NTAU, NTAU1, NT1, NVAR
    Complex (Kind=Kind(0.d0)) :: Z, Z1, Phase_array(N_FL)

    call lg_setup()
    lg_npass = lg_npass + 1
    call lg_g_in(GR)
    do nf_eff = 1, N_FL_eff
       call gpu_stab_reset_r(nf_eff)
    enddo
    NST = 1
    do NTAU = 0, Ltrot - 1
       NTAU1 = NTAU + 1
       call lg_wrap_forces(NTAU1)
       do nf_eff = 1, N_FL_eff
          g_host_ok(Calc_Fl_map(nf_eff)) = .false.       ! G_d moved on: a host reader must fetch it again
       enddo
       if (NTAU1 == Stab_nt(NST)) then
          NT1 = Stab_nt(NST-1)
          call gpu_stab_wrapur(NT1, NTAU1)
          Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
          do nf_eff = 1, N_FL_eff
             nf = Calc_Fl_map(nf_eff)
             call gpu_stab_load_l(NST, nf_eff)               ! udvl = udvst(NST), the storage unchanged
             NVAR = 1
             if (NTAU1 > Ltrot/2) NVAR = 2
             call gpu_stab_cgr(nf_eff, nf, NVAR, GR, Z1, .true., .true.)     ! CGR + Control_PrecisionG, G stays on the device
             call Op_phase(Z1, Op_V, nsigma, nf)
             Phase_array(nf) = Z1
          enddo
          if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
          Z = product(Phase_array)
          Z = Z**N_SUN
          call Control_PrecisionP(Z, Phase)
          Phase = Z
          NST = NST + 1
       endif
       if (NTAU1 >= LOBS_ST .and. NTAU1 <= LOBS_EN .and. Calc_Obser_eq) then
          if (gpu_obs_active()) then
             call ham%Obser(GR, Phase, NTAU1, Delta_t_running)      ! device observables read G_d
          else
             call gpu_g_host(GR)
             if (Symm) then
                call Hop_mod_Symm(GR_Tilde, GR, NTAU1)
                if (reconstruction_needed) call ham%GR_reconstruction(GR_Tilde)
                call ham%Obser(GR_Tilde, Phase, NTAU1, Delta_t_running)
             else
                if (reconstruction_needed) call ham%GR_reconstruction(GR)
                call ham%Obser(GR, Phase, NTAU1, Delta_t_running)
             endif
          endif
       endif
    enddo
    ! D->H of the forces of the pass (16 N_op Ltrot B)
    !$omp target update from(lg_F)
    Forces(1:lg_nop,1:Ltrot) = lg_F
  end subroutine gpu_lang_forces

!--------------------------------------------------------------------
!> Langevin_HMC_Reset_storage on the device: storage from Ltrot down to 0 and G on slice 0 with the phase (CGR NVAR = 1, no
!> precision control, as the CPU routine).  The fields were changed on the host: H->D of nsigma first.
  Subroutine gpu_lang_reset(Phase, GR, Stab_nt, udvst)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
    Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:)
    Integer, intent(in) :: Stab_nt(0:)
    Class (UDV_State), intent(inout) :: udvst(:,:)
    Integer :: NSTM, nf_eff, nf, NST
    Complex (Kind=Kind(0.d0)) :: Z, Phase_array(N_FL)

    call lg_setup()
    lg_nreset = lg_nreset + 1
    nsigma_d = nsigma%f
    !$omp target update to(nsigma_d)
    NSTM = size(Stab_nt) - 1
    do nf_eff = 1, N_FL_eff
       call gpu_stab_reset_l(nf_eff)
       call gpu_stab_reset_st(NSTM, nf_eff)
    enddo
    do NST = NSTM-1, 1, -1
       call gpu_stab_wrapul(Stab_nt(NST+1), Stab_nt(NST))
       do nf_eff = 1, N_FL_eff
          call gpu_stab_store_l(NST, nf_eff)
       enddo
    enddo
    call gpu_stab_wrapul(Stab_nt(1), 0)
    do nf_eff = 1, N_FL_eff
       call gpu_stab_reset_r(nf_eff)
    enddo
    Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
    do nf_eff = 1, N_FL_eff
       nf = Calc_Fl_map(nf_eff)
       call gpu_stab_cgr(nf_eff, nf, 1, GR, Z, .false., .true.)
       call Op_phase(Z, Op_V, nsigma, nf)
       Phase_array(nf) = Z
    enddo
    if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
    Phase = product(Phase_array)
    Phase = Phase**N_SUN
    if (get_Ltau() == 1) then
       ! Tau_m runs on the host after the update (qmc_driver_mod): it reads walker%udvst and GR
       do NST = 1, NSTM
          do nf_eff = 1, N_FL_eff
             call gpu_stab_get_st(NST, nf_eff, udvst(NST, nf_eff))
          enddo
       enddo
       call gpu_g_host(GR)
    endif
  end subroutine gpu_lang_reset

end Module gpu_langevin_mod
#endif

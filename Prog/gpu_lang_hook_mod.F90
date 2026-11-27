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
!> P2.6 (p26-langevin): device hooks of Langevin_HMC_mod.
!> @details
!> Langevin_HMC_mod cannot use the GPU modules directly (gpu_stab_mod -> qmc_walker_mod -> tau_m_mod -> Langevin_HMC_mod would
!> be a module cycle), so it calls the device force pass and storage reset through these pointers.  gpu_langevin_mod sets them
!> (gpu_lang_register, called once the device stabilisation is set up); unset, Langevin_HMC_mod runs its host code.
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_lang_hook_mod

  Use UDV_State_mod, only: UDV_State

  Implicit none
  Private
  Public :: gpu_lang_forces_p, gpu_lang_reset_p, gpu_lang_finalize_p

  Abstract interface
     Subroutine lang_forces_i(Forces, Phase, GR, GR_Tilde, Stab_nt, LOBS_ST, LOBS_EN, Calc_Obser_eq, Delta_t_running)
       Complex (Kind=Kind(0.d0)), intent(inout) :: Forces(:,:), Phase
       Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:), GR_Tilde(:,:,:)
       Integer, intent(in) :: Stab_nt(0:), LOBS_ST, LOBS_EN
       Logical, intent(in) :: Calc_Obser_eq
       Real (Kind=Kind(0.d0)), intent(in) :: Delta_t_running
     end Subroutine lang_forces_i
     Subroutine lang_reset_i(Phase, GR, Stab_nt, udvst)
       Import :: UDV_State
       Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
       Complex (Kind=Kind(0.d0)), intent(inout) :: GR(:,:,:)
       Integer, intent(in) :: Stab_nt(0:)
       Class (UDV_State), intent(inout) :: udvst(:,:)
     end Subroutine lang_reset_i
     Subroutine lang_finalize_i()
     end Subroutine lang_finalize_i
  end interface

  Procedure (lang_forces_i),   Pointer, Save :: gpu_lang_forces_p   => null()
  Procedure (lang_reset_i),    Pointer, Save :: gpu_lang_reset_p    => null()
  Procedure (lang_finalize_i), Pointer, Save :: gpu_lang_finalize_p => null()

end Module gpu_lang_hook_mod
#endif

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
!> @author
!> ALF-project
!
!> @brief
!> Time-slice driver of the low-rank (delayed / submatrix) sequential update.
!> Written with AI assistance on master_ai; not part of the ALF baseline.
!> @details
!> WRAPGRUP and WRAPGRDO (Wrapgr_mod) call Lowrank_wrapgrup / Lowrank_wrapgrdo at their entry.
!> These return .false. if the standard update is selected, and the baseline code then runs
!> unchanged. Otherwise they carry out the whole slice: the same hopping wraps, vertex wraps,
!> proposals and random numbers as the baseline, with Upgrade_lowrank instead of Upgrade2 for
!> eligible vertices (Lowrank_update_mod), and return .true.
!>
!> Runtime parameters, namelist &VAR_LOWRANK in the file parameters (optional; without it the
!> standard update is used):
!> \verbatim
!>  Sequential_scheme = "standard" | "delayed" | "submatrix"
!>  N_delay           = 0 (automatic) or >= 1 (fixed bound on the pending rank)
!> \endverbatim
!> Not supported together with a low-rank scheme (the run stops with a message): sequential MALA
!> updates, and Global_tau moves (N_Global_tau > 0 or N_Global_tau_MALA > 0).
!--------------------------------------------------------------------
module Lowrank_wrapgr_mod

   use runtime_error_mod
   use iso_fortran_env, only: error_unit
   use Lowrank_update_mod
   implicit none
   private

   public :: Lowrank_wrapgrup, Lowrank_wrapgrdo

   Logical :: read_done = .false.
   Character (len=64) :: Sequential_scheme = "standard"
   Integer            :: N_delay = 0

contains

!--------------------------------------------------------------------
!> @brief
!> Reads &VAR_LOWRANK (once, collectively), checks the options and initialises Lowrank_update_mod.
!--------------------------------------------------------------------
   subroutine Lowrank_setup(Nt_sequential_start, Nt_sequential_end, Sequential_MALA, &
        &                   N_Global_tau, N_Global_tau_MALA)

#ifdef MPI
      Use mpi
      Use Hamiltonian_main, only: Group_Comm
#endif
      Implicit none

      Integer, Intent(In) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau, N_Global_tau_MALA
      Logical, Intent(In) :: Sequential_MALA

      Integer :: unit_para, ierr
      Character (len=64) :: file_para
      NAMELIST /VAR_LOWRANK/ Sequential_scheme, N_delay
#ifdef MPI
      Integer :: irank, irank_g, isize_g, igroup, comm, irank_c
#endif

      read_done = .true.
      file_para = "parameters"
#ifdef MPI
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#ifdef PARALLEL_PARAMS
      comm    = Group_Comm
      irank_c = irank_g
      write(file_para,'(A,I0,A)') "Temp_", igroup, "/parameters"
#else
      comm    = MPI_COMM_WORLD
      irank_c = irank
#endif
      if (irank_c == 0) then
#endif
         Open(newunit=unit_para, file=file_para, status='old', action='read', iostat=ierr)
         if (ierr /= 0) then
            Write(error_unit,*) 'Lowrank_setup: unable to open ', trim(file_para), ierr
            CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
         endif
         Read(unit_para, nml=VAR_LOWRANK, iostat=ierr)
         ! ierr < 0: no &VAR_LOWRANK in the file, the defaults (standard update) apply
         if (ierr > 0) then
            Write(error_unit,*) 'Lowrank_setup: error reading &VAR_LOWRANK in ', trim(file_para), ierr
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
         Close(unit_para)
#ifdef MPI
      endif
      CALL MPI_BCAST(Sequential_scheme, 64, MPI_CHARACTER, 0, comm, ierr)
      CALL MPI_BCAST(N_delay, 1, MPI_INTEGER, 0, comm, ierr)
#endif

      call Lowrank_init(Sequential_scheme, N_delay, Nt_sequential_start, Nt_sequential_end)
      if (Lowrank_active() .and. (Sequential_MALA .or. N_Global_tau > 0 .or. N_Global_tau_MALA > 0)) then
         Write(error_unit,*) 'Lowrank_setup: Sequential_scheme = ', trim(Sequential_scheme), &
              & ' is not supported with sequential MALA updates or Global_tau moves; use "standard"'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

   end subroutine Lowrank_setup

!--------------------------------------------------------------------
!> @brief
!> Low-rank version of WRAPGRUP: propagates G from NTAU to NTAU+1 and updates the fields of
!> slice NTAU+1. Same arguments as WRAPGRUP.
!> @return .false. if the standard update is selected (nothing was done), .true. otherwise.
!--------------------------------------------------------------------
   logical function Lowrank_wrapgrup(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, &
        &   N_Global_tau, Sequential_MALA, N_Global_tau_MALA) result(done)

      Use Hamiltonian_main
      Use Hop_mod
      Use Operator_mod
      Use upgrade_mod
      Use Control
      Use Random_wrap
      Implicit none

      COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: GR(:,:,:)
      COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
      INTEGER, INTENT(IN) :: NTAU
      LOGICAL, INTENT(IN) :: Propose_S0, Sequential_MALA
      INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau, N_Global_tau_MALA

      Integer :: nf, nf_eff, N_Type, NTAU1, n
      Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
      Real    (Kind=Kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
      Character (Len=64) :: Mode
      Logical :: Acc, toggle1, LR_n

      if (.not. read_done) call Lowrank_setup(Nt_sequential_start, Nt_sequential_end, Sequential_MALA, &
           &                                  N_Global_tau, N_Global_tau_MALA)
      done = Lowrank_active()
      if (.not. done) return

      NTAU1 = NTAU + 1
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         CALL HOP_MOD_mmthr   (GR(:,:,nf), nf, ntau1)
         CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf, ntau1)
      Enddo
      Do n = Nt_sequential_start, Nt_sequential_end
         LR_n = Lowrank_eligible(n)
         ! Ineligible vertex: apply the pending updates before G0 is wrapped
         if (.not. LR_n) call Lowrank_flush(GR)
         Do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            HS_Field = nsigma%f(n,ntau1)
            N_type = 1
            Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
         enddo
         if (LR_n) call Lowrank_wrap_hook(n, ntau1, HS_Field, 1, +1, GR)
         T0_proposal       = 1.5D0
         T0_Proposal_ratio = 1.D0
         Hs_New            = nsigma%flip(n,ntau1)
         S0_ratio          = ham%S0(n,ntau1, Hs_New)
         if ( Propose_S0 .and. Op_V(n,1)%type == 1 ) then
            T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
            T0_Proposal_ratio = 1.d0/S0_ratio
         Endif
         If ( T0_proposal > ranf_wrap() ) Then
            if (LR_n) then
               Call Upgrade_lowrank(GR,n,ntau1,PHASE,HS_new, S0_ratio,T0_Proposal_ratio, Acc)
            else
               mode = "Final"
               Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
               Call Upgrade2(GR,n,ntau1,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode)
            endif
         else
            toggle1 = .false.
            Call Control_upgrade_eff(toggle1)
         Endif
         do nf_eff = 1,N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            N_type = 2
            Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
         enddo
         if (LR_n) call Lowrank_wrap_hook(n, ntau1, HS_Field, 2, +1, GR)
      Enddo
      ! G is the full Green function on return
      call Lowrank_flush(GR)

   end function Lowrank_wrapgrup

!--------------------------------------------------------------------
!> @brief
!> Low-rank version of WRAPGRDO: updates the fields of slice NTAU and propagates G from NTAU
!> to NTAU-1. Same arguments as WRAPGRDO.
!> @return .false. if the standard update is selected (nothing was done), .true. otherwise.
!--------------------------------------------------------------------
   logical function Lowrank_wrapgrdo(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, &
        &   N_Global_tau, Sequential_MALA, N_Global_tau_MALA) result(done)

      Use Hamiltonian_main
      Use Hop_mod
      Use Operator_mod
      Use upgrade_mod
      Use Control
      Use Random_wrap
      Implicit none

      COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: GR(:,:,:)
      COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
      INTEGER, INTENT(IN) :: NTAU
      LOGICAL, INTENT(IN) :: Propose_S0, Sequential_MALA
      INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau, N_Global_tau_MALA

      Integer :: nf, nf_eff, N_Type, n
      Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
      Real    (Kind=Kind(0.d0)) :: T0_proposal, T0_Proposal_ratio, S0_ratio
      Character (Len=64) :: Mode
      Logical :: Acc, toggle1, LR_n

      if (.not. read_done) call Lowrank_setup(Nt_sequential_start, Nt_sequential_end, Sequential_MALA, &
           &                                  N_Global_tau, N_Global_tau_MALA)
      done = Lowrank_active()
      if (.not. done) return

      Do n = Nt_sequential_end, Nt_sequential_start, -1
         LR_n = Lowrank_eligible(n)
         ! Ineligible vertex: apply the pending updates before G0 is wrapped
         if (.not. LR_n) call Lowrank_flush(GR)
         N_type = 2
         HS_Field = nsigma%f(n,ntau)
         Do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type, ntau)
         enddo
         if (LR_n) call Lowrank_wrap_hook(n, ntau, HS_Field, 2, -1, GR)
         T0_proposal       = 1.5D0
         T0_Proposal_ratio = 1.D0
         Hs_New            = nsigma%flip(n,ntau)
         S0_ratio          = ham%S0(n,ntau, Hs_New)
         if ( Propose_S0 .and. Op_V(n,1)%type == 1 ) then
            T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
            T0_Proposal_ratio = 1.d0/S0_ratio
         Endif
         If ( T0_proposal > ranf_wrap() ) Then
            if (LR_n) then
               Call Upgrade_lowrank(GR,n,ntau,PHASE,HS_new, S0_ratio,T0_Proposal_ratio, Acc)
            else
               mode = "Final"
               Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
               Call Upgrade2(GR,n,ntau,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode)
            endif
         else
            toggle1 = .false.
            Call Control_upgrade_eff(toggle1)
         Endif
         N_type = 1
         HS_Field = nsigma%f(n,ntau)
         Do nf_eff = 1, N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type, ntau)
         enddo
         if (LR_n) call Lowrank_wrap_hook(n, ntau, HS_Field, 1, -1, GR)
      enddo
      ! G must be the full Green function before the hopping propagation and on return
      call Lowrank_flush(GR)
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         Call Hop_mod_mmthl   (GR(:,:,nf), nf, ntau)
         Call Hop_mod_mmthr_m1(GR(:,:,nf), nf, ntau)
      enddo

   end function Lowrank_wrapgrdo

end module Lowrank_wrapgr_mod

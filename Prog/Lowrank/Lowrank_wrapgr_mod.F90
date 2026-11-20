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
!> unchanged. Otherwise they carry out the whole slice and return .true.
!>
!> The two routines are the baseline WRAPGRUP / WRAPGRDO, statement for statement, with these
!> additions only (marked "Low-rank" or calling Lowrank_*):
!>  * per vertex: LR_n = the vertex is eligible (Lowrank_update_mod) and not a sequential MALA vertex;
!>    otherwise the pending updates are flushed first, so that the baseline code sees the full G;
!>  * after each vertex wrap of an eligible vertex: Lowrank_wrap_hook;
!>  * the "Final" Metropolis update of an eligible vertex: Upgrade_lowrank instead of Upgrade2;
!>  * Lowrank_tune_slice_begin / _end around the vertex loop and its flush (timed automatic N_delay = -1);
!>  * a flush at the end of the vertex loop, i.e. before the Global_tau moves (up) and before the
!>    hopping propagation (down), so that G is the full Green function on return.
!> Sequential MALA vertices and the Global_tau moves (N_Global_tau, N_Global_tau_MALA) therefore run
!> exactly as in the baseline. Keep these routines in step with WRAPGRUP / WRAPGRDO when those change.
!>
!> Runtime parameters, namelist &VAR_LOWRANK in the file parameters (optional; without it the
!> standard update is used):
!> \verbatim
!>  Sequential_scheme = "standard" | "delayed" | "submatrix"
!>  N_delay           = 0 (automatic, static rule), -1 (automatic, timed in sweeps 2-3) or >= 1 (fixed bound
!>                      on the pending rank)
!> \endverbatim
!--------------------------------------------------------------------
module Lowrank_wrapgr_mod

   use runtime_error_mod
   use iso_fortran_env, only: error_unit
   use Lowrank_update_mod
   use Hamiltonian_main
   use Hop_mod
   use Operator_mod
   use upgrade_mod
   use Control
   use Random_wrap
   use Langevin_HMC_mod, only: calculate_Force
   implicit none
   private

   public :: Lowrank_wrapgrup, Lowrank_wrapgrdo

   abstract interface
      subroutine Random_update_proc(GR, m, ntau, PHASE, N_Global_tau)
         complex (Kind=Kind(0.d0)), dimension(:,:,:), intent(inout), allocatable :: GR
         integer, intent(inout) :: m
         integer, intent(in)    :: ntau, N_Global_tau
         complex (Kind=Kind(0.d0)), intent(inout) :: PHASE
      end subroutine Random_update_proc
      subroutine Langevin_update_proc(GR, m, ntau, PHASE, N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force)
         complex (Kind=Kind(0.d0)), dimension(:,:,:), intent(inout), allocatable :: GR
         integer, intent(inout) :: m
         integer, intent(in)    :: ntau, N_Global_tau_MALA
         complex (Kind=Kind(0.d0)), intent(inout) :: PHASE
         real (Kind=Kind(0.d0)), intent(in) :: Delta_t_MALA_global_tau, Max_Force
      end subroutine Langevin_update_proc
      subroutine PlaceGR_proc(GR, m, m1, ntau)
         complex (Kind=Kind(0.d0)), dimension(:,:,:), intent(inout), allocatable :: GR
         integer, intent(in) :: m, m1, ntau
      end subroutine PlaceGR_proc
   end interface

   Logical :: read_done = .false.
   Character (len=64) :: Sequential_scheme = "standard"
   Integer            :: N_delay = 0
   !> Copy of G for the sequential MALA accept/reject (as GR_ST in Wrapgr_mod)
   Complex (Kind=Kind(0.d0)), allocatable :: GR_LR_ST(:,:,:)

contains

!--------------------------------------------------------------------
!> @brief
!> Reads &VAR_LOWRANK (once, collectively), checks the options and initialises Lowrank_update_mod.
!--------------------------------------------------------------------
   subroutine Lowrank_setup(Nt_sequential_start, Nt_sequential_end)

#ifdef MPI
      Use mpi
#endif
      Implicit none

      Integer, Intent(In) :: Nt_sequential_start, Nt_sequential_end

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
      if (Lowrank_active()) allocate (GR_LR_ST(Ndim, Ndim, N_FL))

   end subroutine Lowrank_setup

  logical function Lowrank_wrapgrup(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
             &        Sequential_MALA, Delta_t_MALA_sequential, Max_Force_MALA_sequential, &
             &        N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau, &
             &        Wrapgr_Random_update, Wrapgr_Langevin_update, Wrapgr_PlaceGR) result(done)
!--------------------------------------------------------------------
!> @brief
!> Low-rank version of WRAPGRUP (same arguments, plus the three Global_tau routines of Wrapgr_mod).
!> Returns .false. without doing anything if the standard update is selected. As in the baseline:
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
    ! Routines of Wrapgr_mod for the Global_tau moves, passed in by the hook
    procedure(Random_update_proc)   :: Wrapgr_Random_update
    procedure(Langevin_update_proc) :: Wrapgr_Langevin_update
    procedure(PlaceGR_proc)         :: Wrapgr_PlaceGR
    Logical                   :: LR_n

    if (.not. read_done) call Lowrank_setup(Nt_sequential_start, Nt_sequential_end)
    done = Lowrank_active()
    if (.not. done) return
    ! Only used for MALA fields, where they are set before use.
    phase_st  = cmplx(1.d0,0.d0,kind(0.d0))
    nsigma_st = cmplx(0.d0,0.d0,kind(0.d0))

    ! Wrap up, upgrade ntau1.  with B^{1}(tau1)
    NTAU1 = NTAU + 1
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       CALL HOP_MOD_mmthr   (GR(:,:,nf), nf, ntau1 )
       CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf, ntau1 )
    Enddo
    ! Low-rank: timed automatic N_delay brackets the vertex loop and the flush (no-op unless N_delay = -1)
    call Lowrank_tune_slice_begin()
    Do n = Nt_sequential_start,Nt_sequential_end
       L_MALA = Sequential_MALA .and. Op_V(n,1)%type == 3
       ! Low-rank: eligible vertex (not a MALA vertex) or flush of the pending updates
       LR_n = Lowrank_eligible(n) .and. .not. L_MALA
       if (.not. LR_n) call Lowrank_flush(GR)
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          HS_Field =  nsigma%f(n,ntau1)
          N_type = 1
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
       enddo
       if (LR_n) call Lowrank_wrap_hook(n, ntau1, HS_Field, 1, +1, GR)
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       if (L_MALA) then
          nsigma_st = nsigma%f(n,ntau1)
          phase_st = Phase
          GR_LR_ST = Gr
          N_type = 2
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Op_Wrapup(GR_LR_ST(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
          enddo
          force_old = calculate_force(n,ntau1,GR_LR_ST)
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
       if (LR_n) call Lowrank_wrap_hook(n, ntau1, HS_Field, 2, +1, GR)

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
             Gr = GR_LR_ST
             phase = phase_st
          endif

          Call Control_upgrade(acc)
          Call Control_upgrade_eff(acc)

       endif
    Enddo
    ! G is the full Green function for the Global_tau moves and on return
    call Lowrank_flush(GR)
    call Lowrank_tune_slice_end(.false.)

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

  end function Lowrank_wrapgrup

  logical function Lowrank_wrapgrdo(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
             &        Sequential_MALA, Delta_t_MALA_sequential, Max_Force_MALA_sequential, &
             &        N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau, &
             &        Wrapgr_Random_update, Wrapgr_Langevin_update, Wrapgr_PlaceGR) result(done)
!--------------------------------------------------------------------
!> @brief
!> Low-rank version of WRAPGRDO (same arguments, plus the three Global_tau routines of Wrapgr_mod).
!> Returns .false. without doing anything if the standard update is selected. As in the baseline:
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
    ! Routines of Wrapgr_mod for the Global_tau moves, passed in by the hook
    procedure(Random_update_proc)   :: Wrapgr_Random_update
    procedure(Langevin_update_proc) :: Wrapgr_Langevin_update
    procedure(PlaceGR_proc)         :: Wrapgr_PlaceGR
    Logical                   :: LR_n

    if (.not. read_done) call Lowrank_setup(Nt_sequential_start, Nt_sequential_end)
    done = Lowrank_active()
    if (.not. done) return
    ! Only used for MALA fields, where they are set before use.
    phase_st  = cmplx(1.d0,0.d0,kind(0.d0))
    nsigma_st = cmplx(0.d0,0.d0,kind(0.d0))

    m         = Size(OP_V,1)
    If ( N_Global_tau > 0 ) then 
       !Write(6,*) 'Call Ran_up ', m,ntau
       Call Wrapgr_Random_update(GR,m,ntau, PHASE, N_Global_tau )
    Endif

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

    
    ! Low-rank: timed automatic N_delay brackets the vertex loop and the flush (no-op unless N_delay = -1)
    call Lowrank_tune_slice_begin()
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       L_MALA = Sequential_MALA .and. Op_V(n,1)%type == 3
       ! Low-rank: eligible vertex (not a MALA vertex) or flush of the pending updates
       LR_n = Lowrank_eligible(n) .and. .not. L_MALA
       if (.not. LR_n) call Lowrank_flush(GR)
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
       if (LR_n) call Lowrank_wrap_hook(n, ntau, HS_Field, 2, -1, GR)
       !Write(6,*) 'Upgrade : ', ntau,n 
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       if (L_MALA)  then
          GR_LR_ST = Gr
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
                Gr = GR_LR_ST
                phase = phase_st
             endif

             Call Control_upgrade(acc)
             Call Control_upgrade_eff(acc)

          else
             mode = "Final"
             if (LR_n) then
                Call Upgrade_lowrank(GR,n,ntau,PHASE,HS_new, S0_ratio,T0_Proposal_ratio, Acc)
             else
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
       if (LR_n) call Lowrank_wrap_hook(n, ntau, HS_Field, 1, -1, GR)
    enddo
    ! G must be the full Green function before the hopping propagation and on return
    call Lowrank_flush(GR)
    ! NTAU = 1 is the last slice of a sweep: the up sweep is NTAU = 0..Ltrot-1, the down sweep Ltrot..1
    call Lowrank_tune_slice_end(NTAU == 1)
    DO nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       Call Hop_mod_mmthl   (GR(:,:,nf), nf,ntau)
       Call Hop_mod_mmthr_m1(GR(:,:,nf), nf,ntau)
    enddo
    

  end function Lowrank_wrapgrdo

end module Lowrank_wrapgr_mod

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
!> Low-rank (delayed / submatrix) variants of the sequential local update.
!> @details
!> Selected by Sequential_scheme = "standard" | "delayed" | "submatrix" (case-insensitive)
!> and N_delay (maximal number of pending accepted vertices before a flush).
!> A vertex n is eligible for the low-rank path iff, in every flavour nf,
!> Op_V(n,nf)%N_non_zero <= 1 and Op_V(n,nf)%diag. Ineligible vertices are flushed and then
!> updated with the standard Upgrade2. See docs/algorithm-notes.md, sections 3-5.
!>
!> State per computed flavour nf_eff = 1..N_FL_eff (flavour Calc_Fl_map(nf_eff)):
!>  *  delayed  : G = G0 - X Y^T, X(Ndim, N_cols_max, nf_eff), Y(Ndim, N_cols_max, nf_eff)
!>  *  submatrix: site list S(N_delay, nf_eff), slot(Ndim, nf_eff) (site -> position, 0 if absent),
!>                Dvec(N_delay, nf_eff), Ginv(N_delay, N_delay, nf_eff)
!> Columns are counted in vectors (n_cols(nf_eff)), N_cols_max = N_delay*N_P_max, so a rank-N_P
!> extension only changes the column bookkeeping.
!--------------------------------------------------------------------

module Lowrank_update_mod
   use runtime_error_mod
   use iso_fortran_env, only: output_unit, error_unit
   implicit none
   private

   Integer, parameter, public :: SCHEME_STANDARD = 0, SCHEME_DELAYED = 1, SCHEME_SUBMATRIX = 2

   public :: Lowrank_init, Lowrank_finalize, Lowrank_active, Lowrank_eligible, &
        &    Lowrank_pending, Lowrank_assert_flushed

   !> Maximal number of non-zero eigenvalues per eligible vertex and flavour (rank-1 now)
   Integer, parameter :: N_P_max = 1

   ! ---- common state ----------------------------------------------------
   Integer              :: lr_scheme  = SCHEME_STANDARD
   Integer              :: lr_N_delay = 0
   Integer              :: N_cols_max = 0       ! N_delay*N_P_max
   Integer              :: n_pend     = 0       ! pending accepted vertices since the last flush
   Logical, allocatable :: eligible(:)          ! eligible(n), n = 1..size(Op_V,1)
   Integer, allocatable :: n_cols(:)            ! used columns per computed flavour

   ! ---- delayed ---------------------------------------------------------
   Complex (Kind=Kind(0.d0)), allocatable :: X(:,:,:), Y(:,:,:)

   ! ---- submatrix -------------------------------------------------------
   Integer,                   allocatable :: S(:,:), slot(:,:)
   Complex (Kind=Kind(0.d0)), allocatable :: Dvec(:,:), Ginv(:,:,:)

contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Selects the sequential scheme, computes the vertex eligibility from Op_V and allocates
!> the per-flavour state.
!> @details
!> For "standard" nothing is printed and nothing is allocated. For a non-standard scheme
!> the line "Sequential scheme: <NAME>, N_delay = <n>" is written to output_unit and appended
!> to the file info. If vertices are ineligible, one further line gives their number and the
!> first index. Calling Lowrank_init again re-initialises (it calls Lowrank_finalize first).
!> @param [in] scheme_name
!> Character
!> \verbatim
!>  "standard", "delayed" or "submatrix"; case-insensitive. Any other name terminates.
!> \endverbatim
!> @param [in] N_delay
!> Integer
!> \verbatim
!>  Maximal number of pending accepted vertices. Must be >= 1 for a non-standard scheme;
!>  ignored for "standard".
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_init(scheme_name, N_delay)

      Use Hamiltonian_main, only: Op_V, Ndim, N_FL, N_FL_eff
      Use Files_mod, only: str_to_upper
#ifdef MPI
      Use Hamiltonian_main, only: Group_Comm
      Use mpi
#endif
      Implicit none

      Character (len=*), Intent(In) :: scheme_name
      Integer,           Intent(In) :: N_delay

      Character (len=:), allocatable :: name
      Character (len=256)            :: line, line_inel
      Character (len=64)             :: file_info
      Integer :: n, nf, M, n_inel, n_first, unit_info
      Logical :: print_info
#ifdef MPI
      Integer :: ierr, irank, irank_g, isize_g, igroup
#endif

      call Lowrank_finalize()

      name = str_to_upper(trim(adjustl(scheme_name)))
      select case (name)
      case ("STANDARD")
         return
      case ("DELAYED")
         lr_scheme = SCHEME_DELAYED
      case ("SUBMATRIX")
         lr_scheme = SCHEME_SUBMATRIX
      case default
         Write(error_unit,*) 'Lowrank_init: unknown Sequential_scheme "', trim(scheme_name), &
              &              '" (allowed: standard, delayed, submatrix)'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      end select
      if (N_delay < 1) then
         Write(error_unit,'(A,I0,A)') ' Lowrank_init: N_delay = ', N_delay, ' must be >= 1'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      lr_N_delay = N_delay
      N_cols_max = N_delay*N_P_max

      ! Eligibility: diagonal and rank <= 1 in every flavour
      M = size(Op_V, 1)
      allocate (eligible(M))
      n_inel  = 0
      n_first = 0
      do n = 1, M
         eligible(n) = .true.
         do nf = 1, N_FL
            if (Op_V(n,nf)%N_non_zero > 1 .or. .not. Op_V(n,nf)%diag) eligible(n) = .false.
         enddo
         if (.not. eligible(n)) then
            n_inel = n_inel + 1
            if (n_first == 0) n_first = n
         endif
      enddo

      ! Per computed flavour state
      allocate (n_cols(N_FL_eff))
      n_cols = 0
      n_pend = 0
      select case (lr_scheme)
      case (SCHEME_DELAYED)
         allocate (X(Ndim, N_cols_max, N_FL_eff), Y(Ndim, N_cols_max, N_FL_eff))
         X = cmplx(0.d0, 0.d0, kind(0.d0))
         Y = cmplx(0.d0, 0.d0, kind(0.d0))
      case (SCHEME_SUBMATRIX)
         allocate (S(N_delay, N_FL_eff), slot(Ndim, N_FL_eff))
         allocate (Dvec(N_delay, N_FL_eff), Ginv(N_delay, N_delay, N_FL_eff))
         S    = 0
         slot = 0
         Dvec = cmplx(0.d0, 0.d0, kind(0.d0))
         Ginv = cmplx(0.d0, 0.d0, kind(0.d0))
      end select

      ! Report (group rank 0 only)
      write(line,'(A,A,A,I0)') 'Sequential scheme: ', name, ', N_delay = ', N_delay
      if (n_inel > 0) write(line_inel,'(A,A,A,I0,A,I0,A,I0,A)') 'Sequential scheme ', name, ': ', &
           & n_inel, ' of ', M, ' vertices are not diagonal rank-1; they use the standard update (first: n=', &
           & n_first, ')'
      print_info = .true.
      file_info  = "info"
#ifdef MPI
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup     = irank/isize_g
      print_info = (irank_g == 0)
#if defined(TEMPERING)
      write(file_info,'(A,I0,A)') "Temp_", igroup, "/info"
#endif
#endif
      if (print_info) then
         write(output_unit,'(A)') trim(line)
         if (n_inel > 0) write(output_unit,'(A)') trim(line_inel)
         Open(newunit=unit_info, file=file_info, status="unknown", position="append")
         write(unit_info,'(A)') trim(line)
         if (n_inel > 0) write(unit_info,'(A)') trim(line_inel)
         close(unit_info)
      endif

   end subroutine Lowrank_init

!--------------------------------------------------------------------
!> @brief
!> Releases the state and returns to the standard scheme. Safe to call when inactive;
!> Lowrank_init may be called again afterwards.
!--------------------------------------------------------------------
   subroutine Lowrank_finalize()
      Implicit none
      if (allocated(eligible)) deallocate (eligible)
      if (allocated(n_cols))   deallocate (n_cols)
      if (allocated(X))        deallocate (X)
      if (allocated(Y))        deallocate (Y)
      if (allocated(S))        deallocate (S)
      if (allocated(slot))     deallocate (slot)
      if (allocated(Dvec))     deallocate (Dvec)
      if (allocated(Ginv))     deallocate (Ginv)
      lr_scheme  = SCHEME_STANDARD
      lr_N_delay = 0
      N_cols_max = 0
      n_pend     = 0
   end subroutine Lowrank_finalize

!--------------------------------------------------------------------
!> @brief
!> True iff a non-standard (delayed or submatrix) scheme is active.
!--------------------------------------------------------------------
   logical function Lowrank_active()
      Implicit none
      Lowrank_active = (lr_scheme /= SCHEME_STANDARD)
   end function Lowrank_active

!--------------------------------------------------------------------
!> @brief
!> True iff vertex n_op is diagonal rank-1 in every flavour and a low-rank scheme is active.
!> @param [in] n_op
!> Integer
!> \verbatim
!>  Operator (vertex) number, 1..size(Op_V,1)
!> \endverbatim
!--------------------------------------------------------------------
   logical function Lowrank_eligible(n_op)
      Implicit none
      Integer, Intent(In) :: n_op
      if (lr_scheme == SCHEME_STANDARD) then
         Lowrank_eligible = .false.
      else
         Lowrank_eligible = eligible(n_op)
      endif
   end function Lowrank_eligible

!--------------------------------------------------------------------
!> @brief
!> Number of accepted vertices pending since the last flush (0 for the standard scheme).
!--------------------------------------------------------------------
   integer function Lowrank_pending()
      Implicit none
      Lowrank_pending = n_pend
   end function Lowrank_pending

!--------------------------------------------------------------------
!> @brief
!> Terminates if accepted vertices are still pending, i.e. if G is not the full Green function.
!> @param [in] where
!> Character
!> \verbatim
!>  Name of the calling site, written to error_unit.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_assert_flushed(where)
      Implicit none
      Character (len=*), Intent(In) :: where
      if (n_pend > 0) then
         Write(error_unit,'(A,A,A,I0,A)') ' Lowrank_assert_flushed: ', trim(where), ': ', n_pend, &
              &                            ' accepted vertices pending (G not flushed)'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
   end subroutine Lowrank_assert_flushed

end module Lowrank_update_mod

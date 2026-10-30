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

   public :: Lowrank_init, Lowrank_finalize, Lowrank_initialised, Lowrank_active, Lowrank_eligible, &
        &    Lowrank_pending, Lowrank_assert_flushed
   public :: Lowrank_propose, Lowrank_accept, Upgrade_lowrank, Lowrank_wrap_hook, &
        &    Lowrank_flush, Lowrank_current_G

   !> Maximal number of non-zero eigenvalues per eligible vertex and flavour (rank-1 now)
   Integer, parameter :: N_P_max = 1

   ! ---- common state ----------------------------------------------------
   Logical              :: lr_initialised = .false.  ! true after Lowrank_init (any scheme)
   Integer              :: lr_scheme  = SCHEME_STANDARD
   Integer              :: lr_N_delay = 0
   Integer              :: N_cols_max = 0       ! N_delay*N_P_max
   Integer              :: n_pend     = 0       ! pending accepted vertices since the last flush
   Logical, allocatable :: eligible(:)          ! eligible(n), n = 1..size(Op_V,1)
   Integer, allocatable :: n_cols(:)            ! used columns per computed flavour

   ! ---- cache of the last Lowrank_propose (consumed by Lowrank_accept) ---
   Logical                                :: prop_valid = .false.
   Integer                                :: prop_n_op = 0, prop_nt = 0
   Complex (Kind=Kind(0.d0))              :: prop_Hs
   Complex (Kind=Kind(0.d0)), allocatable :: prop_Delta(:), prop_d(:)  ! per computed flavour

   Complex (Kind=Kind(0.d0)), parameter :: cone  = cmplx( 1.d0, 0.d0, kind(0.d0))
   Complex (Kind=Kind(0.d0)), parameter :: cmone = cmplx(-1.d0, 0.d0, kind(0.d0))

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
      lr_initialised = .true.

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
      allocate (prop_Delta(N_FL_eff), prop_d(N_FL_eff))
      prop_Delta = cmplx(0.d0, 0.d0, kind(0.d0))
      prop_d     = cmplx(1.d0, 0.d0, kind(0.d0))
      prop_valid = .false.
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
      if (allocated(prop_Delta)) deallocate (prop_Delta)
      if (allocated(prop_d))     deallocate (prop_d)
      prop_valid = .false.
      lr_initialised = .false.
      lr_scheme  = SCHEME_STANDARD
      lr_N_delay = 0
      N_cols_max = 0
      n_pend     = 0
   end subroutine Lowrank_finalize

!--------------------------------------------------------------------
!> @brief
!> True after Lowrank_init (for every scheme, "standard" included) until Lowrank_finalize.
!> Used by WRAPGRUP/WRAPGRDO to initialise lazily on their first call.
!--------------------------------------------------------------------
   logical function Lowrank_initialised()
      Implicit none
      Lowrank_initialised = lr_initialised
   end function Lowrank_initialised

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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Fermionic ratio of the proposal nsigma%f(n_op,nt) -> Hs_new for an eligible vertex,
!> without random numbers and without changing G.
!> @details
!> The scalar expressions and their order are those of Upgrade2 (rank-1 branch); the only
!> difference is that G_pp is the element of the current G = G0 - X Y^T,
!> G_pp = G0(p,p) - sum_k X(p,k) Y(p,k), p = Op_V(n_op,nf)%P(1). A flavour with
!> N_non_zero = 0 contributes D_mat = 1. Delta_f = e^{Z1 E_1} - 1 and d_f = D_mat are cached
!> per computed flavour for Lowrank_accept. See docs/algorithm-notes.md, sections 1 and 3.
!> @param [in] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0, the flushed part of the Green function. Not modified (intent(inout) only to match
!>  the calling convention of Upgrade2).
!> \endverbatim
!> @param [in] n_op, nt
!> Integer
!> \verbatim
!>  Vertex and time slice. nsigma%f(n_op,nt) holds the old field.
!> \endverbatim
!> @param [in] Hs_new
!> Complex
!> \verbatim
!>  Proposed field.
!> \endverbatim
!> @param [out] Ratiotot
!> Complex
!> \verbatim
!>  (Prod_nf R_nf)^N_SUN * Gama(Hs_new)/Gama(old), as Upgrade2 before the mode branch.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_propose(GR, n_op, nt, Hs_new, Ratiotot)

      Use Hamiltonian_main, only: Op_V, nsigma, ham, Ndim, N_FL, N_FL_eff, Calc_Fl_map, N_SUN, &
           &                      reconstruction_needed
      Use Fields_mod
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In)    :: n_op, nt
      Complex (Kind=Kind(0.d0)), Intent(In)    :: Hs_new
      Complex (Kind=Kind(0.d0)), Intent(Out)   :: Ratiotot

      Type (Fields)             :: nsigma_new
      Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Z1, Z, D_mat, myexp, g_loc, G_pp, Mat11
      Integer                   :: nf, nf_eff, p, j

      if (lr_scheme /= SCHEME_DELAYED) then
         Write(error_unit,*) 'Lowrank_propose: only the delayed scheme is implemented'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      if (.not. Lowrank_eligible(n_op)) then
         Write(error_unit,'(A,I0)') ' Lowrank_propose: vertex not eligible, n_op = ', n_op
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      Call nsigma_new%make(1,1)

      ! Same expressions and order as Upgrade2
      nf = 1
      nsigma_new%f(1,1)  = Hs_new
      nsigma_new%t(1)    = Op_V(n_op,nf)%Type
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         g_loc = Op_V(n_op,nf)%g
         if (Op_V(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
         Z1 = g_loc * ( nsigma_new%Phi(1,1) -  nsigma%Phi(n_op,nt) )
         If (Op_V(n_op,nf)%N_non_zero == 0) then
            D_mat = 1.0d0
            prop_Delta(nf_eff) = cmplx(0.d0, 0.d0, kind(0.d0))
         else
            ! Current element G_pp = G0(p,p) - X(p,1:j) . Y(p,1:j)  (no complex BLAS functions)
            p = Op_V(n_op,nf)%P(1)
            j = n_cols(nf_eff)
            G_pp = GR(p,p,nf)
            if (j > 0) G_pp = G_pp - sum( X(p,1:j,nf_eff) * Y(p,1:j,nf_eff) )
            myexp = exp( Z1* Op_V(n_op,nf)%E(1) )
            Z = myexp - 1.d0
            Mat11 = - Z * G_pp
            Mat11 = myexp + Mat11
            D_mat = Mat11
            prop_Delta(nf_eff) = Z
         endif
         prop_d(nf_eff) = D_mat
         Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
      Enddo

      ! Fill the non-calculated flavours, as Upgrade2
      if (reconstruction_needed) call ham%weight_reconstruction(Ratio)

      Ratiotot = Product(Ratio)
      nf = 1
      Ratiotot = (Ratiotot**dble(N_SUN)) * nsigma_new%Gama(1,1)/nsigma%Gama(n_op,nt)

      Call nsigma_new%clear()

      prop_valid = .true.
      prop_n_op  = n_op
      prop_nt    = nt
      prop_Hs    = Hs_new

   end subroutine Lowrank_propose

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Accepts the proposal evaluated by the preceding Lowrank_propose(GR, n_op, nt, Hs_new, .).
!> @details
!> Per computed flavour with N_non_zero = 1 (p = P(1), j = n_cols):
!>   c = G0(:,p) - X(:,1:j) Y(p,1:j)^T,   r = G0(p,:)^T - Y(:,1:j) X(p,1:j)^T   (ZGEMV)
!>   X(:,j+1) = (Delta/d) c,              Y(:,j+1) = e_p - r
!> i.e. Sherman-Morrison on the current G (docs/algorithm-notes.md, sections 2-3).
!> Then nsigma%f(n_op,nt) = Hs_new; after N_delay pending vertices G is flushed.
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0. Modified only by the flush after N_delay pending vertices.
!> \endverbatim
!> @param [in] n_op, nt, Hs_new
!> \verbatim
!>  Must be those of the preceding Lowrank_propose.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_accept(GR, n_op, nt, Hs_new)

      Use Hamiltonian_main, only: Op_V, nsigma, Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In)    :: n_op, nt
      Complex (Kind=Kind(0.d0)), Intent(In)    :: Hs_new

      Complex (Kind=Kind(0.d0)) :: z
      Integer                   :: nf, nf_eff, p, j

      if (.not. Lowrank_eligible(n_op)) then
         Write(error_unit,'(A,I0)') ' Lowrank_accept: vertex not eligible, n_op = ', n_op
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      if (.not. prop_valid .or. prop_n_op /= n_op .or. prop_nt /= nt .or. prop_Hs /= Hs_new) then
         Write(error_unit,'(A,I0,A,I0)') ' Lowrank_accept: no matching Lowrank_propose for n_op = ', &
              &                          n_op, ', nt = ', nt
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         if (Op_V(n_op,nf)%N_non_zero == 0) cycle
         p = Op_V(n_op,nf)%P(1)
         j = n_cols(nf_eff)
         if (j >= N_cols_max) then
            Write(error_unit,'(A,I0)') ' Lowrank_accept: column overflow, N_cols_max = ', N_cols_max
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
         ! X(:,j+1) <- c = G0(:,p) - X(:,1:j) Y(p,1:j)^T
         call ZCOPY(Ndim, GR(1,p,nf), 1, X(1,j+1,nf_eff), 1)
         ! Y(:,j+1) <- r = G0(p,:)^T - Y(:,1:j) X(p,1:j)^T   (rows: increment Ndim)
         call ZCOPY(Ndim, GR(p,1,nf), Ndim, Y(1,j+1,nf_eff), 1)
         if (j > 0) then
            call ZGEMV('N', Ndim, j, cmone, X(1,1,nf_eff), Ndim, Y(p,1,nf_eff), Ndim, cone, X(1,j+1,nf_eff), 1)
            call ZGEMV('N', Ndim, j, cmone, Y(1,1,nf_eff), Ndim, X(p,1,nf_eff), Ndim, cone, Y(1,j+1,nf_eff), 1)
         endif
         ! X(:,j+1) = (Delta/d) c,  Y(:,j+1) = e_p - r
         z = prop_Delta(nf_eff)/prop_d(nf_eff)
         call ZSCAL(Ndim, z, X(1,j+1,nf_eff), 1)
         Y(:,j+1,nf_eff) = - Y(:,j+1,nf_eff)
         Y(p,j+1,nf_eff) = Y(p,j+1,nf_eff) + 1.d0
         n_cols(nf_eff)  = j + 1
      Enddo
      prop_valid = .false.

      ! Flip the field
      nsigma%f(n_op,nt) = Hs_new

      n_pend = n_pend + 1
      if (n_pend == lr_N_delay) call Lowrank_flush(GR)

   end subroutine Lowrank_accept

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Low-rank counterpart of Upgrade2 in mode "Final" for an eligible vertex.
!> @details
!> Lowrank_propose, then Weight, the single ranf_wrap() draw and the phase update exactly as
!> Upgrade2 (Final, Prev_Ratiotot = 1), Lowrank_accept if accepted, the trace record under
!> ALF_TRACE_UPDATES, and Control_upgrade / Control_upgrade_eff.
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0 (flushed part of G).
!> \endverbatim
!> @param [in] n_op, nt
!> Integer
!> @param [inout] Phase
!> Complex
!> \verbatim
!>  Updated if the move is accepted.
!> \endverbatim
!> @param [in] Hs_new
!> Complex
!> @param [in] S0_ratio, T0_proposal_ratio
!> Real
!> @param [inout] toggle
!> Logical
!> \verbatim
!>  True if the move is accepted.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Upgrade_lowrank(GR, n_op, nt, Phase, Hs_new, S0_ratio, T0_proposal_ratio, toggle)

      Use Hamiltonian_main, only: Ndim, N_FL
      Use Random_wrap
      Use Control
#ifdef ALF_TRACE_UPDATES
      Use Trace_updates_mod
#endif
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In)    :: n_op, nt
      Complex (Kind=Kind(0.d0)), Intent(InOut) :: Phase
      Complex (Kind=Kind(0.d0)), Intent(In)    :: Hs_new
      Real    (Kind=Kind(0.d0)), Intent(In)    :: S0_ratio, T0_proposal_ratio
      Logical,                   Intent(InOut) :: toggle

      Complex (Kind=Kind(0.d0)) :: Ratiotot, Prev_Ratiotot
      Real    (Kind=Kind(0.d0)) :: Weight, ranf_val

      call Lowrank_propose(GR, n_op, nt, Hs_new, Ratiotot)

      ! Upgrade2, mode "Final"
      Prev_Ratiotot = cmplx(1.d0, 0.d0, kind(0.d0))
      Ratiotot = Ratiotot * Prev_Ratiotot
      Weight = S0_ratio * T0_proposal_ratio * abs(  real(Phase * Ratiotot, kind=Kind(0.d0))/real(Phase,kind=Kind(0.d0)) )

      toggle = .false.
      ranf_val = ranf_wrap()           ! exactly one draw, as Upgrade2
      if ( Weight > ranf_val )  Then
         toggle = .true.
         Phase = Phase * Ratiotot/sqrt(Ratiotot*conjg(Ratiotot))
         call Lowrank_accept(GR, n_op, nt, Hs_new)
      endif
      prop_valid = .false.

#ifdef ALF_TRACE_UPDATES
      Call Trace_update_write(nt, n_op, Ratiotot, Weight, ranf_val, toggle)
#endif
      Call Control_upgrade(toggle)
      Call Control_upgrade_eff(toggle)

   end subroutine Upgrade_lowrank

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Carries the pending factors X, Y through the type-1 wrap of an eligible (diagonal) vertex.
!> @details
!> Call once per vertex, after the flavour loop of the type-1 wrap (the routine loops over the
!> computed flavours itself): after Op_Wrapup(.,Op_V(n_op,nf),HS_field,.,1,nt) with dir = +1, and
!> after Op_Wrapdo(.,Op_V(n_op,nf),HS_field,.,1,nt) with dir = -1. Invalidates the cache of
!> Lowrank_propose (the frame changes). The wrap is
!> G <- T G T^{-1} with T = diag(t_i) on rows P(i), t_i = exp(dir g_loc phi(HS_field) E(i));
!> applied to G0 by ALF, to the pending factors here: X(P(i),:) *= t_i, Y(P(i),:) /= t_i
!> (docs/algorithm-notes.md, section 3). Consumes no random numbers. No-op if nothing is pending.
!> @param [in] n_op, nt
!> Integer
!> @param [in] HS_field
!> Complex
!> \verbatim
!>  The field used in the wrap.
!> \endverbatim
!> @param [in] dir
!> Integer
!> \verbatim
!>  +1 after Op_Wrapup type 1, -1 after Op_Wrapdo type 1.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_wrap_hook(n_op, nt, HS_field, dir)

      Use Hamiltonian_main, only: Op_V, N_FL_eff, Calc_Fl_map
      Use Fields_mod
      Implicit none

      Integer,                   Intent(In) :: n_op, nt, dir
      Complex (Kind=Kind(0.d0)), Intent(In) :: HS_field

      Type (Fields)             :: nsigma_single
      Complex (Kind=Kind(0.d0)) :: phi, g_loc, t
      Integer                   :: nf, nf_eff, i, j, q

      prop_valid = .false.
      if (n_pend == 0) return
      if (lr_scheme /= SCHEME_DELAYED) then
         Write(error_unit,*) 'Lowrank_wrap_hook: only the delayed scheme is implemented'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      if (.not. eligible(n_op)) then
         Write(error_unit,'(A,I0)') ' Lowrank_wrap_hook: vertex not eligible, n_op = ', n_op
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      Call nsigma_single%make(1,1)
      Do nf_eff = 1,N_FL_eff
         j = n_cols(nf_eff)
         if (j == 0) cycle
         nf = Calc_Fl_map(nf_eff)
         nsigma_single%f(1,1) = HS_field
         nsigma_single%t(1)   = Op_V(n_op,nf)%type
         phi   = nsigma_single%Phi(1,1)
         g_loc = Op_V(n_op,nf)%g
         if (Op_V(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
         do i = 1,Op_V(n_op,nf)%N
            q = Op_V(n_op,nf)%P(i)
            t = exp( real(dir, kind(0.d0)) * g_loc * phi * Op_V(n_op,nf)%E(i) )
            X(q,1:j,nf_eff) = X(q,1:j,nf_eff) * t
            Y(q,1:j,nf_eff) = Y(q,1:j,nf_eff) / t
         enddo
      Enddo
      Call nsigma_single%clear()

   end subroutine Lowrank_wrap_hook

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Applies the pending updates to G: GR(:,:,nf) <- GR(:,:,nf) - X Y^T (ZGEMM) per computed
!> flavour, and resets the pending counters. No-op when nothing is pending.
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> Complex
!--------------------------------------------------------------------
   subroutine Lowrank_flush(GR)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)

      Integer :: nf, nf_eff, j

      if (n_pend == 0) return
      select case (lr_scheme)
      case (SCHEME_DELAYED)
         Do nf_eff = 1,N_FL_eff
            j = n_cols(nf_eff)
            if (j == 0) cycle
            nf = Calc_Fl_map(nf_eff)
            call ZGEMM('N', 'T', Ndim, Ndim, j, cmone, X(1,1,nf_eff), Ndim, Y(1,1,nf_eff), Ndim, &
                 &     cone, GR(1,1,nf), Ndim)
         Enddo
         n_cols = 0
      case default
         Write(error_unit,*) 'Lowrank_flush: only the delayed scheme is implemented'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      end select
      n_pend     = 0
      prop_valid = .false.

   end subroutine Lowrank_flush

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Materialises the current Green function of flavour nf, Gcur = GR(:,:,nf) - X Y^T,
!> without changing GR or the pending state (tests and debugging).
!> @param [in] GR(Ndim,Ndim,N_FL)
!> Complex
!> @param [in] nf
!> Integer
!> \verbatim
!>  Actual flavour index (third index of GR). A flavour that is not computed is copied.
!> \endverbatim
!> @param [out] Gcur(Ndim,Ndim)
!> Complex
!--------------------------------------------------------------------
   subroutine Lowrank_current_G(GR, nf, Gcur)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In)  :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In)  :: nf
      Complex (Kind=Kind(0.d0)), Intent(Out) :: Gcur(Ndim,Ndim)

      Integer :: nf_eff, nf_e, j

      Gcur = GR(:,:,nf)
      if (n_pend == 0) return
      if (lr_scheme /= SCHEME_DELAYED) then
         Write(error_unit,*) 'Lowrank_current_G: only the delayed scheme is implemented'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      nf_e = 0
      Do nf_eff = 1,N_FL_eff
         if (Calc_Fl_map(nf_eff) == nf) nf_e = nf_eff
      Enddo
      if (nf_e == 0) return
      j = n_cols(nf_e)
      if (j > 0) call ZGEMM('N', 'T', Ndim, Ndim, j, cmone, X(1,1,nf_e), Ndim, Y(1,1,nf_e), Ndim, &
           &                cone, Gcur, Ndim)

   end subroutine Lowrank_current_G

end module Lowrank_update_mod

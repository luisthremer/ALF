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
!> and N_delay (bound on the pending rank per flavour, see below).
!> Eligibility (docs/superpowers/specs/2026-10-31-rankN-design.md, D1), the same for both schemes:
!> vertex n is eligible iff its support max_nf Op_V(n,nf)%N <= N_P_LIM = 8, whatever its structure
!> (diagonal or not, any N_non_zero, field types 1-4, g_t). Ineligible vertices are flushed and then
!> updated with the standard Upgrade2. See docs/algorithm-notes.md, sections 9-12.
!>
!> State per computed flavour nf_eff = 1..N_FL_eff (flavour Calc_Fl_map(nf_eff)):
!>  *  delayed  : G = G0 - X Y^T, X(Ndim, N_cols_max, nf_eff), Y(Ndim, N_cols_max, nf_eff). An accept
!>                adds r = Op_V(n,nf)%N_non_zero columns in flavour nf (none if all Delta_m = 0, D14);
!>                the hook after every wrap G <- L G R of an eligible vertex applies
!>                X(P,:) <- L X(P,:), Y(P,:) <- R^T Y(P,:) (docs/algorithm-notes.md, sections 9, 11).
!>  *  submatrix: M-form (docs/algorithm-notes.md, section 12; spec D3, D4):
!>                G = G0 - G0(:,S) M (E_S^T - G0(S,:)), with the pending sites S(1:k, nf_eff)
!>                (distinct), slot(Ndim, nf_eff) (site -> position in S, 0 if absent) and the dense
!>                M = Msub(1:k, 1:k, nf_eff); k = n_cols(nf_eff). An accept pads S with its active sites
!>                not yet in S (zero rows / columns of M) and applies M <- M + a Delta Mat^{-1} b^T, one
!>                formula for disjoint, overlapping and repeated sites (no flush on a repeated site); none
!>                of this if all Delta_m = 0 in that flavour (D14). The hook after every wrap applies
!>                M <- L_S M R_S (diagonal flavours: row / column scaling); a non-diagonal wrap that only
!>                partially overlaps S first pads the sites P \ S (D8).
!> Capacity (D5): the pending rank p of a flavour is n_cols(nf_eff) (delayed: columns of X;
!> submatrix: k = |S|, padded sites included). An accept adding r' (delayed: r columns; submatrix:
!> active sites not yet in S) flushes first if p > 0 and p + r' > N_delay; every accept, and for the
!> submatrix scheme every hook after its padding and transform, flushes afterwards if p >= N_delay.
!> So the arrays cannot overflow (N_delay = 1: flush after every accept, as for rank 1).
!> Submatrix accuracy guard (spec D15): in a block that contains an accepted non-diagonal update or a
!> non-diagonal hook transform of M, the threshold N_delay is replaced by min(N_delay, K_ND), K_ND = 1,
!> i.e. such a block is flushed right away: the M-form is used for diagonal content only. Reason: with
!> non-diagonal content the M-form's forward error per block reached ~1e-11 (G = G0 - G0(:,S) M (...)
!> cancels), which accumulated in the Kondo statistical run to Precision Green max 2.6e-9 vs 4.0e-11 for
!> standard (same trajectory); K_ND = 2 or 1 restores 4.9e-11, K_ND = 4 and 8 do not (3.2e-9, 2.6e-9).
!> Unit tier: rep_mixed 5.0e-10 without any guard. Delayed has no such loss and is the scheme of choice
!> for models with non-diagonal (bond) vertices.
!> Allocated: delayed N_cols_max = min(max(N_delay, N_P_LIM), size(Op_V,1)*N_P_LIM) (G is flushed at
!> every slice end); submatrix N_cols_max = min(max(N_delay, N_P_LIM), Ndim) + N_P_LIM (headroom for
!> the hook padding).
!>
!> Cross-slice delay (research prototype, opt-in key Lowrank_cross_slice in &VAR_QMC; spec D10,
!> docs/algorithm-notes.md section 12.8): Lowrank_init(..., cross_slice=.true.) keeps the pending
!> state across slice ends. WRAPGRUP/WRAPGRDO then skip their slice-end flush and call
!> Lowrank_hop_up / Lowrank_hop_down right before ALF's hopping wrap of G0, which carry the pending
!> state through that wrap in delayed form:
!>   up  (G <- B_T G B_T^{-1}): X <- B_T X (Hop_mod_mmthr), Y^T <- Y^T B_T^{-1} (Hop_mod_mmthl_m1),
!>   down(G <- B_T^{-1} G B_T): X <- B_T^{-1} X (Hop_mod_mmthr_m1), Y^T <- Y^T B_T (Hop_mod_mmthl),
!> Y^T held in a j x Ndim buffer (YTbuf; the OpT types act on rectangular operands). The
!> representation of the pending state is lr_mode: equal to lr_scheme, except that a submatrix block
!> carried through a hopping wrap is first converted to delayed form,
!>   X(:,1:k) = G0(:,S) M,  Y(:,m) = e_{S_m} - G0(S_m,:)^T  (so that G = G0 - X Y^T exactly),
!> and continues as a delayed block (delayed accepts, hooks, capacity, no K_ND guard) until the next
!> flush, which returns lr_mode to lr_scheme. If that flush is the capacity pre-flush of an accept, the
!> accept is applied in the submatrix form with the submatrix capacity and the K_ND guard (Accept_capacity
!> is re-evaluated after the mode change; debug builds assert the guard, Dbg_check_nd_guard). Hooks cannot
!> change the mode mid-way (their only flush is the final capacity flush). G must be flushed before anything else reads it: main.F90
!> calls Lowrank_flush at every read point when the key is on. Refused (terminate) with
!> N_Global_tau > 0. Lowrank_carry_similarity(GR, Sm, Sminv) is the dense analogue for the unit test.
!> Without the key nothing of this is allocated or called.
!>
!> Debug mode (-DALF_LOWRANK_DEBUG, both schemes): a shadow Gsh(Ndim,Ndim,N_FL) is set to GR at
!> the first accept after a flush (n_pend 0 -> 1), receives the full rank-r Woodbury update
!> (docs/algorithm-notes.md, section 10) on every accept and the general wrap G <- L G R of every
!> hook. After every accept and every hook, max|Lowrank_current_G - Gsh|/max|Gsh| <= TOL_UNIT is
!> checked per computed flavour, and after every flush max|GR - Gsh|/max|Gsh| <= TOL_UNIT; a
!> violation terminates. After every flush the number of checks, the largest relative error so
!> far and (submatrix) the number of padded sites are (re)written to the file
!> lowrank_debug_summary.txt in the run directory. On every accept
!> the matrix Mat = I + (I - G_PP) Delta of the proposal (computed from the low-rank G_PP) is also
!> compared with Mat_sh from the shadow, and the shadow update uses Mat_sh, so a stale G_PP read in
!> Lowrank_propose is caught. Self-tests (never defined by default; they only take effect together
!> with ALF_LOWRANK_DEBUG): -DLOWRANK_DEBUG_SELFTEST skips the X-row scaling of the first diagonal
!> hook (delayed) or the first row transform of M in a hook (submatrix: one M row scaling, or the
!> L product of a non-diagonal wrap); -DLOWRANK_DEBUG_SELFTEST_PROPOSE omits the
!> X.Y correction of G_PP in Lowrank_propose until a perturbed proposal is accepted. Each must
!> make the debug check fire.
!--------------------------------------------------------------------

module Lowrank_update_mod
   use runtime_error_mod
   use iso_fortran_env, only: output_unit, error_unit
   implicit none
   private

   Integer, parameter, public :: SCHEME_STANDARD = 0, SCHEME_DELAYED = 1, SCHEME_SUBMATRIX = 2

   public :: Lowrank_init, Lowrank_finalize, Lowrank_initialised, Lowrank_active, Lowrank_eligible, &
        &    Lowrank_pending, Lowrank_assert_flushed
   public :: Lowrank_max_pending_cols
   public :: Lowrank_propose, Lowrank_accept, Upgrade_lowrank, Lowrank_wrap_hook, &
        &    Lowrank_flush, Lowrank_current_G
   ! Cross-slice delay (spec D10, opt-in)
   public :: Lowrank_cross_slice_active, Lowrank_hop_up, Lowrank_hop_down, Lowrank_carry_similarity, &
        &    Lowrank_read_point
#ifdef ALF_LOWRANK_DEBUG
   public :: Lowrank_debug_check
#endif

   !> Maximal support (Op_V(n,nf)%N) of a vertex eligible for a low-rank scheme (D1)
   Integer, parameter, public :: N_P_LIM = 8

   ! ---- common state ----------------------------------------------------
   Logical              :: lr_initialised = .false.  ! true after Lowrank_init (any scheme)
   Integer              :: lr_scheme  = SCHEME_STANDARD
   !> Representation of the pending state: = lr_scheme, except SCHEME_DELAYED for a submatrix block that
   !> was converted to delayed form at a cross-slice hopping wrap (reset to lr_scheme by every flush)
   Integer              :: lr_mode    = SCHEME_STANDARD
   Logical              :: lr_cross   = .false.  ! cross-slice delay on (low-rank scheme only)
   Integer              :: lr_N_delay = 0
   Integer              :: N_cols_max = 0       ! allocated columns per flavour (capacity, see the header)
   Integer              :: n_pend     = 0       ! pending accepted vertices since the last flush
   Integer              :: lr_M = 0, lr_Ndim = 0, lr_N_FL_eff = 0  ! sizes at init, checked in Lowrank_propose
   Logical, allocatable :: eligible(:)          ! eligible(n), n = 1..size(Op_V,1)
   Integer, allocatable :: n_cols(:)            ! used columns per computed flavour
   Integer              :: last_accept_cols = 0 ! maxval(n_cols) right after the last accept's append (diagnostic)

   ! ---- cache of the last Lowrank_propose (consumed by Lowrank_accept) ---
   Logical                                :: prop_valid = .false.
   Integer                                :: prop_n_op = 0, prop_nt = 0
   Complex (Kind=Kind(0.d0))              :: prop_Hs
   Integer,                   allocatable :: prop_r(:)        ! N_non_zero per computed flavour
   Complex (Kind=Kind(0.d0)), allocatable :: prop_Delta(:,:)   ! Delta_m = e^{Z1 E_m} - 1, (N_P_LIM, N_FL_eff)
   Complex (Kind=Kind(0.d0)), allocatable :: prop_Mat(:,:,:)   ! Mat = I + (I - G_PP) Delta, (N_P_LIM, N_P_LIM, N_FL_eff)
   Complex (Kind=Kind(0.d0)), allocatable :: prop_d(:)         ! det(Mat) as Upgrade2, per computed flavour

   Complex (Kind=Kind(0.d0)), parameter :: cone  = cmplx( 1.d0, 0.d0, kind(0.d0))
   Complex (Kind=Kind(0.d0)), parameter :: cmone = cmplx(-1.d0, 0.d0, kind(0.d0))
   Complex (Kind=Kind(0.d0)), parameter :: czero = cmplx( 0.d0, 0.d0, kind(0.d0))

   ! ---- delayed ---------------------------------------------------------
   Complex (Kind=Kind(0.d0)), allocatable :: X(:,:,:), Y(:,:,:)
   ! workspaces: rows P of X / Y (Xg, Yg) and a product (Tg), (N_P_LIM, N_cols_max); C = G(:,P_a) (Ndim, N_P_LIM)
   Complex (Kind=Kind(0.d0)), allocatable :: Xg(:,:), Yg(:,:), Tg(:,:), Cd(:,:)

   ! ---- submatrix (M-form) ----------------------------------------------
   !> Submatrix: cap on the pending sites k of a block that contains an accepted non-diagonal update
   !> (accuracy guard, spec Sec. 6; see the module header)
   Integer, parameter :: K_ND = 1   ! spec D15: non-diagonal content is never kept in the M-form
   Logical            :: blk_nondiag = .false.   ! the current block contains an accepted non-diagonal update
   ! S(1:k,nf_eff) pending sites, slot(site,nf_eff) position in S (0 if absent), Msub(1:k,1:k,nf_eff) = M
   ! (docs/algorithm-notes.md, section 12.1); all sized N_cols_max = capS
   Integer,                   allocatable :: S(:,:), slot(:,:)
   Complex (Kind=Kind(0.d0)), allocatable :: Msub(:,:,:)
   ! cache of the last proposal per computed flavour (k at proposal time, r active sites):
   ! prop_W(1:k,1:r) = M (J - G0(S,P_a)) (so that G(P_a,P_a) = G0(P_a,P_a) - G0(P_a,S) W and, on the
   ! padded S', a = J - [W; 0]); prop_GPS(1:r,1:k) = G0(P_a,S)
   Complex (Kind=Kind(0.d0)), allocatable :: prop_W(:,:,:), prop_GPS(:,:,:)
   ! workspaces: rows / columns of M on P (Xg, Tg as for delayed; Mc, Tc), a (Aw), a K (Ak), b^T (Bt)
   Complex (Kind=Kind(0.d0)), allocatable :: Mc(:,:), Tc(:,:), Aw(:,:), Ak(:,:), Bt(:,:)
   Complex (Kind=Kind(0.d0)), allocatable :: Cw(:,:), Rw(:,:), Mw(:,:)  ! flush workspace: G0(:,S), E_S^T-G0(S,:), M Rw

   ! ---- cross-slice delay -------------------------------------------------
   ! j x Ndim buffer for Y^T (pointer-remapped to a contiguous (1:j,1:Ndim) array for the Hop_mod routines)
   Complex (Kind=Kind(0.d0)), allocatable, target :: YTbuf(:)

#ifdef ALF_LOWRANK_DEBUG
   ! ---- debug: shadow G with standard rank-1 updates ---------------------
   Real (Kind=Kind(0.d0)), parameter      :: DBG_TOL = 1.d-10     ! TOL_UNIT (CLAUDE.md, Tolerances)
   Complex (Kind=Kind(0.d0)), allocatable :: Gsh(:,:,:)
   Real (Kind=Kind(0.d0))                 :: dbg_max_err  = 0.d0
   Integer (Kind=Kind(0_8))               :: dbg_n_checks = 0
   Integer                                :: dbg_last_nop = 0, dbg_last_nt = 0  ! last accept
   Real (Kind=Kind(0.d0))                 :: dbg_max_derr = 0.d0  ! max relative error of d
#ifdef LOWRANK_DEBUG_SELFTEST_PROPOSE
   Logical                                :: dbg_selftest_prop_done = .false.  ! perturbed proposal accepted
   Logical                                :: dbg_prop_perturbed     = .false.  ! last proposal perturbed
#endif
#ifdef LOWRANK_DEBUG_SELFTEST
   Logical                                :: dbg_selftest_done = .false.
#endif
   Integer (Kind=Kind(0_8))               :: dbg_n_pad_acc  = 0   ! submatrix: sites added to S by accepts
   Integer (Kind=Kind(0_8))               :: dbg_n_pad_hook = 0   ! submatrix: sites padded by hooks (P \ S)
   Integer (Kind=Kind(0_8))               :: dbg_n_carry = 0      ! cross-slice: pending blocks carried through a hopping wrap
   ! Submatrix M-form holds non-diagonal content (an accepted non-diagonal update or a non-diagonal hook
   ! transform), tracked independently of blk_nondiag for the K_ND guard assertion (Dbg_check_nd_guard)
   Logical                                :: dbg_nd_in_M = .false.
   Integer (Kind=Kind(0_8))               :: dbg_n_conv  = 0      ! cross-slice: submatrix blocks converted to delayed form
#endif

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
!>  Bound on the pending rank per flavour (delayed: columns of X; submatrix: pending sites k = |S|,
!>  padded sites included), see the module header. Must be >= 1 for a non-standard scheme; ignored for "standard".
!> \endverbatim
!> @param [in] cross_slice
!> Logical, optional (default .false.)
!> \verbatim
!>  Cross-slice delay (spec D10, module header). Ignored for "standard".
!> \endverbatim
!> @param [in] N_Global_tau
!> Integer, optional (default 0)
!> \verbatim
!>  Number of Global_tau moves per slice. cross_slice with N_Global_tau > 0 and a low-rank scheme
!>  terminates (Wrapgr_Random_update reads G at the WRAPGRDO entry and after the slice-end flush).
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_init(scheme_name, N_delay, cross_slice, N_Global_tau)

      Use Hamiltonian_main, only: Op_V, Ndim, N_FL, N_FL_eff
      Use Files_mod, only: str_to_upper
#ifdef MPI
      Use Hamiltonian_main, only: Group_Comm
      Use mpi
#endif
      Implicit none

      Character (len=*), Intent(In) :: scheme_name
      Integer,           Intent(In) :: N_delay
      Logical, optional, Intent(In) :: cross_slice
      Integer, optional, Intent(In) :: N_Global_tau

      Character (len=:), allocatable :: name
      Character (len=256)            :: line, line_inel
      Character (len=64)             :: file_info
      Integer :: n, nf, M, n_inel, n_first, unit_info
      Logical :: print_info
#ifdef MPI
      Integer :: ierr, irank, irank_g, isize_g, igroup
#endif

      ! A re-init must not silently drop pending state
      if (lr_initialised) call Lowrank_assert_flushed("Lowrank_init")
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
      lr_mode    = lr_scheme
      lr_cross   = .false.
      if (present(cross_slice)) lr_cross = cross_slice
      if (lr_cross .and. present(N_Global_tau)) then
         if (N_Global_tau > 0) then
            Write(error_unit,'(A,A,A,I0,A)') ' Lowrank_init: Lowrank_cross_slice = .T. with Sequential_scheme = ', name, &
                 & ' is not supported with Global_tau moves (N_Global_tau = ', N_Global_tau, &
                 & ' > 0): G must be flushed at every slice end there. Set Lowrank_cross_slice = .F.'
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
      endif
      M = size(Op_V, 1)
      lr_M        = M
      lr_Ndim     = Ndim
      lr_N_FL_eff = N_FL_eff
      ! Capacity (header, D5). Delayed: at most M*N_P_LIM columns can be pending (G is flushed at
      ! every slice end; with the cross-slice delay this cap only bounds the pending rank, since
      ! Lowrank_accept flushes at min(N_delay, N_cols_max)). Submatrix: at most Ndim distinct sites; after an accept k <= max(N_delay, N_P_LIM),
      ! and a hook starts from k <= N_delay - 1 and pads at most N_P_LIM - 1 sites, hence the headroom N_P_LIM.
      ! Computed in 64 bit so that a huge N_delay cannot overflow.
      select case (lr_scheme)
      case (SCHEME_DELAYED)
         N_cols_max = int(min(max(int(N_delay,8), int(N_P_LIM,8)), int(M,8)*N_P_LIM))
      case (SCHEME_SUBMATRIX)
         N_cols_max = int(min(max(int(N_delay,8), int(N_P_LIM,8)), int(Ndim,8))) + N_P_LIM
      end select

      ! Eligibility (D1), both schemes: support <= N_P_LIM in every flavour
      allocate (eligible(M))
      n_inel  = 0
      n_first = 0
      do n = 1, M
         eligible(n) = .true.
         do nf = 1, N_FL
            if (Op_V(n,nf)%N > N_P_LIM) eligible(n) = .false.
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
      allocate (prop_r(N_FL_eff), prop_Delta(N_P_LIM, N_FL_eff), prop_Mat(N_P_LIM, N_P_LIM, N_FL_eff))
      allocate (prop_d(N_FL_eff))
      prop_r     = 0
      prop_Delta = cmplx(0.d0, 0.d0, kind(0.d0))
      prop_Mat   = cmplx(0.d0, 0.d0, kind(0.d0))
      prop_d     = cmplx(1.d0, 0.d0, kind(0.d0))
      prop_valid = .false.
      select case (lr_scheme)
      case (SCHEME_DELAYED)
         allocate (X(Ndim, N_cols_max, N_FL_eff), Y(Ndim, N_cols_max, N_FL_eff))
         allocate (Xg(N_P_LIM, N_cols_max), Yg(N_P_LIM, N_cols_max), Tg(N_P_LIM, N_cols_max))
         allocate (Cd(Ndim, N_P_LIM))
         X = cmplx(0.d0, 0.d0, kind(0.d0))
         Y = cmplx(0.d0, 0.d0, kind(0.d0))
         if (lr_cross) allocate (YTbuf(int(N_cols_max,8)*int(Ndim,8)))
      case (SCHEME_SUBMATRIX)
         allocate (S(N_cols_max, N_FL_eff), slot(Ndim, N_FL_eff))
         allocate (Msub(N_cols_max, N_cols_max, N_FL_eff))
         allocate (prop_W(N_cols_max, N_P_LIM, N_FL_eff), prop_GPS(N_P_LIM, N_cols_max, N_FL_eff))
         allocate (Xg(N_P_LIM, N_cols_max), Tg(N_P_LIM, N_cols_max))
         allocate (Mc(N_cols_max, N_P_LIM), Tc(N_cols_max, N_P_LIM))
         allocate (Aw(N_cols_max, N_P_LIM), Ak(N_cols_max, N_P_LIM), Bt(N_P_LIM, N_cols_max))
         allocate (Cw(Ndim, N_cols_max), Rw(N_cols_max, Ndim), Mw(N_cols_max, Ndim))
         S    = 0
         slot = 0
         Msub = cmplx(0.d0, 0.d0, kind(0.d0))
         if (lr_cross) then
            ! A block converted to delayed form at a slice end (k <= N_cols_max columns) continues as delayed
            allocate (X(Ndim, N_cols_max, N_FL_eff), Y(Ndim, N_cols_max, N_FL_eff))
            allocate (Yg(N_P_LIM, N_cols_max), Cd(Ndim, N_P_LIM))
            allocate (YTbuf(int(N_cols_max,8)*int(Ndim,8)))
            X = cmplx(0.d0, 0.d0, kind(0.d0))
            Y = cmplx(0.d0, 0.d0, kind(0.d0))
         endif
      end select
#ifdef ALF_LOWRANK_DEBUG
      dbg_n_pad_acc  = 0
      dbg_n_pad_hook = 0
      dbg_n_carry    = 0
      dbg_n_conv     = 0
#endif

      ! Report (group rank 0 only)
      write(line,'(A,A,A,I0)') 'Sequential scheme: ', name, ', N_delay = ', N_delay
      if (lr_cross) line = trim(line) // ', cross-slice delay on (Lowrank_cross_slice)'
      if (n_inel > 0) then
         write(line_inel,'(A,A,A,I0,A,I0,A,I0,A,I0,A)') 'Sequential scheme ', name, ': ', n_inel, ' of ', M, &
              & ' vertices have a support larger than ', N_P_LIM, &
              & '; they use the standard update (first: n=', n_first, ')'
      endif
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
      if (lr_initialised) call Lowrank_assert_flushed("Lowrank_finalize")
      if (allocated(eligible)) deallocate (eligible)
      if (allocated(n_cols))   deallocate (n_cols)
      if (allocated(X))        deallocate (X)
      if (allocated(Y))        deallocate (Y)
      if (allocated(Xg))       deallocate (Xg)
      if (allocated(Yg))       deallocate (Yg)
      if (allocated(Tg))       deallocate (Tg)
      if (allocated(Cd))       deallocate (Cd)
      if (allocated(S))        deallocate (S)
      if (allocated(slot))     deallocate (slot)
      if (allocated(Msub))     deallocate (Msub)
      if (allocated(prop_W))   deallocate (prop_W)
      if (allocated(prop_GPS)) deallocate (prop_GPS)
      if (allocated(Mc))       deallocate (Mc)
      if (allocated(Tc))       deallocate (Tc)
      if (allocated(Aw))       deallocate (Aw)
      if (allocated(Ak))       deallocate (Ak)
      if (allocated(Bt))       deallocate (Bt)
      if (allocated(Cw))       deallocate (Cw)
      if (allocated(Rw))       deallocate (Rw)
      if (allocated(Mw))       deallocate (Mw)
      if (allocated(YTbuf))    deallocate (YTbuf)
      if (allocated(prop_r))     deallocate (prop_r)
      if (allocated(prop_Delta)) deallocate (prop_Delta)
      if (allocated(prop_Mat))   deallocate (prop_Mat)
      if (allocated(prop_d))     deallocate (prop_d)
#ifdef ALF_LOWRANK_DEBUG
      if (allocated(Gsh))        deallocate (Gsh)
#endif
      prop_valid = .false.
      lr_initialised = .false.
      lr_scheme  = SCHEME_STANDARD
      lr_mode    = SCHEME_STANDARD
      lr_cross   = .false.
      lr_N_delay = 0
      N_cols_max = 0
      n_pend     = 0
      last_accept_cols = 0
      blk_nondiag = .false.
      lr_M = 0; lr_Ndim = 0; lr_N_FL_eff = 0
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
!> True iff a low-rank scheme is active and vertex n_op is eligible for it (support <= N_P_LIM
!> in every flavour, both schemes).
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
!> Diagnostic (unit tests): the largest pending rank of any computed flavour (delayed: columns of X;
!> submatrix: pending sites) right after the append of the last Lowrank_accept, before its capacity
!> flush. With the capacity rule it never exceeds max(N_delay, r) (r = N_non_zero of that vertex).
!> 0 before the first accept after Lowrank_init.
!--------------------------------------------------------------------
   integer function Lowrank_max_pending_cols()
      Implicit none
      Lowrank_max_pending_cols = last_accept_cols
   end function Lowrank_max_pending_cols

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
!> The scalar expressions and their order are those of Upgrade2: per computed flavour,
!> r = Op_V(n_op,nf)%N_non_zero, Delta_m = e^{Z1 E_m} - 1 and Upgrade2's loop
!> Mat(n,m) = - Delta_m G(P_n,P_m) + delta_nm e^{Z1 E_m}, i.e. Mat = I + (I - G_PP) Delta, with the
!> determinant branches of Upgrade2 (r = 0: 1; r = 1: Mat(1,1); r = 2: the s1/s2 form; r > 2: Det on
!> a copy, since Det factorises in place). The only difference is that G_PP is the block of the
!> current G (docs/algorithm-notes.md, sections 10-12):
!>  *  delayed  : G_PP = G0(P_a,P_a) - X(P_a,1:j) Y(P_a,1:j)^T  (O(r^2 j))
!>  *  submatrix: (section 12.2, any P_a, sites in S or not) W = M (J - G0(S,P_a)) (ZGEMM k x r x k;
!>                J = E_S^T E_{P_a}, columns zero for sites not in S), G_PP = G0(P_a,P_a) - G0(P_a,S) W
!>                (ZGEMM r x r x k), O(k^2 r). W and G0(P_a,S) are cached for Lowrank_accept. No flush:
!>                a site already pending needs nothing special (D3, D4).
!> r, Delta, Mat and d = det(Mat) are cached per computed flavour for Lowrank_accept. See
!> docs/algorithm-notes.md, sections 1, 10-12.
!> @param [in] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0, the flushed part of the Green function.
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
      Use MyMats, only: Det
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In)    :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In)    :: n_op, nt
      Complex (Kind=Kind(0.d0)), Intent(In)    :: Hs_new
      Complex (Kind=Kind(0.d0)), Intent(Out)   :: Ratiotot

      Type (Fields)             :: nsigma_new
      Complex (Kind=Kind(0.d0)) :: Ratio(N_FL), Z1, Z, D_mat, myexp, g_loc, s1, s2
      Complex (Kind=Kind(0.d0)) :: Gb(N_P_LIM,N_P_LIM)
      Complex (Kind=Kind(0.d0)), allocatable :: Matd(:,:)
      Integer                   :: nf, nf_eff, j, m, n, r, op_dim
      Logical                   :: correct_Gpp

      ! The state is sized at Lowrank_init; a change of Op_V, Ndim or N_FL_eff needs a re-init
      if (size(Op_V,1) /= lr_M .or. Ndim /= lr_Ndim .or. N_FL_eff /= lr_N_FL_eff) then
         Write(error_unit,'(A,3(1X,I0),A,3(1X,I0))') &
              & ' Lowrank_propose: (size(Op_V,1), Ndim, N_FL_eff) changed since Lowrank_init:', &
              & size(Op_V,1), Ndim, N_FL_eff, '; at init:', lr_M, lr_Ndim, lr_N_FL_eff
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      if (.not. Lowrank_eligible(n_op)) then
         Write(error_unit,'(A,I0)') ' Lowrank_propose: vertex not eligible, n_op = ', n_op
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      Call nsigma_new%make(1,1)
#if defined(ALF_LOWRANK_DEBUG) && defined(LOWRANK_DEBUG_SELFTEST_PROPOSE)
      dbg_prop_perturbed = .false.
#endif

      ! Workspace of Det (r > 2), sized as Upgrade2's Mat: the largest N_non_zero of the computed flavours
      op_dim = Op_V(n_op,Calc_Fl_map(1))%N_non_zero
      Do nf_eff = 2,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         if (op_dim < Op_V(n_op,nf)%N_non_zero) op_dim = Op_V(n_op,nf)%N_non_zero
      Enddo
      if (op_dim > 2) allocate (Matd(op_dim,op_dim))

      ! Same expressions and order as Upgrade2
      nf = 1
      nsigma_new%f(1,1)  = Hs_new
      nsigma_new%t(1)    = Op_V(n_op,nf)%Type
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         g_loc = Op_V(n_op,nf)%g
         if (Op_V(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
         Z1 = g_loc * ( nsigma_new%Phi(1,1) -  nsigma%Phi(n_op,nt) )
         r = Op_V(n_op,nf)%N_non_zero
         prop_r(nf_eff) = r
         j = n_cols(nf_eff)

         ! G_PP: the r x r block of the current G on the active sites P(1:r)
         select case (lr_mode)
         case (SCHEME_DELAYED)
            ! G_PP = G0(P_a,P_a) - X(P_a,1:j) Y(P_a,1:j)^T  (no complex BLAS functions)
            correct_Gpp = (j > 0)
#if defined(ALF_LOWRANK_DEBUG) && defined(LOWRANK_DEBUG_SELFTEST_PROPOSE)
            ! Self-test of the Mat check: read the stale G0(P_a,P_a) (omit the X.Y correction) until a
            ! perturbed proposal has been accepted
            if (j > 0 .and. r > 0 .and. .not. dbg_selftest_prop_done) then
               dbg_prop_perturbed = .true.
               correct_Gpp = .false.
            endif
#endif
            do m = 1, r
               do n = 1, r
                  Gb(n,m) = GR(Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m), nf)
                  if (correct_Gpp) Gb(n,m) = Gb(n,m) - &
                       & sum( X(Op_V(n_op,nf)%P(n),1:j,nf_eff) * Y(Op_V(n_op,nf)%P(m),1:j,nf_eff) )
               enddo
            enddo
         case (SCHEME_SUBMATRIX)
            call Submatrix_Gpp(GR, n_op, nf, nf_eff, r, Gb)
         end select

         ! Mat = I + (I - G_PP) Delta, Upgrade2's loop
         Do m = 1,r
            myexp = exp( Z1* Op_V(n_op,nf)%E(m) )
            Z = myexp - 1.d0
            prop_Delta(m,nf_eff) = Z
            do n = 1,r
               prop_Mat(n,m,nf_eff) = - Z * Gb(n,m)
            Enddo
            prop_Mat(m,m,nf_eff) = myexp + prop_Mat(m,m,nf_eff)
         Enddo
         ! det(Mat), Upgrade2's branches
         If (r == 0 ) then
            D_mat = 1.0d0
         elseIf (r == 1 ) then
            D_mat = prop_Mat(1,1,nf_eff)
         elseif (r == 2 ) then
            s1 = prop_Mat(1,1,nf_eff)*prop_Mat(2,2,nf_eff)
            s2 = prop_Mat(2,1,nf_eff)*prop_Mat(1,2,nf_eff)
            If (Abs(s1) > Abs(s2)) then
               D_mat = s1*(1.D0 - s2/s1)
            else
               D_mat = s2*(s1/s2 - 1.D0)
            Endif
         else
            ! Det factorises its argument in place: pass a copy
            Matd(1:r,1:r) = prop_Mat(1:r,1:r,nf_eff)
            D_mat = Det(Matd,r)
         endif
         prop_d(nf_eff) = D_mat
         Ratio(nf) =  D_Mat * exp( Z1*Op_V(n_op,nf)%alpha )
      Enddo
      if (allocated(Matd)) deallocate (Matd)

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
!> @brief
!> Submatrix only (docs/algorithm-notes.md, section 12.2): the r x r block Gb = G(P_a,P_a) of the
!> current G of flavour nf (computed index nf_eff, P_a = Op_V(n_op,nf)%P(1:r)),
!>   G_PP = G0(P_a,P_a) - G0(P_a,S) W,   W = M (J - G0(S,P_a))     (O(k^2 r), k = n_cols(nf_eff)),
!> with J = E_S^T E_{P_a} (columns zero for sites not in S). W and G0(P_a,S) are cached in prop_W,
!> prop_GPS for Lowrank_accept. GR is only read.
!--------------------------------------------------------------------
   subroutine Submatrix_Gpp(GR, n_op, nf, nf_eff, r, Gb)

      Use Hamiltonian_main, only: Op_V, Ndim, N_FL
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In)  :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In)  :: n_op, nf, nf_eff, r
      Complex (Kind=Kind(0.d0)), Intent(Out) :: Gb(N_P_LIM,N_P_LIM)

      Integer :: m, n, i, p, j

      j = n_cols(nf_eff)
      do m = 1, r
         do n = 1, r
            Gb(n,m) = GR(Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m), nf)
         enddo
      enddo
      if (j == 0 .or. r == 0) return
      ! Rw(1:k,1:r) = G0(S,P_a) (workspace of the flush), prop_GPS(1:r,1:k) = G0(P_a,S)
      do m = 1, r
         p = Op_V(n_op,nf)%P(m)
         do i = 1, j
            Rw(i,m) = GR(S(i,nf_eff),p,nf)
         enddo
         do i = 1, j
            prop_GPS(m,i,nf_eff) = GR(p,S(i,nf_eff),nf)
         enddo
      enddo
      ! W = - M G0(S,P_a) + M J   (M J: column slot(P_m) of M, if P_m is pending)
      call ZGEMM('N', 'N', j, r, j, cmone, Msub(1,1,nf_eff), N_cols_max, Rw, N_cols_max, &
           &     czero, prop_W(1,1,nf_eff), N_cols_max)
      do m = 1, r
         i = slot(Op_V(n_op,nf)%P(m), nf_eff)
         if (i /= 0) prop_W(1:j,m,nf_eff) = prop_W(1:j,m,nf_eff) + Msub(1:j,i,nf_eff)
      enddo
      ! G_PP = G0(P_a,P_a) - G0(P_a,S) W
      call ZGEMM('N', 'N', r, r, j, cmone, prop_GPS(1,1,nf_eff), N_P_LIM, prop_W(1,1,nf_eff), N_cols_max, &
           &     cone, Gb, N_P_LIM)

   end subroutine Submatrix_Gpp

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Accepts the proposal evaluated by the preceding Lowrank_propose(GR, n_op, nt, Hs_new, .).
!> @details
!>  *  delayed (docs/algorithm-notes.md, sections 10, 11): per computed flavour with r = N_non_zero
!>     and at least one Delta_m /= 0 (otherwise the flavour gets no columns, D14), P_a = P(1:r),
!>     j = n_cols:
!>       C  = G0(:,P_a) - X(:,1:j) Y(P_a,1:j)^T          (ZGEMM, Ndim x r)
!>       Rw = G0(P_a,:) - X(P_a,1:j) Y(:,1:j)^T          (ZGEMM, r x Ndim)
!>       X(:,j+1:j+r) = C Delta Mat^{-1},  Y(:,j+1:j+r) = (E_P^T - Rw)^T
!>     with Mat^{-1} from ZGETRF/ZGETRS on a copy of Mat, i.e. the Woodbury update
!>     G' = G - G(:,P_a) Delta Mat^{-1} (E_P^T - G(P_a,:)) on the current G. r = 1 uses project 1's
!>     in-place sequence (ZCOPY, 2 x ZGEMV, ZSCAL with Delta/d), the same algebra. Capacity (D5): if a
!>     flavour has p = j > 0 pending columns and p + r > N_delay, G is flushed first.
!>  *  submatrix (docs/algorithm-notes.md, section 12.3), per computed flavour with at least one
!>     Delta_m /= 0 (otherwise nothing, D14): pad S with the active sites P_a not yet in S (zero rows
!>     and columns of M; k -> k'), then on S', with J = E_S'^T E_{P_a} and W, G0(P_a,S) from
!>     Lowrank_propose:
!>       a   = J - M~ (J - G0(S',P_a)) = J - [W; 0]                 (k' x r)
!>       b^T = J^T + G0(P_a,S') M~      = J^T + [G0(P_a,S) M, 0]     (r x k', ZGEMM)
!>       M   <- M~ + a (Delta Mat^{-1}) b^T                         (2 x ZGEMM, O(k'^2 r))
!>     with Mat^{-1} from ZGETRF/ZGETRS (r = 1: Delta Mat^{-1} = Delta/d). One formula for disjoint,
!>     overlapping and repeated sites. Capacity (D5): if a flavour has p = k > 0 pending sites and
!>     p + r' > N_delay (r' = active sites not yet in S), G is flushed first (then k = 0: a = J, b = J);
!>     N_delay is replaced by min(N_delay, K_ND) if this vertex or the current block has a non-diagonal
!>     update (accuracy guard, module header).
!> Then nsigma%f(n_op,nt) = Hs_new; G is flushed when a flavour's pending rank reaches N_delay
!> (or, delayed, the capacity N_cols_max), so a column overflow cannot occur.
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0. Modified only by the flushes (capacity).
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
      Complex (Kind=Kind(0.d0)) :: Mi(N_P_LIM,N_P_LIM), K(N_P_LIM,N_P_LIM)
      Integer                   :: nf, nf_eff, p, j, r, m, info, ncap, kp, mode_before
      Logical                   :: nd_here
      Integer                   :: r_add(lr_N_FL_eff), ipiv(N_P_LIM)
#ifdef ALF_LOWRANK_DEBUG
      Complex (Kind=Kind(0.d0)), allocatable :: Csh(:,:), Rsh(:,:)
      Complex (Kind=Kind(0.d0))              :: Mat_sh(N_P_LIM,N_P_LIM), Ksh(N_P_LIM,N_P_LIM)
      Real    (Kind=Kind(0.d0))              :: derr, dscale
      Integer                                :: n
#endif

      if (.not. Lowrank_eligible(n_op)) then
         Write(error_unit,'(A,I0)') ' Lowrank_accept: vertex not eligible, n_op = ', n_op
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      if (.not. prop_valid .or. prop_n_op /= n_op .or. prop_nt /= nt .or. prop_Hs /= Hs_new) then
         Write(error_unit,'(A,I0,A,I0)') ' Lowrank_accept: no matching Lowrank_propose for n_op = ', &
              &                          n_op, ', nt = ', nt
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif

      ! Pending rank added per computed flavour: delayed r columns, submatrix the active sites not yet
      ! in S; none if every Delta_m of the flavour vanishes (D14)
      Do nf_eff = 1,N_FL_eff
         r_add(nf_eff) = prop_r(nf_eff)
         if (r_add(nf_eff) > 0) then
            if (all(prop_Delta(1:r_add(nf_eff),nf_eff) == czero)) r_add(nf_eff) = 0
         endif
         if (lr_mode == SCHEME_SUBMATRIX .and. r_add(nf_eff) > 0) then
            nf = Calc_Fl_map(nf_eff)
            r_add(nf_eff) = count(slot(Op_V(n_op,nf)%P(1:prop_r(nf_eff)), nf_eff) == 0)
         endif
      Enddo
      ! Capacity (D5): flush first if a flavour with pending rank p > 0 would exceed N_delay (delayed:
      ! or the allocated capacity). The proposal cache stays valid: G itself does not change (submatrix:
      ! after a flush k = 0, and the cached W, G0(P_a,S) are not used).
      ! nd_here: this accept updates a non-diagonal flavour (computed in every mode, see below)
      nd_here = .false.
      Do nf_eff = 1,N_FL_eff
         r = prop_r(nf_eff)
         if (r == 0) cycle
         if (all(prop_Delta(1:r,nf_eff) == czero)) cycle      ! D14: no update in this flavour
         if (.not. Op_V(n_op,Calc_Fl_map(nf_eff))%diag) nd_here = .true.
      Enddo
      call Accept_capacity(nd_here, ncap)
      if (any(n_cols > 0 .and. n_cols + r_add > ncap)) then
         mode_before = lr_mode
         call Lowrank_flush(GR)
         ! Cross-slice: the flush may have returned a converted (delayed-form) block to the submatrix form
         ! (lr_mode = lr_scheme), so this accept now takes the submatrix branch: the capacity, including the
         ! accuracy guard K_ND (spec D15), must be that of the mode the accept is applied in. Without the
         ! key lr_mode never changes here and ncap is kept as before.
         if (lr_mode /= mode_before) call Accept_capacity(nd_here, ncap)
      endif
      if (nd_here .and. lr_mode == SCHEME_SUBMATRIX) blk_nondiag = .true.

#ifdef ALF_LOWRANK_DEBUG
      ! Shadow: G before the first accept after a flush is G0; then the rank-r Woodbury update
      ! G' = G - G(:,P_a) Delta Mat^{-1} (E_P^T - G(P_a,:)) (docs/algorithm-notes.md, section 10)
      if (.not. allocated(Gsh)) allocate (Gsh(Ndim,Ndim,N_FL))
      if (n_pend == 0) then
         Do nf_eff = 1,N_FL_eff
            nf = Calc_Fl_map(nf_eff)
            Gsh(:,:,nf) = GR(:,:,nf)
         Enddo
      endif
#ifdef LOWRANK_DEBUG_SELFTEST_PROPOSE
      if (dbg_prop_perturbed) dbg_selftest_prop_done = .true.
#endif
      allocate (Csh(Ndim,N_P_LIM), Rsh(N_P_LIM,Ndim))
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         r  = prop_r(nf_eff)
         if (r == 0) cycle
         ! Mat from the shadow (Upgrade2's loop), independent of the low-rank G_PP read in Lowrank_propose
         do m = 1, r
            do n = 1, r
               Mat_sh(n,m) = - prop_Delta(m,nf_eff) * Gsh(Op_V(n_op,nf)%P(n), Op_V(n_op,nf)%P(m), nf)
            enddo
            Mat_sh(m,m) = 1.d0 + prop_Delta(m,nf_eff) + Mat_sh(m,m)
         enddo
         dscale = max(maxval(abs(Mat_sh(1:r,1:r))), maxval(abs(prop_Mat(1:r,1:r,nf_eff))))
         derr   = maxval(abs(Mat_sh(1:r,1:r) - prop_Mat(1:r,1:r,nf_eff)))
         if (dscale > 0.d0) derr = derr/dscale
         ! .not. (derr <= tol) also catches NaN
         if (.not. (derr <= DBG_TOL)) then
            Write(error_unit,'(A)')                   ' Lowrank_debug: proposal ratio from a stale G (Mat /= Mat_shadow) at Lowrank_accept'
            Write(error_unit,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '   n_op = ', n_op, ', nt = ', nt, ', nf = ', nf, &
                 &                                     ', n_pend = ', n_pend, ', r = ', r
            Write(error_unit,'(A,2ES23.15)')          '   Mat(1,1) (Lowrank_propose) = ', prop_Mat(1,1,nf_eff)
            Write(error_unit,'(A,2ES23.15)')          '   Mat_sh(1,1) (shadow G)     = ', Mat_sh(1,1)
            Write(error_unit,'(A,ES12.4,A,ES12.4)')   '   max|Mat - Mat_sh|/max(|Mat|,|Mat_sh|) = ', derr, &
                 &                                     ' > TOL_UNIT = ', DBG_TOL
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
         if (derr > dbg_max_derr) dbg_max_derr = derr
         if (all(prop_Delta(1:r,nf_eff) == czero)) cycle
         ! Ksh = Delta Mat_sh^{-1}
         Mi(1:r,1:r) = Mat_sh(1:r,1:r)
         call ZGETRF(r, r, Mi, N_P_LIM, ipiv, info)
         Ksh = czero
         do m = 1, r
            Ksh(m,m) = cone
         enddo
         if (info == 0) call ZGETRS('N', r, r, Mi, N_P_LIM, ipiv, Ksh, N_P_LIM, info)
         if (info /= 0) then
            Write(error_unit,'(A,I0)') ' Lowrank_debug: singular Mat_sh at Lowrank_accept, info = ', info
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
         do m = 1, r
            Ksh(m,1:r) = prop_Delta(m,nf_eff) * Ksh(m,1:r)
         enddo
         ! Csh = Gsh(:,P_a),  Rsh = E_P^T - Gsh(P_a,:)  ->  Rsh <- Ksh Rsh,  Gsh <- Gsh - Csh Rsh
         do m = 1, r
            p = Op_V(n_op,nf)%P(m)
            Csh(:,m) = Gsh(:,p,nf)
            Rsh(m,:) = - Gsh(p,:,nf)
            Rsh(m,p) = Rsh(m,p) + 1.d0
         enddo
         Rsh(1:r,:) = matmul(Ksh(1:r,1:r), Rsh(1:r,:))
         call ZGEMM('N', 'N', Ndim, Ndim, r, cmone, Csh, Ndim, Rsh, N_P_LIM, cone, Gsh(1,1,nf), Ndim)
      Enddo
      deallocate (Csh, Rsh)
#endif

      select case (lr_mode)
      case (SCHEME_DELAYED)
         Do nf_eff = 1,N_FL_eff
            r = r_add(nf_eff)
            if (r == 0) cycle
            nf = Calc_Fl_map(nf_eff)
            j  = n_cols(nf_eff)
            if (j + r > N_cols_max) then
               Write(error_unit,'(A,I0)') ' Lowrank_accept: column overflow, N_cols_max = ', N_cols_max
               CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif
            if (r == 1) then
               ! Rank 1: project 1's validated Sherman-Morrison sequence, in place (p = P(1), z = Delta/d)
               p = Op_V(n_op,nf)%P(1)
               z = prop_Delta(1,nf_eff)/prop_d(nf_eff)
               ! X(:,j+1) <- c = G0(:,p) - X(:,1:j) Y(p,1:j)^T
               call ZCOPY(Ndim, GR(1,p,nf), 1, X(1,j+1,nf_eff), 1)
               ! Y(:,j+1) <- r = G0(p,:)^T - Y(:,1:j) X(p,1:j)^T   (rows: increment Ndim)
               call ZCOPY(Ndim, GR(p,1,nf), Ndim, Y(1,j+1,nf_eff), 1)
               if (j > 0) then
                  ! (X / Y appear as the matrix and as the output vector: disjoint columns 1:j and j+1)
                  call ZGEMV('N', Ndim, j, cmone, X(1,1,nf_eff), Ndim, Y(p,1,nf_eff), Ndim, cone, X(1,j+1,nf_eff), 1)
                  call ZGEMV('N', Ndim, j, cmone, Y(1,1,nf_eff), Ndim, X(p,1,nf_eff), Ndim, cone, Y(1,j+1,nf_eff), 1)
               endif
               ! X(:,j+1) = (Delta/d) c,  Y(:,j+1) = e_p - r
               call ZSCAL(Ndim, z, X(1,j+1,nf_eff), 1)
               Y(:,j+1,nf_eff) = - Y(:,j+1,nf_eff)
               Y(p,j+1,nf_eff) = Y(p,j+1,nf_eff) + 1.d0
               n_cols(nf_eff) = j + 1
               cycle
            endif
            ! Rank r >= 2. Cd = C = G0(:,P_a),  Y(:,j+m) = - G0(P_m,:)^T
            do m = 1, r
               p = Op_V(n_op,nf)%P(m)
               call ZCOPY(Ndim, GR(1,p,nf), 1, Cd(1,m), 1)
               Y(:,j+m,nf_eff) = - GR(p,:,nf)
            enddo
            if (j > 0) then
               ! rows P_a of the old columns
               do m = 1, r
                  p = Op_V(n_op,nf)%P(m)
                  Xg(m,1:j) = X(p,1:j,nf_eff)
                  Yg(m,1:j) = Y(p,1:j,nf_eff)
               enddo
               ! C = G0(:,P_a) - X(:,1:j) Y(P_a,1:j)^T,  Y(:,j+1:j+r) = - Rw^T = - G0(P_a,:)^T + Y(:,1:j) X(P_a,1:j)^T
               ! (the second call passes Y as A, columns 1:j, and as C, columns j+1:j+r: disjoint storage, fine
               ! for the external BLAS, which has no explicit interface)
               call ZGEMM('N', 'T', Ndim, r, j, cmone, X(1,1,nf_eff), Ndim, Yg, N_P_LIM, cone, Cd, Ndim)
               call ZGEMM('N', 'T', Ndim, r, j, cone,  Y(1,1,nf_eff), Ndim, Xg, N_P_LIM, cone, Y(1,j+1,nf_eff), Ndim)
            endif
            ! Y(:,j+m) = e_{P_m} - Rw(m,:)^T
            do m = 1, r
               p = Op_V(n_op,nf)%P(m)
               Y(p,j+m,nf_eff) = Y(p,j+m,nf_eff) + 1.d0
            enddo
            ! K = Delta Mat^{-1}
            Mi(1:r,1:r) = prop_Mat(1:r,1:r,nf_eff)
            call ZGETRF(r, r, Mi, N_P_LIM, ipiv, info)
            K = czero
            do m = 1, r
               K(m,m) = cone
            enddo
            if (info == 0) call ZGETRS('N', r, r, Mi, N_P_LIM, ipiv, K, N_P_LIM, info)
            if (info /= 0) then
               Write(error_unit,'(A,I0,A,I0)') ' Lowrank_accept: singular Mat (accepted ratio 0?), n_op = ', &
                    &                          n_op, ', info = ', info
               CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
            endif
            do m = 1, r
               K(m,1:r) = prop_Delta(m,nf_eff) * K(m,1:r)
            enddo
            ! X(:,j+1:j+r) = C Delta Mat^{-1}
            call ZGEMM('N', 'N', Ndim, r, r, cone, Cd, Ndim, K, N_P_LIM, czero, X(1,j+1,nf_eff), Ndim)
            n_cols(nf_eff) = j + r
         Enddo
      case (SCHEME_SUBMATRIX)
         Do nf_eff = 1,N_FL_eff
            r = prop_r(nf_eff)
            if (r == 0) cycle
            if (all(prop_Delta(1:r,nf_eff) == czero)) cycle      ! D14: no padding, no change
            nf = Calc_Fl_map(nf_eff)
            j  = n_cols(nf_eff)                                   ! k (0 after a capacity flush)
#ifdef ALF_LOWRANK_DEBUG
            if (.not. Op_V(n_op,nf)%diag) dbg_nd_in_M = .true.
#endif
            ! Pad S with the active sites not yet pending: zero rows / columns of M (exact)
            kp = j
            do m = 1, r
               p = Op_V(n_op,nf)%P(m)
               if (slot(p,nf_eff) == 0) then
                  call Submatrix_pad(nf_eff, p, kp)
#ifdef ALF_LOWRANK_DEBUG
                  dbg_n_pad_acc = dbg_n_pad_acc + 1
#endif
               endif
            enddo
            ! a = J - [W; 0] (k' x r),  b^T = J^T + [G0(P_a,S) M, 0] (r x k')
            Aw(1:kp,1:r) = czero
            Bt(1:r,1:kp) = czero
            if (j > 0) then
               Aw(1:j,1:r) = - prop_W(1:j,1:r,nf_eff)
               call ZGEMM('N', 'N', r, j, j, cone, prop_GPS(1,1,nf_eff), N_P_LIM, Msub(1,1,nf_eff), N_cols_max, &
                    &     czero, Bt, N_P_LIM)
            endif
            do m = 1, r
               p = slot(Op_V(n_op,nf)%P(m), nf_eff)
               Aw(p,m) = Aw(p,m) + cone
               Bt(m,p) = Bt(m,p) + cone
            enddo
            ! K = Delta Mat^{-1}  (r = 1: Delta/d)
            if (r == 1) then
               K(1,1) = prop_Delta(1,nf_eff)/prop_d(nf_eff)
            else
               Mi(1:r,1:r) = prop_Mat(1:r,1:r,nf_eff)
               call ZGETRF(r, r, Mi, N_P_LIM, ipiv, info)
               K = czero
               do m = 1, r
                  K(m,m) = cone
               enddo
               if (info == 0) call ZGETRS('N', r, r, Mi, N_P_LIM, ipiv, K, N_P_LIM, info)
               if (info /= 0) then
                  Write(error_unit,'(A,I0,A,I0)') ' Lowrank_accept: singular Mat (accepted ratio 0?), n_op = ', &
                       &                          n_op, ', info = ', info
                  CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
               endif
               do m = 1, r
                  K(m,1:r) = prop_Delta(m,nf_eff) * K(m,1:r)
               enddo
            endif
            ! M <- M~ + (a K) b^T
            call ZGEMM('N', 'N', kp, r, r, cone, Aw, N_cols_max, K, N_P_LIM, czero, Ak, N_cols_max)
            call ZGEMM('N', 'N', kp, kp, r, cone, Ak, N_cols_max, Bt, N_P_LIM, cone, Msub(1,1,nf_eff), N_cols_max)
            n_cols(nf_eff) = kp
         Enddo
      end select
      prop_valid = .false.

      ! Flip the field
      nsigma%f(n_op,nt) = Hs_new

      n_pend = n_pend + 1
      last_accept_cols = maxval(n_cols)
#ifdef ALF_LOWRANK_DEBUG
      dbg_last_nop = n_op
      dbg_last_nt  = nt
      call Lowrank_debug_check(GR, n_op, nt, "Lowrank_accept")
#endif
      ! Flush when a flavour's pending rank reaches N_delay (delayed: or the capacity)
      if (maxval(n_cols) >= ncap) call Lowrank_flush(GR)
#ifdef ALF_LOWRANK_DEBUG
      call Dbg_check_nd_guard(n_op, nt, "Lowrank_accept")
#endif

   end subroutine Lowrank_accept

!--------------------------------------------------------------------
!> @brief
!> Capacity (D5) of an accept in the current representation lr_mode: delayed min(N_delay, N_cols_max);
!> submatrix N_delay, or min(N_delay, K_ND) if this accept (nd_here) or the current block (blk_nondiag)
!> has non-diagonal content (accuracy guard, spec D15).
!--------------------------------------------------------------------
   subroutine Accept_capacity(nd_here, ncap)
      Implicit none
      Logical, Intent(In)  :: nd_here
      Integer, Intent(Out) :: ncap
      if (lr_mode == SCHEME_SUBMATRIX) then
         ncap = lr_N_delay
         if (nd_here .or. blk_nondiag) ncap = min(lr_N_delay, K_ND)
      else
         ncap = min(lr_N_delay, N_cols_max)
      endif
   end subroutine Accept_capacity

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
!> Carries the pending state (delayed: X, Y; submatrix: S, M) through a wrap G <- L G R of an
!> eligible vertex.
!> @details
!> Call once per wrap of an eligible vertex, unconditionally, after the flavour loop of the wrap
!> (the routine loops over the computed flavours itself and decides per flavour):
!>   up  : after Op_Wrapup(GR, Op_V(n_op,nf), HS_field, ., 1, nt) with (wkind, dir) = (1, +1) and after
!>         Op_Wrapup(GR, ., HS_field, ., 2, nt) with (2, +1);
!>   down: after Op_Wrapdo(GR, ., HS_field, ., 2, nt) with (2, -1) and after
!>         Op_Wrapdo(GR, ., HS_field, ., 1, nt) (new field) with (1, -1).
!> Invalidates the cache of Lowrank_propose (the frame changes). ALF applies the wrap to G0 (rows P
!> by L, columns P by R; N_P x N_P blocks); with e = exp(g_loc phi(HS_field) E), g_loc = g or g_t(nt)
!> (docs/algorithm-notes.md, section 9):
!>   up type 1: L = diag(e) U^dag, R = U diag(1/e);   up type 2: L = U, R = U^dag;
!>   down type 2: L = U^dag, R = U;                   down type 1: L = U diag(1/e), R = diag(e) U^dag.
!> For a diagonal flavour (U = 1) type 2 is a no-op and type 1 the scaling t_i = exp(dir g_loc phi E(i)).
!> Applied to the pending state:
!>  *  delayed  : X(P,:) <- L X(P,:), Y(P,:) <- R^T Y(P,:) (plain transpose; O(N_P^2 j)); diagonal
!>                flavours: X(P(i),:) *= t_i, Y(P(i),:) /= t_i (docs/algorithm-notes.md, sections 3, 11)
!>  *  submatrix: M <- L_S M R_S with L_S, R_S = L, R restricted to the pending part of P
!>                (docs/algorithm-notes.md, section 12.4): nothing if P and S are disjoint; diagonal
!>                flavours (kind 2: no-op): for m = slot(P(i)) > 0, M(m,1:k) *= t_i, M(1:k,m) /= t_i;
!>                non-diagonal flavours: if P only partially overlaps S, the sites P \ S are padded
!>                first (zero rows / columns of M; no G0 read, so this is exact after ALF's wrap), then
!>                M(P,:) <- L M(P,:), M(:,P) <- M(:,P) R (ZGEMM, O(k N_P^2)). Afterwards G is flushed if a
!>                flavour's pending rank k reaches N_delay (capacity rule, D5; min(N_delay, K_ND) in a
!>                block with a non-diagonal update, module header).
!> Consumes no random numbers. No-op if nothing is pending.
!> @param [in] n_op, nt
!> Integer
!> @param [in] HS_field
!> Complex
!> \verbatim
!>  The field used in the wrap.
!> \endverbatim
!> @param [in] wkind
!> Integer
!> \verbatim
!>  Kind of the wrap = N_type of the wrap just done: 1 or 2.
!> \endverbatim
!> @param [in] dir
!> Integer
!> \verbatim
!>  +1 after Op_Wrapup, -1 after Op_Wrapdo.
!> \endverbatim
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0 after ALF's wrap. Modified only by the capacity flush of the submatrix scheme.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_wrap_hook(n_op, nt, HS_field, wkind, dir, GR)

      Use Hamiltonian_main, only: Op_V, Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Use Fields_mod
      Implicit none

      Integer,                   Intent(In)    :: n_op, nt, wkind, dir
      Complex (Kind=Kind(0.d0)), Intent(In)    :: HS_field
      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)

      Type (Fields)             :: nsigma_single
      Complex (Kind=Kind(0.d0)) :: phi, g_loc, t
      Complex (Kind=Kind(0.d0)) :: L(N_P_LIM,N_P_LIM), Rt(N_P_LIM,N_P_LIM)
      Integer                   :: nf, nf_eff, i, j, q, m, N, nov
      Integer                   :: sl(N_P_LIM)
      Logical                   :: need_scale
#if defined(ALF_LOWRANK_DEBUG) && defined(LOWRANK_DEBUG_SELFTEST)
      Logical                   :: skip_scale
#endif

      prop_valid = .false.
      if (n_pend == 0) return
      if (wkind /= 1 .and. wkind /= 2) then
         Write(error_unit,'(A,I0)') ' Lowrank_wrap_hook: wkind must be 1 or 2, wkind = ', wkind
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      if (.not. eligible(n_op)) then
         Write(error_unit,'(A,I0)') ' Lowrank_wrap_hook: vertex not eligible, n_op = ', n_op
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      ! Is there anything to transform? A flavour with pending state whose wrap is not the identity
      ! (type 1, or a non-diagonal operator); submatrix: and the vertex touches a pending site.
      need_scale = .false.
      Do nf_eff = 1,N_FL_eff
         if (n_cols(nf_eff) == 0) cycle
         nf = Calc_Fl_map(nf_eff)
         select case (lr_mode)
         case (SCHEME_DELAYED)
            if (wkind == 1 .or. .not. Op_V(n_op,nf)%diag) need_scale = .true.
         case (SCHEME_SUBMATRIX)
            if (wkind == 2 .and. Op_V(n_op,nf)%diag) cycle
            do i = 1,Op_V(n_op,nf)%N
               if (slot(Op_V(n_op,nf)%P(i),nf_eff) /= 0) need_scale = .true.
            enddo
         end select
      Enddo
#ifndef ALF_LOWRANK_DEBUG
      if (.not. need_scale) return
#endif

      ! One single-field object per hook call: phi = Phi(HS_field) for the field type of each flavour's
      ! operator (Wrap_phi_g sets the type only: no allocation per flavour)
      Call nsigma_single%make(1,1)
      nsigma_single%f(1,1) = HS_field

      select case (lr_mode)
      case (SCHEME_DELAYED)
         Do nf_eff = 1,N_FL_eff
            j = n_cols(nf_eff)
            if (j == 0) cycle
            nf = Calc_Fl_map(nf_eff)
            N  = Op_V(n_op,nf)%N
            if (Op_V(n_op,nf)%diag) then
               if (wkind == 2) cycle
               call Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)
               do i = 1,N
                  q = Op_V(n_op,nf)%P(i)
                  t = exp( real(dir, kind(0.d0)) * g_loc * phi * Op_V(n_op,nf)%E(i) )
#if defined(ALF_LOWRANK_DEBUG) && defined(LOWRANK_DEBUG_SELFTEST)
                  ! Self-test of the debug check: deliberately skip one X-row scaling (first time only)
                  skip_scale = .not. dbg_selftest_done
                  dbg_selftest_done = .true.
                  if (.not. skip_scale) X(q,1:j,nf_eff) = X(q,1:j,nf_eff) * t
#else
                  X(q,1:j,nf_eff) = X(q,1:j,nf_eff) * t
#endif
                  Y(q,1:j,nf_eff) = Y(q,1:j,nf_eff) / t
               enddo
            else
               call Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)
               call Wrap_LR(n_op, nf, phi, g_loc, wkind, dir, L, Rt)
               ! X(P,1:j) <- L X(P,1:j),  Y(P,1:j) <- R^T Y(P,1:j)
               do i = 1,N
                  Xg(i,1:j) = X(Op_V(n_op,nf)%P(i),1:j,nf_eff)
                  Yg(i,1:j) = Y(Op_V(n_op,nf)%P(i),1:j,nf_eff)
               enddo
               call ZGEMM('N', 'N', N, j, N, cone, L, N_P_LIM, Xg, N_P_LIM, czero, Tg, N_P_LIM)
               do i = 1,N
                  X(Op_V(n_op,nf)%P(i),1:j,nf_eff) = Tg(i,1:j)
               enddo
               call ZGEMM('N', 'N', N, j, N, cone, Rt, N_P_LIM, Yg, N_P_LIM, czero, Tg, N_P_LIM)
               do i = 1,N
                  Y(Op_V(n_op,nf)%P(i),1:j,nf_eff) = Tg(i,1:j)
               enddo
            endif
         Enddo

      case (SCHEME_SUBMATRIX)
         ! M <- L_S M R_S (section 12.4); only if the vertex touches a pending site in some flavour (need_scale)
         Do nf_eff = 1,N_FL_eff
            if (.not. need_scale) exit
            j = n_cols(nf_eff)
            if (j == 0) cycle
            nf = Calc_Fl_map(nf_eff)
            N  = Op_V(n_op,nf)%N
            nov = count(slot(Op_V(n_op,nf)%P(1:N), nf_eff) /= 0)
            if (nov == 0) cycle                                   ! P and S disjoint: nothing to do
            if (Op_V(n_op,nf)%diag) then
               if (wkind == 2) cycle
               call Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)
               do i = 1,N
                  q = Op_V(n_op,nf)%P(i)
                  m = slot(q,nf_eff)
                  if (m == 0) cycle
                  t = exp( real(dir, kind(0.d0)) * g_loc * phi * Op_V(n_op,nf)%E(i) )
#if defined(ALF_LOWRANK_DEBUG) && defined(LOWRANK_DEBUG_SELFTEST)
                  ! Self-test of the debug check: deliberately skip one M row scaling (first row transform)
                  skip_scale = .not. dbg_selftest_done
                  dbg_selftest_done = .true.
                  if (.not. skip_scale) Msub(m,1:j,nf_eff) = Msub(m,1:j,nf_eff) * t
#else
                  Msub(m,1:j,nf_eff) = Msub(m,1:j,nf_eff) * t
#endif
                  Msub(1:j,m,nf_eff) = Msub(1:j,m,nf_eff) / t
               enddo
            else
               ! A non-diagonal transform mixes pending rows/columns: the block counts as non-diagonal for the
               ! accuracy guard K_ND, as after a non-diagonal accept (Task 5 review, minor 1)
               blk_nondiag = .true.
#ifdef ALF_LOWRANK_DEBUG
               dbg_nd_in_M = .true.
#endif
               ! Partial overlap with a non-diagonal wrap: pad P \ S first (zeros only, no G0 read)
               if (nov < N) then
                  do i = 1,N
                     q = Op_V(n_op,nf)%P(i)
                     if (slot(q,nf_eff) == 0) then
                        call Submatrix_pad(nf_eff, q, j)
#ifdef ALF_LOWRANK_DEBUG
                        dbg_n_pad_hook = dbg_n_pad_hook + 1
#endif
                     endif
                  enddo
                  n_cols(nf_eff) = j
               endif
               do i = 1,N
                  sl(i) = slot(Op_V(n_op,nf)%P(i),nf_eff)
               enddo
               call Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)
               call Wrap_LR(n_op, nf, phi, g_loc, wkind, dir, L, Rt)
               ! M(P,1:k) <- L M(P,1:k)
               do i = 1,N
                  Xg(i,1:j) = Msub(sl(i),1:j,nf_eff)
               enddo
               call ZGEMM('N', 'N', N, j, N, cone, L, N_P_LIM, Xg, N_P_LIM, czero, Tg, N_P_LIM)
#if defined(ALF_LOWRANK_DEBUG) && defined(LOWRANK_DEBUG_SELFTEST)
               ! Self-test of the debug check: deliberately skip the row transform (first row transform)
               skip_scale = .not. dbg_selftest_done
               dbg_selftest_done = .true.
               if (skip_scale) Tg(1:N,1:j) = Xg(1:N,1:j)
#endif
               do i = 1,N
                  Msub(sl(i),1:j,nf_eff) = Tg(i,1:j)
               enddo
               ! M(1:k,P) <- M(1:k,P) R = M(1:k,P) (R^T)^T
               do i = 1,N
                  Mc(1:j,i) = Msub(1:j,sl(i),nf_eff)
               enddo
               call ZGEMM('N', 'T', j, N, N, cone, Mc, N_cols_max, Rt, N_P_LIM, czero, Tc, N_cols_max)
               do i = 1,N
                  Msub(1:j,sl(i),nf_eff) = Tc(1:j,i)
               enddo
            endif
         Enddo
      end select

#ifdef ALF_LOWRANK_DEBUG
      ! Shadow: the same wrap G <- L G R for every computed flavour (also those without pending
      ! columns, whose GR was wrapped by ALF as well)
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         N  = Op_V(n_op,nf)%N
         if (Op_V(n_op,nf)%diag) then
            if (wkind == 2) cycle
            call Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)
            do i = 1,N
               q = Op_V(n_op,nf)%P(i)
               t = exp( real(dir, kind(0.d0)) * g_loc * phi * Op_V(n_op,nf)%E(i) )
               Gsh(q,:,nf) = Gsh(q,:,nf) * t
               Gsh(:,q,nf) = Gsh(:,q,nf) / t
            enddo
         else
            call Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)
            call Wrap_LR(n_op, nf, phi, g_loc, wkind, dir, L, Rt)
            Gsh(Op_V(n_op,nf)%P,:,nf) = matmul(L(1:N,1:N), Gsh(Op_V(n_op,nf)%P,:,nf))
            Gsh(:,Op_V(n_op,nf)%P,nf) = matmul(Gsh(:,Op_V(n_op,nf)%P,nf), transpose(Rt(1:N,1:N)))
         endif
      Enddo
#endif
      Call nsigma_single%clear()
      ! Debug check: WRAPGRUP/WRAPGRDO call Lowrank_debug_check right after each hook (under
      ! ALF_LOWRANK_DEBUG); the shadow is wrapped above, before the capacity flush checks GR against it.

      ! Submatrix capacity (D5): flush after the padding and transform if a pending rank reached N_delay
      if (lr_mode == SCHEME_SUBMATRIX) then
         if (blk_nondiag) then
            if (maxval(n_cols) >= min(lr_N_delay, K_ND)) call Lowrank_flush(GR)
         else
            if (maxval(n_cols) >= lr_N_delay) call Lowrank_flush(GR)
         endif
      endif
#ifdef ALF_LOWRANK_DEBUG
      call Dbg_check_nd_guard(n_op, nt, "Lowrank_wrap_hook")
#endif

   end subroutine Lowrank_wrap_hook

!--------------------------------------------------------------------
!> @brief
!> phi = Phi(HS_field) for the field type of Op_V(n_op,nf) and g_loc = g, or g_t(nt) if allocated:
!> the factors of the wraps Op_Wrapup / Op_Wrapdo. nsigma_single (made by the caller, f(1,1) =
!> HS_field) only gets its type set: no allocation.
!--------------------------------------------------------------------
   subroutine Wrap_phi_g(nsigma_single, n_op, nf, nt, phi, g_loc)

      Use Hamiltonian_main, only: Op_V
      Use Fields_mod
      Implicit none

      Type (Fields),             Intent(InOut) :: nsigma_single
      Integer,                   Intent(In)    :: n_op, nf, nt
      Complex (Kind=Kind(0.d0)), Intent(Out)   :: phi, g_loc

      nsigma_single%t(1) = Op_V(n_op,nf)%type
      phi   = nsigma_single%Phi(1,1)
      g_loc = Op_V(n_op,nf)%g
      if (Op_V(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)

   end subroutine Wrap_phi_g

!--------------------------------------------------------------------
!> @brief
!> Wrap matrices of a non-diagonal flavour (docs/algorithm-notes.md, section 9): the wrap
!> (wkind, dir) of Op_V(n_op,nf) with phi, g_loc from Wrap_phi_g is G <- L G R on rows / columns P;
!> returns L and Rt = R^T
!> (N x N blocks of the N_P_LIM x N_P_LIM arguments). As in Op_Wrapup / Op_Wrapdo, for field
!> types 1, 2 without g_t the factors of the zero eigenvalues (index > N_non_zero) are exactly 1
!> (E_Exp in Op_set).
!--------------------------------------------------------------------
   subroutine Wrap_LR(n_op, nf, phi, g_loc, wkind, dir, L, Rt)

      Use Hamiltonian_main, only: Op_V
      Implicit none

      Integer,                   Intent(In)  :: n_op, nf, wkind, dir
      Complex (Kind=Kind(0.d0)), Intent(In)  :: phi, g_loc
      Complex (Kind=Kind(0.d0)), Intent(Out) :: L(N_P_LIM,N_P_LIM), Rt(N_P_LIM,N_P_LIM)

      Complex (Kind=Kind(0.d0)) :: ep(N_P_LIM), em(N_P_LIM)
      Integer                   :: i, k, N

      N = Op_V(n_op,nf)%N
      if (wkind == 1) then
         do i = 1,N
            if (Op_V(n_op,nf)%type < 3 .and. .not. Op_V(n_op,nf)%get_g_t_alloc() .and. &
                 & i > Op_V(n_op,nf)%N_non_zero) then
               ep(i) = cone
               em(i) = cone
            else
               ep(i) = exp( phi*g_loc*Op_V(n_op,nf)%E(i))
               em(i) = exp(-phi*g_loc*Op_V(n_op,nf)%E(i))
            endif
         enddo
      endif
      do i = 1,N
         do k = 1,N
            if (wkind == 1 .and. dir == 1) then
               ! L = diag(e) U^dag,  R^T = (U diag(1/e))^T
               L(i,k)  = ep(i) * conjg(Op_V(n_op,nf)%U(k,i))
               Rt(i,k) = Op_V(n_op,nf)%U(k,i) * em(i)
            elseif (wkind == 2 .and. dir == 1) then
               ! L = U,  R^T = (U^dag)^T = conjg(U)
               L(i,k)  = Op_V(n_op,nf)%U(i,k)
               Rt(i,k) = conjg(Op_V(n_op,nf)%U(i,k))
            elseif (wkind == 2) then
               ! L = U^dag,  R^T = U^T
               L(i,k)  = conjg(Op_V(n_op,nf)%U(k,i))
               Rt(i,k) = Op_V(n_op,nf)%U(k,i)
            else
               ! L = U diag(1/e),  R^T = (diag(e) U^dag)^T = conjg(U) diag(e)
               L(i,k)  = Op_V(n_op,nf)%U(i,k) * em(k)
               Rt(i,k) = conjg(Op_V(n_op,nf)%U(i,k)) * ep(k)
            endif
         enddo
      enddo

   end subroutine Wrap_LR

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Applies the pending updates to G per computed flavour and resets the pending counters.
!> No-op when nothing is pending.
!> @details
!>  *  delayed  : GR(:,:,nf) <- GR(:,:,nf) - X Y^T (ZGEMM)
!>  *  submatrix: C = G0(:,S) (copy), Mw = M (E_S^T - G0(S,:)) (ZGEMM k x Ndim x k),
!>                GR(:,:,nf) <- GR(:,:,nf) - C Mw (ZGEMM Ndim x Ndim x k); slot is cleared
!>                (docs/algorithm-notes.md, section 12.5; k counts padded sites)
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> Complex
!--------------------------------------------------------------------
   subroutine Lowrank_flush(GR)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)

      Integer :: nf, nf_eff, j

      if (n_pend == 0) return
      Do nf_eff = 1,N_FL_eff
         j = n_cols(nf_eff)
         if (j == 0) cycle
         nf = Calc_Fl_map(nf_eff)
         select case (lr_mode)
         case (SCHEME_DELAYED)
            call ZGEMM('N', 'T', Ndim, Ndim, j, cmone, X(1,1,nf_eff), Ndim, Y(1,1,nf_eff), Ndim, &
                 &     cone, GR(1,1,nf), Ndim)
         case (SCHEME_SUBMATRIX)
            ! Cw, Mw from G0 = GR before GR is overwritten
            call Submatrix_factors(GR, nf, nf_eff)
            call ZGEMM('N', 'N', Ndim, Ndim, j, cmone, Cw, Ndim, Mw, N_cols_max, cone, GR(1,1,nf), Ndim)
            slot(S(1:j,nf_eff),nf_eff) = 0
         end select
      Enddo
      n_cols = 0
      blk_nondiag = .false.
      lr_mode = lr_scheme          ! a converted (cross-slice) submatrix block returns to the submatrix form
#ifdef ALF_LOWRANK_DEBUG
      dbg_nd_in_M = .false.
#endif
#ifdef ALF_LOWRANK_DEBUG
      ! n_cols = 0 and n_pend > 0: Lowrank_current_G returns GR, which must equal the shadow
      call Lowrank_debug_check(GR, dbg_last_nop, dbg_last_nt, "Lowrank_flush (n_op, nt of last accept)")
#endif
      n_pend     = 0
      prop_valid = .false.
#ifdef ALF_LOWRANK_DEBUG
      call Lowrank_debug_summary()
#endif

   end subroutine Lowrank_flush

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Materialises the current Green function of flavour nf without changing GR or the pending
!> state (tests and debugging): delayed Gcur = GR(:,:,nf) - X Y^T, submatrix
!> Gcur = GR(:,:,nf) - G0(:,S) M (E_S^T - G0(S,:)) (the algebra of Lowrank_flush).
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
      nf_e = 0
      Do nf_eff = 1,N_FL_eff
         if (Calc_Fl_map(nf_eff) == nf) nf_e = nf_eff
      Enddo
      if (nf_e == 0) return
      j = n_cols(nf_e)
      if (j == 0) return
      select case (lr_mode)
      case (SCHEME_DELAYED)
         call ZGEMM('N', 'T', Ndim, Ndim, j, cmone, X(1,1,nf_e), Ndim, Y(1,1,nf_e), Ndim, &
              &     cone, Gcur, Ndim)
      case (SCHEME_SUBMATRIX)
         call Submatrix_factors(GR, nf, nf_e)
         call ZGEMM('N', 'N', Ndim, Ndim, j, cmone, Cw, Ndim, Mw, N_cols_max, cone, Gcur, Ndim)
      end select

   end subroutine Lowrank_current_G

!--------------------------------------------------------------------
!> @brief
!> Cross-slice delay: called by main.F90 (only when the key Lowrank_cross_slice is on) right before
!> a read of G (GR). Flushes the pending state, so GR is the full Green function. Under
!> ALF_LOWRANK_DEBUG it asserts flushed afterwards and, if expect_flushed, also before the flush (read
!> points at which G must already be flushed, e.g. right after a CGR or at the start of a sweep).
!> No-op when nothing is pending (also before Lowrank_init and for the standard scheme).
!> @param [inout] GR(Ndim,Ndim,N_FL)
!> @param [in] where
!> Character: name of the read point, reported by the assertion.
!> @param [in] expect_flushed
!> Logical: G must already be flushed here (checked under ALF_LOWRANK_DEBUG only).
!--------------------------------------------------------------------
   subroutine Lowrank_read_point(GR, where, expect_flushed)
      Use Hamiltonian_main, only: Ndim, N_FL
      Implicit none
      Complex (Kind=Kind(0.d0)), Intent(InOut) :: GR(Ndim,Ndim,N_FL)
      Character (len=*),         Intent(In)    :: where
      Logical,                   Intent(In)    :: expect_flushed
#ifdef ALF_LOWRANK_DEBUG
      if (expect_flushed) call Lowrank_assert_flushed(where)
#endif
      call Lowrank_flush(GR)
#ifdef ALF_LOWRANK_DEBUG
      call Lowrank_assert_flushed(where)
#endif
   end subroutine Lowrank_read_point

!--------------------------------------------------------------------
!> @brief
!> True iff a low-rank scheme is active with the cross-slice delay (Lowrank_init(..., cross_slice=.true.)).
!--------------------------------------------------------------------
   logical function Lowrank_cross_slice_active()
      Implicit none
      Lowrank_cross_slice_active = lr_cross .and. (lr_scheme /= SCHEME_STANDARD)
   end function Lowrank_cross_slice_active

!--------------------------------------------------------------------
!> @brief
!> Cross-slice only: terminates unless the cross-slice delay is active (X, Y, YTbuf are allocated only then).
!--------------------------------------------------------------------
   subroutine Cross_slice_require(where)
      Implicit none
      Character (len=*), Intent(In) :: where
      if (.not. Lowrank_cross_slice_active()) then
         Write(error_unit,'(A,A,A)') ' ', trim(where), ': requires Lowrank_init(..., cross_slice=.true.) with a low-rank scheme'
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
   end subroutine Cross_slice_require

!--------------------------------------------------------------------
!> @brief
!> Cross-slice only (docs/algorithm-notes.md, section 12.8): converts a pending submatrix block to
!> delayed form, per computed flavour with k = n_cols pending sites,
!>   X(:,1:k) = G0(:,S) M,   Y(:,m) = e_{S_m} - G0(S_m,:)^T   (m = 1..k),
!> so that G = G0 - G0(:,S) M (E_S^T - G0(S,:)) = G0 - X Y^T exactly (O(Ndim k^2)). slot is cleared, n_cols
!> and n_pend are kept, lr_mode becomes SCHEME_DELAYED until the next flush. No-op unless lr_mode is
!> SCHEME_SUBMATRIX. GR (= G0) is only read.
!--------------------------------------------------------------------
   subroutine Lowrank_to_delayed(GR)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)

      Integer :: nf, nf_eff, k, m

      if (lr_mode /= SCHEME_SUBMATRIX) return
      Do nf_eff = 1,N_FL_eff
         k = n_cols(nf_eff)
         if (k == 0) cycle
         nf = Calc_Fl_map(nf_eff)
         ! Cw = G0(:,S),  X(:,1:k) = Cw M
         do m = 1, k
            call ZCOPY(Ndim, GR(1,S(m,nf_eff),nf), 1, Cw(1,m), 1)
         enddo
         call ZGEMM('N', 'N', Ndim, k, k, cone, Cw, Ndim, Msub(1,1,nf_eff), N_cols_max, &
              &     czero, X(1,1,nf_eff), Ndim)
         ! Y(:,m) = e_{S_m} - G0(S_m,:)^T
         do m = 1, k
            Y(:,m,nf_eff) = - GR(S(m,nf_eff),:,nf)
            Y(S(m,nf_eff),m,nf_eff) = Y(S(m,nf_eff),m,nf_eff) + cone
         enddo
         slot(S(1:k,nf_eff),nf_eff) = 0
      Enddo
      lr_mode     = SCHEME_DELAYED
      blk_nondiag = .false.          ! the K_ND guard concerns the M-form only
#ifdef ALF_LOWRANK_DEBUG
      dbg_nd_in_M = .false.
#endif
#ifdef ALF_LOWRANK_DEBUG
      dbg_n_conv = dbg_n_conv + 1
      call Lowrank_debug_check(GR, dbg_last_nop, dbg_last_nt, "Lowrank_to_delayed (n_op, nt of last accept)")
#endif

   end subroutine Lowrank_to_delayed

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Cross-slice delay (spec D10; docs/algorithm-notes.md, section 12.8): carries the pending state
!> through the hopping wrap of WRAPGRUP, G <- B_T G B_T^{-1} (HOP_MOD_mmthr + HOP_MOD_mmthl_m1 on
!> time slice nt), instead of flushing it at the slice end.
!> @details
!> Call right BEFORE ALF's hopping wrap of G0 (a pending submatrix block is first converted to delayed
!> form, which reads G0 in the frame of the pending state). Per computed flavour nf_eff (hopping
!> matrices of the flavour nf = Calc_Fl_map(nf_eff)), j = n_cols:
!>   X(:,1:j) <- B_T X(:,1:j)        (Hop_mod_mmthr on the Ndim x j block),
!>   Y^T      <- Y^T B_T^{-1}        (Hop_mod_mmthl_m1 on a j x Ndim buffer holding Y(:,1:j)^T).
!> Cost O(N_bond j) with the checkerboard decomposition, O(Ndim^2 j) with a dense e^{-Dtau T}.
!> No-op if nothing is pending. Terminates unless the cross-slice delay is active.
!> @param [in] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0 before the hopping wrap (only read, by the submatrix -> delayed conversion).
!> \endverbatim
!> @param [in] nt
!> Integer
!> \verbatim
!>  Time slice of the hopping wrap (ntau1 in WRAPGRUP).
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_hop_up(GR, nt)
      Use Hamiltonian_main, only: Ndim, N_FL
      Implicit none
      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In) :: nt
      call Cross_slice_require("Lowrank_hop_up")
      call Cross_slice_hop(GR, nt, +1)
   end subroutine Lowrank_hop_up

!--------------------------------------------------------------------
!> @brief
!> Cross-slice delay: as Lowrank_hop_up for the hopping wrap at the end of WRAPGRDO,
!> G <- B_T^{-1} G B_T (Hop_mod_mmthl + Hop_mod_mmthr_m1 on time slice nt):
!>   X(:,1:j) <- B_T^{-1} X(:,1:j)   (Hop_mod_mmthr_m1),   Y^T <- Y^T B_T   (Hop_mod_mmthl on the buffer).
!> Call right BEFORE ALF's hopping wrap of G0. GR is only read (conversion).
!--------------------------------------------------------------------
   subroutine Lowrank_hop_down(GR, nt)
      Use Hamiltonian_main, only: Ndim, N_FL
      Implicit none
      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In) :: nt
      call Cross_slice_require("Lowrank_hop_down")
      call Cross_slice_hop(GR, nt, -1)
   end subroutine Lowrank_hop_down

!--------------------------------------------------------------------
!> @brief
!> Lowrank_hop_up (dir = +1) / Lowrank_hop_down (dir = -1), see there. Under ALF_LOWRANK_DEBUG the
!> shadow G receives the same hopping wrap (the caller checks it after ALF's wrap of G0).
!--------------------------------------------------------------------
   subroutine Cross_slice_hop(GR, nt, dir)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Use Hop_mod, only: Hop_mod_mmthr, Hop_mod_mmthr_m1, Hop_mod_mmthl, Hop_mod_mmthl_m1
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In) :: nt, dir

      Complex (Kind=Kind(0.d0)), pointer :: YT(:,:)
      Integer :: nf, nf_eff, j, m

      prop_valid = .false.
      if (n_pend == 0) return
      call Lowrank_to_delayed(GR)
      Do nf_eff = 1,N_FL_eff
         j = n_cols(nf_eff)
         if (j == 0) cycle
         nf = Calc_Fl_map(nf_eff)
         ! X <- B_T X (up) or B_T^{-1} X (down): the Ndim x j block, contiguous
         if (dir == 1) then
            call Hop_mod_mmthr   (X(:,1:j,nf_eff), nf, nt)
         else
            call Hop_mod_mmthr_m1(X(:,1:j,nf_eff), nf, nt)
         endif
         ! Y^T <- Y^T B_T^{-1} (up) or Y^T B_T (down), on the contiguous j x Ndim buffer
         YT(1:j,1:Ndim) => YTbuf(1:int(j,8)*int(Ndim,8))
         do m = 1, j
            YT(m,:) = Y(:,m,nf_eff)
         enddo
         if (dir == 1) then
            call Hop_mod_mmthl_m1(YT, nf, nt)
         else
            call Hop_mod_mmthl   (YT, nf, nt)
         endif
         do m = 1, j
            Y(:,m,nf_eff) = YT(m,:)
         enddo
         nullify (YT)
      Enddo
#ifdef ALF_LOWRANK_DEBUG
      ! Shadow: the same hopping wrap as ALF's on G0 (same calls and order as WRAPGRUP / WRAPGRDO)
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         if (dir == 1) then
            call Hop_mod_mmthr   (Gsh(:,:,nf), nf, nt)
            call Hop_mod_mmthl_m1(Gsh(:,:,nf), nf, nt)
         else
            call Hop_mod_mmthl   (Gsh(:,:,nf), nf, nt)
            call Hop_mod_mmthr_m1(Gsh(:,:,nf), nf, nt)
         endif
      Enddo
      dbg_n_carry = dbg_n_carry + 1
#endif

   end subroutine Cross_slice_hop

!--------------------------------------------------------------------
!> @brief
!> Cross-slice delay, dense analogue of Lowrank_hop_up for tests (unit test 36): carries the pending
!> state through a similarity G <- Sm G Sminv (Sminv = Sm^{-1}, the same for every computed flavour):
!> submatrix blocks are converted to delayed form first (reads GR = G0 before the similarity), then
!>   X(:,1:j) <- Sm X(:,1:j),   Y(:,1:j) <- Sminv^T Y(:,1:j)  (i.e. Y^T <- Y^T Sminv).
!> The caller applies the similarity to G0 afterwards. No-op if nothing is pending. Terminates unless
!> the cross-slice delay is active.
!--------------------------------------------------------------------
   subroutine Lowrank_carry_similarity(GR, Sm, Sminv)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)
      Complex (Kind=Kind(0.d0)), Intent(In) :: Sm(Ndim,Ndim), Sminv(Ndim,Ndim)

      Complex (Kind=Kind(0.d0)), allocatable :: T(:,:)
      Integer :: nf_eff, j
#ifdef ALF_LOWRANK_DEBUG
      Integer :: nf
#endif

      call Cross_slice_require("Lowrank_carry_similarity")
      prop_valid = .false.
      if (n_pend == 0) return
      call Lowrank_to_delayed(GR)
      allocate (T(Ndim, N_cols_max))
      Do nf_eff = 1,N_FL_eff
         j = n_cols(nf_eff)
         if (j == 0) cycle
         call ZGEMM('N', 'N', Ndim, j, Ndim, cone, Sm, Ndim, X(1,1,nf_eff), Ndim, czero, T, Ndim)
         X(:,1:j,nf_eff) = T(:,1:j)
         call ZGEMM('T', 'N', Ndim, j, Ndim, cone, Sminv, Ndim, Y(1,1,nf_eff), Ndim, czero, T, Ndim)
         Y(:,1:j,nf_eff) = T(:,1:j)
      Enddo
      deallocate (T)
#ifdef ALF_LOWRANK_DEBUG
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         Gsh(:,:,nf) = matmul(Sm, matmul(Gsh(:,:,nf), Sminv))
      Enddo
      dbg_n_carry = dbg_n_carry + 1
#endif

   end subroutine Lowrank_carry_similarity

!--------------------------------------------------------------------
!> @brief
!> Submatrix only: fills the workspaces Cw(:,1:k) = G0(:,S) and
!> Mw(1:k,:) = M (E_S^T - G0(S,:)) (via Rw(1:k,:) = E_S^T - G0(S,:), ZGEMM) for
!> flavour nf (computed index nf_eff, k = n_cols(nf_eff) >= 1), so that G = G0 - Cw Mw.
!> GR is only read.
!--------------------------------------------------------------------
   subroutine Submatrix_factors(GR, nf, nf_eff)

      Use Hamiltonian_main, only: Ndim, N_FL
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In) :: nf, nf_eff

      Integer :: k, m, i

      k = n_cols(nf_eff)
      do m = 1, k
         call ZCOPY(Ndim, GR(1,S(m,nf_eff),nf), 1, Cw(1,m), 1)
      enddo
      do i = 1, Ndim
         Rw(1:k,i) = - GR(S(1:k,nf_eff),i,nf)
      enddo
      do m = 1, k
         Rw(m,S(m,nf_eff)) = Rw(m,S(m,nf_eff)) + cone
      enddo
      call ZGEMM('N', 'N', k, Ndim, k, cone, Msub(1,1,nf_eff), N_cols_max, Rw, N_cols_max, &
           &     czero, Mw, N_cols_max)

   end subroutine Submatrix_factors

!--------------------------------------------------------------------
!> @brief
!> Submatrix only: appends site q (not pending) to S of computed flavour nf_eff with a zero row and
!> column of M (exact: G is unchanged; docs/algorithm-notes.md, sections 12.3, 12.4). k is the
!> current number of pending sites on entry and k + 1 on return; n_cols is left to the caller.
!> Terminates on overflow (impossible with the capacity rule, see the module header).
!--------------------------------------------------------------------
   subroutine Submatrix_pad(nf_eff, q, k)

      Implicit none

      Integer, Intent(In)    :: nf_eff, q
      Integer, Intent(InOut) :: k

      if (k >= N_cols_max) then
         Write(error_unit,'(A,I0)') ' Submatrix_pad: overflow of the pending sites, N_cols_max = ', N_cols_max
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      k = k + 1
      S(k,nf_eff)      = q
      slot(q,nf_eff)   = k
      Msub(k,1:k,nf_eff) = czero
      Msub(1:k,k,nf_eff) = czero

   end subroutine Submatrix_pad

#ifdef ALF_LOWRANK_DEBUG
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Debug builds only (-DALF_LOWRANK_DEBUG): checks that the current Green function
!> Lowrank_current_G(GR, nf) equals the shadow Gsh(:,:,nf) for every computed flavour,
!> max|Gcur - Gsh| / max|Gsh| <= TOL_UNIT = 1e-10; terminates otherwise. No-op when nothing is
!> pending (the shadow is then not in use). Records the largest relative error seen.
!> @param [in] GR(Ndim,Ndim,N_FL)
!> Complex
!> \verbatim
!>  G0 (flushed part of G).
!> \endverbatim
!> @param [in] n_op, nt
!> Integer
!> \verbatim
!>  Vertex and time slice, reported on a violation.
!> \endverbatim
!> @param [in] where
!> Character
!> \verbatim
!>  Calling site, reported on a violation.
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Lowrank_debug_check(GR, n_op, nt, where)

      Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map
      Implicit none

      Complex (Kind=Kind(0.d0)), Intent(In) :: GR(Ndim,Ndim,N_FL)
      Integer,                   Intent(In) :: n_op, nt
      Character (len=*),         Intent(In) :: where

      Complex (Kind=Kind(0.d0)), allocatable :: Gcur(:,:)
      Real    (Kind=Kind(0.d0))              :: err, scale
      Integer                                :: nf, nf_eff

      if (n_pend == 0) return
      if (.not. allocated(Gsh)) then
         Write(error_unit,'(A,A)') ' Lowrank_debug_check: shadow G not allocated at ', trim(where)
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
      allocate (Gcur(Ndim,Ndim))
      Do nf_eff = 1,N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         call Lowrank_current_G(GR, nf, Gcur)
         scale = maxval(abs(Gsh(:,:,nf)))
         err   = maxval(abs(Gcur - Gsh(:,:,nf)))
         if (scale > 0.d0) err = err/scale
         dbg_n_checks = dbg_n_checks + 1
         ! .not. (err <= tol) also catches NaN
         if (.not. (err <= DBG_TOL)) then
            Write(error_unit,'(A,A)')               ' Lowrank_debug_check: stale G (low-rank G /= shadow G) at ', &
                 &                                   trim(where)
            Write(error_unit,'(A,I0,A,I0,A,I0,A,I0)') '   n_op = ', n_op, ', nt = ', nt, ', nf = ', nf, &
                 &                                   ', n_pend = ', n_pend
            Write(error_unit,'(A,ES12.4,A,ES12.4)')  '   max|G_lowrank - G_shadow|/max|G_shadow| = ', err, &
                 &                                   ' > TOL_UNIT = ', DBG_TOL
            CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif
         if (err > dbg_max_err) dbg_max_err = err
      Enddo
      deallocate (Gcur)

   end subroutine Lowrank_debug_check

!--------------------------------------------------------------------
!> @brief
!> Debug builds only: (re)writes lowrank_debug_summary.txt in the run directory with the number
!> of shadow-G checks, the largest relative error so far and (submatrix) the numbers of sites added
!> to S by accepts and padded by hooks. Called after every flush, so the file holds the totals of
!> the run when it ends.
!--------------------------------------------------------------------
   subroutine Dbg_check_nd_guard(n_op, nt, where)
!> Debug builds only: the submatrix accuracy guard (spec D15). After an accept or a hook (after its
!> capacity flush), a block whose M-form holds non-diagonal content (dbg_nd_in_M, tracked independently of
!> blk_nondiag) must have been flushed at min(N_delay, K_ND) pending sites, and blk_nondiag must be set.
!> Terminates otherwise.
      Implicit none
      Integer,           Intent(In) :: n_op, nt
      Character (len=*), Intent(In) :: where
      if (lr_mode /= SCHEME_SUBMATRIX .or. .not. dbg_nd_in_M) return
      if (.not. blk_nondiag .or. maxval(n_cols) >= min(lr_N_delay, K_ND)) then
         Write(error_unit,'(A,A)') ' Lowrank_debug: K_ND guard violated (non-diagonal content kept in the M-form) at ', &
              &                    trim(where)
         Write(error_unit,'(A,I0,A,I0,A,I0,A,L1,A,I0)') '   n_op = ', n_op, ', nt = ', nt, ', max pending sites = ', &
              &   maxval(n_cols), ', blk_nondiag = ', blk_nondiag, ', K_ND = ', K_ND
         CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
      endif
   end subroutine Dbg_check_nd_guard

   subroutine Lowrank_debug_summary()
      Implicit none
      Integer :: unit_dbg
      Open(newunit=unit_dbg, file="lowrank_debug_summary.txt", status="replace", action="write")
      Write(unit_dbg,'(A,I0)')    'shadow-G checks         : ', dbg_n_checks
      Write(unit_dbg,'(A,ES12.4)') 'max relative error of G : ', dbg_max_err
      Write(unit_dbg,'(A,ES12.4)') 'max rel. error of Mat(d): ', dbg_max_derr
      Write(unit_dbg,'(A,ES12.4)') 'TOL_UNIT                : ', DBG_TOL
      if (lr_scheme == SCHEME_SUBMATRIX) then
         Write(unit_dbg,'(A,I0)')  'padded sites (accepts)  : ', dbg_n_pad_acc
         Write(unit_dbg,'(A,I0)')  'padded sites (hooks)    : ', dbg_n_pad_hook
      endif
      if (lr_cross) then
         Write(unit_dbg,'(A,I0)')  'cross-slice carries     : ', dbg_n_carry
         Write(unit_dbg,'(A,I0)')  'submatrix -> delayed    : ', dbg_n_conv
      endif
      Close(unit_dbg)
   end subroutine Lowrank_debug_summary
#endif

end module Lowrank_update_mod


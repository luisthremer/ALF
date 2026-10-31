! 36-lowrank-unit: unit-level harness for low-rank (delayed / submatrix) local updates.
!
! Reference engine: the real Upgrade2 in mode "Intermediate" (returns the full ratio
! in Prev_Ratiotot and always accepts), applied to G_ref (accepted proposals) or to a
! copy of G_ref (rejected proposals).
! Candidate engines are dispatched through cand_init / cand_hook / cand_propose /
! cand_accept / cand_flush (select case on the scheme id). Adding a scheme means adding
! one branch to each of these routines plus its state; the driver stays unchanged.
!
! Per vertex the driver mirrors WRAPGRUP / WRAPGRDO (Prog/Wrapgr_mod.F90):
!   up  : pre_wrap -> Op_Wrapup type 1 (old phi) -> hook(1,+1) -> propose -> decide
!         -> Op_Wrapup type 2 -> hook(2,+1)
!   down: pre_wrap -> Op_Wrapdo type 2 -> hook(2,-1) -> propose -> decide
!         -> Op_Wrapdo type 1 (new phi) -> hook(1,-1)
! The hook is called after both wraps of every eligible vertex (kind = N_type of the wrap just done).
! pre_wrap flushes a low-rank candidate before an ineligible vertex; such a vertex then gets no
! hook and the standard Upgrade2 on the candidate's G.
! Eligibility (Lowrank_init): delayed = support max_f Op_V(n,f)%N <= 8; submatrix = diagonal rank-1
! while SUBMATRIX_RANKN = .false. (project 1), the delayed rule once it is .true.; checked per case.
! Case mixed: Mz vertices with interleaved V_SUN bonds (eligible for delayed, not for rank-1 submatrix).
! Case asym: flavour 2's vertex sites are another permutation than flavour 1's (the per-flavour
! pending sets S differ), and every 8th vertex has g = 0 in flavour 2, so Delta = 0 exactly there.
! Case repeat_skip = repeat with a deterministic skip mask: some vertices are wrapped (and hooked) but not
! proposed (no propose/accept, as when the T0/Propose_S0 proposal is rejected before Upgrade), so the
! hook rescale of the pending state stays in use without a following flush.
! Rank-N cases (Mz on-site vertices interleaved unless stated): jz (Predefined_Int_Jz, diagonal N = 2),
! sasb_zz (Predefined_Int_SaSb('z','z'), diagonal N = 4, N_FL = 1 spin-orbital layout), tv_bond (V_SUN ring
! bonds), kondo_j (V_SUN and VJ_SUN between c and f sites), ising1 (Ising_SUN, field type 1), sp_type4
! (type-4 bonds, flip_protocol 2), nd_zero (O = (1 1; 1 1)/2: r = 1 < N_P), asym_rn (flavour 1 V_SUN,
! flavour 2 cycles through nd_zero, V_SUN with g = 0 and a diagonal N = 2 operator), gt_bonds (bonds with
! g_t), diag_e0 (diagonal N = 3 with one E = 0), bonds_skip (tv_bond with the skip mask), bonds_illcond
! (tv_bond's structure at Ndim = 32 with the ill-conditioned similarity S = I + 0.2 R), rep_mixed (rank-1 and bonds revisiting sites within a
! slice), padcap (overlapping bonds only), bigop (diagonal N = 9 vertices: ineligible for both schemes).
! Between slices: flush, compare, then the same dense similarity G <- S G S^{-1} on every G,
! S = I + 0.1 R / sqrt(Ndim) (R complex Gaussian). Cases mz_illcond / bonds_illcond use S = I + 0.1 R /
! S = I + 0.2 R, which drive cond(1+A) to ~1e8: deliberately hard cases with their own tolerance TOL_ILL.
! Cross-slice variant (spec D10, docs/algorithm-notes.md section 12.8; schemes delayed_cs / submat_cs, i.e.
! delayed / submatrix with Lowrank_init(..., cross_slice=.true.)): at a slice end the pending state is NOT
! flushed; the current G (Lowrank_current_G) is compared instead, and Lowrank_carry_similarity(Gc, S, S^{-1})
! carries the pending state through the inter-slice similarity (delayed: X <- S X, Y^T <- Y^T S^{-1}; submatrix:
! first converted to delayed form) before the driver applies the similarity to Gc = G0. Only the last slice
! end flushes. The rows of the variant follow the main table (same tolerances and flags).
!
! nsigma is global and shared by all engines: every engine call sees the OLD field, and
! the driver sets the field after the decision.
!
! The standard path's own drift is measured against an exact recomputation (1+A)^{-1},
! with A tracked explicitly (cf. docs/scripts/check_formulas.py).
!
! Build with -DLOWRANK_UNIT_SELFTEST to perturb one candidate G entry by 1e-8
! (case mz, up, first N_delay, end of slice 1): the harness must then fail.

module lowrank_unit_harness

   use Hamiltonian_main, only: Op_V, nsigma, Ndim, N_FL, N_FL_eff, Calc_Fl_map, N_SUN, &
        &                      reconstruction_needed
   use Operator_mod
   use Fields_mod
   use Predefined_Int
   use upgrade_mod
   use Random_Wrap
   use Lowrank_update_mod, only: Lowrank_init, Lowrank_finalize, Lowrank_propose, Lowrank_accept, &
        &                        Lowrank_wrap_hook, Lowrank_flush, Upgrade_lowrank, Lowrank_eligible, &
        &                        Lowrank_pending, Lowrank_max_pending_cols, Lowrank_current_G, &
        &                        Lowrank_carry_similarity
   implicit none

   integer, parameter :: dp = kind(0.d0)

   ! File to which Lowrank_init appends its scheme line (Temp_<igroup>/info in tempering builds,
   ! the test runs as group 0)
#if defined(MPI) && defined(TEMPERING)
   character(len=*), parameter :: info_path = 'Temp_0/info'
#else
   character(len=*), parameter :: info_path = 'info'
#endif

   ! Tolerances (see CLAUDE.md "Tolerances"): TOL_UNIT for the well-conditioned cases,
   ! TOL_ILL for mz_illcond and bonds_illcond (cond(1+A) ~ 1e8, errors scale ~ cond * eps).
   real(dp), parameter :: TOL_UNIT = 1.0d-10
   real(dp), parameter :: TOL_ILL  = 1.0d-8

   integer, parameter :: N_CASES = 23
   integer, parameter :: CASE_MZ = 1, CASE_DENSITY = 2, CASE_CONTINUOUS = 3, CASE_REPEAT = 4
   integer, parameter :: CASE_MZ_ILLCOND = 5   ! mz with the unnormalised S = I + 0.1 R
   integer, parameter :: CASE_MIXED = 6        ! Mz vertices with interleaved V_SUN bonds
   integer, parameter :: CASE_ASYM = 7         ! flavour-dependent site permutations, Delta = 0 vertices
   integer, parameter :: CASE_REPEAT_SKIP = 8  ! repeat, with unproposed (skipped) vertices
   ! Rank-N cases (see the header)
   integer, parameter :: CASE_JZ = 9, CASE_SASB_ZZ = 10, CASE_TV_BOND = 11, CASE_KONDO_J = 12
   integer, parameter :: CASE_ISING1 = 13, CASE_SP_TYPE4 = 14, CASE_ND_ZERO = 15, CASE_ASYM_RN = 16
   integer, parameter :: CASE_GT_BONDS = 17, CASE_DIAG_E0 = 18, CASE_BONDS_SKIP = 19
   integer, parameter :: CASE_BONDS_ILLCOND = 20, CASE_REP_MIXED = 21, CASE_PADCAP = 22, CASE_BIGOP = 23
   character(len=13), parameter :: case_name(N_CASES) = &
        & [character(len=13) :: 'mz', 'density', 'continuous', 'repeat', 'mz_illcond', 'mixed', 'asym', &
        &  'repeat_skip', 'jz', 'sasb_zz', 'tv_bond', 'kondo_j', 'ising1', 'sp_type4', 'nd_zero', 'asym_rn', &
        &  'gt_bonds', 'diag_e0', 'bonds_skip', 'bonds_illcond', 'rep_mixed', 'padcap', 'bigop']
   integer, parameter :: ASYM_ZERO_PERIOD = 8  ! case asym: vertex n has g = 0 in flavour 2 iff mod(n, 8) == 3
   integer, parameter :: MIXED_BOND_PERIOD = 9 ! case mixed: vertex n is a bond iff mod(n, 9) == 5

   ! Submatrix expectations: .false. = project 1 (eligible iff diagonal rank-1 in every flavour; the
   ! repeated-site flush in Lowrank_propose is asserted); .true. = rank-N submatrix (Task 5: eligibility
   ! as delayed, no repeated-site flush).
   logical, parameter :: SUBMATRIX_RANKN = .true.
   integer, parameter :: N_P_LIM_TEST = 8      ! Lowrank_update_mod's N_P_LIM (support limit)

   ! Candidate schemes
   integer, parameter :: N_SCHEMES = 3
   integer, parameter :: SCHEME_STANDARD = 1, SCHEME_DELAYED = 2, SCHEME_SUBMATRIX = 3
   character(len=10), parameter :: scheme_name(N_SCHEMES) = &
        & [character(len=10) :: 'standard', 'delayed', 'submatrix']
   ! Row labels of the cross-slice variant (delayed, submatrix with cross_slice = .true.)
   character(len=10), parameter :: scheme_name_cs(N_SCHEMES) = &
        & [character(len=10) :: 'standrd_cs', 'delayed_cs', 'submat_cs']

   integer, parameter :: N_NDELAY = 5
   integer, parameter :: N_DELAY_LIST(N_NDELAY) = [1, 2, 7, 32, 1000]

   integer, parameter :: N_SLICES = 3
   real(dp), parameter :: DTAU = 0.1d0, U_HUB = 4.d0
   real(dp), parameter :: ACCEPT_FRACTION = 0.6d0
   real(dp), parameter :: Q_SCALE = 0.3d0     ! Q = I + Q_SCALE * R / sqrt(Ndim)
   real(dp), parameter :: S_SCALE = 0.1d0     ! S = I + S_SCALE * R / sqrt(Ndim)
                                              ! (mz_illcond: S = I + S_SCALE * R)
   ! bonds_illcond: tv_bond's structure at Ndim = 32 with S = I + S_SCALE_BONDS_ILL * R. At Ndim = 16 with the
   ! mz_illcond similarity (0.1) cond(1+A) reached only 7.9e5 (0.2-0.5: <= 2e7); Ndim = 32 with 0.2 gives
   ! cond(1+A) = 6.8e7 (down) / 2.7e8 (up) at the last slice, the regime of mz_illcond (1.9e8).
   real(dp), parameter :: S_SCALE_BONDS_ILL = 0.2d0
   ! Rank-N cases: couplings chosen so that cond(1+A) stays <= ~1e8 (the TOL_UNIT regime of the project-1
   ! cases); stronger real couplings on many overlapping bonds drive cond(1+A) to 1e10-1e16 within 3 slices,
   ! where the double-precision reference itself drifts (measured 2026-10-31: padcap with V = 1: 1.2e8,
   ! rep_mixed with V = 1: 9.6e12, nd_zero with g = 0.5: 1.1e10).
   real(dp), parameter :: V_BOND = 0.25d0     ! V_SUN / VJ_SUN with real coupling g = sqrt(DTAU V)
   real(dp), parameter :: G_ND   = 0.25d0     ! nd_zero operator (1 1; 1 1)/2
   real(dp), parameter :: U_REP  = 2.d0       ! rep_mixed: two Mz vertices per site and slice (U = 4: cond 5.6e8)

   integer :: M_ops                            ! number of vertices per slice

   ! Diagnostics of the last run_case
   integer  :: n_proposals, n_accepted
   real(dp) :: cond_Q_max, cond_S_max
   real(dp) :: cond_1pA_last                   ! max_f cond_2(1+A_f) at the last slice end
   integer  :: n_repeat_flush                  ! proposals at which Lowrank_propose flushed (pending dropped)
   integer  :: n_skipped                       ! vertices wrapped but not proposed (cases repeat_skip, bonds_skip)
   integer  :: n_elig_bad                      ! vertices whose Lowrank_eligible differs from the expectation
   integer  :: n_lowrank_vertices              ! vertices the candidate updates with its low-rank kernel
   integer  :: n_cap_bad                       ! accepts after which the pending rank exceeded max(N_delay, r)
   integer  :: max_pend_cols                   ! largest pending rank seen right after an accept

   ! ---- candidate state -------------------------------------------------
   integer :: cand_scheme = 0, cand_N_delay = 0
   logical :: cand_cross = .false.                  ! cross-slice variant (low-rank schemes only)
   complex(dp), allocatable :: std_Gwork(:,:,:)     ! standard: Upgrade2 on a work copy

contains

   ! =====================================================================
   ! Candidate interface
   ! =====================================================================
   subroutine cand_init(scheme, N_delay, cross_slice)
      integer, intent(in) :: scheme, N_delay
      logical, intent(in) :: cross_slice
      cand_scheme  = scheme
      cand_N_delay = N_delay
      cand_cross   = cross_slice .and. scheme /= SCHEME_STANDARD
      ! Every scheme needs the work copy: low-rank schemes use Upgrade2 at ineligible vertices.
      if (allocated(std_Gwork)) deallocate(std_Gwork)
      allocate(std_Gwork(Ndim, Ndim, N_FL))
      select case (scheme)
      case (SCHEME_STANDARD)
      case (SCHEME_DELAYED)
         ! Works on the driver's Gc (= G0 between flushes); prints the scheme info line.
         call Lowrank_init('delayed', N_delay, cross_slice=cand_cross)
      case (SCHEME_SUBMATRIX)
         call Lowrank_init('submatrix', N_delay, cross_slice=cand_cross)
      case default
         error stop 'cand_init: unknown scheme'
      end select
   end subroutine cand_init

   subroutine cand_finalize()
      if (allocated(std_Gwork)) deallocate(std_Gwork)
      select case (cand_scheme)
      case (SCHEME_DELAYED, SCHEME_SUBMATRIX)
         call Lowrank_finalize()
      end select
   end subroutine cand_finalize

   ! True iff the candidate updates vertex n_op with its low-rank kernel (as WRAPGRUP/WRAPGRDO).
   logical function cand_lowrank(n_op)
      integer, intent(in) :: n_op
      select case (cand_scheme)
      case (SCHEME_DELAYED, SCHEME_SUBMATRIX)
         cand_lowrank = Lowrank_eligible(n_op)
      case default
         cand_lowrank = .false.
      end select
   end function cand_lowrank

   ! Called at the start of every vertex, before any wrap: a low-rank candidate is flushed
   ! before an ineligible vertex (WRAPGRUP/WRAPGRDO: flush, then the standard path).
   subroutine cand_pre_wrap(Gc, n_op)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op
      select case (cand_scheme)
      case (SCHEME_DELAYED, SCHEME_SUBMATRIX)
         if (.not. Lowrank_eligible(n_op)) call Lowrank_flush(Gc)
      end select
   end subroutine cand_pre_wrap

   ! Called after every wrap of vertex n_op: kind = N_type of the wrap just done (1 or 2),
   ! dir = +1 up, -1 down. Unconditional for an eligible vertex (as WRAPGRUP/WRAPGRDO).
   subroutine cand_hook(Gc, n_op, nt, HS_field, kind, dir)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op, nt, kind, dir
      complex(dp), intent(in)    :: HS_field
      if (cand_lowrank(n_op)) call Lowrank_wrap_hook(n_op, nt, HS_field, kind, dir, Gc)
   end subroutine cand_hook

   ! Expected eligibility of vertex n_op for scheme (see the header)
   logical function expected_eligible(scheme, n_op)
      integer, intent(in) :: scheme, n_op
      integer :: nf
      logical :: small, diag_rank1_all
      small = .true.; diag_rank1_all = .true.
      do nf = 1, N_FL
         if (Op_V(n_op,nf)%N > N_P_LIM_TEST) small = .false.
         if (Op_V(n_op,nf)%N_non_zero > 1 .or. .not. Op_V(n_op,nf)%diag) diag_rank1_all = .false.
      enddo
      select case (scheme)
      case (SCHEME_DELAYED)
         expected_eligible = small
      case (SCHEME_SUBMATRIX)
         if (SUBMATRIX_RANKN) then
            expected_eligible = small
         else
            expected_eligible = diag_rank1_all
         endif
      case default
         expected_eligible = .false.
      end select
   end function expected_eligible

   ! Returns the full Upgrade2-style ratio. nsigma%f(n_op,nt) holds the OLD field.
   subroutine cand_propose(Gc, n_op, nt, HS_new, ratio)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op, nt
      complex(dp), intent(in)    :: HS_new
      complex(dp), intent(out)   :: ratio
      complex(dp)        :: phase
      logical            :: toggle
      character(len=64)  :: mode
      integer            :: pend_before
      if (.not. cand_lowrank(n_op)) then
         ! standard scheme, or an ineligible vertex of a low-rank scheme (flushed by cand_pre_wrap)
         std_Gwork = Gc
         ratio = cmplx(1.d0, 0.d0, dp)
         phase = cmplx(1.d0, 0.d0, dp)
         mode  = "Intermediate"
         call Upgrade2(std_Gwork, n_op, nt, phase, HS_new, ratio, 1.d0, 1.d0, toggle, mode)
         return
      endif
      ! delayed or submatrix. A drop of the pending count = a flush inside Lowrank_propose
      ! (project-1 submatrix only: the site is already in the block; never with SUBMATRIX_RANKN).
      pend_before = Lowrank_pending()
      call Lowrank_propose(Gc, n_op, nt, HS_new, ratio)
      if (Lowrank_pending() < pend_before) n_repeat_flush = n_repeat_flush + 1
   end subroutine cand_propose

   ! Called only for accepted proposals, right after cand_propose.
   ! nsigma%f(n_op,nt) holds the OLD field.
   subroutine cand_accept(Gc, n_op, nt, HS_new)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op, nt
      complex(dp), intent(in)    :: HS_new
      integer :: nf, r_last, p_after
      if (.not. cand_lowrank(n_op)) then
         Gc = std_Gwork
      else
         call Lowrank_accept(Gc, n_op, nt, HS_new)   ! sets nsigma%f(n_op,nt); the driver resets it
         ! Capacity rule (D5): with the flush before an accept (p > 0 and p + r > N_delay), the pending
         ! rank right after the append never exceeds max(N_delay, r) (r = this vertex's N_non_zero)
         r_last = 0
         do nf = 1, N_FL
            r_last = max(r_last, Op_V(n_op,nf)%N_non_zero)
         enddo
         p_after = Lowrank_max_pending_cols()
         max_pend_cols = max(max_pend_cols, p_after)
         if (p_after > max(cand_N_delay, r_last)) n_cap_bad = n_cap_bad + 1
      endif
   end subroutine cand_accept

   ! Called at the end of every slice; afterwards Gc must hold the full G.
   subroutine cand_flush(Gc)
      complex(dp), intent(inout) :: Gc(:,:,:)
      select case (cand_scheme)
      case (SCHEME_STANDARD)
         ! no-op
      case (SCHEME_DELAYED, SCHEME_SUBMATRIX)
         call Lowrank_flush(Gc)
      end select
   end subroutine cand_flush

   ! =====================================================================
   ! Linear algebra helpers
   ! =====================================================================
   subroutine invert(A)
      complex(dp), intent(inout) :: A(:,:)
      integer :: n, info, lwork
      integer, allocatable :: ipiv(:)
      complex(dp), allocatable :: work(:)
      n = size(A, 1)
      lwork = 64 * n
      allocate(ipiv(n), work(lwork))
      call ZGETRF(n, n, A, n, ipiv, info)
      if (info /= 0) error stop 'invert: ZGETRF failed'
      call ZGETRI(n, A, n, ipiv, work, lwork, info)
      if (info /= 0) error stop 'invert: ZGETRI failed'
   end subroutine invert

   real(dp) function cond2(A)
      complex(dp), intent(in) :: A(:,:)
      complex(dp), allocatable :: B(:,:), work(:)
      complex(dp) :: dummy(1,1)
      real(dp), allocatable :: s(:), rwork(:)
      integer :: n, info, lwork
      n = size(A, 1)
      lwork = 8 * n
      allocate(B(n,n), work(lwork), s(n), rwork(5*n))
      B = A
      call ZGESVD('N', 'N', n, n, B, n, s, dummy, 1, dummy, 1, work, lwork, rwork, info)
      if (info /= 0) error stop 'cond2: ZGESVD failed'
      cond2 = s(1) / s(n)
   end function cond2

   complex(dp) function crandn()
      real(dp) :: x, y
      x = rang_wrap()
      y = rang_wrap()
      crandn = cmplx(x, y, dp) / sqrt(2.d0)
   end function crandn

   ! Random G0 = (1+A)^{-1}, A = Q diag(e^{x_i}) Q^{-1}, x_i ~ U[-4,4].
   subroutine random_initial_G(G, A)
      complex(dp), intent(out) :: G(:,:), A(:,:)
      complex(dp), allocatable :: Q(:,:), Qinv(:,:)
      integer :: i, j, n
      real(dp) :: x
      n = size(G, 1)
      allocate(Q(n,n), Qinv(n,n))
      do j = 1, n
         do i = 1, n
            Q(i,j) = Q_SCALE * crandn() / sqrt(real(n, dp))
         enddo
         Q(j,j) = Q(j,j) + 1.d0
      enddo
      cond_Q_max = max(cond_Q_max, cond2(Q))
      Qinv = Q
      call invert(Qinv)
      do j = 1, n
         x = 8.d0 * ranf_wrap() - 4.d0
         Q(:,j) = Q(:,j) * exp(x)
      enddo
      A = matmul(Q, Qinv)
      G = A
      do i = 1, n
         G(i,i) = G(i,i) + 1.d0
      enddo
      call invert(G)
   end subroutine random_initial_G

   ! Exact G = (1+A)^{-1}
   subroutine exact_G(A, G)
      complex(dp), intent(in)  :: A(:,:)
      complex(dp), intent(out) :: G(:,:)
      integer :: i
      G = A
      do i = 1, size(A, 1)
         G(i,i) = G(i,i) + 1.d0
      enddo
      call invert(G)
   end subroutine exact_G

   ! =====================================================================
   ! Case setup
   ! =====================================================================
   subroutine clear_case()
      integer :: n, nf
      if (allocated(Op_V)) then
         do nf = 1, size(Op_V, 2)
            do n = 1, size(Op_V, 1)
               call Op_clear(Op_V(n,nf), Op_V(n,nf)%N)
            enddo
         enddo
         deallocate(Op_V)
      endif
      if (allocated(nsigma%f)) call nsigma%clear()
      if (allocated(Calc_Fl_map)) deallocate(Calc_Fl_map)
   end subroutine clear_case

   subroutine random_permutation(perm)
      integer, intent(out) :: perm(:)
      integer :: i, j, tmp
      perm = [(i, i = 1, size(perm))]
      do i = size(perm), 2, -1
         j = 1 + int(real(i, dp) * ranf_wrap())
         if (j > i) j = i
         tmp = perm(i); perm(i) = perm(j); perm(j) = tmp
      enddo
   end subroutine random_permutation

   ! Case mixed: vertex n is a V_SUN bond (rank 2, not diagonal)
   logical function mixed_bond(case_id, n)
      integer, intent(in) :: case_id, n
      mixed_bond = (case_id == CASE_MIXED .and. mod(n, MIXED_BOND_PERIOD) == 5)
   end function mixed_bond

   subroutine setup_case(case_id)
      integer, intent(in) :: case_id
      integer, allocatable :: sites(:), base(:), pairs(:), sites2(:)
      type(Operator) :: Op_tmp
      integer :: n, nf, nt, k, i_mz, i_bond
      real(dp), parameter :: disc(4) = [-2.d0, -1.d0, 1.d0, 2.d0]

      call clear_case()
      if (case_id >= CASE_JZ) then
         call setup_rankn_vertices(case_id)
      else
      select case (case_id)
      case (CASE_MZ, CASE_MZ_ILLCOND, CASE_ASYM)
         Ndim = 64; N_FL = 2; N_SUN = 1; M_ops = Ndim
      case (CASE_DENSITY)
         Ndim = 64; N_FL = 1; N_SUN = 2; M_ops = Ndim
      case (CASE_CONTINUOUS)
         Ndim = 48; N_FL = 2; N_SUN = 1; M_ops = Ndim
      case (CASE_REPEAT, CASE_REPEAT_SKIP)
         Ndim = 32; N_FL = 2; N_SUN = 1; M_ops = 2 * Ndim
      case (CASE_MIXED)
         ! Ndim Mz vertices plus one bond after every 8 of them (n = 5, 14, 23, ...)
         Ndim = 48; N_FL = 2; N_SUN = 1; M_ops = Ndim + Ndim / (MIXED_BOND_PERIOD - 1)
      case default
         error stop 'setup_case: unknown case'
      end select
      N_FL_eff = N_FL
      allocate(Calc_Fl_map(N_FL))
      Calc_Fl_map = [(nf, nf = 1, N_FL)]
      reconstruction_needed = .false.

      ! Vertex -> site: a fixed random permutation (repeat: of the doubled site list)
      allocate(sites(M_ops))
      if (case_id == CASE_REPEAT .or. case_id == CASE_REPEAT_SKIP) then
         allocate(base(M_ops))
         call random_permutation(base)
         do n = 1, M_ops
            sites(n) = mod(base(n) - 1, Ndim) + 1
         enddo
      elseif (case_id == CASE_MIXED) then
         ! sites(1:Ndim): the Mz vertices; pairs: the bond endpoints
         call random_permutation(sites(1:Ndim))
         allocate(pairs(Ndim))
         call random_permutation(pairs)
      else
         call random_permutation(sites)
      endif
      if (case_id == CASE_ASYM) then
         ! flavour 2: another permutation (re-drawn until it differs from flavour 1's at every vertex)
         allocate(sites2(M_ops))
         do
            call random_permutation(sites2)
            if (all(sites2 /= sites)) exit
         enddo
      endif

      allocate(Op_V(M_ops, N_FL))
      i_mz = 0; i_bond = 0
      do n = 1, M_ops
         if (mixed_bond(case_id, n)) then
            i_bond = i_bond + 1
            do nf = 1, N_FL
               call Predefined_Int_V_SUN(Op_V(n,nf), pairs(2*i_bond-1), pairs(2*i_bond), N_SUN, DTAU, 1.d0)
               if (Op_V(n,nf)%diag .or. Op_V(n,nf)%N_non_zero /= 2) &
                    & error stop 'setup_case: expected a non-diagonal rank-2 bond'
            enddo
            cycle
         endif
         i_mz = i_mz + 1
         select case (case_id)
         case (CASE_MZ, CASE_REPEAT, CASE_MZ_ILLCOND, CASE_REPEAT_SKIP)
            call Predefined_Int_U_MZ(Op_V(n,1), Op_V(n,2), sites(n), DTAU, U_HUB)
         case (CASE_ASYM)
            ! flavour 1 from site sites(n), flavour 2 from site sites2(n)
            call Predefined_Int_U_MZ(Op_V(n,1), Op_tmp, sites(n), DTAU, U_HUB)
            call Op_clear(Op_tmp, Op_tmp%N)
            ! every 8th vertex: U = 0 gives g = 0 (E_exp is set from g in Op_set), so Delta = 0 exactly in flavour 2
            if (mod(n, ASYM_ZERO_PERIOD) == 3) then
               call Predefined_Int_U_MZ(Op_tmp, Op_V(n,2), sites2(n), DTAU, 0.d0)
            else
               call Predefined_Int_U_MZ(Op_tmp, Op_V(n,2), sites2(n), DTAU, U_HUB)
            endif
            call Op_clear(Op_tmp, Op_tmp%N)
         case (CASE_MIXED)
            call Predefined_Int_U_MZ(Op_V(n,1), Op_V(n,2), sites(i_mz), DTAU, U_HUB)
         case (CASE_DENSITY)
            ! g = sqrt(-DTAU*U/N_SUN) is purely imaginary for U > 0
            call Predefined_Int_U_SUN(Op_V(n,1), sites(n), N_SUN, DTAU, U_HUB)
         case (CASE_CONTINUOUS)
            call Predefined_Int_U_MZ_continuous_HS(Op_V(n,1), Op_V(n,2), sites(n), DTAU, U_HUB)
         end select
         do nf = 1, N_FL
            if (.not. Op_V(n,nf)%diag .or. Op_V(n,nf)%N /= 1 .or. Op_V(n,nf)%N_non_zero /= 1) &
                 & error stop 'setup_case: expected diagonal rank-1 vertices'
         enddo
      enddo
      if (case_id == CASE_MIXED .and. (i_mz /= Ndim .or. i_bond /= M_ops - Ndim)) &
           & error stop 'setup_case: mixed vertex count'
      if (case_id == CASE_DENSITY) then
         if (abs(aimag(Op_V(1,1)%g)) < 1.d-12) error stop 'setup_case: density g is not complex'
      endif
      endif   ! project-1 cases

      ! Fields per operator type: 1 (Ising) +-1; 2 {+-1, +-2}; 3 continuous real; 4 discrete real part
      ! {+-1, +-2} plus a continuous imaginary part (Spin_Peierls), flip_protocol from the operator (2).
      ! (Op_Wrapup uses E_Exp(nint(real(HS))) for types 1 and 2.)
      call nsigma%make(M_ops, N_SLICES)
      do n = 1, M_ops
         nsigma%t(n) = Op_V(n,1)%type
         nsigma%flip_protocol(n) = Op_V(n,1)%flip_protocol
      enddo
      do nt = 1, N_SLICES
         do n = 1, M_ops
            select case (nsigma%t(n))
            case (1)
               if (ranf_wrap() < 0.5d0) then
                  nsigma%f(n,nt) = cmplx(-1.d0, 0.d0, dp)
               else
                  nsigma%f(n,nt) = cmplx( 1.d0, 0.d0, dp)
               endif
            case (2)
               k = 1 + int(4.d0 * ranf_wrap())
               if (k > 4) k = 4
               nsigma%f(n,nt) = cmplx(disc(k), 0.d0, dp)
            case (4)
               k = 1 + int(4.d0 * ranf_wrap())
               if (k > 4) k = 4
               ! |Im| <= 0.4 initially, one flip (Amplitude 1) moves it by <= 0.5: 1 + Im > 0 (real phi)
               nsigma%f(n,nt) = cmplx(disc(k), 0.8d0 * ranf_wrap() - 0.4d0, dp)
            case default
               nsigma%f(n,nt) = cmplx(3.d0 * ranf_wrap() - 1.5d0, 0.d0, dp)
            end select
         enddo
      enddo
   end subroutine setup_case

   ! Hand-built operator on sites P with matrix O, coupling g, field type optype (type 4: flip_protocol 2),
   ! optional time-dependent coupling g_t (allocated before Op_set, as the Hamiltonians do).
   subroutine make_op(Op, P, O, g, optype, g_t)
      type(Operator), intent(inout)        :: Op
      integer,        intent(in)           :: P(:), optype
      complex(dp),    intent(in)           :: O(:,:), g
      complex(dp),    intent(in), optional :: g_t(:)
      call Op_make(Op, size(P))
      Op%P     = P
      Op%O     = O
      Op%g     = g
      Op%alpha = cmplx(0.d0, 0.d0, dp)
      Op%type  = optype
      if (optype == 4) Op%flip_protocol = 2
      if (present(g_t)) then
         allocate(Op%g_t(size(g_t)))
         Op%g_t = g_t
      endif
      call Op_set(Op)
   end subroutine make_op

   ! Rank-N cases: sizes, Calc_Fl_map and Op_V (fields are set by setup_case). Unless stated otherwise
   ! the vertex list is, per site s = 1..Ndim: an Mz vertex on site perm(s), then the case's bond on
   ! (s, s+1) (ring).
   subroutine setup_rankn_vertices(case_id)
      integer, intent(in) :: case_id
      integer, allocatable :: perm(:)
      integer     :: n, nf, s, s1, s2, Ns, i
      complex(dp) :: Ooff(2,2), Ond(2,2), Odiag2(2,2), O3(3,3), O9(9,9), gt(N_SLICES)
      complex(dp), parameter :: c0 = (0.d0, 0.d0)

      select case (case_id)
      case (CASE_SASB_ZZ)
         ! spin-orbital layout: up orbitals 1..Ns, down orbitals Ns+1..2Ns; N_FL = 1
         Ndim = 16; N_FL = 1; N_SUN = 1; M_ops = 3 * (Ndim / 2)
      case (CASE_KONDO_J)
         ! c sites 1..Ns, f sites Ns+1..2Ns
         Ndim = 16; N_FL = 2; N_SUN = 1; M_ops = 3 * (Ndim / 2)
      case (CASE_REP_MIXED)
         Ndim = 16; N_FL = 2; N_SUN = 1; M_ops = 4 * Ndim
      case (CASE_PADCAP)
         Ndim = 16; N_FL = 2; N_SUN = 1; M_ops = 2 * Ndim
      case (CASE_BIGOP)
         Ndim = 12; N_FL = 2; N_SUN = 1; M_ops = Ndim + 3
      case (CASE_BONDS_ILLCOND)
         ! larger than tv_bond: the ill-conditioned similarity needs Ndim >= 32 (see S_SCALE_BONDS_ILL)
         Ndim = 32; N_FL = 2; N_SUN = 1; M_ops = 2 * Ndim
      case default
         Ndim = 16; N_FL = 2; N_SUN = 1; M_ops = 2 * Ndim
      end select
      N_FL_eff = N_FL
      allocate(Calc_Fl_map(N_FL))
      Calc_Fl_map = [(nf, nf = 1, N_FL)]
      reconstruction_needed = .false.

      allocate(perm(Ndim))
      call random_permutation(perm)
      allocate(Op_V(M_ops, N_FL))

      Ooff = c0; Ooff(1,2) = 1.d0; Ooff(2,1) = 1.d0              ! V_SUN-like bond
      Ond  = cmplx(0.5d0, 0.d0, dp)                               ! (1 1; 1 1)/2: eigenvalues 1, 0
      Odiag2 = c0; Odiag2(1,1) = 1.d0; Odiag2(2,2) = -1.d0       ! diagonal N = 2
      O3 = c0; O3(1,1) = 1.d0; O3(3,3) = -1.d0                    ! diagonal N = 3 with E = 0 in the middle
      O9 = c0
      do i = 1, 9
         O9(i,i) = real(1 - 2 * mod(i, 2), dp) * (1.d0 + 0.1d0 * i)
      enddo

      select case (case_id)
      case (CASE_SASB_ZZ)
         Ns = Ndim / 2
         do s = 1, Ns
            s2 = mod(s, Ns) + 1
            n  = 3 * (s - 1)
            call Predefined_Int_U_SUN(Op_V(n+1,1), s,      N_SUN, DTAU, U_HUB)
            call Predefined_Int_U_SUN(Op_V(n+2,1), s + Ns, N_SUN, DTAU, U_HUB)
            call Predefined_Int_SaSb(Op_V(n+3,1), 'z', 'z', s, s2, s + Ns, s2 + Ns, 1.d0, DTAU, 0.d0)
            if (.not. Op_V(n+3,1)%diag .or. Op_V(n+3,1)%N_non_zero /= 4) error stop 'setup: SaSb zz not diagonal rank 4'
         enddo
      case (CASE_KONDO_J)
         Ns = Ndim / 2
         do s = 1, Ns
            n = 3 * (s - 1)
            call Predefined_Int_U_MZ(Op_V(n+1,1), Op_V(n+1,2), s + Ns, DTAU, U_HUB)
            do nf = 1, N_FL
               call Predefined_Int_V_SUN (Op_V(n+2,nf), s, s + Ns, N_SUN, DTAU, V_BOND)
               call Predefined_Int_VJ_SUN(Op_V(n+3,nf), s, s + Ns, N_SUN, DTAU, V_BOND)
            enddo
         enddo
      case (CASE_REP_MIXED)
         ! per s: Mz(s), V_SUN(s, s+1), Jz(s, s+2), Mz(s+1); the next group starts with Mz(s+1) again
         do s = 1, Ndim
            s1 = mod(s, Ndim) + 1
            s2 = mod(s + 1, Ndim) + 1
            n  = 4 * (s - 1)
            call Predefined_Int_U_MZ(Op_V(n+1,1), Op_V(n+1,2), s, DTAU, U_REP)
            do nf = 1, N_FL
               call Predefined_Int_V_SUN(Op_V(n+2,nf), s, s1, N_SUN, DTAU, V_BOND)
            enddo
            call Predefined_Int_Jz(Op_V(n+3,1), Op_V(n+3,2), s, s2, DTAU, 1.d0)
            call Predefined_Int_U_MZ(Op_V(n+4,1), Op_V(n+4,2), s1, DTAU, U_REP)
         enddo
      case (CASE_PADCAP)
         ! bonds only: (s, s+1) for all s, then (s, s+3): every bond overlaps several earlier ones
         do s = 1, Ndim
            s1 = mod(s, Ndim) + 1
            s2 = mod(s + 2, Ndim) + 1
            do nf = 1, N_FL
               call Predefined_Int_V_SUN(Op_V(s,nf),        s, s1, N_SUN, DTAU, V_BOND)
               call Predefined_Int_V_SUN(Op_V(Ndim + s,nf), s, s2, N_SUN, DTAU, V_BOND)
            enddo
         enddo
      case (CASE_BIGOP)
         ! Mz on perm(1:4), bigop on sites 1..9, Mz on perm(5:8), V_SUN(2,7), Mz on perm(9:12), bigop on 4..12
         n = 0
         do i = 1, Ndim
            n = n + 1
            call Predefined_Int_U_MZ(Op_V(n,1), Op_V(n,2), perm(i), DTAU, U_HUB)
            if (i == 4) then
               n = n + 1
               do nf = 1, N_FL
                  call make_op(Op_V(n,nf), [(s, s = 1, 9)], O9, cmplx(0.2d0 * (3 - 2 * nf), 0.d0, dp), 2)
               enddo
            elseif (i == 8) then
               n = n + 1
               do nf = 1, N_FL
                  call Predefined_Int_V_SUN(Op_V(n,nf), 2, 7, N_SUN, DTAU, V_BOND)
               enddo
            elseif (i == Ndim) then
               n = n + 1
               do nf = 1, N_FL
                  call make_op(Op_V(n,nf), [(s, s = 4, 12)], O9, cmplx(0.2d0, 0.d0, dp), 2)
               enddo
            endif
         enddo
         if (n /= M_ops) error stop 'setup: bigop vertex count'
      case default
         do s = 1, Ndim
            s1 = mod(s, Ndim) + 1
            s2 = mod(s + 1, Ndim) + 1
            n  = 2 * s
            call Predefined_Int_U_MZ(Op_V(n-1,1), Op_V(n-1,2), perm(s), DTAU, U_HUB)
            select case (case_id)
            case (CASE_JZ)
               call Predefined_Int_Jz(Op_V(n,1), Op_V(n,2), s, s1, DTAU, 1.d0)
            case (CASE_TV_BOND, CASE_BONDS_SKIP, CASE_BONDS_ILLCOND)
               ! V = -1: imaginary coupling g (as the tV model with ham_V = -1)
               do nf = 1, N_FL
                  call Predefined_Int_V_SUN(Op_V(n,nf), s, s1, N_SUN, DTAU, -1.d0)
               enddo
            case (CASE_ISING1)
               do nf = 1, N_FL
                  call Predefined_Int_Ising_SUN(Op_V(n,nf), s, s1, DTAU, 2.d0)
               enddo
            case (CASE_SP_TYPE4)
               do nf = 1, N_FL
                  call make_op(Op_V(n,nf), [s, s1], Ooff, cmplx(0.5d0 * sqrt(DTAU), 0.d0, dp), 4)
               enddo
            case (CASE_ND_ZERO)
               call make_op(Op_V(n,1), [s, s1], Ond, cmplx( G_ND, 0.d0, dp), 2)
               call make_op(Op_V(n,2), [s, s1], Ond, cmplx(-G_ND, 0.d0, dp), 2)
               do nf = 1, N_FL
                  if (Op_V(n,nf)%diag .or. Op_V(n,nf)%N_non_zero /= 1) error stop 'setup: nd_zero not rank 1'
               enddo
            case (CASE_ASYM_RN)
               ! flavour 1: V_SUN (non-diagonal, r = 2); flavour 2 cycles through nd_zero on the swapped
               ! pair (r = 1), V_SUN with g = 0 (Delta = 0: no columns) and a diagonal N = 2 operator (r = 2)
               call Predefined_Int_V_SUN(Op_V(n,1), s, s1, N_SUN, DTAU, V_BOND)
               select case (mod(s, 3))
               case (0)
                  call make_op(Op_V(n,2), [s1, s], Ond, cmplx(G_ND, 0.d0, dp), 2)
                  if (Op_V(n,2)%diag .or. Op_V(n,2)%N_non_zero /= 1) error stop 'setup: asym_rn nd_zero'
               case (1)
                  call Predefined_Int_V_SUN(Op_V(n,2), s, s1, N_SUN, DTAU, 0.d0)
                  if (abs(Op_V(n,2)%g) /= 0.d0) error stop 'setup: asym_rn g /= 0'
               case default
                  call make_op(Op_V(n,2), [s, s1], Odiag2, cmplx(0.2d0, 0.d0, dp), 2)
                  if (.not. Op_V(n,2)%diag .or. Op_V(n,2)%N_non_zero /= 2) error stop 'setup: asym_rn diag'
               end select
            case (CASE_GT_BONDS)
               do nf = 1, N_FL
                  gt = [(cmplx((3 - 2 * nf) * sqrt(DTAU * V_BOND) * (0.6d0 + 0.3d0 * i), 0.d0, dp), i = 1, N_SLICES)]
                  call make_op(Op_V(n,nf), [s, s1], Ooff, gt(1), 2, gt)
                  if (.not. Op_V(n,nf)%get_g_t_alloc()) error stop 'setup: g_t not allocated'
               enddo
            case (CASE_DIAG_E0)
               call make_op(Op_V(n,1), [s, s1, s2], O3, cmplx( 0.4d0, 0.d0, dp), 2)
               call make_op(Op_V(n,2), [s, s1, s2], O3, cmplx(-0.4d0, 0.d0, dp), 2)
               do nf = 1, N_FL
                  if (.not. Op_V(n,nf)%diag .or. Op_V(n,nf)%N_non_zero /= 3) error stop 'setup: diag_e0 N_non_zero'
               enddo
            case default
               error stop 'setup_rankn_vertices: unknown case'
            end select
         enddo
      end select
   end subroutine setup_rankn_vertices

   ! =====================================================================
   ! Exact tracker in the site frame. Diagonal rank-1 vertices (vertex frame = site frame):
   ! scalar factors on row/column p. Other vertices (multi-site, non-diagonal): the N x N block
   ! e^{z O} = U e^{z E} U^dagger (Op_exp) on rows/columns P, z with g or g_t(nt).
   ! =====================================================================
   logical function diag_rank1(Op)
      type(Operator), intent(in) :: Op
      diag_rank1 = Op%diag .and. Op%N == 1
   end function diag_rank1

   ! A <- e^{dir g phi O} A e^{-dir g phi O}; diagonal rank-1: T = diag(t at p), t = exp(dir * g * phi * E(1))
   subroutine exact_wrap(A, n_op, nt, dir)
      complex(dp), intent(inout) :: A(:,:,:)
      integer, intent(in) :: n_op, nt, dir
      integer :: nf, p, nb
      complex(dp) :: t, phi, g_loc
      complex(dp), allocatable :: Eb(:,:)
      phi = nsigma%Phi(n_op, nt)
      do nf = 1, N_FL
         g_loc = Op_V(n_op,nf)%g
         if (Op_V(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
         if (diag_rank1(Op_V(n_op,nf))) then
            p = Op_V(n_op,nf)%P(1)
            t = exp(real(dir, dp) * g_loc * phi * Op_V(n_op,nf)%E(1))
            A(p,:,nf) = A(p,:,nf) * t
            A(:,p,nf) = A(:,p,nf) / t
         else
            nb = Op_V(n_op,nf)%N
            allocate(Eb(nb,nb))
            call Op_exp(real(dir, dp) * g_loc * phi, Op_V(n_op,nf), Eb)
            A(Op_V(n_op,nf)%P,:,nf) = matmul(Eb, A(Op_V(n_op,nf)%P,:,nf))
            call Op_exp(-real(dir, dp) * g_loc * phi, Op_V(n_op,nf), Eb)
            A(:,Op_V(n_op,nf)%P,nf) = matmul(A(:,Op_V(n_op,nf)%P,nf), Eb)
            deallocate(Eb)
         endif
      enddo
   end subroutine exact_wrap

   ! Accepted flip: A <- e^{g (phi_new - phi_old) O} A; diagonal rank-1: (1 + Delta e_p e_p^T) A
   ! (g = g_t(nt) if allocated)
   subroutine exact_flip(A, n_op, nt, phi_old, phi_new)
      complex(dp), intent(inout) :: A(:,:,:)
      integer, intent(in) :: n_op, nt
      complex(dp), intent(in) :: phi_old, phi_new
      integer :: nf, p, nb
      complex(dp) :: delta, g_loc
      complex(dp), allocatable :: Eb(:,:)
      do nf = 1, N_FL
         g_loc = Op_V(n_op,nf)%g
         if (Op_V(n_op,nf)%get_g_t_alloc()) g_loc = Op_V(n_op,nf)%g_t(nt)
         if (diag_rank1(Op_V(n_op,nf))) then
            p = Op_V(n_op,nf)%P(1)
            delta = exp(g_loc * (phi_new - phi_old) * Op_V(n_op,nf)%E(1)) - 1.d0
            A(p,:,nf) = A(p,:,nf) * (1.d0 + delta)
         else
            nb = Op_V(n_op,nf)%N
            allocate(Eb(nb,nb))
            call Op_exp(g_loc * (phi_new - phi_old), Op_V(n_op,nf), Eb)
            A(Op_V(n_op,nf)%P,:,nf) = matmul(Eb, A(Op_V(n_op,nf)%P,:,nf))
            deallocate(Eb)
         endif
      enddo
   end subroutine exact_flip

   ! =====================================================================
   ! Driver
   ! =====================================================================
   ! Deterministic skip mask (cases repeat_skip, bonds_skip): vertex n at slice nt is wrapped but not proposed
   logical function skip_vertex(case_id, n, nt)
      integer, intent(in) :: case_id, n, nt
      skip_vertex = ((case_id == CASE_REPEAT_SKIP .or. case_id == CASE_BONDS_SKIP) .and. mod(n + 2 * nt, 3) == 0)
   end function skip_vertex

   ! Ill-conditioned inter-slice similarity S = I + 0.1 R (tolerance TOL_ILL)
   logical function illcond_case(case_id)
      integer, intent(in) :: case_id
      illcond_case = (case_id == CASE_MZ_ILLCOND .or. case_id == CASE_BONDS_ILLCOND)
   end function illcond_case

   subroutine run_case(case_id, scheme, N_delay, direction, max_ratio_err, max_G_err, std_drift, cross_slice)
      integer,          intent(in)  :: case_id, scheme, N_delay
      character(len=*), intent(in)  :: direction
      real(dp),         intent(out) :: max_ratio_err, max_G_err
      real(dp),         intent(out), optional :: std_drift
      logical,          intent(in),  optional :: cross_slice   ! cross-slice variant (see the header)

      complex(dp), allocatable :: Gref(:,:,:), Gc(:,:,:), Gtmp(:,:,:), Aex(:,:,:), Gex(:,:)
      complex(dp), allocatable :: S(:,:), Sinv(:,:)
      integer,     allocatable :: rng_state(:)
      integer :: islice, nt, iv, n, nf, dir, i, j, k_seed
      complex(dp) :: hs_old, hs_new, phi_old, phi_new, R_ref, R_c, phase
      logical :: accept, toggle
      character(len=64) :: mode
      real(dp) :: drift, gmax
      logical  :: cross

      cross = .false.
      if (present(cross_slice)) cross = cross_slice
      if (direction == 'up') then
         dir = +1
      elseif (direction == 'down') then
         dir = -1
      else
         error stop 'run_case: direction must be up or down'
      endif

      ! The configuration depends on the case only (mz_illcond reuses mz's configuration;
      ! bonds_skip reuses tv_bond's; bonds_illcond uses tv_bond's seed at its own Ndim).
      if (case_id == CASE_MZ_ILLCOND) then
         call Ranset([1000 + 17 * CASE_MZ])
      elseif (case_id == CASE_BONDS_SKIP .or. case_id == CASE_BONDS_ILLCOND) then
         call Ranset([1000 + 17 * CASE_TV_BOND])
      else
         call Ranset([1000 + 17 * case_id])
      endif
      call setup_case(case_id)
      n_proposals = 0; n_accepted = 0; n_repeat_flush = 0; n_skipped = 0
      cond_Q_max = 0.d0; cond_S_max = 0.d0; cond_1pA_last = 0.d0

      allocate(Gref(Ndim,Ndim,N_FL), Gc(Ndim,Ndim,N_FL), Gtmp(Ndim,Ndim,N_FL), Aex(Ndim,Ndim,N_FL))
      allocate(Gex(Ndim,Ndim), S(Ndim,Ndim), Sinv(Ndim,Ndim))
      do nf = 1, N_FL
         call random_initial_G(Gref(:,:,nf), Aex(:,:,nf))
      enddo
      Gc = Gref
      call Get_seed_Len(k_seed)
      allocate(rng_state(k_seed))

      call cand_init(scheme, N_delay, cross)
      ! Eligibility as expected for this scheme (standard: no vertex is low-rank)
      n_elig_bad = 0; n_lowrank_vertices = 0; n_cap_bad = 0; max_pend_cols = 0
      do n = 1, M_ops
         if (cand_lowrank(n)) n_lowrank_vertices = n_lowrank_vertices + 1
         if (cand_lowrank(n) .neqv. expected_eligible(scheme, n)) n_elig_bad = n_elig_bad + 1
      enddo
      mode = "Intermediate"
      max_ratio_err = 0.d0; max_G_err = 0.d0; drift = 0.d0

      do islice = 1, N_SLICES
         if (dir == 1) then
            nt = islice
         else
            nt = N_SLICES + 1 - islice
         endif
         do iv = 1, M_ops
            if (dir == 1) then
               n = iv
            else
               n = M_ops + 1 - iv
            endif

            hs_old = nsigma%f(n,nt)
            call cand_pre_wrap(Gc, n)
            if (dir == 1) then
               do nf = 1, N_FL
                  call Op_Wrapup(Gref(:,:,nf), Op_V(n,nf), hs_old, Ndim, 1, nt)
                  call Op_Wrapup(Gc(:,:,nf),   Op_V(n,nf), hs_old, Ndim, 1, nt)
               enddo
               call exact_wrap(Aex, n, nt, +1)
               call cand_hook(Gc, n, nt, hs_old, 1, +1)
            else
               do nf = 1, N_FL
                  call Op_Wrapdo(Gref(:,:,nf), Op_V(n,nf), hs_old, Ndim, 2, nt)
                  call Op_Wrapdo(Gc(:,:,nf),   Op_V(n,nf), hs_old, Ndim, 2, nt)
               enddo
               call cand_hook(Gc, n, nt, hs_old, 2, -1)
            endif

            ! ---- skipped vertex: no propose, no accept (nothing between the two wraps) ----
            if (skip_vertex(case_id, n, nt)) then
               n_skipped = n_skipped + 1
            else
            ! ---- proposal (drawn once, shared) and fixed decision ----
            hs_new = nsigma%flip(n,nt)
            accept = ranf_wrap() < ACCEPT_FRACTION
            n_proposals = n_proposals + 1
            if (accept) n_accepted = n_accepted + 1
            ! Engines may consume random numbers (Upgrade2 does); keep the proposal
            ! stream independent of the engines.
            call random_seed(get=rng_state)

            ! ---- reference: Upgrade2 Intermediate on G_ref or a copy ----
            nsigma%f(n,nt) = hs_old
            R_ref = cmplx(1.d0, 0.d0, dp)
            phase = cmplx(1.d0, 0.d0, dp)
            if (accept) then
               call Upgrade2(Gref, n, nt, phase, hs_new, R_ref, 1.d0, 1.d0, toggle, mode)
            else
               Gtmp = Gref
               call Upgrade2(Gtmp, n, nt, phase, hs_new, R_ref, 1.d0, 1.d0, toggle, mode)
            endif
            if (.not. toggle) error stop 'run_case: Upgrade2 Intermediate did not accept'

            ! ---- candidate ----
            nsigma%f(n,nt) = hs_old
            call cand_propose(Gc, n, nt, hs_new, R_c)
            nsigma%f(n,nt) = hs_old
            if (accept) then
               call cand_accept(Gc, n, nt, hs_new)
               nsigma%f(n,nt) = hs_old
            endif
            max_ratio_err = max(max_ratio_err, abs(R_c - R_ref) / abs(R_ref))

            ! ---- decision ----
            if (accept) then
               phi_old = nsigma%Phi(n,nt)
               nsigma%f(n,nt) = hs_new
               phi_new = nsigma%Phi(n,nt)
               call exact_flip(Aex, n, nt, phi_old, phi_new)
            else
               nsigma%f(n,nt) = hs_old
            endif
            call random_seed(put=rng_state)
            endif   ! not skipped

            if (dir == 1) then
               do nf = 1, N_FL
                  call Op_Wrapup(Gref(:,:,nf), Op_V(n,nf), hs_old, Ndim, 2, nt)
                  call Op_Wrapup(Gc(:,:,nf),   Op_V(n,nf), hs_old, Ndim, 2, nt)
               enddo
               call cand_hook(Gc, n, nt, hs_old, 2, +1)
            else
               hs_new = nsigma%f(n,nt)       ! post-decision field
               do nf = 1, N_FL
                  call Op_Wrapdo(Gref(:,:,nf), Op_V(n,nf), hs_new, Ndim, 1, nt)
                  call Op_Wrapdo(Gc(:,:,nf),   Op_V(n,nf), hs_new, Ndim, 1, nt)
               enddo
               call exact_wrap(Aex, n, nt, -1)
               call cand_hook(Gc, n, nt, hs_new, 1, -1)
            endif
         enddo

         ! ---- slice end: flush, compare (cross-slice variant: compare the current G, flush only at the end) ----
         gmax = maxval(abs(Gref))
         if (cand_cross .and. islice < N_SLICES) then
            do nf = 1, N_FL
               call Lowrank_current_G(Gc, nf, Gtmp(:,:,nf))
            enddo
            max_G_err = max(max_G_err, maxval(abs(Gtmp - Gref)) / gmax)
         else
            call cand_flush(Gc)
#ifdef LOWRANK_UNIT_SELFTEST
            if (case_id == CASE_MZ .and. dir == 1 .and. N_delay == N_DELAY_LIST(1) .and. islice == 1) &
                 & Gc(1,1,1) = Gc(1,1,1) + 1.0d-8
#endif
            max_G_err = max(max_G_err, maxval(abs(Gc - Gref)) / gmax)
         endif
         do nf = 1, N_FL
            call exact_G(Aex(:,:,nf), Gex)
            drift = max(drift, maxval(abs(Gref(:,:,nf) - Gex)) / gmax)
            if (islice == N_SLICES) then
               Gtmp(:,:,1) = Aex(:,:,nf)
               do i = 1, Ndim
                  Gtmp(i,i,1) = Gtmp(i,i,1) + 1.d0
               enddo
               cond_1pA_last = max(cond_1pA_last, cond2(Gtmp(:,:,1)))
            endif
         enddo

         ! ---- between slices: the same dense similarity on every G (and on A) ----
         if (islice < N_SLICES) then
            do j = 1, Ndim
               do i = 1, Ndim
                  if (case_id == CASE_BONDS_ILLCOND) then
                     S(i,j) = S_SCALE_BONDS_ILL * crandn()
                  else
                     S(i,j) = S_SCALE * crandn()
                     if (.not. illcond_case(case_id)) S(i,j) = S(i,j) / sqrt(real(Ndim, dp))
                  endif
               enddo
               S(j,j) = S(j,j) + 1.d0
            enddo
            cond_S_max = max(cond_S_max, cond2(S))
            Sinv = S
            call invert(Sinv)
            ! Cross-slice: carry the pending state through G <- S G S^{-1} (reads Gc = G0 before it is transformed)
            if (cand_cross) call Lowrank_carry_similarity(Gc, S, Sinv)
            do nf = 1, N_FL
               Gref(:,:,nf) = matmul(S, matmul(Gref(:,:,nf), Sinv))
               Gc(:,:,nf)   = matmul(S, matmul(Gc(:,:,nf),   Sinv))
               Aex(:,:,nf)  = matmul(S, matmul(Aex(:,:,nf),  Sinv))
            enddo
         endif
      enddo

      call cand_finalize()
      if (present(std_drift)) std_drift = drift
   end subroutine run_case

end module lowrank_unit_harness


! Eligibility section: Lowrank_init on hand-built vertex lists, independent of the
! candidate-dispatch table above. Delayed: a vertex is eligible iff its support
! max_f Op_V(n,f)%N <= 8, whatever its structure. Submatrix: iff Op_V(n,nf)%N_non_zero <= 1 and
! Op_V(n,nf)%diag in every flavour while SUBMATRIX_RANKN = .false., the delayed rule otherwise.
! A vertex with support > 8 (bigop: diagonal N = 9) is ineligible for both schemes, and the
! delayed scheme reports it with exactly one info line.
module lowrank_eligibility_section
   use lowrank_unit_harness, only: dp, clear_case, info_path, SUBMATRIX_RANKN
   use Hamiltonian_main, only: Op_V, Ndim, N_FL, N_FL_eff, Calc_Fl_map
   use Operator_mod
   use Predefined_Int
   use Lowrank_update_mod, only: Lowrank_init, Lowrank_finalize, Lowrank_active, &
        &                        Lowrank_eligible, Lowrank_pending
   implicit none

   real(dp), parameter :: ELIG_DTAU = 0.1d0, ELIG_U = 4.d0

contains

   subroutine elig_check(ok, label, n_fail)
      logical,          intent(in)    :: ok
      character(len=*), intent(in)    :: label
      integer,          intent(inout) :: n_fail
      if (ok) then
         write(6,'(A,A)') '  eligibility: ok   ', label
      else
         write(6,'(A,A)') '  eligibility: FAIL ', label
         n_fail = n_fail + 1
      endif
   end subroutine elig_check

   ! Non-diagonal rank-1 operator on sites (i,j): O = (1 1; 1 1)/2, eigenvalues 1 and 0.
   subroutine make_nondiag_rank1(Op, i, j)
      type(Operator), intent(inout) :: Op
      integer,        intent(in)    :: i, j
      call Op_make(Op, 2)
      Op%P = [i, j]
      Op%O = cmplx(0.5d0, 0.d0, dp)
      Op%g     = cmplx(sqrt(ELIG_DTAU), 0.d0, dp)
      Op%alpha = cmplx(0.d0, 0.d0, dp)
      Op%type  = 2
      call Op_set(Op)
   end subroutine make_nondiag_rank1

   ! Diagonal rank-1 operator on site i (O = 1).
   subroutine make_diag_rank1(Op, i)
      type(Operator), intent(inout) :: Op
      integer,        intent(in)    :: i
      call Op_make(Op, 1)
      Op%P(1)   = i
      Op%O(1,1) = cmplx(1.d0, 0.d0, dp)
      Op%g      = cmplx(sqrt(ELIG_DTAU), 0.d0, dp)
      Op%alpha  = cmplx(0.d0, 0.d0, dp)
      Op%type   = 2
      call Op_set(Op)
   end subroutine make_diag_rank1

   ! Diagonal operator on the sites P, O = diag(+1, -1, +1, ...): support size(P), N_non_zero = size(P).
   subroutine make_diag_n(Op, P)
      type(Operator), intent(inout) :: Op
      integer,        intent(in)    :: P(:)
      integer :: i
      call Op_make(Op, size(P))
      Op%P = P
      do i = 1, size(P)
         Op%O(i,i) = cmplx(real(1 - 2 * mod(i + 1, 2), dp), 0.d0, dp)
      enddo
      Op%g     = cmplx(sqrt(ELIG_DTAU), 0.d0, dp)
      Op%alpha = cmplx(0.d0, 0.d0, dp)
      Op%type  = 2
      call Op_set(Op)
   end subroutine make_diag_n

   ! Number of lines of the file info equal (after trimming) to want.
   integer function info_count_line(want)
      character(len=*), intent(in) :: want
      character(len=256) :: line
      integer :: u, ios
      logical :: ex
      info_count_line = 0
      inquire(file=info_path, exist=ex)
      if (.not. ex) return
      open(newunit=u, file=info_path, status='old', action='read')
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (trim(adjustl(line)) == want) info_count_line = info_count_line + 1
      enddo
      close(u)
   end function info_count_line

   ! True iff the file 'info' contains a line equal (after trimming) to want.
   logical function info_has_line(want)
      character(len=*), intent(in) :: want
      info_has_line = (info_count_line(want) > 0)
   end function info_has_line

   subroutine run_eligibility_section(n_fail)
      integer, intent(inout) :: n_fail
      integer :: nf, u, n, c0
      logical :: info_existed, all_elig, sub_ok
      character(len=*), parameter :: big_line = 'Sequential scheme DELAYED: 1 of 4 vertices have a support '// &
           & 'larger than 8; they use the standard update (first: n=2)'
      character(len=*), parameter :: big_line_sub = 'Sequential scheme SUBMATRIX: 1 of 4 vertices have a support '// &
           & 'larger than 8; they use the standard update (first: n=2)'

      write(6,'(A)') '36-lowrank-unit: eligibility section'
      call clear_case()
      Ndim = 8; N_FL = 2; N_FL_eff = 2
      allocate(Calc_Fl_map(N_FL))
      Calc_Fl_map = [1, 2]

      ! Vertices: 1, 2 = Mz (diagonal rank-1); 3 = V_SUN bond (non-diagonal rank 2);
      ! 4 = diagonal rank-1 in flavour 1 but non-diagonal rank-1 in flavour 2; 5 = Jz (diagonal rank 2).
      allocate(Op_V(5, N_FL))
      call Predefined_Int_U_MZ(Op_V(1,1), Op_V(1,2), 1, ELIG_DTAU, ELIG_U)
      call Predefined_Int_U_MZ(Op_V(2,1), Op_V(2,2), 2, ELIG_DTAU, ELIG_U)
      do nf = 1, N_FL
         call Predefined_Int_V_SUN(Op_V(3,nf), 3, 4, 1, ELIG_DTAU, 1.d0)
      enddo
      call make_diag_rank1(Op_V(4,1), 7)
      call make_nondiag_rank1(Op_V(4,2), 5, 6)
      call Predefined_Int_Jz(Op_V(5,1), Op_V(5,2), 7, 8, ELIG_DTAU, 1.d0)

      ! The hand-built operators must have the intended structure.
      call elig_check(Op_V(3,1)%N_non_zero == 2 .and. .not. Op_V(3,1)%diag, &
           & 'V_SUN bond: N_non_zero = 2, not diagonal', n_fail)
      call elig_check(Op_V(4,1)%N_non_zero == 1 .and. Op_V(4,1)%diag, &
           & 'hand-built diagonal rank-1: N_non_zero = 1, diagonal', n_fail)
      call elig_check(Op_V(4,2)%N_non_zero == 1 .and. .not. Op_V(4,2)%diag, &
           & 'hand-built (1 1;1 1)/2: N_non_zero = 1, not diagonal', n_fail)
      call elig_check(Op_V(5,1)%N_non_zero == 2 .and. Op_V(5,1)%diag, &
           & 'Jz bond: N_non_zero = 2, diagonal', n_fail)

      inquire(file=info_path, exist=info_existed)

      ! Delayed, mixed-case name: every vertex has support <= 8
      call Lowrank_init('Delayed', 7)
      call elig_check(Lowrank_active(), "Lowrank_init('Delayed', 7): active", n_fail)
      call elig_check(Lowrank_pending() == 0, "Lowrank_init('Delayed', 7): pending = 0", n_fail)
      call elig_check(Lowrank_eligible(1) .and. Lowrank_eligible(2), 'delayed: Mz vertices 1, 2 eligible', n_fail)
      call elig_check(Lowrank_eligible(3), 'delayed: V_SUN bond vertex 3 eligible', n_fail)
      call elig_check(Lowrank_eligible(4), 'delayed: vertex 4 (non-diagonal rank-1 in flavour 2) eligible', n_fail)
      call elig_check(Lowrank_eligible(5), 'delayed: Jz bond vertex 5 eligible', n_fail)
      call elig_check(info_has_line('Sequential scheme: DELAYED, N_delay = 7'), &
           & "info contains 'Sequential scheme: DELAYED, N_delay = 7'", n_fail)
      call Lowrank_finalize()
      call elig_check(.not. Lowrank_active(), 'Lowrank_finalize: inactive', n_fail)

      ! Re-init after finalize (re-init path), submatrix with a mixed-case name
      call Lowrank_init('Submatrix', 3)
      call elig_check(Lowrank_active() .and. Lowrank_pending() == 0, &
           & "Lowrank_init('Submatrix', 3) after finalize: active, pending = 0", n_fail)
      if (SUBMATRIX_RANKN) then
         all_elig = .true.
         do n = 1, 5
            all_elig = all_elig .and. Lowrank_eligible(n)
         enddo
         call elig_check(all_elig, 'submatrix (rank-N): every vertex eligible', n_fail)
      else
         call elig_check(Lowrank_eligible(1) .and. Lowrank_eligible(2) .and. .not. Lowrank_eligible(3) &
              & .and. .not. Lowrank_eligible(4) .and. .not. Lowrank_eligible(5), &
              & 'submatrix (rank-1): only the diagonal rank-1 vertices 1, 2 eligible', n_fail)
      endif
      call elig_check(info_has_line('Sequential scheme: SUBMATRIX, N_delay = 3'), &
           & "info contains 'Sequential scheme: SUBMATRIX, N_delay = 3'", n_fail)
      call Lowrank_finalize()

      ! Standard: nothing active, nothing pending
      call Lowrank_init('standard', 32)
      call elig_check(.not. Lowrank_active(), "Lowrank_init('standard', 32): not active", n_fail)
      call elig_check(Lowrank_pending() == 0, "Lowrank_init('standard', 32): pending = 0", n_fail)
      call Lowrank_finalize()

      ! bigop: Ndim = 12; 1 = Mz, 2 = diagonal N = 9 on sites 1..9, 3 = V_SUN bond, 4 = Mz.
      ! (Its agreement with the reference is checked by the driver, case bigop.)
      call clear_case()
      Ndim = 12; N_FL = 2; N_FL_eff = 2
      allocate(Calc_Fl_map(N_FL))
      Calc_Fl_map = [1, 2]
      allocate(Op_V(4, N_FL))
      call Predefined_Int_U_MZ(Op_V(1,1), Op_V(1,2), 1, ELIG_DTAU, ELIG_U)
      do nf = 1, N_FL
         call make_diag_n(Op_V(2,nf), [(n, n = 1, 9)])
         call Predefined_Int_V_SUN(Op_V(3,nf), 10, 11, 1, ELIG_DTAU, 1.d0)
      enddo
      call Predefined_Int_U_MZ(Op_V(4,1), Op_V(4,2), 12, ELIG_DTAU, ELIG_U)
      call elig_check(Op_V(2,1)%N == 9 .and. Op_V(2,1)%diag .and. Op_V(2,1)%N_non_zero == 9, &
           & 'bigop: diagonal, N = N_non_zero = 9', n_fail)
      c0 = info_count_line(big_line)
      call Lowrank_init('delayed', 5)
      call elig_check(Lowrank_eligible(1) .and. .not. Lowrank_eligible(2) .and. Lowrank_eligible(3) &
           & .and. Lowrank_eligible(4), 'delayed: bigop vertex 2 ineligible, vertices 1, 3, 4 eligible', n_fail)
      call elig_check(info_count_line(big_line) == c0 + 1, &
           & "info contains '"//big_line//"' exactly once", n_fail)
      call Lowrank_finalize()
      c0 = info_count_line(big_line_sub)
      call Lowrank_init('submatrix', 5)
      if (SUBMATRIX_RANKN) then
         sub_ok = Lowrank_eligible(1) .and. .not. Lowrank_eligible(2) .and. Lowrank_eligible(3) &
              & .and. Lowrank_eligible(4) .and. info_count_line(big_line_sub) == c0 + 1
      else
         sub_ok = Lowrank_eligible(1) .and. .not. Lowrank_eligible(2) .and. .not. Lowrank_eligible(3) &
              & .and. Lowrank_eligible(4)
      endif
      call elig_check(sub_ok, 'submatrix: bigop vertex 2 ineligible (rank-N: reported once in info)', n_fail)
      call Lowrank_finalize()

      ! Remove the info file if this section created it.
      if (.not. info_existed) then
         open(newunit=u, file=info_path, status='old')
         close(u, status='delete')
      endif
      call clear_case()
   end subroutine run_eligibility_section

end module lowrank_eligibility_section


! Kernel-equivalence section: from a flushed state with the same RNG state, N_PROP_EQ
! proposals (case mz, mode "Final", S0_ratio = T0_proposal_ratio = 1, proposals drawn once
! beforehand) through the real Upgrade2 and through Upgrade_lowrank (delayed and submatrix),
! each on its own G copy and its own copy of the global nsigma. The proposals hit random
! vertices, so sites repeat within a block (rank-N submatrix: no flush, one M-form accept formula). Asserted: identical toggle sequence and
! fields, |Delta Phase| <= PHASE_TOL, the same next ranf_wrap() value (same RNG consumption),
! and G after the final flush (and Lowrank_current_G before it) within TOL_UNIT.
module lowrank_kernel_equiv_section
   use lowrank_unit_harness, only: dp, TOL_UNIT, N_SLICES, N_DELAY_LIST, N_NDELAY, CASE_MZ, &
        &                          M_ops, setup_case, random_initial_G, clear_case
   use Hamiltonian_main, only: nsigma, Ndim, N_FL
   use upgrade_mod
   use Random_Wrap
   use Lowrank_update_mod, only: Lowrank_init, Lowrank_finalize, Upgrade_lowrank, Lowrank_flush, &
        &                        Lowrank_pending, Lowrank_current_G
   implicit none

   integer,  parameter :: N_PROP_EQ = 200
   real(dp), parameter :: PHASE_TOL = 1.0d-12

contains

   subroutine run_kernel_equivalence(n_fail)
      integer, intent(inout) :: n_fail

      complex(dp), allocatable :: G0(:,:,:), Ga(:,:,:), Gb(:,:,:), Aex(:,:,:), Gcur(:,:)
      complex(dp), allocatable :: f_save(:,:), f_a(:,:)
      integer,     allocatable :: rng_state(:)
      integer     :: n_eq(N_PROP_EQ), nt_eq(N_PROP_EQ)
      complex(dp) :: hs_eq(N_PROP_EQ)
      logical     :: tog_a(N_PROP_EQ), tog_b(N_PROP_EQ)
      integer     :: id, k, nf, k_seed, n_acc, n_mis, n_pend_end, isch
      character(len=10), parameter :: eq_scheme(2) = [character(len=10) :: 'delayed', 'submatrix']
      complex(dp) :: phase0, phase_a, phase_b, prev
      real(dp)    :: next_a, next_b, dphase, gerr_cur, gerr, gmax
      character(len=64) :: mode
      character(len=4)  :: flag
      logical     :: same_f

      write(6,'(A,I0,A)') '36-lowrank-unit: kernel equivalence, Upgrade2(Final) vs Upgrade_lowrank, case mz, ', &
           & N_PROP_EQ, ' proposals'
      write(6,'(A10,1X,A7,1X,A5,1X,A8,1X,A8,1X,A10,1X,A22,1X,A22,1X,A11,1X,A11,1X,A6,1X,A4)') 'scheme', 'N_delay', &
           & 'acc', 'tog_mis', 'fld_same', '|dPhase|', 'next ranf (Upgrade2)', 'next ranf (lowrank)', 'G_cur_err', &
           & 'G_err', 'n_pend', ''
      mode = "Final"
      call Get_seed_Len(k_seed)
      allocate(rng_state(k_seed))
      do isch = 1, size(eq_scheme)
      do id = 1, N_NDELAY
         call Ranset([1000 + 17 * CASE_MZ])
         call setup_case(CASE_MZ)
         allocate(G0(Ndim,Ndim,N_FL), Ga(Ndim,Ndim,N_FL), Gb(Ndim,Ndim,N_FL), Aex(Ndim,Ndim,N_FL), Gcur(Ndim,Ndim))
         do nf = 1, N_FL
            call random_initial_G(G0(:,:,nf), Aex(:,:,nf))
         enddo
         ! Proposals, drawn once (relative to the initial configuration)
         do k = 1, N_PROP_EQ
            n_eq(k)  = min(M_ops,    1 + int(real(M_ops, dp)    * ranf_wrap()))
            nt_eq(k) = min(N_SLICES, 1 + int(real(N_SLICES, dp) * ranf_wrap()))
            hs_eq(k) = nsigma%flip(n_eq(k), nt_eq(k))
         enddo
         phase0 = cmplx(cos(0.3d0), sin(0.3d0), dp)
         f_save = nsigma%f
         call random_seed(get=rng_state)

         ! ---- pass A: Upgrade2, mode Final ----
         Ga = G0
         phase_a = phase0
         do k = 1, N_PROP_EQ
            prev = cmplx(1.d0, 0.d0, dp)
            call Upgrade2(Ga, n_eq(k), nt_eq(k), phase_a, hs_eq(k), prev, 1.d0, 1.d0, tog_a(k), mode)
         enddo
         next_a = ranf_wrap()
         f_a = nsigma%f

         ! ---- pass B: Upgrade_lowrank, same RNG state and fields ----
         nsigma%f = f_save
         call random_seed(put=rng_state)
         call Lowrank_init(trim(eq_scheme(isch)), N_DELAY_LIST(id))
         Gb = G0
         phase_b = phase0
         do k = 1, N_PROP_EQ
            call Upgrade_lowrank(Gb, n_eq(k), nt_eq(k), phase_b, hs_eq(k), 1.d0, 1.d0, tog_b(k))
         enddo
         next_b = ranf_wrap()
         n_pend_end = Lowrank_pending()
         gmax = maxval(abs(Ga))
         gerr_cur = 0.d0
         do nf = 1, N_FL
            call Lowrank_current_G(Gb, nf, Gcur)
            gerr_cur = max(gerr_cur, maxval(abs(Gcur - Ga(:,:,nf))) / gmax)
         enddo
         call Lowrank_flush(Gb)
         gerr = maxval(abs(Gb - Ga)) / gmax
         same_f = all(nsigma%f == f_a)
         call Lowrank_finalize()
         nsigma%f = f_save

         n_acc  = count(tog_a)
         n_mis  = count(tog_a .neqv. tog_b)
         dphase = abs(phase_a - phase_b)
         flag = ''
         if (n_mis /= 0 .or. .not. same_f .or. .not. (dphase <= PHASE_TOL) .or. next_a /= next_b &
              & .or. .not. (gerr <= TOL_UNIT) .or. .not. (gerr_cur <= TOL_UNIT)) then
            flag = 'FAIL'
            n_fail = n_fail + 1
         endif
         write(6,'(A10,1X,I7,1X,I5,1X,I8,1X,L8,1X,ES10.3,1X,F22.17,1X,F22.17,1X,ES11.4,1X,ES11.4,1X,I6,1X,A4)') &
              & eq_scheme(isch), N_DELAY_LIST(id), n_acc, n_mis, same_f, dphase, next_a, next_b, gerr_cur, gerr, &
              & n_pend_end, flag
         deallocate(G0, Ga, Gb, Aex, Gcur)
      enddo
      enddo
      call clear_case()
   end subroutine run_kernel_equivalence

end module lowrank_kernel_equiv_section


program lowrank_unit
   use lowrank_unit_harness
   use lowrank_eligibility_section, only: run_eligibility_section
   use lowrank_kernel_equiv_section, only: run_kernel_equivalence
#ifdef MPI
   use mpi
   use Hamiltonian_main, only: Group_Comm
#endif
   implicit none
   integer :: ic, is, id, idir, n_fail, ics, is_first
   character(len=10) :: row_scheme
   real(dp) :: rerr, gerr, drift, max_drift, tol, drift_rel
   character(len=4), parameter :: dirs(2) = [character(len=4) :: 'up', 'down']
   character(len=4) :: flag
   logical :: info_existed
   integer :: u_info
#ifdef MPI
   integer :: ierr
#endif
#if defined(MPI) && defined(TEMPERING)
   logical :: temp_dir_existed
#endif

#ifdef MPI
   ! Lowrank_init queries the communicators: a single rank, group 0
   call MPI_Init(ierr)
   Group_Comm = MPI_COMM_WORLD
#endif
   call Fields_init(1.d0)
   n_fail = 0
#if defined(MPI) && defined(TEMPERING)
   ! Lowrank_init appends to Temp_0/info
   inquire(file='Temp_0', exist=temp_dir_existed)
   if (.not. temp_dir_existed) call execute_command_line('mkdir -p Temp_0')
#endif
   ! Lowrank_init (delayed candidate, eligibility section) appends to the file info_path.
   inquire(file=info_path, exist=info_existed)
   max_drift = 0.d0

   write(6,'(A,ES9.2,A,ES9.2,A)') '36-lowrank-unit: TOL_UNIT = ', TOL_UNIT, ', TOL_ILL = ', TOL_ILL, &
        & ' (mz_illcond, bonds_illcond; the fixed bound drift <= TOL_UNIT/2 is not asserted there)'
   write(6,'(A,L1)') '  SUBMATRIX_RANKN = ', SUBMATRIX_RANKN
   write(6,'(A)') '  std_drift/ce = std_drift / (cond(1+A) * epsilon(1d0)); asserted <= 10 for every case'
   write(6,'(A)') '  rflush = flushes inside Lowrank_propose (site already pending); asserted > 0 iff submatrix,'// &
        & ' case repeat, repeat_skip or rep_mixed, N_delay > 1 (only while SUBMATRIX_RANKN = .false.;'// &
        & ' repeat/repeat_skip with N_delay = 2 not asserted)'
   write(6,'(A)') '  skip = vertices wrapped but not proposed; asserted > 0 iff case repeat_skip or bonds_skip'
   write(6,'(A)') '  lrv = vertices updated by the low-rank kernel; Lowrank_eligible asserted per vertex (flag ELG)'
   write(6,'(A)') '  pmax = largest pending rank right after an accept; asserted <= max(N_delay, r) per accept (flag CAP)'
   write(6,'(A13,1X,A10,1X,A4,1X,A7,1X,A13,1X,A13,1X,A13,1X,A10,1X,A12,1X,A7,1X,A6,1X,A6,1X,A6,1X,A5,1X,A4,1X,A4,1X,A4)') &
        & 'case', 'scheme', 'dir', 'N_delay', 'max_ratio_err', 'max_G_err', 'std_drift', &
        & 'cond(1+A)', 'std_drift/ce', 'acc', 'condQ', 'condS', 'rflush', 'skip', 'lrv', 'pmax', ''
   ! ics = 0: the main table; ics = 1: the cross-slice variant (delayed, submatrix only)
   do ics = 0, 1
   if (ics == 1) write(6,'(A)') '36-lowrank-unit: cross-slice variant (pending state carried through the '// &
        & 'inter-slice similarity by Lowrank_carry_similarity; flushed only at the last slice end)'
   is_first = 1
   if (ics == 1) is_first = SCHEME_DELAYED
   do ic = 1, N_CASES
      do idir = 1, 2
         do is = is_first, N_SCHEMES
            do id = 1, N_NDELAY
               call run_case(ic, is, N_DELAY_LIST(id), trim(dirs(idir)), rerr, gerr, drift, cross_slice=(ics == 1))
               row_scheme = scheme_name(is)
               if (ics == 1) row_scheme = scheme_name_cs(is)
               flag = ''
               tol = TOL_UNIT
               if (illcond_case(ic)) tol = TOL_ILL
               if (rerr > tol .or. gerr > tol .or. .not. (rerr == rerr .and. gerr == gerr)) then
                  flag = 'FAIL'
                  n_fail = n_fail + 1
               endif
               ! TOL_UNIT must stay >= 2x the standard path's own drift vs exact. Not asserted
               ! for the ill-conditioned cases: the double-precision "exact" inverse is unreliable
               ! at cond ~1e8. Every rank-N case meets it (no exemption; couplings: V_BOND, G_ND, U_REP).
               if (.not. illcond_case(ic)) then
                  if (drift > 0.5d0 * TOL_UNIT .or. drift /= drift) then
                     flag = 'DRFT'
                     n_fail = n_fail + 1
                  endif
                  max_drift = max(max_drift, drift)
               endif
               ! Sanity: the standard path must sit at the double-precision floor,
               ! drift <= 10 * cond(1+A) * eps (cond at the last slice = the largest one).
               drift_rel = drift / (cond_1pA_last * epsilon(1.d0))
               if (drift_rel > 10.d0 .or. drift_rel /= drift_rel) then
                  flag = 'COND'
                  n_fail = n_fail + 1
               endif
               ! Repeated-site flush (project-1 submatrix only): only submatrix flushes in Lowrank_propose,
               ! and only when a site recurs within a block: cases repeat, repeat_skip (random order of the
               ! doubled site list; with N_delay = 2 no site of this configuration recurs within a block of
               ! two accepts, so N_delay = 2 is not asserted there) and rep_mixed (the last Mz vertex of a
               ! group and the first of the next share the site, so N_delay >= 2 suffices); N_delay = 1
               ! never leaves a block open. Not asserted for the rank-N submatrix (SUBMATRIX_RANKN = .true.).
               if (.not. SUBMATRIX_RANKN .and. .not. (is == SCHEME_SUBMATRIX .and. N_DELAY_LIST(id) == 2 .and. &
                    & (ic == CASE_REPEAT .or. ic == CASE_REPEAT_SKIP))) then
                  if ((n_repeat_flush > 0) .neqv. (is == SCHEME_SUBMATRIX .and. N_DELAY_LIST(id) > 1 .and. &
                       & (ic == CASE_REPEAT .or. ic == CASE_REPEAT_SKIP .or. ic == CASE_REP_MIXED))) then
                     flag = 'RFL'
                     n_fail = n_fail + 1
                  endif
               endif
               ! The skip mask must act exactly in cases repeat_skip, bonds_skip
               if ((n_skipped > 0) .neqv. (ic == CASE_REPEAT_SKIP .or. ic == CASE_BONDS_SKIP)) then
                  flag = 'SKP'
                  n_fail = n_fail + 1
               endif
               ! Lowrank_eligible must match the expected eligibility rule at every vertex
               if (n_elig_bad /= 0) then
                  flag = 'ELG'
                  n_fail = n_fail + 1
               endif
               ! Capacity rule: pending rank after every accept <= max(N_delay, r of that vertex)
               if (n_cap_bad /= 0) then
                  flag = 'CAP'
                  n_fail = n_fail + 1
               endif
               write(6,'(A13,1X,A10,1X,A4,1X,I7,1X,ES13.4,1X,ES13.4,1X,ES13.4,1X,ES10.2,1X,ES12.2,1X,F7.3,1X,F6.2,1X,F6.2,1X,I6,1X,I5,1X,I4,1X,I4,1X,A4)') &
                    & case_name(ic), row_scheme, dirs(idir), N_DELAY_LIST(id), rerr, gerr, drift, &
                    & cond_1pA_last, drift_rel, real(n_accepted, dp) / real(n_proposals, dp), cond_Q_max, cond_S_max, &
                    & n_repeat_flush, n_skipped, n_lowrank_vertices, max_pend_cols, flag
            enddo
         enddo
      enddo
   enddo
   enddo   ! ics
   write(6,'(A,ES10.3)') 'max std_drift_vs_exact (well-conditioned cases) = ', max_drift
   call run_eligibility_section(n_fail)
   call run_kernel_equivalence(n_fail)
   if (.not. info_existed) then
      inquire(file=info_path, exist=info_existed)
      if (info_existed) then
         open(newunit=u_info, file=info_path, status='old')
         close(u_info, status='delete')
      endif
   endif
#if defined(MPI) && defined(TEMPERING)
   if (.not. temp_dir_existed) call execute_command_line('rmdir Temp_0')
#endif
   if (n_fail > 0) then
      write(6,'(A,I0,A)') '36-lowrank-unit: ', n_fail, ' row(s) FAILED'
#ifdef MPI
      call MPI_Finalize(ierr)
#endif
      stop 1
   endif
   write(6,'(A)') '36-lowrank-unit: all rows passed'
#ifdef MPI
   call MPI_Finalize(ierr)
#endif
end program lowrank_unit

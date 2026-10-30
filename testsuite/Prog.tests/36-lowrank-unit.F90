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
!   up  : Op_Wrapup type 1 (old phi) -> hook(+1) -> propose -> decide -> Op_Wrapup type 2
!   down: Op_Wrapdo type 2 -> propose -> decide -> Op_Wrapdo type 1 (new phi) -> hook(-1)
! Between slices: flush, compare, then the same dense similarity G <- S G S^{-1} on every G,
! S = I + 0.1 R / sqrt(Ndim) (R complex Gaussian). Case mz_illcond = mz with S = I + 0.1 R,
! which drives cond(1+A) to ~1e8: a deliberately hard case with its own tolerance TOL_ILL.
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
        &                        Lowrank_wrap_hook, Lowrank_flush, Upgrade_lowrank
   implicit none

   integer, parameter :: dp = kind(0.d0)

   ! Tolerances (see CLAUDE.md "Tolerances"): TOL_UNIT for the well-conditioned cases,
   ! TOL_ILL for mz_illcond (cond(1+A) ~ 1e8, errors scale ~ cond * eps).
   real(dp), parameter :: TOL_UNIT = 1.0d-10
   real(dp), parameter :: TOL_ILL  = 1.0d-8

   integer, parameter :: N_CASES = 5
   integer, parameter :: CASE_MZ = 1, CASE_DENSITY = 2, CASE_CONTINUOUS = 3, CASE_REPEAT = 4
   integer, parameter :: CASE_MZ_ILLCOND = 5   ! mz with the unnormalised S = I + 0.1 R
   character(len=10), parameter :: case_name(N_CASES) = &
        & [character(len=10) :: 'mz', 'density', 'continuous', 'repeat', 'mz_illcond']

   ! Candidate schemes. Task 11 adds SCHEME_SUBMATRIX here.
   integer, parameter :: N_SCHEMES = 2
   integer, parameter :: SCHEME_STANDARD = 1, SCHEME_DELAYED = 2
   character(len=10), parameter :: scheme_name(N_SCHEMES) = [character(len=10) :: 'standard', 'delayed']

   integer, parameter :: N_NDELAY = 4
   integer, parameter :: N_DELAY_LIST(N_NDELAY) = [1, 7, 32, 1000]

   integer, parameter :: N_SLICES = 3
   real(dp), parameter :: DTAU = 0.1d0, U_HUB = 4.d0
   real(dp), parameter :: ACCEPT_FRACTION = 0.6d0
   real(dp), parameter :: Q_SCALE = 0.3d0     ! Q = I + Q_SCALE * R / sqrt(Ndim)
   real(dp), parameter :: S_SCALE = 0.1d0     ! S = I + S_SCALE * R / sqrt(Ndim)
                                              ! (mz_illcond: S = I + S_SCALE * R)

   integer :: M_ops                            ! number of vertices per slice

   ! Diagnostics of the last run_case
   integer  :: n_proposals, n_accepted
   real(dp) :: cond_Q_max, cond_S_max
   real(dp) :: cond_1pA_last                   ! max_f cond_2(1+A_f) at the last slice end

   ! ---- candidate state -------------------------------------------------
   integer :: cand_scheme = 0, cand_N_delay = 0
   complex(dp), allocatable :: std_Gwork(:,:,:)     ! standard: Upgrade2 on a work copy

contains

   ! =====================================================================
   ! Candidate interface
   ! =====================================================================
   subroutine cand_init(scheme, N_delay)
      integer, intent(in) :: scheme, N_delay
      cand_scheme  = scheme
      cand_N_delay = N_delay
      select case (scheme)
      case (SCHEME_STANDARD)
         if (allocated(std_Gwork)) deallocate(std_Gwork)
         allocate(std_Gwork(Ndim, Ndim, N_FL))
      case (SCHEME_DELAYED)
         ! Works on the driver's Gc (= G0 between flushes); prints the scheme info line.
         call Lowrank_init('delayed', N_delay)
      case default
         error stop 'cand_init: unknown scheme'
      end select
   end subroutine cand_init

   subroutine cand_finalize()
      select case (cand_scheme)
      case (SCHEME_STANDARD)
         if (allocated(std_Gwork)) deallocate(std_Gwork)
      case (SCHEME_DELAYED)
         call Lowrank_finalize()
      end select
   end subroutine cand_finalize

   ! Called after every type-1 wrap of vertex n_op (dir=+1 up, -1 down).
   subroutine cand_hook(Gc, n_op, nt, HS_field, dir)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op, nt, dir
      complex(dp), intent(in)    :: HS_field
      select case (cand_scheme)
      case (SCHEME_STANDARD)
         ! no-op
      case (SCHEME_DELAYED)
         call Lowrank_wrap_hook(n_op, nt, HS_field, dir)
      end select
   end subroutine cand_hook

   ! Returns the full Upgrade2-style ratio. nsigma%f(n_op,nt) holds the OLD field.
   subroutine cand_propose(Gc, n_op, nt, HS_new, ratio)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op, nt
      complex(dp), intent(in)    :: HS_new
      complex(dp), intent(out)   :: ratio
      complex(dp)        :: phase
      logical            :: toggle
      character(len=64)  :: mode
      select case (cand_scheme)
      case (SCHEME_STANDARD)
         std_Gwork = Gc
         ratio = cmplx(1.d0, 0.d0, dp)
         phase = cmplx(1.d0, 0.d0, dp)
         mode  = "Intermediate"
         call Upgrade2(std_Gwork, n_op, nt, phase, HS_new, ratio, 1.d0, 1.d0, toggle, mode)
      case (SCHEME_DELAYED)
         call Lowrank_propose(Gc, n_op, nt, HS_new, ratio)
      end select
   end subroutine cand_propose

   ! Called only for accepted proposals, right after cand_propose.
   ! nsigma%f(n_op,nt) holds the OLD field.
   subroutine cand_accept(Gc, n_op, nt, HS_new)
      complex(dp), intent(inout) :: Gc(:,:,:)
      integer,     intent(in)    :: n_op, nt
      complex(dp), intent(in)    :: HS_new
      select case (cand_scheme)
      case (SCHEME_STANDARD)
         Gc = std_Gwork
      case (SCHEME_DELAYED)
         call Lowrank_accept(Gc, n_op, nt, HS_new)   ! sets nsigma%f(n_op,nt); the driver resets it
      end select
   end subroutine cand_accept

   ! Called at the end of every slice; afterwards Gc must hold the full G.
   subroutine cand_flush(Gc)
      complex(dp), intent(inout) :: Gc(:,:,:)
      select case (cand_scheme)
      case (SCHEME_STANDARD)
         ! no-op
      case (SCHEME_DELAYED)
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

   subroutine setup_case(case_id)
      integer, intent(in) :: case_id
      integer, allocatable :: sites(:), base(:)
      integer :: n, nf, nt, k
      real(dp), parameter :: disc(4) = [-2.d0, -1.d0, 1.d0, 2.d0]

      call clear_case()
      select case (case_id)
      case (CASE_MZ, CASE_MZ_ILLCOND)
         Ndim = 64; N_FL = 2; N_SUN = 1; M_ops = Ndim
      case (CASE_DENSITY)
         Ndim = 64; N_FL = 1; N_SUN = 2; M_ops = Ndim
      case (CASE_CONTINUOUS)
         Ndim = 48; N_FL = 2; N_SUN = 1; M_ops = Ndim
      case (CASE_REPEAT)
         Ndim = 32; N_FL = 2; N_SUN = 1; M_ops = 2 * Ndim
      case default
         error stop 'setup_case: unknown case'
      end select
      N_FL_eff = N_FL
      allocate(Calc_Fl_map(N_FL))
      Calc_Fl_map = [(nf, nf = 1, N_FL)]
      reconstruction_needed = .false.

      ! Vertex -> site: a fixed random permutation (repeat: of the doubled site list)
      allocate(sites(M_ops))
      if (case_id == CASE_REPEAT) then
         allocate(base(M_ops))
         call random_permutation(base)
         do n = 1, M_ops
            sites(n) = mod(base(n) - 1, Ndim) + 1
         enddo
      else
         call random_permutation(sites)
      endif

      allocate(Op_V(M_ops, N_FL))
      do n = 1, M_ops
         select case (case_id)
         case (CASE_MZ, CASE_REPEAT, CASE_MZ_ILLCOND)
            call Predefined_Int_U_MZ(Op_V(n,1), Op_V(n,2), sites(n), DTAU, U_HUB)
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
      if (case_id == CASE_DENSITY) then
         if (abs(aimag(Op_V(1,1)%g)) < 1.d-12) error stop 'setup_case: density g is not complex'
      endif

      call nsigma%make(M_ops, N_SLICES)
      do n = 1, M_ops
         nsigma%t(n) = Op_V(n,1)%type
      enddo
      do nt = 1, N_SLICES
         do n = 1, M_ops
            if (nsigma%t(n) == 2) then
               k = 1 + int(4.d0 * ranf_wrap())
               if (k > 4) k = 4
               nsigma%f(n,nt) = cmplx(disc(k), 0.d0, dp)
            else
               nsigma%f(n,nt) = cmplx(3.d0 * ranf_wrap() - 1.5d0, 0.d0, dp)
            endif
         enddo
      enddo
   end subroutine setup_case

   ! =====================================================================
   ! Exact tracker (diagonal rank-1 vertices; vertex frame = site frame)
   ! =====================================================================
   ! A <- T A T^{-1}, T = diag(t at p), t = exp(dir * g * phi * E(1))
   subroutine exact_wrap(A, n_op, nt, dir)
      complex(dp), intent(inout) :: A(:,:,:)
      integer, intent(in) :: n_op, nt, dir
      integer :: nf, p
      complex(dp) :: t, phi
      phi = nsigma%Phi(n_op, nt)
      do nf = 1, N_FL
         p = Op_V(n_op,nf)%P(1)
         t = exp(real(dir, dp) * Op_V(n_op,nf)%g * phi * Op_V(n_op,nf)%E(1))
         A(p,:,nf) = A(p,:,nf) * t
         A(:,p,nf) = A(:,p,nf) / t
      enddo
   end subroutine exact_wrap

   ! Accepted flip: A <- (1 + Delta e_p e_p^T) A
   subroutine exact_flip(A, n_op, phi_old, phi_new)
      complex(dp), intent(inout) :: A(:,:,:)
      integer, intent(in) :: n_op
      complex(dp), intent(in) :: phi_old, phi_new
      integer :: nf, p
      complex(dp) :: delta
      do nf = 1, N_FL
         p = Op_V(n_op,nf)%P(1)
         delta = exp(Op_V(n_op,nf)%g * (phi_new - phi_old) * Op_V(n_op,nf)%E(1)) - 1.d0
         A(p,:,nf) = A(p,:,nf) * (1.d0 + delta)
      enddo
   end subroutine exact_flip

   ! =====================================================================
   ! Driver
   ! =====================================================================
   subroutine run_case(case_id, scheme, N_delay, direction, max_ratio_err, max_G_err, std_drift)
      integer,          intent(in)  :: case_id, scheme, N_delay
      character(len=*), intent(in)  :: direction
      real(dp),         intent(out) :: max_ratio_err, max_G_err
      real(dp),         intent(out), optional :: std_drift

      complex(dp), allocatable :: Gref(:,:,:), Gc(:,:,:), Gtmp(:,:,:), Aex(:,:,:), Gex(:,:)
      complex(dp), allocatable :: S(:,:), Sinv(:,:)
      integer,     allocatable :: rng_state(:)
      integer :: islice, nt, iv, n, nf, dir, i, j, k_seed
      complex(dp) :: hs_old, hs_new, phi_old, phi_new, R_ref, R_c, phase
      logical :: accept, toggle
      character(len=64) :: mode
      real(dp) :: drift, gmax

      if (direction == 'up') then
         dir = +1
      elseif (direction == 'down') then
         dir = -1
      else
         error stop 'run_case: direction must be up or down'
      endif

      ! The configuration depends on the case only (mz_illcond reuses mz's configuration).
      if (case_id == CASE_MZ_ILLCOND) then
         call Ranset([1000 + 17 * CASE_MZ])
      else
         call Ranset([1000 + 17 * case_id])
      endif
      call setup_case(case_id)
      n_proposals = 0; n_accepted = 0
      cond_Q_max = 0.d0; cond_S_max = 0.d0; cond_1pA_last = 0.d0

      allocate(Gref(Ndim,Ndim,N_FL), Gc(Ndim,Ndim,N_FL), Gtmp(Ndim,Ndim,N_FL), Aex(Ndim,Ndim,N_FL))
      allocate(Gex(Ndim,Ndim), S(Ndim,Ndim), Sinv(Ndim,Ndim))
      do nf = 1, N_FL
         call random_initial_G(Gref(:,:,nf), Aex(:,:,nf))
      enddo
      Gc = Gref
      call Get_seed_Len(k_seed)
      allocate(rng_state(k_seed))

      call cand_init(scheme, N_delay)
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
            if (dir == 1) then
               do nf = 1, N_FL
                  call Op_Wrapup(Gref(:,:,nf), Op_V(n,nf), hs_old, Ndim, 1, nt)
                  call Op_Wrapup(Gc(:,:,nf),   Op_V(n,nf), hs_old, Ndim, 1, nt)
               enddo
               call exact_wrap(Aex, n, nt, +1)
               call cand_hook(Gc, n, nt, hs_old, +1)
            else
               do nf = 1, N_FL
                  call Op_Wrapdo(Gref(:,:,nf), Op_V(n,nf), hs_old, Ndim, 2, nt)
                  call Op_Wrapdo(Gc(:,:,nf),   Op_V(n,nf), hs_old, Ndim, 2, nt)
               enddo
            endif

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
               call exact_flip(Aex, n, phi_old, phi_new)
            else
               nsigma%f(n,nt) = hs_old
            endif
            call random_seed(put=rng_state)

            if (dir == 1) then
               do nf = 1, N_FL
                  call Op_Wrapup(Gref(:,:,nf), Op_V(n,nf), hs_old, Ndim, 2, nt)
                  call Op_Wrapup(Gc(:,:,nf),   Op_V(n,nf), hs_old, Ndim, 2, nt)
               enddo
            else
               hs_new = nsigma%f(n,nt)       ! post-decision field
               do nf = 1, N_FL
                  call Op_Wrapdo(Gref(:,:,nf), Op_V(n,nf), hs_new, Ndim, 1, nt)
                  call Op_Wrapdo(Gc(:,:,nf),   Op_V(n,nf), hs_new, Ndim, 1, nt)
               enddo
               call exact_wrap(Aex, n, nt, -1)
               call cand_hook(Gc, n, nt, hs_new, -1)
            endif
         enddo

         ! ---- slice end: flush, compare ----
         call cand_flush(Gc)
#ifdef LOWRANK_UNIT_SELFTEST
         if (case_id == CASE_MZ .and. dir == 1 .and. N_delay == N_DELAY_LIST(1) .and. islice == 1) &
              & Gc(1,1,1) = Gc(1,1,1) + 1.0d-8
#endif
         gmax = maxval(abs(Gref))
         max_G_err = max(max_G_err, maxval(abs(Gc - Gref)) / gmax)
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
                  S(i,j) = S_SCALE * crandn()
                  if (case_id /= CASE_MZ_ILLCOND) S(i,j) = S(i,j) / sqrt(real(Ndim, dp))
               enddo
               S(j,j) = S(j,j) + 1.d0
            enddo
            cond_S_max = max(cond_S_max, cond2(S))
            Sinv = S
            call invert(Sinv)
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


! Eligibility section: Lowrank_init on a mixed vertex list, independent of the
! candidate-dispatch table above. A vertex is eligible iff, in every flavour,
! Op_V(n,nf)%N_non_zero <= 1 and Op_V(n,nf)%diag.
module lowrank_eligibility_section
   use lowrank_unit_harness, only: dp, clear_case
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

   ! True iff the file 'info' contains a line equal (after trimming) to want.
   logical function info_has_line(want)
      character(len=*), intent(in) :: want
      character(len=256) :: line
      integer :: u, ios
      logical :: ex
      info_has_line = .false.
      inquire(file='info', exist=ex)
      if (.not. ex) return
      open(newunit=u, file='info', status='old', action='read')
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (trim(adjustl(line)) == want) info_has_line = .true.
      enddo
      close(u)
   end function info_has_line

   subroutine run_eligibility_section(n_fail)
      integer, intent(inout) :: n_fail
      integer :: nf, u
      logical :: info_existed

      write(6,'(A)') '36-lowrank-unit: eligibility section'
      call clear_case()
      Ndim = 8; N_FL = 2; N_FL_eff = 2
      allocate(Calc_Fl_map(N_FL))
      Calc_Fl_map = [1, 2]

      ! Vertices: 1, 2 = Mz (eligible); 3 = V_SUN bond, rank 2 (ineligible);
      ! 4 = diagonal rank-1 in flavour 1 but non-diagonal rank-1 in flavour 2 (ineligible).
      allocate(Op_V(4, N_FL))
      call Predefined_Int_U_MZ(Op_V(1,1), Op_V(1,2), 1, ELIG_DTAU, ELIG_U)
      call Predefined_Int_U_MZ(Op_V(2,1), Op_V(2,2), 2, ELIG_DTAU, ELIG_U)
      do nf = 1, N_FL
         call Predefined_Int_V_SUN(Op_V(3,nf), 3, 4, 1, ELIG_DTAU, 1.d0)
      enddo
      call make_diag_rank1(Op_V(4,1), 7)
      call make_nondiag_rank1(Op_V(4,2), 5, 6)

      ! The hand-built operators must have the intended structure.
      call elig_check(Op_V(3,1)%N_non_zero == 2 .and. .not. Op_V(3,1)%diag, &
           & 'V_SUN bond: N_non_zero = 2, not diagonal', n_fail)
      call elig_check(Op_V(4,1)%N_non_zero == 1 .and. Op_V(4,1)%diag, &
           & 'hand-built diagonal rank-1: N_non_zero = 1, diagonal', n_fail)
      call elig_check(Op_V(4,2)%N_non_zero == 1 .and. .not. Op_V(4,2)%diag, &
           & 'hand-built (1 1;1 1)/2: N_non_zero = 1, not diagonal', n_fail)

      inquire(file='info', exist=info_existed)

      ! Delayed, mixed-case name
      call Lowrank_init('Delayed', 7)
      call elig_check(Lowrank_active(), "Lowrank_init('Delayed', 7): active", n_fail)
      call elig_check(Lowrank_pending() == 0, "Lowrank_init('Delayed', 7): pending = 0", n_fail)
      call elig_check(Lowrank_eligible(1) .and. Lowrank_eligible(2), 'Mz vertices 1, 2 eligible', n_fail)
      call elig_check(.not. Lowrank_eligible(3), 'V_SUN bond vertex 3 ineligible', n_fail)
      call elig_check(.not. Lowrank_eligible(4), 'vertex 4 (non-diagonal rank-1 in flavour 2) ineligible', n_fail)
      call elig_check(info_has_line('Sequential scheme: DELAYED, N_delay = 7'), &
           & "info contains 'Sequential scheme: DELAYED, N_delay = 7'", n_fail)
      call Lowrank_finalize()
      call elig_check(.not. Lowrank_active(), 'Lowrank_finalize: inactive', n_fail)

      ! Re-init after finalize (re-init path). Lowrank_init('submatrix') terminates until Phase 4 implements the
      ! scheme (it cannot be tested in-process); Phase 4 must restore a submatrix re-init check here.
      call Lowrank_init('delayed', 3)
      call elig_check(Lowrank_active() .and. Lowrank_pending() == 0, &
           & "Lowrank_init('delayed', 3) after finalize: active, pending = 0", n_fail)
      call elig_check(Lowrank_eligible(1) .and. Lowrank_eligible(2) .and. .not. Lowrank_eligible(3) &
           & .and. .not. Lowrank_eligible(4), 're-init: same eligibility pattern', n_fail)
      call Lowrank_finalize()

      ! Standard: nothing active, nothing pending
      call Lowrank_init('standard', 32)
      call elig_check(.not. Lowrank_active(), "Lowrank_init('standard', 32): not active", n_fail)
      call elig_check(Lowrank_pending() == 0, "Lowrank_init('standard', 32): pending = 0", n_fail)
      call Lowrank_finalize()

      ! Remove the info file if this section created it.
      if (.not. info_existed) then
         open(newunit=u, file='info', status='old')
         close(u, status='delete')
      endif
      call clear_case()
   end subroutine run_eligibility_section

end module lowrank_eligibility_section


! Kernel-equivalence section: from a flushed state with the same RNG state, N_PROP_EQ
! proposals (case mz, mode "Final", S0_ratio = T0_proposal_ratio = 1, proposals drawn once
! beforehand) through the real Upgrade2 and through Upgrade_lowrank (delayed), each on its
! own G copy and its own copy of the global nsigma. Asserted: identical toggle sequence and
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
      integer     :: id, k, nf, k_seed, n_acc, n_mis, n_pend_end
      complex(dp) :: phase0, phase_a, phase_b, prev
      real(dp)    :: next_a, next_b, dphase, gerr_cur, gerr, gmax
      character(len=64) :: mode
      character(len=4)  :: flag
      logical     :: same_f

      write(6,'(A,I0,A)') '36-lowrank-unit: kernel equivalence, Upgrade2(Final) vs Upgrade_lowrank (delayed), case mz, ', &
           & N_PROP_EQ, ' proposals'
      write(6,'(A7,1X,A5,1X,A8,1X,A8,1X,A10,1X,A22,1X,A22,1X,A11,1X,A11,1X,A6,1X,A4)') 'N_delay', 'acc', 'tog_mis', &
           & 'fld_same', '|dPhase|', 'next ranf (Upgrade2)', 'next ranf (lowrank)', 'G_cur_err', 'G_err', 'n_pend', ''
      mode = "Final"
      call Get_seed_Len(k_seed)
      allocate(rng_state(k_seed))
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

         ! ---- pass B: Upgrade_lowrank (delayed), same RNG state and fields ----
         nsigma%f = f_save
         call random_seed(put=rng_state)
         call Lowrank_init('delayed', N_DELAY_LIST(id))
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
         write(6,'(I7,1X,I5,1X,I8,1X,L8,1X,ES10.3,1X,F22.17,1X,F22.17,1X,ES11.4,1X,ES11.4,1X,I6,1X,A4)') &
              & N_DELAY_LIST(id), n_acc, n_mis, same_f, dphase, next_a, next_b, gerr_cur, gerr, n_pend_end, flag
         deallocate(G0, Ga, Gb, Aex, Gcur)
      enddo
      call clear_case()
   end subroutine run_kernel_equivalence

end module lowrank_kernel_equiv_section


program lowrank_unit
   use lowrank_unit_harness
   use lowrank_eligibility_section, only: run_eligibility_section
   use lowrank_kernel_equiv_section, only: run_kernel_equivalence
   implicit none
   integer :: ic, is, id, idir, n_fail
   real(dp) :: rerr, gerr, drift, max_drift, tol, drift_rel
   character(len=4), parameter :: dirs(2) = [character(len=4) :: 'up', 'down']
   character(len=4) :: flag
   logical :: info_existed
   integer :: u_info

   call Fields_init(1.d0)
   n_fail = 0
   ! Lowrank_init (delayed candidate, eligibility section) appends to the file 'info'.
   inquire(file='info', exist=info_existed)
   max_drift = 0.d0

   write(6,'(A,ES9.2,A,ES9.2,A)') '36-lowrank-unit: TOL_UNIT = ', TOL_UNIT, ', TOL_ILL = ', TOL_ILL, &
        & ' (mz_illcond; the fixed bound drift <= TOL_UNIT/2 is not asserted there)'
   write(6,'(A)') '  std_drift/ce = std_drift / (cond(1+A) * epsilon(1d0)); asserted <= 10 for every case'
   write(6,'(A10,1X,A10,1X,A4,1X,A7,1X,A13,1X,A13,1X,A13,1X,A10,1X,A12,1X,A7,1X,A6,1X,A6,1X,A4)') &
        & 'case', 'scheme', 'dir', 'N_delay', 'max_ratio_err', 'max_G_err', 'std_drift', &
        & 'cond(1+A)', 'std_drift/ce', 'acc', 'condQ', 'condS', ''
   do ic = 1, N_CASES
      do idir = 1, 2
         do is = 1, N_SCHEMES
            do id = 1, N_NDELAY
               call run_case(ic, is, N_DELAY_LIST(id), trim(dirs(idir)), rerr, gerr, drift)
               flag = ''
               tol = TOL_UNIT
               if (ic == CASE_MZ_ILLCOND) tol = TOL_ILL
               if (rerr > tol .or. gerr > tol .or. .not. (rerr == rerr .and. gerr == gerr)) then
                  flag = 'FAIL'
                  n_fail = n_fail + 1
               endif
               ! TOL_UNIT must stay >= 2x the standard path's own drift vs exact. Not asserted
               ! for mz_illcond: the double-precision "exact" inverse is unreliable at cond ~1e8.
               if (ic /= CASE_MZ_ILLCOND) then
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
               write(6,'(A10,1X,A10,1X,A4,1X,I7,1X,ES13.4,1X,ES13.4,1X,ES13.4,1X,ES10.2,1X,ES12.2,1X,F7.3,1X,F6.2,1X,F6.2,1X,A4)') &
                    & case_name(ic), scheme_name(is), dirs(idir), N_DELAY_LIST(id), rerr, gerr, drift, &
                    & cond_1pA_last, drift_rel, real(n_accepted, dp) / real(n_proposals, dp), cond_Q_max, cond_S_max, flag
            enddo
         enddo
      enddo
   enddo
   write(6,'(A,ES10.3)') 'max std_drift_vs_exact (well-conditioned cases) = ', max_drift
   call run_eligibility_section(n_fail)
   call run_kernel_equivalence(n_fail)
   if (.not. info_existed) then
      inquire(file='info', exist=info_existed)
      if (info_existed) then
         open(newunit=u_info, file='info', status='old')
         close(u_info, status='delete')
      endif
   endif
   if (n_fail > 0) then
      write(6,'(A,I0,A)') '36-lowrank-unit: ', n_fail, ' row(s) FAILED'
      stop 1
   endif
   write(6,'(A)') '36-lowrank-unit: all rows passed'
end program lowrank_unit

! Test of the device stabilisation (P2.4 step 1d), registered only when the build has -DALF_GPU.
! Usage: 52-gpu-stab <mode> [Nwrap]; one ctest per mode (see CMakeLists.txt), every mode in its own directory because
! Ham_set and the mirror write ./info.  The model is Hubbard on a 4x4 square lattice, beta = 2, Dtau = 0.1 (Ltrot = 20).
!
!   selftest [Nwrap]  Setup_gpu with ALF_GPU_SELFTEST=3 in the environment: the self-test of gpu_stab_mod runs the
!                     left and right chains (WRAPUL/WRAPUR + decompose) and CGR NVAR = 1, 2 at every stabilisation point
!                     against the upstream host routines, and the precision reduction against COMPARE.  With
!                     ALF_GPU_FAULT_STAB=perm|phase|nan it must fail (negative controls).
!   phase             determinant phase of the LU factors (no model): large N where the raw product underflows
!                     (n = 3200, diagonal moduli exp(-0.6) / exp(+0.6); n = 40, moduli 1e-20 / 1e+20), the normalised phase
!                     of k_det_phase against the analytic one, the selection of gpu_lu_phase_select bitwise against
!                     LU_PHASE_SAFE, a case where the raw product is usable (kept bit for bit), a zero pivot, and a negative
!                     control (the raw-product loop cannot give a unit phase).
!   refdet            determinant phase from the QR reflectors (no model; the formula of the tracked det(U), k_det_qr): random
!                     complex A of order 4, 50, 200, 1000 against the phase of the LU (LU_PHASE_SAFE), a random unitary matrix,
!                     columns with tau = 0, a complex last reflector, and a negative control (one reflector factor dropped).
!   refdet3200        the same formula at n = 3200 (diagonal A with moduli exp(-0.6) / exp(+0.6): the raw product of the R_ii
!                     is not representable), modulus 1 and the analytic phase.
!   untracked         Setup_gpu, then slot udvl of flavour 1 is marked untracked and one CGR is called.  With ALF_GPU_DEBUG=1 and
!                     ALF_GPU_CGR_LU=0 or 2 it must stop with "untracked U in CGR"; with ALF_GPU_CGR_LU=1 (LU) it must pass.
! This program should exit successfully (return code 0) in the modes that are expected to pass.
Program Test_gpu_stab
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: Clear_Ham_state, Ltrot
#ifdef MPI
  Use Hamiltonian_main, only: Set_Group_Comm
#endif
  Use QMC_runtime_var, only: set_QMC_runtime_default_var, set_ham_name, set_Nwrap, set_Ltau, &
       &                     set_default_global_tau_sampling, set_NSweep
  Use qmc_setup_mod, only: Setup_model, Setup_fields, Setup_update_schemes, Setup_gpu, Finalize_gpu
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc
  Use gpu_state_mod, only: gpu_active
  Use gpu_stab_mod, only: gpu_stab_on, gpu_lu_phase_select, k_det_phase, k_det_qr, gpu_stab_cgr, gpu_stab_untrack
  Use Hamiltonian_main, only: Ndim, N_FL, Calc_Fl_map
  Use linalg_mod, only: la_zgeqrf, la_zgetrf
  Use MyMats, only: LU_PHASE_SAFE
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit, int64
  Implicit none

  Type (qmc_walker_t) :: walker
  Character (len=64) :: mode, arg
  Integer :: n_fail, N_op, nwrap_arg, ios
  Real (Kind=Kind(0.d0)), parameter :: pi = 3.14159265358979323846d0
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Call Set_Group_Comm(MPI_COMM_WORLD)
#endif
  n_fail = 0
  Call get_command_argument(1, mode)
  Call get_command_argument(2, arg)

  Select case (trim(mode))
  Case ('selftest')
     nwrap_arg = 5
     If (len_trim(arg) > 0) Read(arg,*,iostat=ios) nwrap_arg
     Call setup_model_for_run(nwrap_arg, N_op)
     Call Walker_layout(walker, Ltrot, nwrap_arg)
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')           ! sets up the stabilisation and runs its self-test (ALF_GPU_SELFTEST=3); a failure ends the run
     If (.not. gpu_active) then
        Write(error_unit,*) 'GPU path inactive'
        n_fail = n_fail + 1
     Else If (.not. gpu_stab_on) then
        Write(error_unit,*) 'GPU stabilisation inactive'
        n_fail = n_fail + 1
     Endif
     Call Finalize_gpu()
     Call Walker_dealloc(walker)
     Call Clear_Ham_state()
  Case ('phase')
     Call phase_tests()
  Case ('refdet')
     Call refdet_tests()
  Case ('refdet3200')
     Call refdet_3200()
  Case ('untracked')
     Call setup_model_for_run(5, N_op)
     Call Walker_layout(walker, Ltrot, 5)
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')
     If (.not. gpu_stab_on) then
        Write(error_unit,*) 'GPU stabilisation inactive'
        n_fail = n_fail + 1
     Else
        Call untracked_cgr()
     Endif
     Call Finalize_gpu()
     Call Walker_dealloc(walker)
     Call Clear_Ham_state()
  Case default
     Write(error_unit,*) 'usage: 52-gpu-stab selftest [Nwrap] | phase | refdet | refdet3200 | untracked'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  End Select

  If (n_fail > 0) then
     Write(error_unit,*) 'gpu-stab test (', trim(mode), '): ', n_fail, ' failure(s)'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) 'gpu-stab test passed: ', trim(mode)

contains

  !> Hubbard model, parameters and seeds written here, model and fields set up as main does.
  Subroutine setup_model_for_run(nwrap, N_op)
    Integer, intent(in) :: nwrap
    Integer, intent(out) :: N_op
    Integer :: u, i
    Open(newunit=u, file='parameters', status='replace', action='write')
    Write(u,'(A)') '&VAR_lattice'
    Write(u,'(A)') 'L1 = 4'
    Write(u,'(A)') 'L2 = 4'
    Write(u,'(A)') 'Lattice_type = "Square"'
    Write(u,'(A)') 'Model = "Hubbard"'
    Write(u,'(A)') '/'
    Write(u,'(A)') '&VAR_Model_Generic'
    Write(u,'(A)') 'Checkerboard = .T.'
    Write(u,'(A)') 'Symm = .T.'
    Write(u,'(A)') 'N_SUN = 2'
    Write(u,'(A)') 'N_FL = 1'
    Write(u,'(A)') 'Dtau = 0.1d0'
    Write(u,'(A)') 'Beta = 2.d0'
    Write(u,'(A)') 'Projector = .F.'
    Write(u,'(A)') '/'
    Write(u,'(A)') '&VAR_Hubbard'
    Write(u,'(A)') 'Mz = .T.'
    Write(u,'(A)') 'Continuous = .F.'
    Write(u,'(A)') 'ham_T = 1.d0'
    Write(u,'(A)') 'ham_chem = 0.d0'
    Write(u,'(A)') 'ham_U = 4.d0'
    Write(u,'(A)') '/'
    Close(u)
    Open(newunit=u, file='seeds', status='replace', action='write')
    Do i = 1, 12
       Write(u,'(I10)') 1000003 * i + 17
    Enddo
    Close(u)
    Call set_QMC_runtime_default_var()
    Call set_ham_name('Hubbard')
    Call set_Nwrap(nwrap)
    Call set_NSweep(1)          ! the start-up check (KNOWN_ISSUES H3) needs at least one sweep
    Call set_Ltau(0)
    Call Setup_model(N_op)
    Call Setup_fields(N_op)
    Call Setup_update_schemes(N_op, 'info')
  end Subroutine setup_model_for_run

  !> Determinant phase tests (no model, no mirror)
  Subroutine phase_tests()
    Call phase_case('n = 3200, moduli exp(-0.6) / exp(+0.6)', 3200, exp(-0.6d0), exp(0.6d0), .true.)
    Call phase_case('n = 40, moduli 1e-20 / 1e+20', 40, 1.d-20, 1.d+20, .true.)
    Call phase_case('n = 12, moduli 0.7 / 1.3 (raw product usable)', 12, 0.7d0, 1.3d0, .false.)
    Call phase_zero_pivot()
  end Subroutine phase_tests

  !> Diagonal LU array of order n: |A(i,i)| = lo for i <= n/2 and hi otherwise, random phases, a few pivot swaps.
  Subroutine phase_case(label, n, lo, hi, underflow)
    Character (len=*), intent(in) :: label
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(in) :: lo, hi
    Logical, intent(in) :: underflow
    Complex (Kind=Kind(0.d0)), Allocatable :: A(:,:)
    Integer, Allocatable :: ipiv(:)
    Complex (Kind=Kind(0.d0)) :: cs(8), raw, p, sel, pref, ph0, analytic, det
    Real (Kind=Kind(0.d0)) :: turns, t, md
    Integer(int64) :: seed
    Integer :: i, nsw
    Logical :: ok_a, ok_b, ok_c, ok_d

    allocate(A(n,n), ipiv(n))
    A = (0.d0, 0.d0)
    seed = 20261126_int64
    turns = 0.d0
    nsw = 0
    do i = 1, n
       ipiv(i) = i
    enddo
    do i = 3, n - 4, max(n/5, 1)             ! a few pivot swaps
       ipiv(i) = i + 3
       nsw = nsw + 1
    enddo
    do i = 1, n
       t = rnd01(seed)
       if (i <= n/2) then
          md = lo
       else
          md = hi
       endif
       A(i,i) = md * cmplx(cos(2.d0*pi*t), sin(2.d0*pi*t), kind(0.d0))
       turns = turns + t
       turns = turns - floor(turns)         ! keeps the accumulated angle small: rounding error ~ n eps, not n eps * 2 pi n
    enddo
    turns = turns + 0.5d0*dble(nsw)          ! each pivot swap flips the sign
    turns = turns - floor(turns)
    analytic = cmplx(cos(2.d0*pi*turns), sin(2.d0*pi*turns), kind(0.d0))

    cs = (0.d0, 0.d0)
    call k_det_phase(n, A, n, ipiv, cs)
    p = cs(4); raw = cs(5)
    sel = gpu_lu_phase_select(raw, p)
    ph0 = conjg(sel); ph0 = ph0/abs(ph0)
    ! p_ref: the normalised product as LU_PHASE_SAFE forms it (a raw value 0 is never usable, so it returns that product)
    pref = LU_PHASE_SAFE((0.d0, 0.d0), A, ipiv, n)

    ! (a) the raw product has underflowed (or is usable, in the last case)
    if (underflow) then
       ok_a = .not. (abs(raw) > 1.d-250 .and. abs(raw) < 1.d250)
    else
       ok_a = (abs(raw) > 1.d-250 .and. abs(raw) < 1.d250)
    endif
    ! (b) the normalised phase has modulus 1 and equals the analytic phase and the LU_PHASE_SAFE product (to rounding: the
    !     order of the rounding of the two loops is the compiler's, so these two are not compared bit for bit)
    ok_b = abs(abs(p) - 1.d0) < 1.d-12 .and. abs(p - analytic) < 1.d-11 .and. abs(p - pref) < 1.d-12
    ! (c) the selection rule is the one of LU_PHASE_SAFE, bit for bit, for the raw product of the kernel and for raw values that
    !     are usable with the right phase (kept), with a phase off by 1e-6 (replaced), out of range, NaN and Inf; the phase
    !     used by CGR has modulus 1; with a usable raw product the raw value is kept
    ok_c = abs(abs(ph0) - 1.d0) < 1.d-12
    ok_c = ok_c .and. same_select(raw, pref, A, ipiv, n)
    ok_c = ok_c .and. same_select(1.d3*pref*cmplx(cos(1.d-10), sin(1.d-10), kind(0.d0)), pref, A, ipiv, n)
    ok_c = ok_c .and. same_select(1.d3*pref*cmplx(cos(1.d-6), sin(1.d-6), kind(0.d0)), pref, A, ipiv, n)
    ok_c = ok_c .and. same_select(1.d-300*pref, pref, A, ipiv, n)
    ok_c = ok_c .and. same_select(1.d+300*pref, pref, A, ipiv, n)
    ok_c = ok_c .and. same_select(transfer([9221120237041090560_int64, 0_int64], pref), pref, A, ipiv, n)   ! NaN
    ok_c = ok_c .and. same_select(transfer([9218868437227405312_int64, 0_int64], pref), pref, A, ipiv, n)   ! Inf
    if (.not. underflow) ok_c = ok_c .and. (real(sel) == real(raw) .and. aimag(sel) == aimag(raw))
    ! (d) negative control: the raw-product loop (upstream det_C before the fix) cannot give a phase of modulus 1 where (a) holds
    det = (1.d0, 0.d0)
    do i = 1, n
       det = det*A(i,i)
    enddo
    if (underflow) then
       ok_d = .not. (abs(det) > 1.d-250 .and. abs(det) < 1.d250)    ! det/abs(det) would be 0/0
    else
       ok_d = .true.
    endif
    write(*,'(A,A,A,L1,A,L1,A,L1,A,L1)') ' phase test ', label, ': raw-product state as expected ', ok_a, ', phase = analytic ', ok_b, &
         & ', selection rule = LU_PHASE_SAFE bitwise ', ok_c, ', negative control seen ', ok_d
    write(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') '    |p| - 1 = ', abs(p) - 1.d0, ', |p - analytic| = ', abs(p - analytic), &
         & ', |p - p_ref| = ', abs(p - pref)
    if (.not. (ok_a .and. ok_b .and. ok_c .and. ok_d)) n_fail = n_fail + 1
    deallocate(A, ipiv)
  end Subroutine phase_case

  !> gpu_lu_phase_select(z_raw, p_ref) equals LU_PHASE_SAFE(z_raw, A, ipiv, n) bit for bit (also NaN results)
  Logical function same_select(z_raw, pref, A, ipiv, n)
    Complex (Kind=Kind(0.d0)), intent(in) :: z_raw, pref
    Integer, intent(in) :: n, ipiv(:)
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:)
    Complex (Kind=Kind(0.d0)) :: s, l
    Real (Kind=Kind(0.d0)) :: v(4)
    s = gpu_lu_phase_select(z_raw, pref)
    l = LU_PHASE_SAFE(z_raw, A, ipiv, n)
    v = [real(s), aimag(s), real(l), aimag(l)]
    same_select = (transfer(v(1), 0_int64) == transfer(v(3), 0_int64)) .and. (transfer(v(2), 0_int64) == transfer(v(4), 0_int64))
  end function same_select

  !> A zero pivot: the kernel marks it (p = 0), the selection returns the raw product, as LU_PHASE_SAFE
  Subroutine phase_zero_pivot()
    Complex (Kind=Kind(0.d0)) :: A(6,6), cs(8), sel, lps
    Integer :: ipiv(6), i
    Logical :: ok
    A = (0.d0, 0.d0)
    do i = 1, 6
       A(i,i) = cmplx(0.5d0*i, 0.25d0*i, kind(0.d0))
       ipiv(i) = i
    enddo
    A(3,3) = (0.d0, 0.d0)
    cs = (0.d0, 0.d0)
    call k_det_phase(6, A, 6, ipiv, cs)
    sel = gpu_lu_phase_select(cs(5), cs(4))
    lps = LU_PHASE_SAFE(cs(5), A, ipiv, 6)
    ok = (abs(cs(4)) == 0.d0) .and. (real(sel) == real(lps) .and. aimag(sel) == aimag(lps)) .and. &
         & (real(sel) == real(cs(5)) .and. aimag(sel) == aimag(cs(5)))
    write(*,'(A,L1)') ' phase test zero pivot: marked, raw product returned as LU_PHASE_SAFE: ', ok
    if (.not. ok) n_fail = n_fail + 1
  end Subroutine phase_zero_pivot

  !> CGR on the identity states with udvl marked untracked: stops with "untracked U in CGR" if the phase is tracked (modes 0, 2)
  Subroutine untracked_cgr()
    Complex (Kind=Kind(0.d0)), Allocatable :: G(:,:,:)
    Complex (Kind=Kind(0.d0)) :: ph
    Allocate(G(Ndim,Ndim,N_FL))
    G = (0.d0, 0.d0)
    Call gpu_stab_untrack(1, 1)
    Call gpu_stab_cgr(1, Calc_Fl_map(1), 1, G, ph, .false.)
    ! reached only in mode 1 (LU phase): the phase is that of the identity states
    If (abs(ph - (1.d0, 0.d0)) > 1.d-12) n_fail = n_fail + 1
    Deallocate(G)
  end Subroutine untracked_cgr

  !> phase of det(A) from the QR of a copy: the tracked-det formula (k_det_qr); a = copy of A on entry, overwritten
  Complex (Kind=Kind(0.d0)) function qr_phase(n, a, tau_out)
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: a(n,n)
    Complex (Kind=Kind(0.d0)), intent(out) :: tau_out(n)
    Complex (Kind=Kind(0.d0)) :: wq(1), d(1)
    Complex (Kind=Kind(0.d0)), Allocatable :: work(:)
    Integer :: info
    Call la_zgeqrf(n, n, a, n, tau_out, wq, -1, info)
    Allocate(work(max(1, int(real(wq(1))))))
    Call la_zgeqrf(n, n, a, n, tau_out, work, size(work), info)
    d(1) = (0.d0, 0.d0)
    Call k_det_qr(n, a, n, tau_out, d)
    qr_phase = d(1)
    Deallocate(work)
  end function qr_phase

  !> phase of det(A) from the LU (LU_PHASE_SAFE normalised product, A is overwritten)
  Complex (Kind=Kind(0.d0)) function lu_phase(n, a)
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: a(n,n)
    Integer :: ipiv(n), info
    Call la_zgetrf(n, n, a, n, ipiv, info)
    lu_phase = LU_PHASE_SAFE((0.d0, 0.d0), a, ipiv, n)
  end function lu_phase

  !> det of the reflector with scalar tau, as the formula says: -tau/conj(tau) (1 for tau = 0)
  Complex (Kind=Kind(0.d0)) function hdet(tau)
    Complex (Kind=Kind(0.d0)), intent(in) :: tau
    If (abs(tau) == 0.d0) then
       hdet = (1.d0, 0.d0)
    else
       hdet = -tau/conjg(tau)
    endif
  end function hdet

  !> Random matrix, entries re + i im uniform in [-1,1)
  Subroutine rnd_matrix(n, a, seed)
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(out) :: a(n,n)
    Integer(int64), intent(inout) :: seed
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: re, im
    Do j = 1, n
       Do i = 1, n
          re = 2.d0*rnd01(seed) - 1.d0; im = 2.d0*rnd01(seed) - 1.d0
          a(i,j) = cmplx(re, im, kind(0.d0))
       Enddo
    Enddo
  end Subroutine rnd_matrix

  !> One line and one verdict per case: |d_qr - p_lu| (or the analytic phase) against tol, |d_qr| - 1 within 1e-12
  Subroutine report_det(label, n, d, pref, tol)
    Character (len=*), intent(in) :: label
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: d, pref
    Real (Kind=Kind(0.d0)), intent(in) :: tol
    Real (Kind=Kind(0.d0)) :: dif, md
    Logical :: ok
    dif = abs(d - pref); md = abs(abs(d) - 1.d0)
    ok = (dif <= tol) .and. (md <= 1.d-12)
    write(*,'(A,A,A,I0,A,ES10.3,A,ES10.3,A,ES8.1,A,L1)') ' refdet ', label, ' n = ', n, ': |phase(QR reflectors) - phase(ref)| = ', dif, &
         & ', ||d| - 1| = ', md, ' (tol ', tol, ') ok = ', ok
    If (.not. ok) n_fail = n_fail + 1
  end Subroutine report_det

  !> The reflector formula of the tracked det(U) against the LU phase and analytic cases (no model)
  Subroutine refdet_tests()
    Integer, parameter :: ns(4) = [4, 50, 200, 1000]
    Complex (Kind=Kind(0.d0)), Allocatable :: a(:,:), b(:,:), tau(:)
    Complex (Kind=Kind(0.d0)) :: d, pl, dnl
    Real (Kind=Kind(0.d0)) :: dif, dneg
    Integer(int64) :: seed
    Integer :: k, n, i
    Logical :: ok
    seed = 20261127_int64
    ! random complex A against the LU phase.  The phase of det A from either factorisation is accurate to about
    ! N eps cond(A); the tolerance is 1e-13 N (cond of a random square matrix is O(N)).
    Do k = 1, size(ns)
       n = ns(k)
       Allocate(a(n,n), b(n,n), tau(n))
       Call rnd_matrix(n, a, seed)
       b = a
       pl = lu_phase(n, b)
       b = a
       d = qr_phase(n, b, tau)
       Call report_det('random complex A', n, d, pl, 1.d-13*n)
       ! negative control: dropping the determinant of the last reflector (tau(n) /= 0 for a complex last column) must be seen
       dnl = d/hdet(tau(n))
       dneg = abs(dnl - pl)
       ok = (abs(tau(n)) > 0.d0) .and. (dneg > 1.d-3)
       write(*,'(A,I0,A,ES10.3,A,L1)') ' refdet negative control (last reflector factor dropped) n = ', n, ': difference ', dneg, ' seen = ', ok
       If (.not. ok) n_fail = n_fail + 1
       Deallocate(a, b, tau)
    Enddo
    ! random unitary matrix (Q of a random matrix), against the LU phase (cond = 1: tolerance 1e-12)
    Do k = 2, 3
       n = ns(k)
       Allocate(a(n,n), b(n,n), tau(n))
       Call rnd_matrix(n, a, seed)
       Call unitary_from_qr(n, a)
       b = a
       pl = lu_phase(n, b)
       b = a
       d = qr_phase(n, b, tau)
       Call report_det('random unitary', n, d, pl, 1.d-12)
       Deallocate(a, b, tau)
    Enddo
    ! columns with tau = 0: real diagonal, a few negative entries (every reflector is the identity, R_ii = A_ii)
    n = 50
    Allocate(a(n,n), tau(n))
    a = (0.d0, 0.d0)
    Do i = 1, n
       a(i,i) = cmplx(1.d0 + 0.01d0*i, 0.d0, kind(0.d0))
       If (mod(i,7) == 0) a(i,i) = -a(i,i)
    Enddo
    d = qr_phase(n, a, tau)
    ok = all(abs(tau) == 0.d0)
    write(*,'(A,L1)') ' refdet tau = 0 columns: all reflector scalars are 0: ', ok
    If (.not. ok) n_fail = n_fail + 1
    Call report_det('real diagonal, tau = 0, 7 negative entries', n, d, cmplx(real((-1)**(n/7), kind(0.d0)), 0.d0, kind(0.d0)), 1.d-14)
    Deallocate(a, tau)
    ! complex last reflector: leading block real positive diagonal, last entry exp(i 0.7): the only nontrivial reflector is the last
    n = 3
    Allocate(a(n,n), tau(n))
    a = (0.d0, 0.d0)
    a(1,1) = (2.d0, 0.d0); a(2,2) = (3.d0, 0.d0); a(3,3) = cmplx(cos(0.7d0), sin(0.7d0), kind(0.d0))
    d = qr_phase(n, a, tau)
    ok = (abs(tau(1)) == 0.d0 .and. abs(tau(2)) == 0.d0 .and. abs(tau(3)) > 0.d0)
    write(*,'(A,L1)') ' refdet complex last reflector: tau(1:2) = 0, tau(3) /= 0: ', ok
    If (.not. ok) n_fail = n_fail + 1
    Call report_det('complex last reflector', n, d, cmplx(cos(0.7d0), sin(0.7d0), kind(0.d0)), 1.d-14)
    ! negative control on the same matrix: without det H_3 the phase is not exp(i 0.7)
    dneg = abs(d/hdet(tau(3)) - cmplx(cos(0.7d0), sin(0.7d0), kind(0.d0)))
    ok = dneg > 1.d-3
    write(*,'(A,ES10.3,A,L1)') ' refdet negative control (complex last reflector factor dropped): difference ', dneg, ' seen = ', ok
    If (.not. ok) n_fail = n_fail + 1
    Deallocate(a, tau)
  end Subroutine refdet_tests

  !> a (random) -> its unitary factor Q by ZGEQRF + back-substitution of the reflectors on the unit matrix (host, test only)
  Subroutine unitary_from_qr(n, a)
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(inout) :: a(n,n)
    Complex (Kind=Kind(0.d0)), Allocatable :: tau(:), work(:), q(:,:), v(:)
    Complex (Kind=Kind(0.d0)) :: wq(1)
    Integer :: info, i, j
    Complex (Kind=Kind(0.d0)) :: s
    Allocate(tau(n), q(n,n), v(n))
    Call la_zgeqrf(n, n, a, n, tau, wq, -1, info)
    Allocate(work(max(1, int(real(wq(1))))))
    Call la_zgeqrf(n, n, a, n, tau, work, size(work), info)
    ! Q = H_1 ... H_n applied to the unit matrix, H_i = I - tau_i v v^H, applied from i = n down to 1
    q = (0.d0, 0.d0)
    Do i = 1, n
       q(i,i) = (1.d0, 0.d0)
    Enddo
    Do i = n, 1, -1
       v = (0.d0, 0.d0)
       v(i) = (1.d0, 0.d0)
       v(i+1:n) = a(i+1:n,i)
       Do j = 1, n
          s = dot_product(v(i:n), q(i:n,j))
          q(i:n,j) = q(i:n,j) - tau(i)*s*v(i:n)
       Enddo
    Enddo
    a = q
    Deallocate(tau, work, q, v)
  end Subroutine unitary_from_qr

  !> n = 3200, diagonal A: |A(i,i)| = exp(-0.6) for i <= n/2, exp(+0.6) otherwise, random phases.  The raw product of the R_ii
  !> underflows (1e-417 after n/2 factors); the tracked phase has modulus 1 within 1e-12 and equals the analytic phase.
  Subroutine refdet_3200()
    Integer, parameter :: n = 3200
    Complex (Kind=Kind(0.d0)), Allocatable :: a(:,:), tau(:)
    Complex (Kind=Kind(0.d0)) :: d, analytic
    Real (Kind=Kind(0.d0)) :: turns, t, md, lg
    Integer(int64) :: seed
    Integer :: i
    Logical :: ok
    Allocate(a(n,n), tau(n))
    a = (0.d0, 0.d0)
    seed = 20261127_int64
    turns = 0.d0
    Do i = 1, n
       t = rnd01(seed)
       If (i <= n/2) then
          md = exp(-0.6d0)
       else
          md = exp(0.6d0)
       Endif
       a(i,i) = md * cmplx(cos(2.d0*pi*t), sin(2.d0*pi*t), kind(0.d0))
       turns = turns + t
       turns = turns - floor(turns)
    Enddo
    analytic = cmplx(cos(2.d0*pi*turns), sin(2.d0*pi*turns), kind(0.d0))
    d = qr_phase(n, a, tau)
    ! the raw product of the first n/2 diagonal entries of R is exp(-960) = 1e-417, below the smallest double (4.9e-324): a raw
    ! product of the R_ii cannot be represented (decided in the log domain: the order of a compiled product loop is not fixed)
    lg = 0.d0
    Do i = 1, n/2
       lg = lg + log(abs(a(i,i)))
    Enddo
    ok = lg < log(tiny(1.d0))
    write(*,'(A,L1,A,F9.2,A,F9.2,A)') ' refdet n = 3200: the raw product of the first n/2 R_ii is not representable: ', ok, &
         & ' (log = ', lg, ', log(tiny) = ', log(tiny(1.d0)), ')'
    If (.not. ok) n_fail = n_fail + 1
    Call report_det('diagonal 3200, moduli exp(-0.6)/exp(0.6)', n, d, analytic, 1.d-11)
    Deallocate(a, tau)
  end Subroutine refdet_3200

  !> uniform number in [0,1)
  Real (Kind=Kind(0.d0)) function rnd01(s)
    Integer(int64), intent(inout) :: s
    s = mod(1103515245_int64*s + 12345_int64, 2147483648_int64)
    rnd01 = real(s, kind(0.d0))/2147483648.d0
  end function rnd01

end Program Test_gpu_stab

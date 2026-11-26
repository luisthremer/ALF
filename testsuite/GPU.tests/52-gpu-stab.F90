! Test of the device stabilisation (P2.4 step 1d), registered only when the build has -DALF_GPU.
! Usage: 52-gpu-stab <mode> [Nwrap]; one ctest per mode (see CMakeLists.txt), every mode in its own directory because
! Ham_set and the mirror write ./info.  The model is Hubbard on a 4x4 square lattice, beta = 2, Dtau = 0.1 (Ltrot = 20).
!
!   selftest [Nwrap]  stabilisation set-up with ALF_GPU_SELFTEST=3 in the environment: the self-test of gpu_stab_mod runs the
!                     left and right chains (WRAPUL/WRAPUR + decompose) and CGR NVAR = 1, 2 at every stabilisation point
!                     against the upstream host routines, and the precision reduction against COMPARE.  With
!                     ALF_GPU_FAULT_STAB=perm|phase|nan it must fail (negative controls).
!   phase             determinant phase of the LU factors (no model): large N where the raw product underflows
!                     (n = 3200, diagonal moduli exp(-0.6) / exp(+0.6); n = 40, moduli 1e-20 / 1e+20), the normalised phase
!                     of k_det_phase against the analytic one, the selection of gpu_lu_phase_select bitwise against
!                     LU_PHASE_SAFE, a case where the raw product is usable (kept bit for bit), a zero pivot, and a negative
!                     control (the raw-product loop cannot give a unit phase).
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
  Use gpu_stab_mod, only: gpu_stab_on, gpu_stab_setup, gpu_stab_hook, gpu_stab_finalize, gpu_lu_phase_select, k_det_phase
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
     Call Setup_gpu(walker, 'info')
     If (.not. gpu_active) then
        Write(error_unit,*) 'GPU path inactive'
        n_fail = n_fail + 1
     Else
        Call gpu_stab_setup(walker)
        If (.not. gpu_stab_on) then
           Write(error_unit,*) 'GPU stabilisation inactive'
           n_fail = n_fail + 1
        Else
           Call gpu_stab_hook()            ! runs the self-test (ALF_GPU_SELFTEST=3); a failure ends the run
           Call gpu_stab_finalize()
        Endif
     Endif
     Call Finalize_gpu()
     Call Walker_dealloc(walker)
     Call Clear_Ham_state()
  Case ('phase')
     Call phase_tests()
  Case default
     Write(error_unit,*) 'usage: 52-gpu-stab selftest [Nwrap] | phase'
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

  !> uniform number in [0,1)
  Real (Kind=Kind(0.d0)) function rnd01(s)
    Integer(int64), intent(inout) :: s
    s = mod(1103515245_int64*s + 12345_int64, 2147483648_int64)
    rnd01 = real(s, kind(0.d0))/2147483648.d0
  end function rnd01

end Program Test_gpu_stab

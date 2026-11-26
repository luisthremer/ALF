! Unit tests of gpu_state_mod (P2.4 step 1b), registered only when the build has -DALF_GPU.
! Usage: 50-gpu-state <mode> [argument]; one ctest per mode (see CMakeLists.txt), every mode in its own directory
! because Ham_set and the mirror write ./info.  The model is Hubbard on a 4x4 square lattice, beta = 2, Dtau = 0.1.
!
!   tables [dense]    mirror set-up on the model; with ALF_GPU_SELFTEST=1 in the environment the set-up runs the table
!                     self-test (flat hopping and Op_V tables against Hop_mod_* and Op_mmultL/R/Wrapup/Wrapdo, walker
!                     pack/unpack, device kernel).  'dense': Checkerboard = .F. (one 16 x 16 hopping element).
!   roundtrip         bitwise round trip walker -> mirror -> walker of GR, udvst, udvl and udvr (D, or L with STABLOG)
!   inactive <text>   the set-up must leave the GPU path inactive, with a reason that contains <text>
!   refusals          gpu_supported refuses each runtime-detectable unsupported option, with the right reason
!   envreader         environment reader: ALF_GPU_DISABLE=1 is seen by both readers, ALF_GPU_FOO=1 only by /proc
!   plain             set-up only (used for the cases where the set-up has to stop the run)
! This program should exit successfully (return code 0) in the modes that are expected to pass.
Program Test_gpu_state
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: Clear_Ham_state, Set_model_for_testing, Setup_flavor_map, Ndim, N_FL, N_FL_eff, Ltrot, &
       &                      Calc_Fl_map
#ifdef MPI
  Use Hamiltonian_main, only: Set_Group_Comm
#endif
  Use Operator_mod, only: Operator
  Use QMC_runtime_var, only: set_QMC_runtime_default_var, set_ham_name, set_Nwrap, set_Ltau, get_Nwrap, &
       &                     set_Langevin, set_HMC, set_Global_moves, set_Global_tau_moves, set_N_Global_tau, &
       &                     set_Propose_S0, set_Sequential_MALA, set_Global_MALA_moves, set_Global_tau_MALA_moves, &
       &                     set_N_Global_tau_MALA, set_default_global_tau_sampling, set_NSweep
  Use qmc_setup_mod, only: Setup_model, Setup_fields, Setup_update_schemes, Setup_gpu, Finalize_gpu
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc
  Use UDV_State_mod, only: UDV_State
  Use gpu_state_mod, only: gpu_active, gpu_reason, gpu_supported, gpu_put_G, gpu_get_G, gpu_put_udv, gpu_get_udv, &
       &                   gpu_env_get, gpu_env_is, gpu_env_reader_name, rand_mat, rnd, same_bits_c, same_bits
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit, int64
  Implicit none

  Type (qmc_walker_t) :: walker
  Character (len=64) :: mode, arg
  Integer :: n_fail, N_op
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Call Set_Group_Comm(MPI_COMM_WORLD)
#endif
  n_fail = 0
  Call get_command_argument(1, mode)
  Call get_command_argument(2, arg)

  Select case (trim(mode))
  Case ('tables')
     Call setup_model_for_run(trim(arg) /= 'dense', N_op, .true.)
     Call Walker_layout(walker, Ltrot, get_Nwrap())
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')
     Call check(gpu_active, 'mirror active')
     Call check(info_has('ACTIVE (mirror; wraps, updates and stabilisation on the device path)'), 'info: ACTIVE line')
     Call check(info_has('GPU build: backend'), 'info: backend line')
     Call check(info_has('GPU debug env reader:'), 'info: env reader line')
     Call Finalize_gpu()
     Call check(.not. gpu_active, 'inactive after finalize')
     Call Walker_dealloc(walker)
  Case ('roundtrip')
     Call setup_model_for_run(.true., N_op, .true.)
     Call Walker_layout(walker, Ltrot, get_Nwrap())
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')
     Call check(gpu_active, 'mirror active')
     Call test_roundtrip()
     Call Finalize_gpu()
     Call Walker_dealloc(walker)
  Case ('inactive')
     Call setup_model_for_run(.true., N_op, .true.)
     Call Walker_layout(walker, Ltrot, get_Nwrap())
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')
     Call check(.not. gpu_active, 'mirror inactive')
     Call check(index(gpu_reason, trim(arg)) > 0, 'gpu_reason contains: '//trim(arg))
     Call check(info_has('INACTIVE, reason: '), 'info: INACTIVE line')
     Call check(info_has(trim(arg)), 'info: reason text')
     Call check(.not. info_has('ACTIVE (mirror built'), 'info: no ACTIVE line')
     Call Finalize_gpu()                           ! no-op on the inactive path
     Call Walker_dealloc(walker)
  Case ('plain')
     Call setup_model_for_run(.true., N_op, .true.)
     Call Walker_layout(walker, Ltrot, get_Nwrap())
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')
     Write(*,'(A,L1)') ' gpu_active = ', gpu_active
     Call Finalize_gpu()
     Call Walker_dealloc(walker)
  Case ('refusals')
     Call setup_model_for_run(.true., N_op, .false.)
     Call test_refusals()
  Case ('envreader')
     Call test_envreader()
  Case default
     Write(error_unit,*) 'usage: 50-gpu-state tables [dense] | roundtrip | inactive <text> | refusals | envreader | plain'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  End Select

  ! (envreader has no model; refusals ends with a toy projector model that has no trial wave functions to clear)
  If (trim(mode) /= 'envreader' .and. trim(mode) /= 'refusals') Call Clear_Ham_state()
  If (n_fail > 0) then
     Write(error_unit,*) 'gpu-state test (', trim(mode), '): ', n_fail, ' failure(s)'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) 'gpu-state test passed: ', trim(mode)

contains

  !> Hubbard model, parameters and seeds written here, model and fields set up as main does.  lock_schemes:
  !> Setup_update_schemes (locks the run-time settings) or only the default global-tau sampling (setters stay usable).
  Subroutine setup_model_for_run(checkerboard, N_op, lock_schemes)
    Logical, intent(in) :: checkerboard, lock_schemes
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
    If (checkerboard) then
       Write(u,'(A)') 'Checkerboard = .T.'
    else
       Write(u,'(A)') 'Checkerboard = .F.'
    endif
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
    Call set_Nwrap(5)
    Call set_NSweep(1)          ! the start-up check (KNOWN_ISSUES H3) needs at least one sweep
    Call set_Ltau(1)
    Call Setup_model(N_op)
    Call Setup_fields(N_op)
    If (lock_schemes) then
       Call Setup_update_schemes(N_op, 'info')
    else
       Call set_default_global_tau_sampling(N_op)
    endif
  end Subroutine setup_model_for_run

  !> True if the info file has a line containing text.
  Logical Function info_has(text)
    Character (len=*), intent(in) :: text
    Character (len=1024) :: line
    Integer :: u, ios
    info_has = .false.
    Open(newunit=u, file='info', status='old', action='read', iostat=ios)
    If (ios /= 0) return
    Do
       Read(u,'(A)',iostat=ios) line
       If (ios /= 0) exit
       If (index(line, text) > 0) then
          info_has = .true.
          exit
       endif
    Enddo
    Close(u)
  end Function info_has

  !> Fills a walker with a reproducible pattern (same seed, same bits).
  Subroutine fill(w, seed)
    Type (qmc_walker_t), intent(inout) :: w
    Integer(int64), intent(in) :: seed
    Integer(int64) :: s
    Integer :: nf, nfe, n
    s = seed
    Do nf = 1, N_FL
       Call rand_mat(w%GR(:,:,nf), s)
    Enddo
    Do nfe = 1, N_FL_eff
       Do n = 1, w%NSTM
          Call fill_udv(w%udvst(n,nfe), s)
       Enddo
       Call fill_udv(w%udvl(nfe), s)
       Call fill_udv(w%udvr(nfe), s)
    Enddo
  end Subroutine fill

  Subroutine fill_udv(u, s)
    Type (UDV_State), intent(inout) :: u
    Integer(int64), intent(inout) :: s
    Integer :: i
    Call rand_mat(u%U, s)
    Call rand_mat(u%V, s)
#if defined(STABLOG)
    Do i = 1, size(u%L)
       u%L(i) = rnd(s)
    Enddo
#else
    Do i = 1, size(u%D)
       u%D(i) = cmplx(rnd(s), rnd(s), kind(0.d0))
    Enddo
#endif
  end Subroutine fill_udv

  !> Number of arrays (U, V, D or L) of two UDV states that differ in at least one bit.
  Integer Function n_diff(a, b)
    Type (UDV_State), intent(in) :: a, b
    n_diff = 0
    If (.not. same_bits_c(size(a%U), a%U, b%U)) n_diff = n_diff + 1
    If (.not. same_bits_c(size(a%V), a%V, b%V)) n_diff = n_diff + 1
#if defined(STABLOG)
    If (.not. same_bits(size(a%L), a%L, b%L)) n_diff = n_diff + 1
#else
    If (.not. same_bits_c(size(a%D), a%D, b%D)) n_diff = n_diff + 1
#endif
  end Function n_diff

  !> walker -> mirror -> walker, bit for bit.  The reference walker has the same shapes as the real one.
  Subroutine test_roundtrip()
    Type (qmc_walker_t) :: ref
    Integer :: nf, nfe, n, nbad_g, nbad_udv
    Allocate(ref%GR(Ndim,Ndim,N_FL), ref%udvl(N_FL_eff), ref%udvr(N_FL_eff), ref%udvst(walker%NSTM,N_FL_eff))
    ref%NSTM = walker%NSTM
    Do nfe = 1, N_FL_eff
       Call ref%udvl(nfe)%alloc(Ndim)
       Call ref%udvr(nfe)%alloc(Ndim)
       Do n = 1, walker%NSTM
          Call ref%udvst(n,nfe)%alloc(Ndim)
       Enddo
    Enddo
    Call fill(ref, 20251126_int64)
    Call fill(walker, 20251126_int64)
    Call gpu_put_G(walker)
    Call gpu_put_udv(walker)
    Call fill(walker, 7_int64)                       ! overwrite the walker with different data
    Call check(.not. same_bits_c(size(ref%GR(:,:,1)), ref%GR(:,:,1), walker%GR(:,:,1)), 'poison differs from pattern')
    Call gpu_get_G(walker)
    Call gpu_get_udv(walker)
    nbad_g = 0; nbad_udv = 0
    Do nfe = 1, N_FL_eff
       nf = Calc_Fl_map(nfe)
       If (.not. same_bits_c(size(ref%GR(:,:,nf)), ref%GR(:,:,nf), walker%GR(:,:,nf))) nbad_g = nbad_g + 1
       Do n = 1, walker%NSTM
          nbad_udv = nbad_udv + n_diff(walker%udvst(n,nfe), ref%udvst(n,nfe))
       Enddo
       nbad_udv = nbad_udv + n_diff(walker%udvl(nfe), ref%udvl(nfe)) + n_diff(walker%udvr(nfe), ref%udvr(nfe))
    Enddo
    Call check(nbad_g == 0, 'GR round trip is bitwise')
    Call check(nbad_udv == 0, 'udvst, udvl, udvr round trip is bitwise')
    Write(*,'(A,I0,A,I0,A,I0,A,I0)') ' round trip: Ndim ', Ndim, ', N_FL_eff ', N_FL_eff, ', NSTM ', walker%NSTM, &
         & ', mismatching arrays: GR ', nbad_g
    Write(*,'(A,I0)') ' round trip: mismatching UDV arrays ', nbad_udv
    Do nfe = 1, N_FL_eff
       Call ref%udvl(nfe)%dealloc
       Call ref%udvr(nfe)%dealloc
       Do n = 1, walker%NSTM
          Call ref%udvst(n,nfe)%dealloc
       Enddo
    Enddo
    Deallocate(ref%GR, ref%udvl, ref%udvr, ref%udvst)
  end Subroutine test_roundtrip

  !> One unsupported option after the other (the set-up leaves the settings unlocked): gpu_supported must say no, with
  !> the reason of that option; with the option off again it must say yes.  Builds that refuse everything (tempering
  !> without PP, STAB1/2/3) refuse with their own reason.
  Subroutine test_refusals()
    Character (len=256) :: reason
    Character (len=64) :: build_reason
    Logical :: ok
    Type (Operator), dimension(:,:), allocatable :: Op_T, Op_V
    build_reason = ''
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
    build_reason = 'tempering'
#elif defined(STAB1) || defined(STAB2) || defined(STAB3)
    build_reason = 'STAB1/STAB2/STAB3'
#endif
    ok = gpu_supported(reason)
    If (len_trim(build_reason) == 0) then
       Call check(ok, 'baseline supported: '//trim(reason))
    else
       Call check(.not. ok .and. index(reason, trim(build_reason)) > 0, 'build-wide refusal: '//trim(reason))
    endif
    Call refuse('Langevin', 'Langevin/HMC')
    Call refuse('HMC', 'Langevin/HMC')
    Call refuse('Global_moves', 'global moves')
    Call refuse('Global_tau_moves', 'global moves')
    Call refuse('N_Global_tau', 'global moves')
    Call refuse('Propose_S0', 'Propose_S0')
    Call refuse('Sequential_MALA', 'MALA')
    Call refuse('Global_MALA_moves', 'MALA')
    Call refuse('Global_tau_MALA_moves', 'MALA')
    Call refuse('N_Global_tau_MALA', 'MALA')
    ! projector code: a toy model flagged as projector (the refusal comes before any operator is looked at)
    Call Clear_Ham_state()
    Allocate(Op_T(0,2), Op_V(0,2))
    Call Set_model_for_testing(4, 2, 1, 10, 0, .true., .false., Op_T, Op_V)
    Call Setup_flavor_map()
    ok = gpu_supported(reason)
    If (len_trim(build_reason) == 0) then
       Call check(.not. ok .and. index(reason, 'projector') > 0, 'projector refused: '//trim(reason))
    else
       Call check(.not. ok .and. index(reason, trim(build_reason)) > 0, 'projector (build-wide): '//trim(reason))
    endif
  end Subroutine test_refusals

  Subroutine refuse(option, expected)
    Character (len=*), intent(in) :: option, expected
    Character (len=256) :: reason
    Character (len=64) :: build_reason
    Logical :: ok
    build_reason = expected
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
    build_reason = 'tempering'
#elif defined(STAB1) || defined(STAB2) || defined(STAB3)
    build_reason = 'STAB1/STAB2/STAB3'
#endif
    Call switch(option, .true.)
    ok = gpu_supported(reason)
    Call check(.not. ok, 'refused: '//option)
    Call check(index(reason, trim(build_reason)) > 0, 'reason of '//option//': '//trim(reason))
    Write(*,'(A,A,A,A)') ' refusal ', option, ': ', trim(reason)
    Call switch(option, .false.)
  end Subroutine refuse

  Subroutine switch(option, on)
    Character (len=*), intent(in) :: option
    Logical, intent(in) :: on
    Integer :: n
    n = 0
    If (on) n = 1
    Select case (option)
    Case ('Langevin');              Call set_Langevin(on)
    Case ('HMC');                   Call set_HMC(on)
    Case ('Global_moves');          Call set_Global_moves(on)
    Case ('Global_tau_moves');      Call set_Global_tau_moves(on)
    Case ('N_Global_tau');          Call set_N_Global_tau(n)
    Case ('Propose_S0');            Call set_Propose_S0(on)
    Case ('Sequential_MALA');       Call set_Sequential_MALA(on)
    Case ('Global_MALA_moves');     Call set_Global_MALA_moves(on)
    Case ('Global_tau_MALA_moves'); Call set_Global_tau_MALA_moves(on)
    Case ('N_Global_tau_MALA');     Call set_N_Global_tau_MALA(n)
    Case default
       Call check(.false., 'unknown option '//option)
    End Select
  end Subroutine switch

  !> The reader that ran is named; ALF_GPU_DISABLE=1 (on the allow-list) is seen by both readers; ALF_GPU_FOO=1 (not on
  !> the allow-list) only by the /proc reader.
  Subroutine test_envreader()
    Character (len=256) :: v
    Character (len=96) :: reader
    Logical :: foo
    reader = gpu_env_reader_name()
    Write(*,*) 'reader: ', trim(reader)
    Call check(gpu_env_is('ALF_GPU_DISABLE', '1'), 'ALF_GPU_DISABLE=1 seen')
    foo = gpu_env_get('ALF_GPU_FOO', v)
    If (reader(1:5) == '/proc') then
       Call check(foo .and. trim(v) == '1', 'ALF_GPU_FOO=1 seen by the /proc reader')
    else
       Call check(reader(1:10) == 'allow-list', 'reader name')
       Call check(.not. foo, 'ALF_GPU_FOO not seen by the allow-list reader')
    endif
  end Subroutine test_envreader

  Subroutine check(ok, label)
    Logical, intent(in) :: ok
    Character (len=*), intent(in) :: label
    If (.not. ok) then
       Write(error_unit,*) 'FAIL ', label
       n_fail = n_fail + 1
    endif
  end Subroutine check

End Program Test_gpu_state

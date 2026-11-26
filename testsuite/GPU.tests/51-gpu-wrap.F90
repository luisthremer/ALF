! Test of the device wrap kernels (P2.4 step 1c), registered only when the build has -DALF_GPU.
! Usage: 51-gpu-wrap <mode>; one ctest per mode (see CMakeLists.txt), every mode in its own directory because Ham_set
! and the mirror write ./info.  The model is Hubbard on a 4x4 square lattice, beta = 2, Dtau = 0.1.
!
!   selftest [dense]  mirror set-up with ALF_GPU_SELFTEST=2 in the environment: Setup_gpu runs the gpu_wrap_mod self-test
!                     (hopping kernels, Op_V kernels, slice products and WRAPGRUP/DO sequences against the host routines
!                     Hop_mod_mm*, Op_mmult*, Op_Wrap*).  With ALF_GPU_FAULT=<fault> it must fail (negative control).
!                     'dense': Checkerboard = .F. (one 16 x 16 hopping element, la_zgemm / la_zgemm_cr paths).
! This program should exit successfully (return code 0) in the modes that are expected to pass.
Program Test_gpu_wrap
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: Clear_Ham_state, Ltrot
#ifdef MPI
  Use Hamiltonian_main, only: Set_Group_Comm
#endif
  Use QMC_runtime_var, only: set_QMC_runtime_default_var, set_ham_name, set_Nwrap, set_Ltau, get_Nwrap, &
       &                     set_default_global_tau_sampling, set_NSweep
  Use qmc_setup_mod, only: Setup_model, Setup_fields, Setup_update_schemes, Setup_gpu, Finalize_gpu
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc
  Use gpu_state_mod, only: gpu_active
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit
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
  Case ('selftest')
     Call setup_model_for_run(trim(arg) /= 'dense', N_op, .true.)
     Call Walker_layout(walker, Ltrot, get_Nwrap())
     Call Walker_alloc(walker)
     Call Setup_gpu(walker, 'info')           ! runs the wrap self-test (ALF_GPU_SELFTEST=2); a failure ends the run
     If (.not. gpu_active) then
        Write(error_unit,*) 'GPU path inactive'
        n_fail = n_fail + 1
     Endif
     Call Finalize_gpu()
     Call Walker_dealloc(walker)
  Case default
     Write(error_unit,*) 'usage: 51-gpu-wrap selftest [dense]'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  End Select

  Call Clear_Ham_state()
  If (n_fail > 0) then
     Write(error_unit,*) 'gpu-wrap test (', trim(mode), '): ', n_fail, ' failure(s)'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) 'gpu-wrap test passed: ', trim(mode)

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

end Program Test_gpu_wrap

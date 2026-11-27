! Test of the device tau pass (P2.4 step 1e), registered only when the build has -DALF_GPU.
! Usage: 54-gpu-tau selftest <Nwrap> <mz|su2>; one ctest per case, every run in its own directory.
!
!   selftest   Hubbard 4x4 (Ltau = 1), Symm and Checkerboard on, beta = 2, Dtau = 0.1: initial stacks and Green function
!              (Walker_init_stack / Walker_init_GR), handed to the device (Setup_gpu_stack), then gpu_tau_selftest: one pass of
!              the device tau pass against the upstream PropR / PropRm1, WRAPUR and CGR2_2 (matrices 1e-10 per slice and at every
!              stabilisation point, precision reductions 1e-12).  With ALF_GPU_FAULT_TAU = perm|prop|prec|blocks it must fail.
! This program should exit successfully (return code 0) in the modes that are expected to pass.
Program Test_gpu_tau
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: ham, Clear_Ham_state, Ltrot
#ifdef MPI
  Use Hamiltonian_main, only: Set_Group_Comm
#endif
  Use QMC_runtime_var, only: set_QMC_runtime_default_var, set_ham_name, set_Nwrap, set_Ltau, set_NSweep, set_NBin
  Use qmc_setup_mod, only: Setup_model, Setup_fields, Setup_update_schemes, Setup_gpu, Setup_gpu_stack, Finalize_gpu
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc, Walker_init_stack, Walker_init_GR
  Use gpu_state_mod, only: gpu_active
  Use gpu_stab_mod, only: gpu_stab_on
  Use gpu_tau_mod, only: gpu_tau_on, gpu_tau_setup, gpu_tau_finalize, gpu_tau_selftest
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit
  Implicit none

  Type (qmc_walker_t) :: walker
  Character (len=64) :: mode, arg, model
  Integer :: nwrap, N_op, ios, u, i
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Call Set_Group_Comm(MPI_COMM_WORLD)
#endif
  Call get_command_argument(1, mode)
  Call get_command_argument(2, arg)
  Call get_command_argument(3, model)
  If (trim(mode) /= 'selftest') then
     Write(error_unit,*) 'usage: 54-gpu-tau selftest <Nwrap> <mz|su2>'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  Endif
  nwrap = 5
  If (len_trim(arg) > 0) Read(arg,*,iostat=ios) nwrap
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
  If (trim(model) == 'mz') then
     Write(u,'(A)') 'N_FL = 2'
  else
     Write(u,'(A)') 'N_FL = 1'
  endif
  Write(u,'(A)') 'Dtau = 0.1d0'
  Write(u,'(A)') 'Beta = 2.d0'
  Write(u,'(A)') 'Projector = .F.'
  Write(u,'(A)') '/'
  Write(u,'(A)') '&VAR_Hubbard'
  If (trim(model) == 'mz') then
     Write(u,'(A)') 'Mz = .T.'
  else
     Write(u,'(A)') 'Mz = .F.'
  endif
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
  Call set_NBin(1)            ! the start-up check of P2.8 needs at least one bin
  Call set_Ltau(1)
  Call Setup_model(N_op)
  Call Setup_fields(N_op)
  Call ham%Alloc_obs(1)
  Call Setup_update_schemes(N_op, 'info')
  Call Walker_layout(walker, Ltrot, nwrap)
  Call Walker_alloc(walker)
  Call Setup_gpu(walker, 'info')
  If (.not. gpu_active .or. .not. gpu_stab_on) then
     Write(error_unit,*) 'GPU path or device stabilisation inactive'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  Endif
  If (.not. gpu_tau_on) Call gpu_tau_setup(walker)      ! Setup_gpu does it from step 1e hooks on
  Call ham%Init_obs(1)
  Call Walker_init_stack(walker)
  Call Walker_init_GR(walker)
  Call Setup_gpu_stack(walker)
  Call gpu_tau_selftest(walker)
  Call gpu_tau_finalize()
  Call Finalize_gpu()
  Call Walker_dealloc(walker)
  Call Clear_Ham_state()
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) 'gpu-tau test passed: ', trim(mode), ' ', nwrap, ' ', trim(model)
end Program Test_gpu_tau

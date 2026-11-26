! Unit test for start-up checks outside QMC_runtime_var. Usage: 43-prog-checks-fail <case>
! Each case must terminate with the message that ctest matches (PASS_REGULAR_EXPRESSION).
! KNOWN_ISSUES.md: H4 (confin field 0), M1 (Ltrot from Dtau = 0), M3 (Nwrap < 1).
Program Test_prog_checks_fail

  Use Fields_mod, only: Fields
  Use Hamiltonian_main, only: Set_model_for_testing, Validate_Ham_Variables
  Use Operator_mod, only: Operator
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout
#ifdef MPI
  Use mpi
#endif
  Implicit none
  Character (len=32) :: name
  Type (qmc_walker_t) :: walker
  Type (Fields) :: nsigma
  Type (Operator), dimension(:,:), allocatable :: Op_T, Op_V
  Integer :: unit_seeds
  Integer :: group_comm
  Complex (Kind=Kind(0.d0)), allocatable :: f0(:,:)
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  group_comm = MPI_COMM_WORLD
#else
  group_comm = 0
#endif
  Call get_command_argument(1, name)

  Select Case (trim(name))
  Case ("field0")   ! a start field with value 0 (type 1: +-1) must be rejected, as if read from confin_0
     Open(newunit=unit_seeds, file='seeds', status='replace', action='write')
     write(unit_seeds,*) 1234
     Close(unit_seeds)
     Call nsigma%make(2, 3)
     nsigma%t = 1
     allocate(f0(2,3));  f0 = (0.d0, 0.d0)
     Call nsigma%in(group_comm, Initial_field=f0)
  Case ("ltrot")    ! Dtau = 0 gives Ltrot = nint(Inf) = huge
     allocate(Op_T(0,1), Op_V(0,1))
     Call Set_model_for_testing(4, 1, 1, huge(1), 0, .false., .false., Op_T, Op_V)
     Call Validate_Ham_Variables()
  Case ("nwrap")
     Call Walker_layout(walker, 10, 0)
  Case Default
     write(*,*) "ERROR: unknown case ", trim(name)
  end Select

  write(*,*) "ERROR: invalid value did not terminate"
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
end Program Test_prog_checks_fail

! Precision Green warning (Prog/control_mod.F90, Control_Print; KNOWN_ISSUES.md, M6).
! Two Green functions that differ in one element by the deviation given as argument ("large": 0.5,
! "small": 2**-20 ~ 1e-6) are compared by Control_PrecisionG, then Control_Print writes the summary.
!   large -> above the warning level 1e-3: "Warning: Precision Green Max" on stderr and in info
!   small -> no warning
Program precision_green_warning_test
  Use Control, only: control_init, Control_PrecisionG, Control_Print
#ifdef MPI
  Use mpi
#endif
  Implicit none
  Integer, parameter :: Ndim = 3
  Complex (Kind=Kind(0.d0)) :: A(Ndim,Ndim), B(Ndim,Ndim)
  Character (len=64) :: scheme
  Character (len=32) :: name
  Integer :: i, j, Group_Comm
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Group_Comm = MPI_COMM_WORLD
#else
  Group_Comm = 0
#endif
  Call get_command_argument(1, name)
  Call control_init(Group_Comm)
  do j = 1, Ndim
     do i = 1, Ndim
        A(i,j) = cmplx(1.d0/dble(i+j), 0.25d0*dble(i-j), kind(0.d0))
     enddo
  enddo
  B = A
  if (trim(name) == "large") then
     B(2,3) = B(2,3) + cmplx(0.5d0, 0.d0, kind(0.d0))
  else
     B(2,3) = B(2,3) + cmplx(2.d0**(-20), 0.d0, kind(0.d0))
  endif
  Call Control_PrecisionG(A, B, Ndim)
  scheme = "None"
  Call Control_Print(Group_Comm, scheme, .false., .false., .false.)
  write(*,*) 'PASS: summary written'
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
end Program precision_green_warning_test

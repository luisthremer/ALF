! Unit test for Control_Check_Phase: a phase of unit modulus passes, a NaN phase must terminate
! the program with the phase message (ctest PASS_REGULAR_EXPRESSION; KNOWN_ISSUES.md,
! "Determinant phase underflow").
Program Test_check_phase_fail

  Use Control, only: Control_Check_Phase
  Use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
#ifdef MPI
  Use mpi
#endif
  Implicit none
  Complex (Kind=Kind(0.d0)) :: Z
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
#endif

  Call Control_Check_Phase(cmplx(0.6d0, 0.8d0, kind(0.d0)))   ! |Z| = 1: passes
  write(*,*) "unit-modulus phase accepted"

  Z = cmplx(ieee_value(1.d0, ieee_quiet_nan), 0.d0, kind(0.d0))  ! NaN without an FP exception
  Call Control_Check_Phase(Z)                                   ! must terminate

  write(*,*) "ERROR: NaN phase did not terminate"
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
end Program Test_check_phase_fail

! ALF_UNSTABLE_THRESHOLD (Prog/control_mod.F90, control_read_threshold): abort threshold of the G
! precision check Control_PrecisionG (upstream: 10).
! The program reads the environment variable as control_init does and then compares two Green
! functions that differ by exactly 0.5 in one element (max |A - B| = 0.5, exact in binary).
!   unset            -> 0.5 < 10    -> passes, exit 0                       (38-unstable-threshold-default)
!   = 1.0            -> 0.5 < 1     -> passes, exit 0                       (38-unstable-threshold-1.0)
!   = 0.25           -> 0.5 > 0.25  -> "exceeding the threshold", abort     (38-unstable-threshold-0.25*)
!   = abc / -1 / 0   -> not a positive number -> alf_exit(2), exit status 2 (38-unstable-threshold-bad-*)
Program unstable_threshold_test
  Use Control, only: control_read_threshold, Control_PrecisionG
  Implicit none
  Integer, parameter :: Ndim = 3
  Complex (Kind=Kind(0.d0)) :: A(Ndim,Ndim), B(Ndim,Ndim)
  Integer :: i, j
  call control_read_threshold()
  do j = 1, Ndim
     do i = 1, Ndim
        A(i,j) = cmplx(1.d0/dble(i+j), 0.25d0*dble(i-j), kind(0.d0))
     enddo
  enddo
  B = A
  B(2,3) = B(2,3) + cmplx(0.5d0, 0.d0, kind(0.d0))
  call Control_PrecisionG(A, B, Ndim)
  write(*,*) 'PASS: deviation 0.5 below the threshold'
end Program unstable_threshold_test

! Self-test of the test helpers: a NaN / Inf measurement must be reported as a failure, finite values must pass.
program t_checks
  use la_test_util
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  implicit none
  real(8) :: nan, inf
  nan = ieee_value(1.d0, ieee_quiet_nan)
  inf = ieee_value(1.d0, ieee_positive_inf)
  call test_begin('t_checks')
  call check('finite zero', 0.d0, 1.d-12)
  if (nfail /= 0) error stop 2
  call check('NaN (must be counted as failure)', nan, 1.d-12)
  if (nfail /= 1) error stop 3
  call check_inv('NaN inv (must be counted as failure)', nan, 1.d-12)
  if (nfail /= 2) error stop 4
  call check('Inf (must be counted as failure)', inf, 1.d-12)
  if (nfail /= 3) error stop 5
  call check_inv('Inf inv (must be counted as failure)', inf, 1.d-12)
  if (nfail /= 4) error stop 6
  call check_inv('large finite above tol', 1.d-3, 1.d-12)
  if (nfail /= 5) error stop 7
  block
    complex(8) :: X(4,3), Y(4,3)
    X = (1.d0, 0.d0); Y = X
    call check('maxdiff finite', maxdiff(4, 3, X, Y), 0.d0)
    if (nfail /= 5) error stop 8
    X(2,2) = cmplx(nan, 0.d0, kind=8)
    call check('maxdiff with one NaN (must fail)', maxdiff(4, 3, X, Y), 1.d-12)
    if (nfail /= 6) error stop 9
    X(2,2) = (1.d0, 0.d0); X(3,1) = cmplx(inf, 0.d0, kind=8)
    call check('maxdiff with one Inf (must fail)', maxdiff(4, 3, X, Y), 1.d-12)
    if (nfail /= 7) error stop 10
    X(3,1) = (1.d0, 0.d0); X(1,3) = cmplx(nan, 0.d0, kind=8)      ! NaN outside the 2x2 block at (2,2)
    call check('outside_diff with one NaN (must fail)', outside_diff(4, 3, X, Y, 2, 2, 2, 2), 0.d0)
    if (nfail /= 8) error stop 11
    X(1,3) = (1.d0, 0.d0); X(4,1) = cmplx(inf, 0.d0, kind=8)
    call check('outside_diff with one Inf (must fail)', outside_diff(4, 3, X, Y, 2, 2, 2, 2), 0.d0)
    if (nfail /= 9) error stop 12
    X(4,1) = (1.d0, 0.d0); X(2,2) = cmplx(nan, 0.d0, kind=8)      ! NaN inside the block is not "outside"
    call check('outside_diff, NaN inside block (must pass)', outside_diff(4, 3, X, Y, 2, 2, 2, 2), 0.d0)
    if (nfail /= 9) error stop 13
    call check('fmax NaN (must fail)', fmax(0.d0, nan), 1.d-12)
    if (nfail /= 10) error stop 14
  end block
  write(*,*) '  (the 10 FAIL lines above are the intended detections)'
  nfail = 0; maxval_seen = 0.d0
  call test_end()
end program t_checks

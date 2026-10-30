! Helpers for the alf_la unit tests (ticket 202): reproducible random data, error norms, device mapping, reporting.
! HOST backend: every check must be bit-identical (tolerance 0). CUDA backend: tolerance as in the ticket.
! A NaN/Inf measurement is always a failure (ieee_is_finite). This file is compiled WITHOUT finite-math assumptions
! (see CMakeLists.txt: -fno-finite-math-only / -fp-model=precise / -Kieee) because the test flags contain -ffast-math etc.
module la_test_util
  use alf_la_mod, only: la_backend_name
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  integer, parameter :: sizes(5) = [1, 7, 64, 144, 400]
  real(8), parameter :: TOL_CUDA = 1.d-12
  integer(8), save :: rng = 20261027_8
  integer, save :: nfail = 0, nchk = 0
  character(len=32), save :: tname = ''
  real(8), save :: maxval_seen = 0.d0
contains

  logical function is_host()
    is_host = (trim(la_backend_name()) == 'host')
  end function

  subroutine test_begin(name)
    character(len=*), intent(in) :: name
    tname = name
    nfail = 0; nchk = 0; maxval_seen = 0.d0
    write(*,'(3A)') 'TEST ', trim(name), '  backend=' // trim(la_backend_name())
  end subroutine

  ! Uniform in (-0.5, 0.5); 31-bit LCG in int64: identical on all compilers.
  real(8) function rnd()
    rng = mod(rng * 1103515245_8 + 12345_8, 2147483648_8)
    rnd = dble(rng) / 2147483648.d0 - 0.5d0
  end function

  subroutine fill_z(ld, nc, A)
    integer, intent(in) :: ld, nc
    complex(8), intent(out) :: A(ld, nc)
    integer :: i, j
    real(8) :: re, im
    do j = 1, nc
      do i = 1, ld
        re = rnd(); im = rnd()
        A(i,j) = cmplx(re, im, kind=8)
      end do
    end do
  end subroutine

  ! Frobenius norm of the m x n block of A (leading dimension ld)
  real(8) function nrmf(m, n, A, ld)
    integer, intent(in) :: m, n, ld
    complex(8), intent(in) :: A(ld, *)
    integer :: i, j
    real(8) :: s
    s = 0.d0
    do j = 1, n
      do i = 1, m
        s = s + dble(A(i,j))**2 + aimag(A(i,j))**2
      end do
    end do
    nrmf = sqrt(s)
  end function

  ! || X - Y ||_F / || Y ||_F over an m x n block (0 if both are zero)
  real(8) function relfro(m, n, X, ldx, Y, ldy)
    integer, intent(in) :: m, n, ldx, ldy
    complex(8), intent(in) :: X(ldx, *), Y(ldy, *)
    integer :: i, j
    real(8) :: s, t
    s = 0.d0; t = 0.d0
    do j = 1, n
      do i = 1, m
        s = s + abs(X(i,j) - Y(i,j))**2
        t = t + abs(Y(i,j))**2
      end do
    end do
    if (t == 0.d0) then
      relfro = sqrt(s)
    else
      relfro = sqrt(s / t)
    end if
  end function

  ! max(a,b) that returns NaN if a or b is NaN/Inf (Fortran MAX may drop NaN)
  real(8) function fmax(a, b)
    real(8), intent(in) :: a, b
    if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) then
      fmax = ieee_value(1.d0, ieee_quiet_nan)
    else
      fmax = max(a, b)
    end if
  end function

  ! max |x(i) - y(i)| over rows r0..r1 of two ld x nc arrays, explicit loops (no temporaries). NaN as soon as any
  ! difference is NaN or Inf (Fortran MAXVAL ignores NaN, so it must not be used for error measures).
  real(8) function maxdiff_rows(ld, nc, r0, r1, X, Y)
    integer, intent(in) :: ld, nc, r0, r1
    complex(8), intent(in) :: X(ld, nc), Y(ld, nc)
    integer :: i, j
    real(8) :: v
    maxdiff_rows = 0.d0
    do j = 1, nc
      do i = r0, r1
        v = abs(X(i,j) - Y(i,j))
        if (.not. ieee_is_finite(v)) then
          maxdiff_rows = ieee_value(1.d0, ieee_quiet_nan)
          return
        end if
        maxdiff_rows = max(maxdiff_rows, v)
      end do
    end do
  end function

  ! flat vectors
  real(8) function maxabs_diff(n, x, y)
    integer, intent(in) :: n
    complex(8), intent(in) :: x(n), y(n)
    maxabs_diff = maxdiff_rows(n, 1, 1, n, x, y)
  end function

  ! largest |X - Y| over the whole ld x nc array (bit-identity checks; NaN-aware)
  real(8) function maxdiff(ld, nc, X, Y)
    integer, intent(in) :: ld, nc
    complex(8), intent(in) :: X(ld, nc), Y(ld, nc)
    maxdiff = maxdiff_rows(ld, nc, 1, ld, X, Y)
  end function

  ! largest |X - Y| over the part of the ld x nc array OUTSIDE the m x n block at (oi,oj) (must be exactly 0)
  real(8) function outside_diff(ld, nc, X, Y, oi, oj, m, n)
    integer, intent(in) :: ld, nc, oi, oj, m, n
    complex(8), intent(in) :: X(ld, nc), Y(ld, nc)
    integer :: i, j
    real(8) :: v
    outside_diff = 0.d0
    do j = 1, nc
      do i = 1, ld
        if (i >= oi .and. i <= oi+m-1 .and. j >= oj .and. j <= oj+n-1) cycle
        v = abs(X(i,j) - Y(i,j))
        if (.not. ieee_is_finite(v)) then
          outside_diff = ieee_value(1.d0, ieee_quiet_nan)
          return
        end if
        outside_diff = max(outside_diff, v)
      end do
    end do
  end function

  ! Record a check: HOST needs val == 0 (bit-identical), CUDA needs val <= tol.
  subroutine check(label, val, tol)
    character(len=*), intent(in) :: label
    real(8), intent(in) :: val, tol
    real(8) :: t
    logical :: ok
    t = tol
    if (is_host()) t = 0.d0
    ok = ieee_is_finite(val) .and. (val <= t)
    nchk = nchk + 1
    if (ieee_is_finite(val)) maxval_seen = max(maxval_seen, val)
    if (.not. ok) then
      nfail = nfail + 1
      write(*,'(A,A,A,ES10.3,A,ES9.2)') '  FAIL ', trim(label), ' value ', val, ' > tol ', t
    end if
  end subroutine

  ! Invariant check: tolerance applies to both backends (a residual is never bit-identical to anything)
  subroutine check_inv(label, val, tol)
    character(len=*), intent(in) :: label
    real(8), intent(in) :: val, tol
    nchk = nchk + 1
    if (ieee_is_finite(val)) maxval_seen = max(maxval_seen, val)
    if (.not. (ieee_is_finite(val) .and. val <= tol)) then
      nfail = nfail + 1
      write(*,'(A,A,A,ES10.3,A,ES9.2)') '  FAIL ', trim(label), ' value ', val, ' > tol ', tol
    end if
  end subroutine

  ! Exact (integer) comparison counted as a check (value = number of differing entries)
  subroutine check_int(label, n, a, b)
    character(len=*), intent(in) :: label
    integer, intent(in) :: n
    integer(4), intent(in) :: a(n), b(n)
    call check(label, dble(count(a /= b)), 0.d0)
  end subroutine

  subroutine test_end()
    write(*,'(A,A,A,I0,A,I0,A,ES10.3)') '  ', trim(tname), ': checks ', nchk, ' failed ', nfail, &
         '  max measured ', maxval_seen
    if (nfail > 0) then
      write(*,'(A,A)') 'FAILED ', trim(tname)
      error stop 1
    end if
    write(*,'(A,A)') 'PASSED ', trim(tname)
  end subroutine

  ! ---- device mapping helpers (no-ops for the HOST backend) -------------------------------------
  subroutine dev_enter_z(a, n)
    integer, intent(in) :: n
    complex(8) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target enter data map(to: a)
#endif
  end subroutine
  subroutine dev_enter_i(a, n)
    integer, intent(in) :: n
    integer(4) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target enter data map(to: a)
#endif
  end subroutine
  subroutine dev_exit_z(a, n)
    integer, intent(in) :: n
    complex(8) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target exit data map(delete: a)
#endif
  end subroutine
  subroutine dev_exit_i(a, n)
    integer, intent(in) :: n
    integer(4) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target exit data map(delete: a)
#endif
  end subroutine
  subroutine dev_from_z(a, n)
    integer, intent(in) :: n
    complex(8) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target update from(a)
#endif
  end subroutine
  subroutine dev_from_i(a, n)
    integer, intent(in) :: n
    integer(4) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target update from(a)
#endif
  end subroutine
  subroutine dev_enter_d(a, n)
    integer, intent(in) :: n
    real(8) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target enter data map(to: a)
#endif
  end subroutine
  subroutine dev_exit_d(a, n)
    integer, intent(in) :: n
    real(8) :: a(n)
#if defined(ALF_LA_CUDA) || defined(ALF_LA_ONEMKL)
    !$omp target exit data map(delete: a)
#endif
  end subroutine
end module la_test_util

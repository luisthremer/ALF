! alf_la unit test: la_geru vs direct ZGERU (ticket 202). Cases: square, rectangular with lda > m and
! strided vectors (incx = 2, incy = 3), sub-block at (3,2).
program t_geru
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgeru
  integer :: in, cs
  call test_begin('t_geru')
  call la_init()
  do in = 1, size(sizes)
    do cs = 1, 3
      call run(sizes(in), cs)
    end do
  end do
  call la_finalize()
  call test_end()
contains
  subroutine run(N, cs)
    integer, intent(in) :: N, cs
    integer :: m, nn, oi, oj, lda, incx, incy
    complex(8), allocatable :: A(:,:), Aref(:,:), x(:), y(:)
    complex(8) :: alpha
    real(8) :: e, e2
    character(len=64) :: lab
    m = N; nn = N; oi = 1; oj = 1; lda = m; incx = 1; incy = 1
    if (cs >= 2) then
      m = N + 3; nn = N + 5; lda = m + 5; incx = 2; incy = 3
    end if
    if (cs == 3) then
      oi = 3; oj = 2; lda = m + oi - 1 + 4
    end if
    allocate(A(lda, nn+oj), Aref(lda, nn+oj), x(1+(m-1)*incx), y(1+(nn-1)*incy))
    call fill_z(lda, nn+oj, A); call fill_z(size(x), 1, x); call fill_z(size(y), 1, y)
    Aref = A
    alpha = cmplx(0.6d0, -0.3d0, kind=8)
    call zgeru(m, nn, alpha, x, incx, y, incy, Aref(oi,oj), lda)
    call dev_enter_z(A, size(A)); call dev_enter_z(x, size(x)); call dev_enter_z(y, size(y))
    call la_geru(m, nn, alpha, x, incx, y, incy, A(oi,oj), lda)
    call la_sync()
    call dev_from_z(A, size(A))
    call dev_exit_z(A, size(A)); call dev_exit_z(x, size(x)); call dev_exit_z(y, size(y))
    e = relfro(m, nn, A(oi,oj), lda, Aref(oi,oj), lda)
    A(oi:oi+m-1, oj:oj+nn-1) = (0.d0, 0.d0); Aref(oi:oi+m-1, oj:oj+nn-1) = (0.d0, 0.d0)
    e2 = maxdiff(lda, nn+oj, A, Aref)
    write(lab,'(A,I0,A,I0)') 'geru N=', N, ' case', cs
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine
end program t_geru

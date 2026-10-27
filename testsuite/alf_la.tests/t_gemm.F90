! alf_la unit test: la_gemm, la_gemm_cr, la_gemm_batched vs direct host BLAS (ticket 202).
! Cases per size N: 1 square, ld = rows; 2 rectangular (m=N+3, n=N+5, k=N+2) with lda > rows;
! 3 rectangular sub-block starting at element (3,2) of a larger array (checks that nothing outside the block is written).
program t_gemm
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgemm, dgemm
  integer :: in, cs, ia, ib
  integer :: ops(2,6)
  call test_begin('t_gemm')
  call la_init()
  ops = reshape([LA_N,LA_N, LA_C,LA_N, LA_N,LA_C, LA_N,LA_T, LA_T,LA_N, LA_C,LA_C], [2,6])
  do in = 1, size(sizes)
    do cs = 1, 3
      do ia = 1, 6
        call run_gemm(sizes(in), cs, ops(1,ia), ops(2,ia))
      end do
      call run_cr(sizes(in), cs)
    end do
    do ib = 1, 2
      call run_batched(sizes(in), merge(1, 32, ib == 1), 0)
      call run_batched(sizes(in), merge(1, 32, ib == 1), 3)      ! lda > m, strides with gaps
    end do
  end do
  call la_finalize()
  call test_end()

contains

  character function ch(op)
    integer, intent(in) :: op
    select case (op)
    case (LA_N); ch = 'N'
    case (LA_T); ch = 'T'
    case default; ch = 'C'
    end select
  end function

  subroutine dims(N, cs, m, nn, k, oi, oj, pad)
    integer, intent(in) :: N, cs
    integer, intent(out) :: m, nn, k, oi, oj, pad
    m = N; nn = N; k = N; oi = 1; oj = 1; pad = 0
    if (cs >= 2) then
      m = N + 3; nn = N + 5; k = N + 2; pad = 5
    end if
    if (cs == 3) then
      oi = 3; oj = 2; pad = 4
    end if
  end subroutine

  subroutine zero_block(ld, nc, X, oi, oj, r, c)
    integer, intent(in) :: ld, nc, oi, oj, r, c
    complex(8), intent(inout) :: X(ld, nc)
    X(oi:oi+r-1, oj:oj+c-1) = (0.d0, 0.d0)
  end subroutine

  subroutine run_gemm(N, cs, opa, opb)
    integer, intent(in) :: N, cs, opa, opb
    integer :: m, nn, k, oi, oj, pad, ra, ca, rb, cb, lda, ldb, ldc
    complex(8), allocatable :: A(:,:), B(:,:), C(:,:), Cref(:,:)
    complex(8) :: alpha, beta
    real(8) :: e, e2
    character(len=64) :: lab
    call dims(N, cs, m, nn, k, oi, oj, pad)
    ra = merge(m, k, opa == LA_N); ca = merge(k, m, opa == LA_N)
    rb = merge(k, nn, opb == LA_N); cb = merge(nn, k, opb == LA_N)
    lda = ra + oi - 1 + pad; ldb = rb + oi - 1 + pad; ldc = m + oi - 1 + pad
    allocate(A(lda, ca+oj), B(ldb, cb+oj), C(ldc, nn+oj), Cref(ldc, nn+oj))
    call fill_z(lda, ca+oj, A); call fill_z(ldb, cb+oj, B); call fill_z(ldc, nn+oj, C)
    Cref = C
    alpha = cmplx(0.6d0, -0.3d0, kind=8); beta = cmplx(0.2d0, 0.7d0, kind=8)
    call zgemm(ch(opa), ch(opb), m, nn, k, alpha, A(oi,oj), lda, B(oi,oj), ldb, beta, Cref(oi,oj), ldc)
    call dev_enter_z(A, size(A)); call dev_enter_z(B, size(B)); call dev_enter_z(C, size(C))
    call la_gemm(opa, opb, m, nn, k, alpha, A(oi,oj), lda, B(oi,oj), ldb, beta, C(oi,oj), ldc)
    call la_sync()
    call dev_from_z(C, size(C))
    call dev_exit_z(A, size(A)); call dev_exit_z(B, size(B)); call dev_exit_z(C, size(C))
    e = relfro(m, nn, C(oi,oj), ldc, Cref(oi,oj), ldc)
    call zero_block(ldc, nn+oj, C, oi, oj, m, nn); call zero_block(ldc, nn+oj, Cref, oi, oj, m, nn)
    e2 = maxdiff(ldc, nn+oj, C, Cref)
    write(lab,'(A,I0,A,I0,A,A,A)') 'gemm N=', N, ' case', cs, ' ', ch(opa)//ch(opb), ''
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine

  ! C(complex m x n) = A(complex m x k) * B(real k x n)
  subroutine run_cr(N, cs)
    integer, intent(in) :: N, cs
    integer :: m, nn, k, oi, oj, pad, lda, ldb, ldc, i, j
    complex(8), allocatable :: A(:,:), C(:,:), Cref(:,:)
    real(8), allocatable :: B(:,:), Ar(:,:), Cr(:,:)
    real(8) :: e, e2
    character(len=64) :: lab
    call dims(N, cs, m, nn, k, oi, oj, pad)
    lda = m + oi - 1 + pad; ldb = k + oi - 1 + pad; ldc = m + oi - 1 + pad
    allocate(A(lda, k+oj), B(ldb, nn+oj), C(ldc, nn+oj), Cref(ldc, nn+oj), Ar(2*lda, k+oj), Cr(2*ldc, nn+oj))
    call fill_z(lda, k+oj, A); call fill_z(ldc, nn+oj, C)
    do j = 1, nn + oj
      do i = 1, ldb
        B(i,j) = rnd()
      end do
    end do
    do j = 1, k + oj
      do i = 1, lda
        Ar(2*i-1, j) = dble(A(i,j)); Ar(2*i, j) = aimag(A(i,j))
      end do
    end do
    Cr = 0.d0
    do j = 1, nn + oj
      do i = 1, ldc
        Cr(2*i-1, j) = dble(C(i,j)); Cr(2*i, j) = aimag(C(i,j))
      end do
    end do
    ! direct host reference: DGEMM on the interleaved real view (beta = 0 as in la_gemm_cr)
    call dgemm('N', 'N', 2*m, nn, k, 1.d0, Ar(2*oi-1, oj), 2*lda, B(oi,oj), ldb, 0.d0, Cr(2*oi-1, oj), 2*ldc)
    do j = 1, nn + oj
      do i = 1, ldc
        Cref(i,j) = cmplx(Cr(2*i-1, j), Cr(2*i, j), kind=8)
      end do
    end do
    call dev_enter_z(A, size(A)); call dev_enter_d(B, size(B)); call dev_enter_z(C, size(C))
    call la_gemm_cr(m, nn, k, A(oi,oj), lda, B(oi,oj), ldb, C(oi,oj), ldc)
    call la_sync()
    call dev_from_z(C, size(C))
    call dev_exit_z(A, size(A)); call dev_exit_d(B, size(B)); call dev_exit_z(C, size(C))
    e = relfro(m, nn, C(oi,oj), ldc, Cref(oi,oj), ldc)
    call zero_block(ldc, nn+oj, C, oi, oj, m, nn); call zero_block(ldc, nn+oj, Cref, oi, oj, m, nn)
    e2 = maxdiff(ldc, nn+oj, C, Cref)
    write(lab,'(A,I0,A,I0)') 'gemm_cr N=', N, ' case', cs
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine

  subroutine run_batched(N, nb, pad)
    integer, intent(in) :: N, nb, pad
    integer :: m, nn, k, ld, s, ibat, j
    complex(8), allocatable :: A(:), B(:), C(:), Cref(:)
    complex(8) :: alpha, beta
    real(8) :: e, e2
    character(len=64) :: lab
    m = N; nn = N; k = N; ld = N + pad; s = ld*N + 2*pad
    allocate(A(s*nb), B(s*nb), C(s*nb), Cref(s*nb))
    call fill_z(s*nb, 1, A); call fill_z(s*nb, 1, B); call fill_z(s*nb, 1, C)
    Cref = C
    alpha = cmplx(0.6d0, -0.3d0, kind=8); beta = cmplx(0.2d0, 0.7d0, kind=8)
    do ibat = 0, nb - 1
      call zgemm('N', 'C', m, nn, k, alpha, A(1+ibat*s), ld, B(1+ibat*s), ld, beta, Cref(1+ibat*s), ld)
    end do
    call dev_enter_z(A, size(A)); call dev_enter_z(B, size(B)); call dev_enter_z(C, size(C))
    call la_gemm_batched(LA_N, LA_C, m, nn, k, alpha, A, ld, s, B, ld, s, beta, C, ld, s, nb)
    call la_sync()
    call dev_from_z(C, size(C))
    call dev_exit_z(A, size(A)); call dev_exit_z(B, size(B)); call dev_exit_z(C, size(C))
    e = relfro(s*nb, 1, C, s*nb, Cref, s*nb)
    ! padding rows (i > N inside each column) and the gaps between the matrices must be untouched
    e2 = 0.d0
    do ibat = 0, nb - 1
      do j = 1, nn
        C(1+ibat*s+(j-1)*ld : ibat*s+(j-1)*ld+N) = (0.d0, 0.d0)
        Cref(1+ibat*s+(j-1)*ld : ibat*s+(j-1)*ld+N) = (0.d0, 0.d0)
      end do
    end do
    e2 = maxabs_diff(size(C), C, Cref)
    write(lab,'(A,I0,A,I0,A,I0)') 'gemm_batched N=', N, ' k=', nb, ' pad=', pad
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine

end program t_gemm

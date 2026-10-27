! alf_la unit test: la_getrf (square, rectangular, lda > m, sub-block) and la_getrs (N, T, C) vs direct
! ZGETRF/ZGETRS (ticket 202). getrf cases: 1 square, 2 tall (m=N+3), 3 tall sub-block, 4 wide (n=N+4).
! Test matrices: random + boosted diagonal (condition number ~ 10), rows then permuted at random so that the pivots are
! nontrivial but unambiguous: LU factors, ipiv and the solutions are then comparable between cuSOLVER and LAPACK
! (a generic random matrix would put cond * eps ~ 1e-12 into the solution and make the 1e-12 check meaningless).
program t_lu
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgetrf, zgetrs
  integer :: in, cs, it
  call test_begin('t_lu')
  call la_init()
  do in = 1, size(sizes)
    do cs = 1, 3
      call run_getrf(sizes(in), cs)
      do it = 0, 2
        call run_solve(sizes(in), cs, it)
      end do
    end do
    call run_getrf(sizes(in), 4)          ! wide: m < n

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

  ! m x n block at (oi,oj) of A(ld, .)
  subroutine make_matrix(ld, nc, A, oi, oj, m, n)
    integer, intent(in) :: ld, nc, oi, oj, m, n
    complex(8), intent(out) :: A(ld, nc)
    complex(8) :: t(n)
    integer :: i, p
    call fill_z(ld, nc, A)
    do i = 1, min(m, n)
      A(oi+i-1, oj+i-1) = A(oi+i-1, oj+i-1) + cmplx(1.5d0*sqrt(dble(n)) + 1.d0, 0.d0, kind=8)
    end do
    do i = m, 2, -1                              ! random row permutation inside the block
      p = 1 + int((rnd() + 0.5d0) * dble(i))
      p = min(max(p, 1), i)
      t = A(oi+i-1, oj:oj+n-1); A(oi+i-1, oj:oj+n-1) = A(oi+p-1, oj:oj+n-1); A(oi+p-1, oj:oj+n-1) = t
    end do
  end subroutine

  subroutine run_getrf(N, cs)
    integer, intent(in) :: N, cs
    integer :: m, nn, oi, oj, lda
    complex(8), allocatable :: A(:,:), Aref(:,:)
    integer(4), allocatable :: ip(:), ipref(:)
    integer :: info
    real(8) :: e, e2
    character(len=64) :: lab
    m = N; nn = N; oi = 1; oj = 1; lda = N
    if (cs >= 2 .and. cs <= 3) then
      m = N + 3; lda = m + 5
    end if
    if (cs == 3) then
      oi = 3; oj = 2; lda = m + oi - 1 + 4
    end if
    if (cs == 4) then                      ! wide matrix, m < n
      nn = N + 4; lda = m + 5
    end if
    allocate(A(lda, nn+oj), Aref(lda, nn+oj), ip(min(m,nn)), ipref(min(m,nn)))
    call make_matrix(lda, nn+oj, A, oi, oj, m, nn)
    Aref = A
    call zgetrf(m, nn, Aref(oi,oj), lda, ipref, info)
    ip = -7
    call dev_enter_z(A, size(A)); call dev_enter_i(ip, size(ip))
    call la_getrf(m, nn, A(oi,oj), lda, ip)
    call la_sync()
    call dev_from_z(A, size(A)); call dev_from_i(ip, size(ip))
    call dev_exit_z(A, size(A)); call dev_exit_i(ip, size(ip))
    write(lab,'(A,I0,A,I0)') 'getrf N=', N, ' case', cs
    call check_int(trim(lab)//' ipiv', size(ip), ip, ipref)
    e = relfro(m, nn, A(oi,oj), lda, Aref(oi,oj), lda)
    A(oi:oi+m-1, oj:oj+nn-1) = (0.d0, 0.d0); Aref(oi:oi+m-1, oj:oj+nn-1) = (0.d0, 0.d0)
    e2 = maxdiff(lda, nn+oj, A, Aref)
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine

  subroutine run_solve(N, cs, it)
    integer, intent(in) :: N, cs, it
    integer :: nrhs, oi, oj, lda, ldb, info, trans
    complex(8), allocatable :: A(:,:), Aref(:,:), B(:,:), Bref(:,:)
    integer(4), allocatable :: ip(:), ipref(:)
    real(8) :: e, e2
    character(len=64) :: lab
    trans = it
    nrhs = N; oi = 1; oj = 1; lda = N; ldb = N
    if (cs >= 2) then
      nrhs = N + 2; lda = N + 5; ldb = N + 6
    end if
    if (cs == 3) then
      oi = 3; oj = 2; lda = N + oi - 1 + 4; ldb = N + oi - 1 + 4
    end if
    allocate(A(lda, N+oj), Aref(lda, N+oj), B(ldb, nrhs+oj), Bref(ldb, nrhs+oj), ip(N), ipref(N))
    call make_matrix(lda, N+oj, A, oi, oj, N, N)
    call fill_z(ldb, nrhs+oj, B)
    Aref = A; Bref = B
    call zgetrf(N, N, Aref(oi,oj), lda, ipref, info)
    call zgetrs(ch(trans), N, nrhs, Aref(oi,oj), lda, ipref, Bref(oi,oj), ldb, info)
    ip = -7
    call dev_enter_z(A, size(A)); call dev_enter_z(B, size(B)); call dev_enter_i(ip, size(ip))
    call la_getrf(N, N, A(oi,oj), lda, ip)
    call la_getrs(trans, N, nrhs, A(oi,oj), lda, ip, B(oi,oj), ldb)
    call la_sync()
    call dev_from_z(B, size(B))
    call dev_exit_z(A, size(A)); call dev_exit_z(B, size(B)); call dev_exit_i(ip, size(ip))
    write(lab,'(A,I0,A,I0,A,A)') 'getrs N=', N, ' case', cs, ' ', ch(trans)
    e = relfro(N, nrhs, B(oi,oj), ldb, Bref(oi,oj), ldb)
    B(oi:oi+N-1, oj:oj+nrhs-1) = (0.d0, 0.d0); Bref(oi:oi+N-1, oj:oj+nrhs-1) = (0.d0, 0.d0)
    e2 = maxdiff(ldb, nrhs+oj, B, Bref)
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine
end program t_lu

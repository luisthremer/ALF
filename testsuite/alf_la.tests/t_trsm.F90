! alf_la unit test: la_trsm vs direct ZTRSM (ticket 202). The three combinations used by CGR/CGR2_2
! (L,U,N,N), (R,U,C,N), (L,U,C,N) plus (L,L,N,U), (R,L,T,N) and (R,U,N,N).
program t_trsm
  use alf_la_mod
  use la_test_util
  implicit none
  external :: ztrsm
  character(len=1), parameter :: sd(6) = ['L','R','L','L','R','R'], up(6) = ['U','U','U','L','L','U'], &
       dg(6) = ['N','N','N','U','N','N']
  integer, parameter :: tr(6) = [LA_N, LA_C, LA_C, LA_N, LA_T, LA_N]
  integer :: in, cs, ic
  call test_begin('t_trsm')
  call la_init()
  do in = 1, size(sizes)
    do cs = 1, 3
      do ic = 1, 6
        call run(sizes(in), cs, ic)
      end do
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
  subroutine run(N, cs, ic)
    integer, intent(in) :: N, cs, ic
    integer :: m, nn, oi, oj, na, lda, ldb, i
    complex(8), allocatable :: A(:,:), B(:,:), Bref(:,:)
    complex(8) :: alpha
    real(8) :: e, e2
    character(len=64) :: lab
    m = N; nn = N; oi = 1; oj = 1
    if (cs >= 2) then
      m = N + 3; nn = N + 5
    end if
    if (cs == 3) then
      oi = 3; oj = 2
    end if
    na = merge(m, nn, sd(ic) == 'L')
    lda = na + oi - 1 + merge(0, 5, cs == 1); ldb = m + oi - 1 + merge(0, 5, cs == 1)
    allocate(A(lda, na+oj), B(ldb, nn+oj), Bref(ldb, nn+oj))
    call fill_z(lda, na+oj, A); call fill_z(ldb, nn+oj, B)
    do i = 1, na                                  ! well-conditioned triangle: |diagonal| >= 1.5
      A(oi+i-1, oj+i-1) = A(oi+i-1, oj+i-1) + cmplx(2.d0, 0.d0, kind=8)
    end do
    Bref = B
    alpha = cmplx(0.7d0, -0.3d0, kind=8)
    call ztrsm(sd(ic), up(ic), ch(tr(ic)), dg(ic), m, nn, alpha, A(oi,oj), lda, Bref(oi,oj), ldb)
    call dev_enter_z(A, size(A)); call dev_enter_z(B, size(B))
    call la_trsm(sd(ic), up(ic), tr(ic), dg(ic), m, nn, alpha, A(oi,oj), lda, B(oi,oj), ldb)
    call la_sync()
    call dev_from_z(B, size(B))
    call dev_exit_z(A, size(A)); call dev_exit_z(B, size(B))
    e = relfro(m, nn, B(oi,oj), ldb, Bref(oi,oj), ldb)
    B(oi:oi+m-1, oj:oj+nn-1) = (0.d0, 0.d0); Bref(oi:oi+m-1, oj:oj+nn-1) = (0.d0, 0.d0)
    e2 = maxdiff(ldb, nn+oj, B, Bref)
    write(lab,'(A,I0,A,I0,A,4A)') 'trsm N=', N, ' case', cs, ' ', sd(ic), up(ic), ch(tr(ic)), dg(ic)
    call check(trim(lab), e, TOL_CUDA)
    call check(trim(lab)//' outside', e2, 0.d0)
  end subroutine
end program t_trsm

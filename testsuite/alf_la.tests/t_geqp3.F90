! alf_la unit test: la_geqp3_host (STABLOG QRCP on the host, ticket 108/202) vs direct ZGEQP3.
! Both backends must be bit-identical (CUDA: the host factors are copied back to the device; the test clears the
! host copy, copies the DEVICE data back and compares).
program t_geqp3
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgeqp3
  integer :: in, cs
  call test_begin('t_geqp3')
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
    integer :: m, n2, oi, oj, lda, info, lw
    complex(8), allocatable :: A(:,:), Aref(:,:), A0(:,:), tau(:), tauref(:), w(:)
    real(8), allocatable :: rw(:)
    integer(4), allocatable :: ip(:), ipref(:)
    complex(8) :: wq(1)
    character(len=64) :: lab
    n2 = N; m = N; oi = 1; oj = 1; lda = N
    if (cs >= 2) then
      m = N + 9; lda = m + 5
    end if
    if (cs == 3) then
      oi = 3; oj = 2; lda = m + oi - 1 + 4
    end if
    allocate(A(lda, n2+oj), Aref(lda, n2+oj), A0(lda, n2+oj), tau(n2), tauref(n2), ip(n2), ipref(n2), rw(2*n2))
    call fill_z(lda, n2+oj, A)
    Aref = A; A0 = A
    ipref = 0
    call zgeqp3(m, n2, Aref(oi,oj), lda, ipref, tauref, wq, -1, rw, info)
    lw = int(dble(wq(1))); allocate(w(lw))
    call zgeqp3(m, n2, Aref(oi,oj), lda, ipref, tauref, w, lw, rw, info)
    ip = 5                                          ! must be ignored (all columns free)
    call dev_enter_z(A, size(A)); call dev_enter_z(tau, size(tau)); call dev_enter_i(ip, size(ip))
    call la_geqp3_host(m, n2, A(oi,oj), lda, ip, tau)
    call la_sync()
    if (.not. is_host()) then                       ! clear the host copies: the compare then sees device data only
      A(oi:oi+m-1, oj:oj+n2-1) = (-1.d0, -1.d0); tau = (-1.d0, -1.d0); ip = -1
    end if
    call dev_from_z(A, size(A)); call dev_from_z(tau, size(tau)); call dev_from_i(ip, size(ip))
    call dev_exit_z(A, size(A)); call dev_exit_z(tau, size(tau)); call dev_exit_i(ip, size(ip))
    write(lab,'(A,I0,A,I0)') 'geqp3_host N=', N, ' case', cs
    call check(trim(lab)//' A', relfro(m, n2, A(oi,oj), lda, Aref(oi,oj), lda), 0.d0)
    call check(trim(lab)//' tau', maxabs_diff(n2, tau, tauref), 0.d0)
    call check_int(trim(lab)//' ipvt', n2, ip, ipref)
    call check(trim(lab)//' outside', outside_diff(lda, n2+oj, A, A0, oi, oj, m, n2), 0.d0)
  end subroutine
end program t_geqp3

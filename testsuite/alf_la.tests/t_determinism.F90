! alf_la unit test (ticket 212): la_qrpp is bitwise reproducible from call to call (same binary, same input).
! Fixed inputs, including matrices with EXACTLY equal column norms (duplicated columns) and with columns that are row
! permutations of each other (equal norms up to the summation order): the pivot order, A (R and the reflectors) and tau of
! call 1 are compared bit for bit with those of calls 2..NCALL.  A device reduction with an unspecified combination order
! (the former k_colnorms) flips the pivots of near-tied columns from call to call.
! Test-only switch ALF_LA_TEST_PLANT_ULP=1 (positive control): one input element is changed by 1 ulp in one of the
! calls, so the test must fail (` FAILED t_determinism`, error stop 1).  Runs under the HOST backend too (and passes).
program t_determinism
  use alf_la_mod
  use la_test_util
  implicit none
  integer, parameter :: NCALL = 1000
  integer, parameter :: ms(4) = [7, 36, 144, 150], ns(4) = [7, 36, 144, 144]
  logical :: plant
  integer :: ic, st
  character(len=16) :: env
  env = ''
  call get_environment_variable('ALF_LA_TEST_PLANT_ULP', env, status=st)
  plant = (st == 0 .and. len_trim(env) > 0 .and. trim(env) /= '0')
  call test_begin('t_determinism')
  call la_init()
  do ic = 1, size(ms)
    call run(ms(ic), ns(ic), plant)
  end do
  call la_finalize()
  call test_end()
contains
  ! number of entries whose bit patterns differ
  integer function nbits_diff(nel, X, Y)
    integer, intent(in) :: nel
    complex(8), intent(in) :: X(nel), Y(nel)
    integer :: i
    nbits_diff = 0
    do i = 1, nel
      if (transfer(dble(X(i)), 0_8) /= transfer(dble(Y(i)), 0_8) .or. &
          transfer(aimag(X(i)), 0_8) /= transfer(aimag(Y(i)), 0_8)) nbits_diff = nbits_diff + 1
    end do
  end function

  subroutine run(m, n, plant)
    integer, intent(in) :: m, n
    logical, intent(in) :: plant
    complex(8), allocatable :: A0(:,:), A(:,:), A1(:,:), tau(:), tau1(:)
    integer(4), allocatable :: ip(:), ip1(:)
    integer :: icall, j, i, k, ndA, ndT, ndP, nbadA, nbadT, nbadP, kplant
    character(len=64) :: lab
    allocate(A0(m,n), A(m,n), A1(m,n), tau(n), tau1(n), ip(n), ip1(n))
    call fill_z(m, n, A0)
    ! exact ties: duplicated columns; near ties: a column that is a row permutation of another one
    if (n >= 12) then
      A0(:,5) = A0(:,2); A0(:,9) = A0(:,2)
      A0(:,6) = A0(:,3); A0(:,10) = A0(:,3)
      do i = 1, m
        A0(i,11) = A0(m+1-i,4)
        A0(i,12) = A0(mod(i+3,m)+1,4)
      end do
    else if (n >= 4) then
      A0(:,3) = A0(:,1); A0(:,4) = A0(:,1)
    end if
    kplant = NCALL / 2
    nbadA = 0; nbadT = 0; nbadP = 0
    do icall = 1, NCALL
      A = A0
      if (plant .and. icall == kplant) then      ! positive control: one element differs by 1 ulp in this call
        k = max(n/2, 1)
        A(max(m/2,1), k) = cmplx(nearest(dble(A(max(m/2,1), k)), 1.d0), aimag(A(max(m/2,1), k)), kind=8)
      end if
      tau = (0.d0, 0.d0); ip = 0
      call dev_enter_z(A, m*n); call dev_enter_z(tau, n); call dev_enter_i(ip, n)
      call la_qrpp(m, n, A, m, ip, tau)
      call la_sync()
      call dev_from_z(A, m*n); call dev_from_z(tau, n); call dev_from_i(ip, n)
      call dev_exit_z(A, m*n); call dev_exit_z(tau, n); call dev_exit_i(ip, n)
      if (icall == 1) then
        A1 = A; tau1 = tau; ip1 = ip
      else
        nbadA = nbadA + merge(1, 0, nbits_diff(m*n, A, A1) > 0)
        nbadT = nbadT + merge(1, 0, nbits_diff(n, tau, tau1) > 0)
        nbadP = nbadP + merge(1, 0, count(ip /= ip1) > 0)
      end if
    end do
    write(lab,'(A,I0,A,I0,A)') 'determinism ', m, 'x', n, ' calls with differing pivots'
    call check(trim(lab), dble(nbadP), 0.d0)
    write(lab,'(A,I0,A,I0,A)') 'determinism ', m, 'x', n, ' calls with differing factors A'
    call check(trim(lab), dble(nbadA), 0.d0)
    write(lab,'(A,I0,A,I0,A)') 'determinism ', m, 'x', n, ' calls with differing tau'
    call check(trim(lab), dble(nbadT), 0.d0)
    ndA = nbadA; ndT = nbadT; ndP = nbadP
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  ', m, 'x', n, ': differing calls pivots ', ndP, ' A ', ndA, ' tau ', ndT
  end subroutine
end program t_determinism

module bench_util
  use omp_lib
  use bench_la
  implicit none
  integer, parameter :: MAXREC = 8000
  character(len=40) :: r_sec(MAXREC), r_var(MAXREC)
  character(len=8)  :: r_unit(MAXREC)
  integer :: r_n(MAXREC), nrec = 0
  real(8) :: r_us(MAXREC), r_rate(MAXREC)
  integer :: nfail = 0, nchk = 0
  abstract interface
    subroutine op_i()
    end subroutine
  end interface
contains

  real(8) function now()
    now = omp_get_wtime()
  end function

  ! Warm-up call, then 21 timed repetitions (11 if one call takes > 50 ms, 5 if > 0.5 s); median in
  ! microseconds.
  ! Each repetition is bracketed by device synchronisation.
  subroutine timeit(op, us)
    procedure(op_i) :: op
    real(8), intent(out) :: us
    real(8) :: t0, t(21)
    integer :: r, nr
    call la_sync(); t0 = now(); call op(); call la_sync()
    t(1) = now() - t0
    nr = 21
    if (t(1) > 0.05d0) nr = 11
    if (t(1) > 0.5d0) nr = 5
    do r = 1, nr
      call la_sync(); t0 = now()
      call op()
      call la_sync()
      t(r) = now() - t0
    end do
    us = 1d6 * median(t, nr)
  end subroutine

  real(8) function median(t, n)
    real(8), intent(in) :: t(*)
    integer, intent(in) :: n
    real(8) :: s(n), x
    integer :: i, j
    s = t(1:n)
    do i = 2, n
      x = s(i); j = i - 1
      do while (j >= 1)
        if (s(j) <= x) exit
        s(j+1) = s(j); j = j - 1
      end do
      s(j+1) = x
    end do
    median = s((n+1)/2)
  end function

  subroutine add_rec(sec, var, n, us, rate, unit)
    character(len=*), intent(in) :: sec, var, unit
    integer, intent(in) :: n
    real(8), intent(in) :: us, rate
    nrec = nrec + 1
    if (nrec > MAXREC) error stop 'too many records'
    r_sec(nrec) = sec; r_var(nrec) = var; r_n(nrec) = n
    r_us(nrec) = us; r_rate(nrec) = rate; r_unit(nrec) = unit
  end subroutine

  ! PASS/FAIL line; fail if val > tol or val is NaN
  subroutine check(name, n, val, tol)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(8), intent(in) :: val, tol
    character(len=4) :: s
    nchk = nchk + 1
    if (val <= tol) then
      s = 'PASS'
    else
      s = 'FAIL'; nfail = nfail + 1
    end if
    write(*,'(a,1x,a56,1x,a,i5,a,es10.3,a,es9.2)') s, name, 'N=', n, '  value=', val, '  tol=', tol
  end subroutine

  subroutine report(name, n, val)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(8), intent(in) :: val
    write(*,'(a,1x,a56,1x,a,i5,a,es10.3)') 'INFO', name, 'N=', n, '  value=', val
  end subroutine

  real(8) function relfro(A, B, n)
    integer, intent(in) :: n
    complex(8), intent(in) :: A(n), B(n)
    real(8) :: d, s
    integer :: i
    d = 0d0; s = 0d0
    do i = 1, n
      d = d + abs(A(i) - B(i))**2
      s = s + abs(B(i))**2
    end do
    relfro = sqrt(d / max(s, tiny(1d0)))
  end function

  subroutine rand_c(A, n, scale)
    integer, intent(in) :: n
    complex(8), intent(out) :: A(n)
    real(8), intent(in) :: scale
    real(8), allocatable :: x(:), y(:)     ! heap, not automatic: 2*N*N reals overflow ifx's 8 MB stack at N=576
    allocate(x(n), y(n))
    call random_number(x); call random_number(y)
    A = scale * cmplx(x - 0.5d0, y - 0.5d0, kind=8)
  end subroutine

  subroutine print_tables(u)
    integer, intent(in) :: u
    character(len=40) :: secs(64), vars(32)
    integer :: ns, nv, nn, ns_i, i, j, k, ns2, ns_n(40)
    integer :: ns_list(64)
    integer :: nlist(32), tab
    logical :: found
    real(8) :: v
    character(len=8) :: unit
    ns = 0
    do i = 1, nrec
      found = .false.
      do j = 1, ns
        if (secs(j) == r_sec(i)) found = .true.
      end do
      if (.not. found) then; ns = ns + 1; secs(ns) = r_sec(i); end if
    end do
    do ns_i = 1, ns
      nv = 0; nn = 0
      unit = ''
      do i = 1, nrec
        if (r_sec(i) /= secs(ns_i)) cycle
        found = .false.
        do j = 1, nv
          if (vars(j) == r_var(i)) found = .true.
        end do
        if (.not. found) then; nv = nv + 1; vars(nv) = r_var(i); end if
        found = .false.
        do j = 1, nn
          if (nlist(j) == r_n(i)) found = .true.
        end do
        if (.not. found) then; nn = nn + 1; nlist(nn) = r_n(i); end if
        if (r_rate(i) > 0d0) unit = r_unit(i)
      end do
      do tab = 1, 2
        if (tab == 1) then
          write(u,'(/,a,a,a)') '== ', trim(secs(ns_i)), ' : median time [us]'
        else
          if (len_trim(unit) == 0) cycle
          write(u,'(/,a,a,a,a,a)') '== ', trim(secs(ns_i)), ' : ', trim(unit), ' (effective)'
        end if
        write(u,'(a6)', advance='no') 'N'
        do j = 1, nv
          write(u,'(1x,a18)', advance='no') adjustr(vars(j)(1:18))
        end do
        write(u,*)
        do k = 1, nn
          write(u,'(i6)', advance='no') nlist(k)
          do j = 1, nv
            v = -1d0
            do i = 1, nrec
              if (r_sec(i) == secs(ns_i) .and. r_var(i) == vars(j) .and. r_n(i) == nlist(k)) then
                if (tab == 1) then; v = r_us(i); else; v = r_rate(i); end if
              end if
            end do
            if (v < 0d0) then
              write(u,'(1x,a18)', advance='no') '-'
            else if (v >= 1d5) then
              write(u,'(1x,f18.0)', advance='no') v
            else
              write(u,'(1x,f18.2)', advance='no') v
            end if
          end do
          write(u,*)
        end do
      end do
    end do
  end subroutine

end module bench_util

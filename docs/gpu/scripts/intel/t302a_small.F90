! 302a development: cost of tiny solver calls (t_order does 8292 8x8 geqrf per pass)
program t302a_small
  use alf_la_mod
  implicit none
  complex(8) :: A(8,8), tau(8)
  integer :: i, n
  integer(8) :: t0, t1, rate
  real(8) :: r = 0.d0
  call system_clock(count_rate=rate)
  n = 2000
  call la_init()
  A = (1.d0, 0.1d0); do i = 1, 8; A(i,i) = A(i,i) + 5.d0*i; end do; tau = 0.d0
  !$omp target enter data map(to: A, tau)
  call la_geqrf(8, 8, A, 8, tau)
  call la_sync()
  call system_clock(t0)
  do i = 1, n
    call la_geqrf(8, 8, A, 8, tau)
  end do
  call la_sync()
  call system_clock(t1)
  print '(A,F10.1,A)', 'la_geqrf 8x8: ', 1.d6*dble(t1-t0)/rate/n, ' us/call'
  call system_clock(t0)
  do i = 1, n
    r = r + la_geqrf_lwork(8, 8, A, 8)
  end do
  call system_clock(t1)
  print '(A,F10.1,A)', 'la_geqrf_lwork 8x8: ', 1.d6*dble(t1-t0)/rate/n, ' us/call', r
  call system_clock(t0)
  do i = 1, n
    call la_gemm(LA_N, LA_N, 8, 8, 8, (1.d0,0.d0), A, 8, A, 8, (0.d0,0.d0), tau, 8)
  end do
  call la_sync()
  call system_clock(t1)
  print '(A,F10.1,A)', 'la_gemm 8x8x8 (writes into tau region): ', 1.d6*dble(t1-t0)/rate/n, ' us/call'
  !$omp target exit data map(delete: A, tau)
  call la_finalize()
end program

! Intel batched study, Step 1(c): launch latency of an OpenMP target region from ifx (the path of the ALF omp kernels).
! sync = target region (implicitly waits); nowait = target nowait + one taskwait at the end.
program ib_omp
  use omp_lib
  implicit none
  integer, parameter :: nrep = 20000, n = 64
  integer :: a(n), i, r
  real(8) :: t0, ts, ta
  a = 0
  !$omp target enter data map(to: a)
  do r = 1, 200
    !$omp target teams distribute parallel do
    do i = 1, n
      a(i) = i + r
    end do
  end do
  t0 = omp_get_wtime()
  do r = 1, nrep
    !$omp target teams distribute parallel do
    do i = 1, n
      a(i) = i + r
    end do
  end do
  ts = (omp_get_wtime() - t0) / nrep
  t0 = omp_get_wtime()
  !$omp parallel
  !$omp single
  do r = 1, nrep
    !$omp target teams distribute parallel do nowait
    do i = 1, n
      a(i) = i + r
    end do
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel
  ta = (omp_get_wtime() - t0) / nrep
  !$omp target exit data map(from: a)
  write(*, '(a,f8.2,a,f8.2,a,i0)') 'omp sync_us=', ts * 1d6, ' nowait_us=', ta * 1d6, ' check=', a(1)
end program ib_omp

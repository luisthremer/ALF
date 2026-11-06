program fma
  use omp_lib
  implicit none
  integer, parameter :: n = 2**20, nit = 2000
  real(8), allocatable :: a8(:)
  real(4), allocatable :: a4(:)
  real(8) :: t0, t1, x8
  real(4) :: x4
  integer :: i, k, rep
  allocate(a8(n), a4(n)); a8 = 1.0d0; a4 = 1.0
  !$omp target enter data map(to: a8, a4)
  do rep = 1, 3
    !$omp target teams distribute parallel do private(k,x8)
    do i = 1, n
      x8 = a8(i)
      do k = 1, 16
        x8 = x8 * 0.999999d0 + 1.0d-9
      enddo
      a8(i) = x8
    enddo
  enddo
  t0 = omp_get_wtime()
  do rep = 1, 5
    !$omp target teams distribute parallel do private(k,x8)
    do i = 1, n
      x8 = a8(i)
      do k = 1, nit
        x8 = x8 * 0.999999d0 + 1.0d-9
      enddo
      a8(i) = x8
    enddo
  enddo
  t1 = omp_get_wtime()
  print '(A,F10.1,A)', 'FP64 FMA throughput: ', 2.d0*n*nit*5/(t1-t0)/1.d9, ' GFlop/s'
  do rep = 1, 3
    !$omp target teams distribute parallel do private(k,x4)
    do i = 1, n
      x4 = a4(i)
      do k = 1, 16
        x4 = x4 * 0.999999 + 1.0e-9
      enddo
      a4(i) = x4
    enddo
  enddo
  t0 = omp_get_wtime()
  do rep = 1, 5
    !$omp target teams distribute parallel do private(k,x4)
    do i = 1, n
      x4 = a4(i)
      do k = 1, nit
        x4 = x4 * 0.999999 + 1.0e-9
      enddo
      a4(i) = x4
    enddo
  enddo
  t1 = omp_get_wtime()
  print '(A,F10.1,A)', 'FP32 FMA throughput: ', 2.d0*n*nit*5/(t1-t0)/1.d9, ' GFlop/s'
  !$omp target exit data map(from: a8, a4)
  print *, a8(1) > 0.d0, a4(1) > 0.0
end program

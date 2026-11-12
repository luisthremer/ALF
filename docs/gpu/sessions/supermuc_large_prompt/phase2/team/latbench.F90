include "mkl_blas_omp_offload.f90"
! Latency microbenchmark for OpenMP offload on PVC (ifx): cost per operation of the primitives ALF's sweep is made of.
program latbench
  use omp_lib
  use iso_c_binding
  use onemkl_blas_omp_offload_lp64, only: zgemm
  implicit none
  integer, parameter :: NREP = 2000
  integer :: i, rc, n, k, kk
  real(8) :: t0, t1, dummy
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), s(:)
  complex(8) :: one, zero
  real(8), allocatable :: x(:)
  integer :: ndev
  ndev = omp_get_num_devices()
  call omp_set_default_device(0)
  one = (1.d0,0.d0); zero = (0.d0,0.d0)
  n = 1024; k = 32
  allocate(A(n,n), B(n,k), C(n,k), s(64), x(n)); A = one; B = one; C = zero; s = one; x = 1.d0
  !$omp target enter data map(to: A, B, C, s, x)
  ! warm-up
  do i = 1, 50
    !$omp target
    dummy = 0.d0
    !$omp end target
  end do
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target
    dummy = 0.d0
    !$omp end target
  end do
  t1 = omp_get_wtime(); call rep('empty target region (sync)', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target nowait
    dummy = 0.d0
    !$omp end target
  end do
  !$omp taskwait
  t1 = omp_get_wtime(); call rep('empty target nowait + 1 taskwait', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target update to(s(1:1))
  end do
  t1 = omp_get_wtime(); call rep('target update to, 16 B', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target update to(s(1:64))
  end do
  t1 = omp_get_wtime(); call rep('target update to, 1 KB', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target update from(s(1:1))
  end do
  t1 = omp_get_wtime(); call rep('target update from, 16 B', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target update from(B(1:n,1:k))
  end do
  t1 = omp_get_wtime(); call rep('target update from, 512 KB', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target teams distribute parallel do
    do rc = 1, n
      x(rc) = x(rc) + 1.d0
    end do
  end do
  t1 = omp_get_wtime(); call rep('small kernel (1024 elems)', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp target teams distribute parallel do collapse(2)
    do rc = 1, 1024
      do kk = 1, 1024
        A(kk,rc) = A(kk,rc) * 1.0000001d0
      end do
    end do
  end do
  t1 = omp_get_wtime(); call rep('1M-element kernel (16 MB rw)', t1-t0)
  n = 1024
  do i = 1, 20
    !$omp dispatch
    call zgemm('N','N',n,k,k,one,A,n,B,n,zero,C,n)
  end do
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp dispatch
    call zgemm('N','N',n,k,k,one,A,n,B,n,zero,C,n)
  end do
  t1 = omp_get_wtime(); call rep('dispatch zgemm 1024x32x32 (incl. implicit sync)', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp dispatch
    call zgemm('N','T',n,n,k,one,B,n,B,n,zero,A,n)
  end do
  t1 = omp_get_wtime(); call rep('dispatch zgemm 1024x1024x32', t1-t0)
  t0 = omp_get_wtime()
  do i = 1, NREP
    !$omp dispatch
    call zgemm('N','N',n,n,n,one,A,n,A,n,zero,A,n)
  end do
  t1 = omp_get_wtime(); call rep('dispatch zgemm 1024^3', t1-t0)
contains
  subroutine rep(name, t)
    character(*), intent(in) :: name
    real(8), intent(in) :: t
    write(*,'(A,T58,F9.2,A)') name, 1.d6*t/NREP, ' us/op'
  end subroutine
end program

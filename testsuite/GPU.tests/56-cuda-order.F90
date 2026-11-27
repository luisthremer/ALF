!> @brief Canary of the CUDA backend (GPU=nvidia): cuBLAS/cuSOLVER calls of linalg_mod and the caller's synchronous OpenMP
!> target regions must run in order on one stream (linalg_cuda_mod header).
!> @details Ported from gpu-port testsuite/alf_la.tests/t_order.F90. RAW (region writes A, library reads it), WAR (library
!> reads A, region overwrites it), WAW (library writes C, region overwrites it), 24 repetitions each with N = 1024, no
!> host-blocking read in between except the reductions that check the result. Then 2*4096+100 small device QRs with one
!> la_cuda_sync at the end (info slot read-back). Must be re-run for every compiler/driver version.
!> Negative control: built with -DALF_LA_TEST_BADSTREAM the backend uses a created stream that is not ordered against the
!> target regions, and this test must FAIL.
!> Device kernels of the canary (module procedures: nvfortran does not offload internal procedures of a program)
module cuda_order_kernels
  implicit none
  integer, parameter :: N = 1024
  complex(8), parameter :: z0 = (0.d0, 0.d0)
contains
  subroutine set_diag(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
    !$omp target teams distribute parallel do collapse(2)
    do jj = 1, N
       do ii = 1, N
          if (ii == jj) then                    ! no merge on complex: nvfortran drops such a kernel (CDSELECT)
             M(ii,jj) = cmplx(v, 0.d0, kind=8)
          else
             M(ii,jj) = z0
          endif
       enddo
    enddo
  end subroutine set_diag
  subroutine fill_const(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
    !$omp target teams distribute parallel do collapse(2)
    do jj = 1, N
       do ii = 1, N
          M(ii,jj) = cmplx(v, 0.d0, kind=8)
       enddo
    enddo
  end subroutine fill_const
  real(8) function sum_diag(M)
    complex(8), intent(in) :: M(N,N)
    integer :: ii
    real(8) :: s
    s = 0.d0
    !$omp target teams distribute parallel do reduction(+: s)
    do ii = 1, N
       s = s + dble(M(ii,ii))
    enddo
    sum_diag = s
  end function sum_diag
  real(8) function sum_all(M)
    complex(8), intent(in) :: M(N,N)
    integer :: ii, jj
    real(8) :: s
    s = 0.d0
    !$omp target teams distribute parallel do collapse(2) reduction(+: s)
    do jj = 1, N
       do ii = 1, N
          s = s + dble(M(ii,jj))
       enddo
    enddo
    sum_all = s
  end function sum_all
end module cuda_order_kernels

program cuda_order
  use linalg_mod, only: la_zgemm, la_zgeqrf
  use linalg_cuda_mod, only: la_cuda_sync, la_cuda_finalize
  use cuda_order_kernels, only: N, set_diag, fill_const, sum_diag, sum_all
  implicit none
  integer, parameter :: NREP = 24, NSLOT = 4096
  complex(8), parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)
  integer :: i, rep, info, nfail
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), C1(:,:), small(:,:), tau(:), w(:)
  real(8) :: acc(NREP, 2), bsum, worst

  nfail = 0
  allocate(A(N,N), B(N,N), C(N,N), C1(N,N))
  do i = 1, N
     do rep = 1, N
        B(rep,i) = cmplx(sin(0.37d0*rep + 1.1d0*i), cos(0.53d0*rep - 0.7d0*i), kind=8)
     enddo
  enddo
  A = z0; C = z0; C1 = z0
  bsum = 0.d0
  do i = 1, N
     bsum = bsum + dble(B(i,i))
  enddo
  !$omp target enter data map(to: A, B, C, C1)
  ! RAW
  do rep = 1, NREP
     call set_diag(A, dble(rep))
     call la_zgemm('N', 'N', N, N, N, z1, A, N, B, N, z0, C, N)
     acc(rep,1) = sum_diag(C)
  enddo
  worst = 0.d0
  do rep = 1, NREP
     worst = max(worst, abs(acc(rep,1) - dble(rep)*bsum)/abs(dble(rep)*bsum))
  enddo
  call check('RAW', worst)
  ! WAR
  do rep = 1, NREP
     call set_diag(A, 1.d0)
     call la_zgemm('N', 'N', N, N, N, z1, A, N, B, N, z0, C1, N)
     call set_diag(A, 2.d0)
     call la_zgemm('N', 'N', N, N, N, z1, A, N, B, N, z0, C, N)
     acc(rep,1) = sum_diag(C1)
     acc(rep,2) = sum_diag(C)
  enddo
  worst = 0.d0
  do rep = 1, NREP
     worst = max(worst, abs(acc(rep,1) - bsum)/abs(bsum), abs(acc(rep,2) - 2.d0*bsum)/abs(2.d0*bsum))
  enddo
  call check('WAR', worst)
  ! WAW
  call set_diag(A, 3.d0)
  do rep = 1, NREP
     call la_zgemm('N', 'N', N, N, N, z1, A, N, B, N, z0, C, N)
     call fill_const(C, 7.d0)
     acc(rep,1) = sum_all(C)
  enddo
  worst = 0.d0
  do rep = 1, NREP
     worst = max(worst, abs(acc(rep,1) - 7.d0*dble(N)*dble(N))/(7.d0*dble(N)*dble(N)))
  enddo
  call check('WAW', worst)
  call la_cuda_sync()
  !$omp target exit data map(delete: A, B, C, C1)
  ! info slot traffic: more QRs than slots, one sync at the end
  allocate(small(8,8), tau(8), w(1))
  small = B(1:8,1:8); tau = z0
  !$omp target enter data map(to: small, tau)
  do rep = 1, 2*NSLOT + 100
     call la_zgeqrf(8, 8, small, 8, tau, w, 1, info)
  enddo
  call la_cuda_sync()
  !$omp target exit data map(delete: small, tau)
  write(*,'(A)') 'info slot traffic: ok'
  call la_cuda_finalize()
  if (nfail > 0) then
     write(*,'(A,I0,A)') '56-cuda-order: FAIL (', nfail, ' hazards)'
     error stop 1
  endif
  write(*,'(A)') '56-cuda-order: PASS'
contains
  subroutine check(what, w)
    character(len=*), intent(in) :: what
    real(8), intent(in) :: w
    write(*,'(A,A,ES10.3)') what, ': worst relative error ', w
    if (.not. (w <= 1.d-12)) nfail = nfail + 1
  end subroutine check
end program cuda_order

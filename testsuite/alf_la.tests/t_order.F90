! alf_la unit test: ordering between asynchronous alf_la library calls and the caller's own target regions
! (canary for the stream decision in alf_la_mod.F90; re-run for every compiler / driver version).
! The host waits for the caller's synchronous target regions but NEVER for library calls; besides the region waits the
! only implicit device->host transfers are the reduction scalars (sum_diag / sum_all), which come after each hazard chain.
! No la_sync or other blocking call sits between the operations of a hazard chain:
!   RAW : A = r I (region) ; C = A B (library) ; sum diag C (region, reduction)
!   WAR : C1 = A B (library, A = I) ; A = 2 I (region, right away) ; C2 = A B ; sums of diag C1 and diag C2
!         (if the region overwrote A before the library call read it, C1 would be 2 B)
!   WAW : C = A B (library) ; C = 7 (region, right away) ; sum C   (if the library call finished after the
!         region, C would not be 7)
! Also: more than LA_NSLOT solver calls and la_init / la_finalize twice.
!
! NEGATIVE CONTROL (not registered in ctest): build the module and this test with -DALF_LA_TEST_BADSTREAM
! (la_init then uses a cudaStreamCreate stream) -- t_order must fail on RAW, WAR and WAW. Recipe:
!   docs/gpu/scripts/nividia/t202_negctl.sh
program t_order
  use alf_la_mod
  use la_test_util
  implicit none
  integer, parameter :: N = 1024, NREP = 24
  integer :: i, pass, rep, ns
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), C1(:,:), small(:,:), tau(:)
  real(8) :: acc(NREP, 2), bsum, worst, expect
  call test_begin('t_order')
  do pass = 1, 2
    call la_init()
    allocate(A(N,N), B(N,N), C(N,N), C1(N,N))
    call fill_z(N, N, B)
    A = (0.d0, 0.d0); C = (0.d0, 0.d0); C1 = (0.d0, 0.d0); acc = 0.d0
    bsum = 0.d0
    do i = 1, N
      bsum = bsum + dble(B(i,i))
    end do
    call dev_enter_z(A, N*N); call dev_enter_z(B, N*N); call dev_enter_z(C, N*N); call dev_enter_z(C1, N*N)
    ! ---- RAW
    do rep = 1, NREP
      call set_diag(A, dble(rep))
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
      call sum_diag(C, rep, 1)
    end do
    worst = 0.d0
    do rep = 1, NREP
      expect = dble(rep) * bsum
      worst = fmax(worst, abs(acc(rep,1) - expect) / abs(expect))
    end do
    call check('ordering RAW pass ' // achar(48+pass), worst, 1.d-12)
    ! ---- WAR
    acc = 0.d0
    do rep = 1, NREP
      call set_diag(A, 1.d0)
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C1, N)
      call set_diag(A, 2.d0)
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
      call sum_diag(C1, rep, 1)
      call sum_diag(C, rep, 2)
    end do
    worst = 0.d0
    do rep = 1, NREP
      worst = fmax(worst, abs(acc(rep,1) - bsum) / abs(bsum))
      worst = fmax(worst, abs(acc(rep,2) - 2.d0*bsum) / abs(2.d0*bsum))
    end do
    call check('ordering WAR pass ' // achar(48+pass), worst, 1.d-12)
    ! ---- WAW
    acc = 0.d0
    call set_diag(A, 3.d0)
    do rep = 1, NREP
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
      call fill_const(C, 7.d0)
      call sum_all(C, rep)
    end do
    worst = 0.d0
    do rep = 1, NREP
      worst = fmax(worst, abs(acc(rep,1) - 7.d0*dble(N)*dble(N)) / (7.d0*dble(N)*dble(N)))
    end do
    call check('ordering WAW pass ' // achar(48+pass), worst, 1.d-12)
    call la_sync()
    call dev_exit_z(A, N*N); call dev_exit_z(B, N*N); call dev_exit_z(C, N*N); call dev_exit_z(C1, N*N)
    deallocate(A, B, C, C1)
    ! info-slot traffic: 2*LA_NSLOT+100 tiny QR factorizations, one la_sync at the end
    ns = 2*LA_NSLOT + 100
    allocate(small(8,8), tau(8))
    call fill_z(8, 8, small); tau = (0.d0, 0.d0)
    call dev_enter_z(small, 64); call dev_enter_z(tau, 8)
    do rep = 1, ns
      call la_geqrf(8, 8, small, 8, tau)
    end do
    call la_sync()
    call dev_exit_z(small, 64); call dev_exit_z(tau, 8)
    deallocate(small, tau)
    call check('info slot traffic', 0.d0, 0.d0)
    call la_finalize()
  end do
  call test_end()
contains
  subroutine set_diag(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
#ifdef ALF_LA_CUDA
    !$omp target teams distribute parallel do collapse(2)
#endif
    do jj = 1, N
      do ii = 1, N
        if (ii == jj) then
          M(ii,jj) = cmplx(v, 0.d0, kind=8)
        else
          M(ii,jj) = (0.d0, 0.d0)
        end if
      end do
    end do
  end subroutine
  subroutine fill_const(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
#ifdef ALF_LA_CUDA
    !$omp target teams distribute parallel do collapse(2)
#endif
    do jj = 1, N
      do ii = 1, N
        M(ii,jj) = cmplx(v, 0.d0, kind=8)
      end do
    end do
  end subroutine
  subroutine sum_diag(M, r, k)
    complex(8), intent(in) :: M(N,N)
    integer, intent(in) :: r, k
    integer :: ii
    real(8) :: s
    s = 0.d0
#ifdef ALF_LA_CUDA
    !$omp target teams distribute parallel do reduction(+: s)
#endif
    do ii = 1, N
      s = s + dble(M(ii,ii))
    end do
    acc(r,k) = s
  end subroutine
  subroutine sum_all(M, r)
    complex(8), intent(in) :: M(N,N)
    integer, intent(in) :: r
    integer :: ii, jj
    real(8) :: s
    s = 0.d0
#ifdef ALF_LA_CUDA
    !$omp target teams distribute parallel do collapse(2) reduction(+: s)
#endif
    do jj = 1, N
      do ii = 1, N
        s = s + dble(M(ii,jj))
      end do
    end do
    acc(r,1) = s
  end subroutine
end program t_order

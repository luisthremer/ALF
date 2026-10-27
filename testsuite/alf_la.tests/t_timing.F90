! alf_la timing table (ticket 202; not a ctest). Per-call time in microseconds, library calls pipelined with one la_sync
! at the end ("pipe") and with la_sync after every call ("sync"). Compare with docs/gpu/BENCHMARKS.md (ticket 101).
! In-place routines (geqrf, getrf, ...) are re-applied to their own output: the work is data-independent.
program t_timing
  use alf_la_mod
  use la_test_util
  implicit none
  integer, parameter :: ns = 2, nsz(ns) = [400, 1024]
  integer :: is, N
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), x(:), y(:), tau(:), Q(:,:), A2(:,:)
  real(8), allocatable :: Br(:,:)
  integer(4), allocatable :: ip(:), ipq(:)
  complex(8) :: one, zero, a1
  one = (1.d0, 0.d0); zero = (0.d0, 0.d0); a1 = (0.5d0, 0.1d0)
  write(*,'(A,A)') 'backend: ', trim(la_backend_name())
  write(*,'(A)') 'routine                          N     pipe_us    sync_us'
  call la_init()
  do is = 1, ns
    N = nsz(is)
    call setup()
    call bench('gemm NN (ZGEMM)', 20, 1)
    call bench('2x gemm (dense B_T wrap)', 20, 2)
    call bench('gemm_cr (DGEMM real view)', 20, 3)
    call bench('geru', 200, 4)
    call bench('trsm L,U,N,N', 20, 5)
    call bench('geqrf', 5, 6)
    call bench('qrpp + ungqr', 5, 7)
    call bench('ormqr L,C', 5, 8)
    call bench('getrf', 5, 9)
    call bench('getrs N', 5, 10)
    call bench('geqp3_host (+transfers)', 3, 11)
    call teardown()
  end do
  call la_finalize()
contains
  real(8) function now()
    integer(8) :: c, r
    call system_clock(c, r)
    now = dble(c) / dble(r) * 1.d6
  end function

  subroutine setup()
    integer :: i
    allocate(A(N,N), B(N,N), C(N,N), x(N), y(N), tau(N), Q(N,N), A2(N,N), Br(N,N), ip(N), ipq(N))
    call fill_z(N, N, A); call fill_z(N, N, B); call fill_z(N, N, C); call fill_z(N, 1, x); call fill_z(N, 1, y)
    call fill_z(N, N, Q); call fill_z(N, N, A2)
    do i = 1, N
      A(i,i) = A(i,i) + cmplx(2.d0 * sqrt(dble(N)), 0.d0, kind=8)
      B(i,i) = B(i,i) + cmplx(2.d0 * sqrt(dble(N)), 0.d0, kind=8)   ! triangular solves / LU stay tame
      Br(:, i) = dble(B(:, i))
    end do
    tau = zero; ip = 0; ipq = 0
    call dev_enter_z(A, N*N); call dev_enter_z(B, N*N); call dev_enter_z(C, N*N); call dev_enter_z(x, N)
    call dev_enter_z(y, N); call dev_enter_z(tau, N); call dev_enter_z(Q, N*N); call dev_enter_d(Br, N*N)
    call dev_enter_i(ip, N); call dev_enter_i(ipq, N); call dev_enter_z(A2, N*N)
  end subroutine

  subroutine teardown()
    call dev_exit_z(A, N*N); call dev_exit_z(B, N*N); call dev_exit_z(C, N*N); call dev_exit_z(x, N)
    call dev_exit_z(y, N); call dev_exit_z(tau, N); call dev_exit_z(Q, N*N); call dev_exit_d(Br, N*N)
    call dev_exit_i(ip, N); call dev_exit_i(ipq, N); call dev_exit_z(A2, N*N)
    deallocate(A, B, C, x, y, tau, Q, A2, Br, ip, ipq)
  end subroutine

  subroutine one_call(k)
    integer, intent(in) :: k
    select case (k)
    case (1); call la_gemm(LA_N, LA_N, N, N, N, one, A, N, B, N, zero, C, N)
    case (2)
      call la_gemm(LA_N, LA_N, N, N, N, one, A, N, B, N, zero, C, N)
      call la_gemm(LA_N, LA_N, N, N, N, one, C, N, B, N, zero, Q, N)
    case (3); call la_gemm_cr(N, N, N, A, N, Br, N, C, N)
    case (4); call la_geru(N, N, a1, x, 1, y, 1, C, N)
    case (5); call la_trsm('L', 'U', LA_N, 'N', N, N, one, B, N, C, N)
    case (6); call la_geqrf(N, N, A2, N, tau)
    case (7)
      call la_qrpp(N, N, A2, N, ipq, tau)
      call la_ungqr(N, N, N, A2, N, tau)
    case (8); call la_ormqr('L', LA_C, N, N, N, A2, N, tau, C, N)
    case (9); call la_getrf(N, N, B, N, ip)
    case (10); call la_getrs(LA_N, N, N, B, N, ip, C, N)
    case (11); call la_geqp3_host(N, N, A2, N, ipq, tau)
    end select
  end subroutine

  subroutine bench(name, reps, k)
    character(len=*), intent(in) :: name
    integer, intent(in) :: reps, k
    integer :: r
    real(8) :: t0, tp, ts
    call one_call(k); call la_sync()                 ! warm-up (workspace growth, library init)
    t0 = now()
    do r = 1, reps
      call one_call(k)
    end do
    call la_sync()
    tp = (now() - t0) / dble(reps)
    t0 = now()
    do r = 1, reps
      call one_call(k)
      call la_sync()
    end do
    ts = (now() - t0) / dble(reps)
    write(*,'(A30,I6,2F12.1)') name, N, tp, ts
  end subroutine
end program t_timing

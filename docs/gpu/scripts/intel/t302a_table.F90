! 302a per-routine table (acceptance 3, measured; development tool, not a test).
! For every public alf_la routine and N = 144, 400, 1024: wall time per call of the alf_la call on the device (la_sync
! after every call), CPU time / wall time of the process during those calls (a pure device call keeps this near 0-1
! cores, host work shows more), and the same operation as a direct host MKL call (host-pinned reference).
! With argument "trace" every routine runs once per size between "SECTION" markers, for an MKL_VERBOSE=1 trace.
program t302a_table
  use alf_la_mod
  implicit none
  integer, parameter :: ns = 3, nsz(ns) = [144, 400, 1024]
  integer, parameter :: nr = 14
  character(len=28), parameter :: names(nr) = [character(len=28) :: 'la_gemm NN', 'la_gemm ZGEMM NC', 'la_gemm_cr', &
       'la_gemm_batched (nb=4)', 'la_geru', 'la_trsm L,U,N,N', 'la_geqrf', 'la_ungqr', 'la_ormqr L,C', 'la_getrf', &
       'la_getrs N', 'la_qrpp', 'la_geqp3_host', 'la_geqrf_lwork+la_ormqr_lwork']
  integer :: is, N, k
  logical :: trace
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), x(:), y(:), tau(:), A2(:,:), Ah(:,:), Ch(:,:), Bb(:), Ab(:), Cb(:)
  real(8), allocatable :: Br(:,:), Cr(:,:), rw(:)
  complex(8), allocatable :: hw(:)
  integer(4), allocatable :: ip(:), ipq(:)
  complex(8) :: one, zero, a1
  character(len=16) :: arg
  external :: zgemm, dgemm, zgeru, ztrsm, zgeqrf, zungqr, zunmqr, zgetrf, zgetrs, zgeqp3
  one = (1.d0, 0.d0); zero = (0.d0, 0.d0); a1 = (0.5d0, 0.1d0)
  call get_command_argument(1, arg)
  trace = (trim(arg) == 'trace')
  call la_init()
  do is = 1, ns
    N = nsz(is)
    call setup()
    do k = 1, nr - 1
      call bench(k)
    end do
    call teardown()
  end do
  call la_finalize()
contains
  real(8) function now()
    integer(8) :: c, r
    call system_clock(c, r)
    now = dble(c) / dble(r)
  end function

  subroutine fill(M)
    complex(8), intent(out) :: M(:,:)
    integer :: ii, jj
    integer(8), save :: r = 20261030_8
    real(8) :: a, b
    do jj = 1, size(M,2)
      do ii = 1, size(M,1)
        r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); a = dble(r) / 2147483648.d0 - 0.5d0
        r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); b = dble(r) / 2147483648.d0 - 0.5d0
        M(ii,jj) = cmplx(a, b, kind=8)
      end do
    end do
  end subroutine

  subroutine setup()
    integer :: i
    allocate(A(N,N), B(N,N), C(N,N), x(N), y(N), tau(N), A2(N,N), Ah(N,N), Ch(N,N), Br(N,N), Cr(N,N), ip(N), ipq(N))
    allocate(Ab(4*N*N), Bb(4*N*N), Cb(4*N*N), rw(2*N))
    call fill(A); call fill(B); call fill(C); call fill(A2)
    x = A(:,1); y = B(:,1)
    do i = 0, 4*N-1
      Ab(1+i*N:(i+1)*N) = A2(:, 1 + mod(i, N)); Bb(1+i*N:(i+1)*N) = B(:, 1 + mod(i+1, N)); Cb(1+i*N:(i+1)*N) = C(:, 1 + mod(i+2, N))
    end do
    do i = 1, N
      A(i,i) = A(i,i) + cmplx(2.d0 * sqrt(dble(N)), 0.d0, kind=8)
      B(i,i) = B(i,i) + cmplx(2.d0 * sqrt(dble(N)), 0.d0, kind=8)
    end do
    Br = dble(B); Cr = 0.d0; tau = zero; ip = 0; ipq = 0
    !$omp target enter data map(to: A, B, C, x, y, tau, A2, Br, ip, ipq, Ab, Bb, Cb)
  end subroutine

  subroutine teardown()
    !$omp target exit data map(delete: A, B, C, x, y, tau, A2, Br, ip, ipq, Ab, Bb, Cb)
    deallocate(A, B, C, x, y, tau, A2, Ah, Ch, Br, Cr, ip, ipq, Ab, Bb, Cb, rw)
  end subroutine

  ! one alf_la call (device) or one host-pinned call (host=.true.), operation k
  subroutine one_call(k, host)
    integer, intent(in) :: k
    logical, intent(in) :: host
    integer :: info, lw
    complex(8) :: wq(1)
    if (.not. host) then
      select case (k)
      case (1); call la_gemm(LA_N, LA_N, N, N, N, one, A, N, B, N, zero, C, N)
      case (2); call la_gemm(LA_N, LA_C, N, N, N, one, A, N, B, N, zero, C, N)
      case (3); call la_gemm_cr(N, N, N, A, N, Br, N, C, N)
      case (4); call la_gemm_batched(LA_N, LA_N, N/2, N/2, N/2, one, Ab, N, N*N/2, Bb, N, N*N/2, zero, Cb, N, N*N/2, 4)
      case (5); call la_geru(N, N, a1, x, 1, y, 1, C, N)
      case (6); call la_trsm('L', 'U', LA_N, 'N', N, N, one, B, N, C, N)
      case (7); call la_geqrf(N, N, A2, N, tau)
      case (8); call la_ungqr(N, N, N, A2, N, tau)
      case (9); call la_ormqr('L', LA_C, N, N, N, A2, N, tau, C, N)
      case (10); call la_getrf(N, N, B, N, ip)
      case (11); call la_getrs(LA_N, N, N, B, N, ip, C, N)
      case (12); call la_qrpp(N, N, A2, N, ipq, tau)
      case (13); call la_geqp3_host(N, N, A2, N, ipq, tau)
      end select
    else
      select case (k)
      case (1); call zgemm('N', 'N', N, N, N, one, Ah, N, Ch, N, zero, Ch, N)
      case (2); call zgemm('N', 'C', N, N, N, one, Ah, N, Ch, N, zero, Ch, N)
      case (3); call dgemm('N', 'N', 2*N, N, N, 1.d0, Ah, 2*N, Br, N, 0.d0, Ch, 2*N)
      case (4)
        do lw = 0, 3
          call zgemm('N', 'N', N/2, N/2, N/2, one, Ab(1+lw*N*N/2), N, Bb(1+lw*N*N/2), N, zero, Cb(1+lw*N*N/2), N)
        end do
      case (5); call zgeru(N, N, a1, x, 1, y, 1, Ch, N)
      case (6); call ztrsm('L', 'U', 'N', 'N', N, N, one, Ah, N, Ch, N)
      case (7); call zgeqrf(N, N, Ah, N, tau, hw, size(hw), info)
      case (8); call zungqr(N, N, N, Ah, N, tau, hw, size(hw), info)
      case (9); call zunmqr('L', 'C', N, N, N, Ah, N, tau, Ch, N, hw, size(hw), info)
      case (10); call zgetrf(N, N, Ah, N, ip, info)
      case (11); call zgetrs('N', N, N, Ah, N, ip, Ch, N, info)
      case (12); call zgeqrf(N, N, Ah, N, tau, hw, size(hw), info)
      case (13); call zgeqp3(N, N, Ah, N, ipq, tau, hw, size(hw), rw, info)
      end select
    end if
  end subroutine

  subroutine bench(k)
    integer, intent(in) :: k
    integer :: r, reps, info
    real(8) :: t0, tg, th, c0, c1, cg
    complex(8) :: wq(1)
    reps = max(1, min(50, int(2.d9 / (8.d0*dble(N)**3))))
    if (k == 8 .or. k == 7 .or. k == 9 .or. k == 10 .or. k == 11 .or. k == 12 .or. k == 13) reps = max(1, min(20, reps))
    if (allocated(hw)) deallocate(hw)
    allocate(hw(max(1, 64*N)))
    Ah = A2; Ch = C
    if (trace) then
      write(*,'(A,I0,A,A)') 'SECTION N=', N, ' ', trim(names(k))
      call one_call(k, .false.); call la_sync()
      write(*,'(A,I0,A,A)') 'SECTION_END N=', N, ' ', trim(names(k))
      return
    end if
    call one_call(k, .false.); call la_sync()                   ! warm-up: JIT, scratch growth
    call cpu_time(c0); t0 = now()
    do r = 1, reps
      call one_call(k, .false.)
      call la_sync()
    end do
    tg = (now() - t0) / reps; call cpu_time(c1); cg = (c1 - c0) / reps
    call one_call(k, .true.)
    call cpu_time(c0); t0 = now()
    do r = 1, reps
      call one_call(k, .true.)
    end do
    th = (now() - t0) / reps
    write(*,'(A,A30,I6,3(1X,ES12.4))') 'TABLE ', names(k), N, tg*1.d3, cg*1.d3, th*1.d3   ! ms per call: device wall, process CPU, host-pinned wall
  end subroutine
end program t302a_table

! timing copy of t_order (development only): wall time of each phase. Optional argument: NREPS of the RAW loop (default 24).
program tprof
  use alf_la_mod
  use la_test_util, only: fill_z, dev_enter_z, dev_exit_z
  implicit none
  integer, parameter :: N = 1024
  integer :: i, pass, rep, ns, nrep
  character(len=16) :: arg
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), C1(:,:), small(:,:), tau(:)
  real(8) :: acc(100, 2), t0, t1, tstart
  nrep = 24
  if (command_argument_count() > 0) then
    call get_command_argument(1, arg); read(arg,*) nrep
  end if
  tstart = wt()
  do pass = 1, 2
    t0 = wt(); call la_init(); call tp('la_init')
    allocate(A(N,N), B(N,N), C(N,N), C1(N,N))
    call fill_z(N, N, B)
    A = (0.d0, 0.d0); C = (0.d0, 0.d0); C1 = (0.d0, 0.d0); acc = 0.d0
    t0 = wt(); call dev_enter_z(A, N*N); call dev_enter_z(B, N*N); call dev_enter_z(C, N*N); call dev_enter_z(C1, N*N); call tp('enter data')
    t0 = wt()
    call set_diag(A, 1.d0); call tp('first set_diag (JIT)')
    t0 = wt()
    call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N); call la_sync(); call tp('first gemm+sync (JIT)')
    t0 = wt()
    call sum_diag(C, 1, 1); call tp('first sum_diag (JIT)')
    t0 = wt()
    call fill_const(C, 7.d0); call tp('first fill_const (JIT)')
    t0 = wt()
    call sum_all(C, 1); call tp('first sum_all (JIT)')
    t0 = wt()
    do rep = 1, nrep
      call set_diag(A, dble(rep))
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
      call sum_diag(C, rep, 1)
    end do
    call tp('RAW loop')
    t0 = wt()
    do rep = 1, nrep
      call set_diag(A, 1.d0)
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C1, N)
      call set_diag(A, 2.d0)
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
      call sum_diag(C1, rep, 1)
      call sum_diag(C, rep, 2)
    end do
    call tp('WAR loop')
    t0 = wt()
    call set_diag(A, 3.d0)
    do rep = 1, nrep
      call la_gemm(LA_N, LA_N, N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
      call fill_const(C, 7.d0)
      call sum_all(C, rep)
    end do
    call tp('WAW loop')
    t0 = wt(); call la_sync(); call tp('la_sync')
    t0 = wt()
    call dev_exit_z(A, N*N); call dev_exit_z(B, N*N); call dev_exit_z(C, N*N); call dev_exit_z(C1, N*N)
    deallocate(A, B, C, C1); call tp('exit data')
    ns = 2*LA_NSLOT + 100
    allocate(small(8,8), tau(8))
    call fill_z(8, 8, small); tau = (0.d0, 0.d0)
    call dev_enter_z(small, 64); call dev_enter_z(tau, 8)
    t0 = wt()
    call la_geqrf(8, 8, small, 8, tau); call la_sync(); call tp('first geqrf+sync')
    t0 = wt()
    do rep = 1, ns
      call la_geqrf(8, 8, small, 8, tau)
    end do
    call la_sync(); call tp('geqrf traffic')
    call dev_exit_z(small, 64); call dev_exit_z(tau, 8)
    deallocate(small, tau)
    t0 = wt(); call la_finalize(); call tp('la_finalize')
  end do
  t0 = tstart; call tp('TOTAL in main')
contains
  real(8) function wt()
    integer(8) :: c, r
    call system_clock(c, r); wt = dble(c)/dble(r)
  end function
  subroutine tp(s)
    character(*), intent(in) :: s
    write(*,'(A24,F9.3)') s, wt() - t0
  end subroutine
  subroutine set_diag(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
    !$omp target teams distribute parallel do collapse(2)
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
    !$omp target teams distribute parallel do collapse(2)
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
    !$omp target teams distribute parallel do reduction(+: s)
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
    !$omp target teams distribute parallel do collapse(2) reduction(+: s)
    do jj = 1, N
      do ii = 1, N
        s = s + dble(M(ii,jj))
      end do
    end do
    acc(r,1) = s
  end subroutine
end program tprof

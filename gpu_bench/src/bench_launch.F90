! Ticket 102: launch overhead of per-site work, OpenMP target (-DMODEL_OMP) vs OpenACC (-DMODEL_ACC).
! Standalone program; reuses bench_util / bench_ref / bench_la (cuda build flags) unchanged.
! All kernels act on device-resident arrays (mapped once with the model's own "enter data").
!
! Decision model of the single-site update (stand-in for ALF's accept test, upgrade_mod.F90:222):
!   r = |1 + G(i,i)|, accept  <=>  u(ks) < r/(1+r)   (u = pre-generated uniform number, accept rate ~ 50%)
!   on accept: G += x_ks y_ks^T (rank-1, N x N complex).
! Site index i = mod(ks-1, N) + 1 as in 101.

#define RATIO(g) sqrt((1d0+real(g,8))**2 + aimag(g)**2)

module launch_kern
  implicit none
contains

  ! ---- empty-ish kernel: one scalar update, a single thread ----
  subroutine k_empty(d)
    real(8) :: d(1)
#ifdef MODEL_ACC
    !$acc serial present(d)
    d(1) = d(1) + 1d0
    !$acc end serial
#else
    !$omp target
    d(1) = d(1) + 1d0
    !$omp end target
#endif
  end subroutine

  subroutine k_empty_a(d)
    real(8) :: d(1)
#ifdef MODEL_ACC
    !$acc serial present(d) async(1)
    d(1) = d(1) + 1d0
    !$acc end serial
#else
    !$omp target nowait depend(inout: d)
    d(1) = d(1) + 1d0
    !$omp end target
#endif
  end subroutine

#ifndef MODEL_ACC
  ! informational: OpenMP nowait without the depend clause (ordering then not guaranteed by the standard)
  subroutine k_empty_nd(d)
    real(8) :: d(1)
    !$omp target nowait
    d(1) = d(1) + 1d0
    !$omp end target
  end subroutine
#endif

  ! ---- rank-1 update G += x y^T (y has stride incy) ----
  subroutine k_ger(G, N, x, y, incy)
    integer, intent(in) :: N, incy
    complex(8) :: G(N,N), x(N), y(incy*(N-1)+1)
    integer :: k, l
#ifdef MODEL_ACC
    !$acc parallel loop collapse(2) present(G, x, y)
#else
    !$omp target teams distribute parallel do collapse(2)
#endif
    do l = 1, N
      do k = 1, N
        G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
  end subroutine

  subroutine k_ger_a(G, N, x, y, incy)
    integer, intent(in) :: N, incy
    complex(8) :: G(N,N), x(N), y(incy*(N-1)+1)
    integer :: k, l
#ifdef MODEL_ACC
    !$acc parallel loop collapse(2) present(G, x, y) async(1)
#else
    !$omp target teams distribute parallel do collapse(2) nowait depend(inout: G)
#endif
    do l = 1, N
      do k = 1, N
        G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
  end subroutine

  ! ---- Op_V wrap of one site: row i * ea, column i / ea ----
  subroutine k_opv(G, N, i, ea)
    integer, intent(in) :: N, i
    complex(8) :: G(N,N)
    real(8), intent(in) :: ea
    integer :: m
#ifdef MODEL_ACC
    !$acc parallel loop present(G)
#else
    !$omp target teams distribute parallel do
#endif
    do m = 1, N
      if (m /= i) then
        G(i,m) = G(i,m) * ea
        G(m,i) = G(m,i) / ea
      end if
    end do
  end subroutine

  subroutine k_opv_a(G, N, i, ea)
    integer, intent(in) :: N, i
    complex(8) :: G(N,N)
    real(8), intent(in) :: ea
    integer :: m
#ifdef MODEL_ACC
    !$acc parallel loop present(G) async(1)
#else
    !$omp target teams distribute parallel do nowait depend(inout: G)
#endif
    do m = 1, N
      if (m /= i) then
        G(i,m) = G(i,m) * ea
        G(m,i) = G(m,i) / ea
      end if
    end do
  end subroutine

  ! ---- device decision, kernel 1 of 2: one thread reads G(i,i), compares with u(ks), stores the flag ----
  subroutine k_dec(G, N, i, u, flag, ks)
    integer, intent(in) :: N, i, ks
    complex(8) :: G(N,N)
    real(8) :: u(*)
    integer :: flag(*)
    real(8) :: r
#ifdef MODEL_ACC
    !$acc serial present(G, u, flag) async(1)
#else
    !$omp target nowait depend(inout: G)
#endif
    r = RATIO(G(i,i))
    if (u(ks) < r/(1d0+r)) then
      flag(ks) = 1
    else
      flag(ks) = 0
    end if
#ifdef MODEL_ACC
    !$acc end serial
#else
    !$omp end target
#endif
  end subroutine

  ! ---- device decision, kernel 2 of 2: conditional rank-1 update ----
  subroutine k_ger_c(G, N, x, y, incy, flag, ks)
    integer, intent(in) :: N, incy, ks
    complex(8) :: G(N,N), x(N), y(incy*(N-1)+1)
    integer :: flag(*)
    integer :: k, l
#ifdef MODEL_ACC
    !$acc parallel loop collapse(2) present(G, x, y, flag) async(1)
#else
    !$omp target teams distribute parallel do collapse(2) nowait depend(inout: G)
#endif
    do l = 1, N
      do k = 1, N
        if (flag(ks) /= 0) G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
  end subroutine

  ! ---- device decision fused: ONE gang/team; decision, barrier, conditional rank-1 update ----
  subroutine k_fused(G, N, i, x, y, incy, u, flag, ks)
    integer, intent(in) :: N, i, incy, ks
    complex(8) :: G(N,N), x(N), y(incy*(N-1)+1)
    real(8) :: u(*)
    integer :: flag(*)
    integer :: k, l
    real(8) :: r
    logical :: take
#ifdef MODEL_ACC
    !$acc parallel num_gangs(1) vector_length(512) present(G, x, y, u, flag) async(1) private(r, take)
    r = RATIO(G(i,i))
    take = u(ks) < r/(1d0+r)       ! gang-redundant scalar code, completes before the vector loop starts
    if (take) then
      flag(ks) = 1
    else
      flag(ks) = 0
    end if
    !$acc loop vector collapse(2)
    do l = 1, N
      do k = 1, N
        if (take) G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
    !$acc end parallel
#else
    !$omp target teams num_teams(1) thread_limit(512) nowait depend(inout: G) private(r, take)
    !$omp parallel private(r, take)
    r = RATIO(G(i,i))
    take = u(ks) < r/(1d0+r)
    !$omp barrier
    !$omp single
    if (take) then
      flag(ks) = 1
    else
      flag(ks) = 0
    end if
    !$omp end single
    !$omp do collapse(2)
    do l = 1, N
      do k = 1, N
        if (take) G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
    !$omp end do
    !$omp end parallel
    !$omp end target teams
#endif
  end subroutine

end module launch_kern

module launch_sec
  use iso_c_binding
  use bench_util
  use bench_ref
  use launch_kern
  implicit none
  integer, parameter :: NU = 256, NE = 1000
  complex(8), parameter :: cone = (1d0, 0d0)
  type(c_ptr), save :: hb = c_null_ptr
  complex(8), allocatable, target, save :: cab(:)
  interface
    integer(c_int) function cublasCreate_v2(h) bind(C, name='cublasCreate_v2')
      import :: c_ptr, c_int
      type(c_ptr) :: h
    end function
    integer(c_int) function cublasZgeru_v2(h, m, n, alpha, x, incx, y, incy, A, lda) bind(C, name='cublasZgeru_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, x, y, A
      integer(c_int), value :: m, n, incx, incy, lda
    end function
  end interface
  real(8) :: sink = 0d0
contains

  ! cuBLAS ZGERU on device-resident data, device pointers from the model's own host_data / use_device_addr
  subroutine geru_dev(N, x, incx, y, incy, A)
    integer, intent(in) :: N, incx, incy
    complex(8), target :: x(*), y(*), A(N,*)
    integer :: st
    cab(1) = cone     ! cuBLAS wants a 16-byte aligned host alpha (see README gotchas)
#ifdef MODEL_ACC
    !$acc host_data use_device(x, y, A)
#else
    !$omp target data use_device_addr(x, y, A)
#endif
    st = cublasZgeru_v2(hb, N, N, c_loc(cab(1)), c_loc(x(1)), incx, c_loc(y(1)), incy, c_loc(A(1,1)), N)
#ifdef MODEL_ACC
    !$acc end host_data
#else
    !$omp end target data
#endif
    if (st /= 0) error stop 'cublasZgeru failed'
  end subroutine

  subroutine run_n(N)
    integer, intent(in) :: N
    complex(8), allocatable :: X(:,:), Y(:,:), G0(:,:), Gd(:,:), Gr(:,:), Gr1(:,:)
    real(8), allocatable :: rndh(:), dmy(:), ea(:), a(:)
    integer, allocatable :: flagH(:), flagR(:), flagD(:)
    integer :: ks, i, nmis
    real(8) :: us, r
    logical :: take
    character(len=*), parameter :: s1 = 'launch latency [us/launch]'
    character(len=*), parameter :: s2 = 'rank-1 update [us/update]'
    character(len=*), parameter :: s3 = 'single-site path [us/site]'
    character(len=*), parameter :: s4 = 'Op_V per-site launches [us/site]'

    allocate(X(N,NU), Y(NU,N), G0(N,N), Gd(N,N), Gr(N,N), Gr1(N,N), rndh(NU), dmy(1), ea(N), a(N))
    allocate(flagH(NU), flagR(NU), flagD(NU))
    call rand_c(G0, N*N, 1d0)
    call rand_c(X, N*NU, 0.1d0); call rand_c(Y, N*NU, 0.1d0)
    call random_number(rndh)
    call random_number(a); ea = exp(a - 0.5d0)
    dmy = 0d0; flagH = 0; flagD = 0; Gd = G0
#ifdef MODEL_ACC
    !$acc enter data copyin(X, Y, Gd, rndh, dmy, ea, flagD)
#else
    !$omp target enter data map(to: X, Y, Gd, rndh, dmy, ea, flagD)
#endif
    ! ---- host references ----
    Gr1 = G0                              ! all 256 updates accepted
    do ks = 1, NU
      call zgeru(N, N, cone, X(1,ks), 1, Y(ks,1), NU, Gr1, N)
    end do
    Gr = G0                               ! with decision
    do ks = 1, NU
      i = mod(ks-1, N) + 1
      r = RATIO(Gr(i,i)); take = rndh(ks) < r/(1d0+r)
      flagR(ks) = merge(1, 0, take)
      if (take) call zgeru(N, N, cone, X(1,ks), 1, Y(ks,1), NU, Gr, N)
    end do
    call report('accepted sites of 256 (reference decision)', N, dble(sum(flagR)))

    ! ---- correctness ----
    call verify_g('rank-1 x256 kernel sync vs host zgeru', r1_sync, Gr1)
    call verify_g('rank-1 x256 kernel async vs host zgeru', r1_async, Gr1)
    call verify_g('rank-1 x256 cublas zgeru vs host zgeru', r1_cublas, Gr1)
    flagH = -1
    call verify_g('host decision + rank-1 vs host', dec_host, Gr)
    nmis = count(flagH /= flagR)
    call check('host decision accept pattern mismatches', N, dble(nmis), 0.5d0)
    call verify_g('device decision 2 kernels vs host', dec_dev2, Gr)
    call get_flags(); nmis = count(flagD /= flagR)
    call check('device decision 2 kernels: accept pattern mismatches', N, dble(nmis), 0.5d0)
    call verify_g('device decision 1 gang vs host', dec_dev1, Gr)
    call get_flags(); nmis = count(flagD /= flagR)
    call check('device decision 1 gang: accept pattern mismatches', N, dble(nmis), 0.5d0)
    Gr = G0; call ref_opv(Gr, N, ea)
    call verify_g('Op_V per-site launches sync vs host', opv_sync, Gr)
    call verify_g('Op_V per-site launches async vs host', opv_async, Gr)

    ! ---- timing ----
    call timeit(emp_sync, us);  call add_rec(s1, 'empty sync', N, us/NE, 0d0, '')
    call timeit(emp_async, us); call add_rec(s1, 'empty async', N, us/NE, 0d0, '')
#ifndef MODEL_ACC
    call timeit(emp_nodep, us); call add_rec(s1, 'empty async nodep', N, us/NE, 0d0, '')
#endif
    call timeit(r1_sync, us);   call add_rec(s2, 'kernel sync', N, us/NU, 0d0, '')
    call timeit(r1_async, us);  call add_rec(s2, 'kernel async', N, us/NU, 0d0, '')
    call timeit(r1_cublas, us); call add_rec(s2, 'cublas zgeru', N, us/NU, 0d0, '')
    call timeit(dec_host, us);  call add_rec(s3, 'host decision', N, us/NU, 0d0, '')
    call timeit(dec_dev2, us);  call add_rec(s3, 'dev decision 2 kern', N, us/NU, 0d0, '')
    call timeit(dec_dev1, us);  call add_rec(s3, 'dev decision 1 gang', N, us/NU, 0d0, '')
    call timeit(opv_sync, us);  call add_rec(s4, 'sync', N, us/N, 0d0, '')
    call timeit(opv_async, us); call add_rec(s4, 'async', N, us/N, 0d0, '')
#ifdef MODEL_ACC
    !$acc exit data delete(X, Y, Gd, rndh, dmy, ea, flagD)
#else
    !$omp target exit data map(delete: X, Y, Gd, rndh, dmy, ea, flagD)
#endif
  contains
    subroutine restore()      ! verification only: put the start matrix back on the device
      Gd = G0
#ifdef MODEL_ACC
      !$acc update device(Gd)
#else
      !$omp target update to(Gd)
#endif
    end subroutine
    subroutine verify_g(name, op, Gref)
      character(len=*), intent(in) :: name
      procedure(op_i) :: op
      complex(8), intent(in) :: Gref(N,N)
      call restore()
      call op(); call la_sync()
#ifdef MODEL_ACC
      !$acc update host(Gd)      ! verification only
#else
      !$omp target update from(Gd)       ! verification only
#endif
      call check(name, N, relfro(Gd, Gref, N*N), 1d-12)
    end subroutine
    subroutine get_flags()
#ifdef MODEL_ACC
      !$acc update host(flagD)   ! verification only
#else
      !$omp target update from(flagD)    ! verification only
#endif
    end subroutine
    subroutine wait_all()
#ifdef MODEL_ACC
      !$acc wait(1)
#else
      !$omp taskwait
#endif
    end subroutine

    subroutine emp_sync()
      integer :: k
      do k = 1, NE
        call k_empty(dmy)
      end do
    end subroutine
    subroutine emp_async()
      integer :: k
      do k = 1, NE
        call k_empty_a(dmy)
      end do
      call wait_all()
    end subroutine
#ifndef MODEL_ACC
    subroutine emp_nodep()
      integer :: k
      do k = 1, NE
        call k_empty_nd(dmy)
      end do
      call wait_all()
    end subroutine
#endif
    subroutine r1_sync()
      integer :: k
      do k = 1, NU
        call k_ger(Gd, N, X(1,k), Y(k,1), NU)
      end do
    end subroutine
    subroutine r1_async()
      integer :: k
      do k = 1, NU
        call k_ger_a(Gd, N, X(1,k), Y(k,1), NU)
      end do
      call wait_all()
    end subroutine
    subroutine r1_cublas()
      integer :: k
      do k = 1, NU
        call geru_dev(N, X(1,k), 1, Y(k,1), NU, Gd)
      end do
    end subroutine
    subroutine dec_host()
      integer :: k, i
      complex(8) :: g
      real(8) :: r
      logical :: take
      do k = 1, NU
        i = mod(k-1, N) + 1
        ! unavoidable here: the host takes the decision, so G(i,i) must come to the host ...
#ifdef MODEL_ACC
        !$acc update host(Gd(i:i,i:i))
#else
        !$omp target update from(Gd(i:i,i:i))
#endif
        g = Gd(i,i)
        r = RATIO(g); take = rndh(k) < r/(1d0+r)
        flagH(k) = merge(1, 0, take)
        sink = sink + r
        if (take) call k_ger(Gd, N, X(1,k), Y(k,1), NU)   ! ... and the decision is "sent" as a launch (no data back)
      end do
    end subroutine
    subroutine dec_dev2()
      integer :: k, i
      do k = 1, NU
        i = mod(k-1, N) + 1
        call k_dec(Gd, N, i, rndh, flagD, k)
        call k_ger_c(Gd, N, X(1,k), Y(k,1), NU, flagD, k)
      end do
      call wait_all()
    end subroutine
    subroutine dec_dev1()
      integer :: k, i
      do k = 1, NU
        i = mod(k-1, N) + 1
        call k_fused(Gd, N, i, X(1,k), Y(k,1), NU, rndh, flagD, k)
      end do
      call wait_all()
    end subroutine
    subroutine opv_sync()
      integer :: i
      do i = 1, N
        call k_opv(Gd, N, i, ea(i))
      end do
    end subroutine
    subroutine opv_async()
      integer :: i
      do i = 1, N
        call k_opv_a(Gd, N, i, ea(i))
      end do
      call wait_all()
    end subroutine
  end subroutine

end module launch_sec

program bench_launch
  use launch_sec
#ifdef MODEL_ACC
  use openacc
#endif
  implicit none
  integer, parameter :: LL = 40
  integer :: Ns(16), nN, i, k, nc, seed(64), ns_seed, stat
  character(len=256) :: arg
  nN = 5; Ns(1:5) = [64, 144, 256, 400, 1024]
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    nc = 1
    do i = 1, len_trim(arg)
      if (arg(i:i) == ',') nc = nc + 1
    end do
    nN = nc
    read(arg,*) Ns(1:nN)
  end if
  call random_seed(size=ns_seed)
  seed = 0
  do i = 1, ns_seed
    seed(i) = 12345 + 7*i
  end do
  call random_seed(put=seed(1:ns_seed))
  call la_init()
  stat = cublasCreate_v2(hb)
  if (stat /= 0) error stop 'cublasCreate failed'
  allocate(cab(2))
#ifdef MODEL_ACC
  write(*,'(a,i0,a)') 'model: OpenACC (-acc=gpu), devices: ', acc_get_num_devices(acc_device_nvidia), ' (nvidia)'
#else
  write(*,'(a,i0,a)') 'model: OpenMP target (-mp=gpu), devices: ', omp_get_num_devices()
#endif
  do k = 1, nN
    write(*,'(/,a,i0)') '--- N = ', Ns(k)
    call run_n(Ns(k))
  end do
  write(*,'(/,a,i0,a,i0,a)') 'checks: ', nchk, ' run, ', nfail, ' FAILED'
  call print_tables(6)
  ! projection: single-site sweep = 2*L*N sites, L = 40, at the measured per-site cost (accept rate ~ 50%)
  write(*,'(/,a,i0,a)') '== projected s/sweep for a single-site sweep (2*L*N sites, L = ', LL, ')'
  write(*,'(a6,1x,a24,1x,a14)') 'N', 'variant', 's/sweep'
  do i = 1, nrec
    if (r_sec(i) == 'single-site path [us/site]' .or. r_var(i) == 'kernel async') then
      if (r_sec(i) == 'single-site path [us/site]' .or. r_sec(i) == 'rank-1 update [us/update]') &
        write(*,'(i6,1x,a24,1x,f14.4)') r_n(i), trim(r_var(i)), 2d0*LL*r_n(i)*r_us(i)*1d-6
    end if
  end do
  if (nfail > 0) stop 1
end program

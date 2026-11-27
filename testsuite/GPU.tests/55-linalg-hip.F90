!> @brief
!> GPU=amd (ALF_LA_HIP): the device path of linalg_mod against the host path, same wrappers, same inputs.
!> @details Mode "compare" (default): for la_zgemm (N,N / C,N / N,C), la_zgemm_cr, la_ztrsm (L,U,N,N / R,U,C,N / L,U,C,N),
!> la_zgeqrf, la_zunmqr (L,C / L,N / R,N), la_zqrpp, la_zgetrf + la_zgetri the result on device-resident (mapped) arrays is
!> compared with the host call on unmapped copies: max |dev - host| / max(1, max |host|) <= TOL; the pivots of la_zqrpp and
!> la_zgetrf must be equal; a NaN in the matrix must give info = -101 from the device la_zqrpp.  Run once with the own HIP QR
!> and once with ALF_HIP_QR=0 (hipSOLVER).  Mode "mixed": la_zgemm with A on the device and B, C on the host must stop the run
!> (ctest WILL_FAIL).  Identity input at n = 200 and 2400 (all tau = 0; rocSOLVER's blocked zgeqrf gives NaN there,
!> KNOWN_ISSUES R1): the device result must equal the host one.  Prints "linalg HIP test PASSED" or "... FAILED".
program linalg_hip_test
  use linalg_mod, only: la_zgemm, la_zgemm_cr, la_ztrsm, la_zgeqrf, la_zunmqr, la_zqrpp, la_zgetrf, la_zgetri
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  integer, parameter :: dp = kind(0.d0)
  real(dp), parameter :: TOL = 1.d-12
  integer, parameter :: n = 48, m = 64
  complex(dp), allocatable :: A(:,:), B(:,:), C(:,:), Ch(:,:), T(:,:), tau(:), tauh(:), wk(:), Th(:,:), Ah(:,:)
  real(dp), allocatable :: Br(:,:), re(:,:), im(:,:)
  integer, allocatable :: ip(:), iph(:)
  character(len=16) :: mode
  integer :: info, nfail, lw
  complex(dp) :: wq(1)

  mode = 'compare'
  if (command_argument_count() >= 1) call get_command_argument(1, mode)
  nfail = 0
  allocate(A(m,m), B(m,m), C(m,m), Ch(m,m), T(n,n), Th(n,n), Ah(n,n), tau(n), tauh(n), Br(m,m), re(m,m), im(m,m))
  allocate(ip(n), iph(n))
  call random_seed(put=[(17*lw, lw=1,64)])
  call random_number(re); call random_number(im)
  A = cmplx(re - 0.5d0, im - 0.5d0, dp)
  call random_number(re); call random_number(im)
  B = cmplx(re - 0.5d0, im - 0.5d0, dp)
  call random_number(Br)
  Br = Br - 0.5d0
  call random_number(re); call random_number(im)
  C = cmplx(re - 0.5d0, im - 0.5d0, dp)
  call la_zgeqrf(n, n, T, n, tau, wq, -1, info)       ! host LWORK query (never on the device)
  lw = max(int(dble(wq(1))), 4*m)
  call la_zunmqr('L', 'C', n, n, n, T, n, tau, Th, n, wq, -1, info)
  lw = max(lw, int(dble(wq(1))))
  allocate(wk(lw))

  if (trim(mode) == 'mixed') then
     !$omp target enter data map(to: A)
     call la_zgemm('N', 'N', m, m, m, (1.d0,0.d0), A, m, B, m, (0.d0,0.d0), C, m)
     nfail = 1                                         ! not stopped: the run ends normally, ctest WILL_FAIL sees a pass
  else
     call compare_all()
  endif
  if (nfail == 0) then
     write(*,'(A)') 'linalg HIP test PASSED'
  else
     write(*,'(A,I0,A)') 'linalg HIP test FAILED (', nfail, ' checks)'
  endif

contains

  subroutine compare_all()

  ! ---- GEMM ----
  call gemm_case('N', 'N'); call gemm_case('C', 'N'); call gemm_case('N', 'C')
  Ch = C
  call la_zgemm_cr(m, m, m, A, m, Br, m, Ch, m)
  block
    complex(dp) :: Cd(m,m)
    Cd = C
    !$omp target enter data map(to: A, Br, Cd)
    call la_zgemm_cr(m, m, m, A, m, Br, m, Cd, m)
    !$omp target exit data map(from: Cd) map(delete: A, Br)
    call check('la_zgemm_cr', Cd, Ch)
  end block

  ! ---- triangular solves (well-conditioned upper triangle: diagonal + 4) ----
  Th = A(1:n,1:n)
  do lw = 1, n
     Th(lw,lw) = Th(lw,lw) + 4.d0
  enddo
  call trsm_case('L', 'N'); call trsm_case('R', 'C'); call trsm_case('L', 'C')

  ! ---- QR and Q-apply ----
  Ah = A(1:n,1:n)
  T = Ah
  call la_zgeqrf(n, n, T, n, tauh, wk, size(wk), info)
  call need0('host zgeqrf', info)
  block
    complex(dp) :: Td(n,n), tdv(n)
    Td = Ah
    !$omp target enter data map(to: Td) map(alloc: tdv)
    call la_zgeqrf(n, n, Td, n, tdv, wk, size(wk), info)
    call need0('device la_zgeqrf', info)
    !$omp target exit data map(from: Td, tdv)
    call check('la_zgeqrf factors', Td, T)
    call check('la_zgeqrf tau', reshape(tdv, [n,1]), reshape(tauh, [n,1]))
  end block
  call unmqr_case('L', 'C'); call unmqr_case('L', 'N'); call unmqr_case('R', 'N')

  ! ---- pre-pivoted QR ----
  T = Ah
  call la_zqrpp(n, n, T, n, iph, tauh, wk, size(wk), info)
  call need0('host la_zqrpp', info)
  block
    complex(dp) :: Td(n,n), tdv(n)
    integer :: ipd(n)
    Td = Ah
    !$omp target enter data map(to: Td) map(alloc: tdv, ipd)
    call la_zqrpp(n, n, Td, n, ipd, tdv, wk, size(wk), info)
    call need0('device la_zqrpp', info)
    !$omp target exit data map(from: Td, tdv) map(delete: ipd)
    if (any(ipd /= iph)) call fail('la_zqrpp: pivots differ from the host')
    call check('la_zqrpp factors', Td, T)
    ! negative input: a NaN must give info = -101
    Td = Ah
    Td(3,5) = cmplx(ieee_value(1.d0, ieee_quiet_nan), 0.d0, dp)
    !$omp target enter data map(to: Td) map(alloc: tdv, ipd)
    call la_zqrpp(n, n, Td, n, ipd, tdv, wk, size(wk), info)
    !$omp target exit data map(delete: Td, tdv, ipd)
    if (info /= -101) call fail('la_zqrpp: NaN not reported (info /= -101)')
  end block

  ! ---- LU and inverse ----
  T = Th
  call la_zgetrf(n, n, T, n, iph, info)
  call need0('host zgetrf', info)
  block
    complex(dp) :: Td(n,n)
    integer :: ipd(n)
    Td = Th
    !$omp target enter data map(to: Td) map(alloc: ipd)
    call la_zgetrf(n, n, Td, n, ipd, info)
    call need0('device la_zgetrf', info)
    !$omp target update from(Td, ipd)
    if (any(ipd /= iph)) call fail('la_zgetrf: pivots differ from the host')
    call check('la_zgetrf factors', Td, T)
    call la_zgetri(n, T, n, iph, wk, size(wk), info)
    call need0('host zgetri', info)
    call la_zgetri(n, Td, n, ipd, wk, size(wk), info)
    call need0('device la_zgetri', info)
    !$omp target exit data map(from: Td) map(delete: ipd)
    call check('la_zgetri inverse', Td, T)
  end block
  ! identity input (every reflector has tau = 0): n = 200 and 2400 (> 2304: hipSOLVER only); rocSOLVER's blocked zgeqrf
  ! returns NaN there (KNOWN_ISSUES R1), linalg_hip_mod must give the exact host result (A = identity, tau = 0)
  call ident_case(200); call ident_case(2400)
  end subroutine compare_all

  subroutine ident_case(ni)
    integer, intent(in) :: ni
    complex(dp), allocatable :: Id(:,:), Ih(:,:), ti(:), th(:)
    integer :: i, inf
    allocate(Id(ni,ni), Ih(ni,ni), ti(ni), th(ni))
    Ih = (0.d0, 0.d0)
    do i = 1, ni
       Ih(i,i) = (1.d0, 0.d0)
    enddo
    Id = Ih
    call la_zgeqrf(ni, ni, Ih, ni, th, wk, size(wk), inf)
    call need0('host zgeqrf (identity)', inf)
    !$omp target enter data map(to: Id) map(alloc: ti)
    call la_zgeqrf(ni, ni, Id, ni, ti, wk, size(wk), inf)
    call need0('device la_zgeqrf (identity)', inf)
    !$omp target exit data map(from: Id, ti)
    if (any(Id /= Id)) call fail('la_zgeqrf (identity): NaN in the result')
    call check('la_zgeqrf identity factors', Id, Ih)
    call check('la_zgeqrf identity tau', reshape(ti, [ni,1]), reshape(th, [ni,1]))
  end subroutine ident_case

  subroutine gemm_case(ta, tb)
    character(len=1), intent(in) :: ta, tb
    complex(dp) :: Cd(m,m)
    Ch = C
    call la_zgemm(ta, tb, m, m, m, (0.7d0,-0.2d0), A, m, B, m, (0.3d0,0.1d0), Ch, m)
    Cd = C
    !$omp target enter data map(to: A, B, Cd)
    call la_zgemm(ta, tb, m, m, m, (0.7d0,-0.2d0), A, m, B, m, (0.3d0,0.1d0), Cd, m)
    !$omp target exit data map(from: Cd) map(delete: A, B)
    call check('la_zgemm '//ta//tb, Cd, Ch)
  end subroutine gemm_case

  subroutine trsm_case(side, op)
    character(len=1), intent(in) :: side, op
    complex(dp) :: Bh(n,n), Bd(n,n)
    Bh = B(1:n,1:n)
    call la_ztrsm(side, 'U', op, 'N', n, n, (1.d0,0.d0), Th, n, Bh, n)
    Bd = B(1:n,1:n)
    !$omp target enter data map(to: Th, Bd)
    call la_ztrsm(side, 'U', op, 'N', n, n, (1.d0,0.d0), Th, n, Bd, n)
    !$omp target exit data map(from: Bd) map(delete: Th)
    call check('la_ztrsm '//side//op, Bd, Bh)
  end subroutine trsm_case

  ! T, tauh hold the host QR of Ah (zgeqrf)
  subroutine unmqr_case(side, op)
    character(len=1), intent(in) :: side, op
    complex(dp) :: Bh(n,n), Bd(n,n)
    integer :: inf
    Bh = B(1:n,1:n)
    call la_zunmqr(side, op, n, n, n, T, n, tauh, Bh, n, wk, size(wk), inf)
    call need0('host zunmqr', inf)
    Bd = B(1:n,1:n)
    !$omp target enter data map(to: T, tauh, Bd)
    call la_zunmqr(side, op, n, n, n, T, n, tauh, Bd, n, wk, size(wk), inf)
    call need0('device la_zunmqr', inf)
    !$omp target exit data map(from: Bd) map(delete: T, tauh)
    call check('la_zunmqr '//side//op, Bd, Bh)
  end subroutine unmqr_case

  subroutine check(what, X, Y)
    character(len=*), intent(in) :: what
    complex(dp), intent(in) :: X(:,:), Y(:,:)
    real(dp) :: d
    d = maxval(abs(X - Y)) / max(1.d0, maxval(abs(Y)))
    write(*,'(A,A,ES10.3)') what, ': max rel diff device vs host ', d
    if (.not. (d <= TOL)) call fail(what)
  end subroutine check

  subroutine need0(what, inf)
    character(len=*), intent(in) :: what
    integer, intent(in) :: inf
    if (inf /= 0) call fail(what//': info /= 0')
  end subroutine need0

  subroutine fail(what)
    character(len=*), intent(in) :: what
    write(*,'(2A)') 'FAIL: ', what
    nfail = nfail + 1
  end subroutine fail

end program linalg_hip_test

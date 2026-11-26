! Unit test for linalg_mod: every la_<name> must give bitwise the same result
! as the direct call of <name>.
program test_linalg
  use linalg_mod
  implicit none
  integer, parameter :: dp = kind(0.d0), n = 7, m = 5
  complex(dp), parameter :: alpha = (0.7d0, -0.3d0), beta = (0.2d0, 0.5d0)
  complex(dp) :: A(n,n), B(n,n), C1(n,n), C2(n,n), H(n,n), x(n), y1(n), y2(n)
  complex(dp) :: tau1(n), tau2(n), work1(4*n*n), work2(4*n*n)
  real(dp)    :: rwork1(5*n), rwork2(5*n), s1(n), s2(n)
  complex(dp) :: U1(n,n), U2(n,n), V1(n,n), V2(n,n), w1(n), w2(n)
  integer     :: ipiv1(n), ipiv2(n), jpvt1(n), jpvt2(n), info1, info2, n_fail
  real(dp)    :: Br(n,m), Cr2(2*n,m)
  complex(dp) :: Cc(n,m), Cd(n,m)
  real(dp), external :: dznrm2
  external :: dgemm, zgemm, zhemm, ztrsm, ztrmm, zgemv, zgeru, zscal, zcopy, zlacpy, &
       &      zgetrf, zgetrs, zgetri, zgeqrf, zgeqp3, zungqr, zunmqr, &
       &      zlaset, zlapmr, zlapmt, zgesvd, zheev, zgeev

  n_fail = 0
  call fill(A); call fill(B); call fill(H); H = H + conjg(transpose(H)); call fillv(x)

  C1 = B; C2 = B
  call zgemm('N', 'C', n, m, n, alpha, A, n, H, n, beta, C1, n)
  call la_zgemm('N', 'C', n, m, n, alpha, A, n, H, n, beta, C2, n)
  call check('zgemm', all(C1 == C2))

  ! la_zgemm_cr: bitwise the DGEMM on the interleaved real view; 1e-14 vs la_zgemm with cmplx(B)
  call random_number(Br)
  Cc = (0.d0, 0.d0); Cd = (0.d0, 0.d0)
  call la_zgemm_cr(n, m, n, A, n, Br, n, Cc, n)
  call dgemm('N', 'N', 2*n, m, n, 1.d0, transfer(A, 0.d0, 2*n*n), 2*n, Br, n, 0.d0, Cr2, 2*n)
  call check('zgemm_cr bitwise', all(Cc == cmplx(Cr2(1::2,:), Cr2(2::2,:), kind=dp)))
  call la_zgemm('N', 'N', n, m, n, (1.d0, 0.d0), A, n, cmplx(Br, 0.d0, kind=dp), n, (0.d0, 0.d0), Cd, n)
  call check('zgemm_cr vs zgemm', maxval(abs(Cc - Cd)) < 1.d-14)

  ! la_dznrm2: bitwise the direct dznrm2 (a column and a strided row)
  call check('dznrm2', la_dznrm2(n, A(1,3), 1) == dznrm2(n, A(1,3), 1) .and. la_dznrm2(m, A(2,1), n) == dznrm2(m, A(2,1), n))

  call test_zqrpp()
#if defined(ALF_LA_ONEMKL)
  call test_device()
#endif

  C1 = B; C2 = B
  call zhemm('L', 'U', n, n, alpha, H, n, A, n, beta, C1, n)
  call la_zhemm('L', 'U', n, n, alpha, H, n, A, n, beta, C2, n)
  call check('zhemm', all(C1 == C2))

  C1 = B; C2 = B
  call ztrsm('L', 'U', 'N', 'N', n, n, alpha, A + 5.d0, n, C1, n)
  call la_ztrsm('L', 'U', 'N', 'N', n, n, alpha, A + 5.d0, n, C2, n)
  call check('ztrsm', all(C1 == C2))

  C1 = B; C2 = B
  call ztrmm('R', 'U', 'C', 'N', n, n, alpha, A, n, C1, n)
  call la_ztrmm('R', 'U', 'C', 'N', n, n, alpha, A, n, C2, n)
  call check('ztrmm', all(C1 == C2))

  y1 = x(n:1:-1); y2 = y1
  call zgemv('T', n, m, alpha, A, n, x, 1, beta, y1, 1)
  call la_zgemv('T', n, m, alpha, A, n, x, 1, beta, y2, 1)
  call check('zgemv', all(y1 == y2))

  C1 = B; C2 = B
  call zgeru(n, n, alpha, x, 1, y1, 1, C1, n)
  call la_zgeru(n, n, alpha, x, 1, y1, 1, C2, n)
  call check('zgeru', all(C1 == C2))

  C1 = B; C2 = B          ! a row, stride n
  call zscal(n, alpha, C1(2,1), n)
  call la_zscal(n, alpha, C2(2,1), n)
  call check('zscal', all(C1 == C2))

  y1 = 0; y2 = 0
  call zcopy(n, A(1,3), 1, y1, 1)
  call la_zcopy(n, A(1,3), 1, y2, 1)
  call check('zcopy', all(y1 == y2))

  C1 = 0; C2 = 0
  call zlacpy('U', n, n, A, n, C1, n)
  call la_zlacpy('U', n, n, A, n, C2, n)
  call check('zlacpy', all(C1 == C2))

  C1 = A; C2 = A
  call zgetrf(n, n, C1, n, ipiv1, info1)
  call la_zgetrf(n, n, C2, n, ipiv2, info2)
  call check('zgetrf', all(C1 == C2) .and. all(ipiv1 == ipiv2) .and. info1 == info2)

  y1 = B(:,1); y2 = y1    ! C1 = C2 = LU factors from here on
  call zgetrs('C', n, 1, C1, n, ipiv1, y1, n, info1)
  call la_zgetrs('C', n, 1, C2, n, ipiv2, y2, n, info2)
  call check('zgetrs', all(y1 == y2) .and. info1 == info2)

  call zgetri(n, C1, n, ipiv1, work1, size(work1), info1)
  call la_zgetri(n, C2, n, ipiv2, work2, size(work2), info2)
  call check('zgetri', all(C1 == C2) .and. info1 == info2)

  C1 = A; C2 = A
  call zgeqrf(n, n, C1, n, tau1, work1, size(work1), info1)
  call la_zgeqrf(n, n, C2, n, tau2, work2, size(work2), info2)
  call check('zgeqrf', all(C1 == C2) .and. all(tau1 == tau2) .and. info1 == info2)

  H = C1                  ! QR factors (and tau1) of zgeqrf, used by zunmqr and zungqr
  call fill(B)
  C1 = B; C2 = B
  call zunmqr('R', 'N', n, n, n, H, n, tau1, C1, n, work1, size(work1), info1)
  call la_zunmqr('R', 'N', n, n, n, H, n, tau1, C2, n, work2, size(work2), info2)
  call check('zunmqr', all(C1 == C2) .and. info1 == info2)

  C1 = H; C2 = H
  call zungqr(n, n, n, C1, n, tau1, work1, size(work1), info1)
  call la_zungqr(n, n, n, C2, n, tau1, work2, size(work2), info2)
  call check('zungqr', all(C1 == C2) .and. info1 == info2)

  C1 = A; C2 = A; jpvt1 = 0; jpvt2 = 0
  call zgeqp3(n, n, C1, n, jpvt1, tau1, work1, size(work1), rwork1, info1)
  call la_zgeqp3(n, n, C2, n, jpvt2, tau2, work2, size(work2), rwork2, info2)
  call check('zgeqp3', all(C1 == C2) .and. all(tau1 == tau2) .and. all(jpvt1 == jpvt2) .and. info1 == info2)

  C1 = B; C2 = B
  call zlaset('L', n, m, alpha, beta, C1, n)
  call la_zlaset('L', n, m, alpha, beta, C2, n)
  call check('zlaset', all(C1 == C2))

  ipiv1 = [(n + 1 - info1, info1 = 1, n)]; ipiv2 = ipiv1      ! a permutation
  C1 = A; C2 = A
  call zlapmr(.true., n, m, C1, n, ipiv1)
  call la_zlapmr(.true., n, m, C2, n, ipiv2)
  call check('zlapmr', all(C1 == C2) .and. all(ipiv1 == ipiv2))

  C1 = A; C2 = A
  call zlapmt(.false., m, n, C1, n, ipiv1)
  call la_zlapmt(.false., m, n, C2, n, ipiv2)
  call check('zlapmt', all(C1 == C2) .and. all(ipiv1 == ipiv2))

  C1 = A; C2 = A
  call zgesvd('A', 'A', n, n, C1, n, s1, U1, n, V1, n, work1, size(work1), rwork1, info1)
  call la_zgesvd('A', 'A', n, n, C2, n, s2, U2, n, V2, n, work2, size(work2), rwork2, info2)
  call check('zgesvd', all(s1 == s2) .and. all(U1 == U2) .and. all(V1 == V2) .and. info1 == info2)

  C1 = H; C2 = H
  call zheev('V', 'U', n, C1, n, s1, work1, size(work1), rwork1, info1)
  call la_zheev('V', 'U', n, C2, n, s2, work2, size(work2), rwork2, info2)
  call check('zheev', all(C1 == C2) .and. all(s1 == s2) .and. info1 == info2)

  C1 = A; C2 = A
  call zgeev('V', 'V', n, C1, n, w1, U1, n, V1, n, work1, size(work1), rwork1, info1)
  call la_zgeev('V', 'V', n, C2, n, w2, U2, n, V2, n, work2, size(work2), rwork2, info2)
  call check('zgeev', all(w1 == w2) .and. all(U1 == U2) .and. all(V1 == V2) .and. info1 == info2)

  if (n_fail > 0) then
     write(*,*) 'linalg_mod test: ', n_fail, ' failure(s)'
     error stop 1
  endif
  write(*,*) 'linalg_mod test: all wrappers bitwise equal to the direct calls'

contains

  !> la_zqrpp: bitwise against a hand-written reference (norms, stable selection by descending norm with the lowest index
  !> first among equal norms, ZGEQRF); permutation sorted by norm; A(:,ipvt) = Q R; NaN/Inf -> info = -101 and A unchanged;
  !> equal columns keep their order; workspace query.
  subroutine test_zqrpp()
    complex(dp) :: A0(n,n), A1(n,n), A2(n,n), Q(n,n), R(n,n), t1(n), t2(n), wq(1)
    real(dp) :: nr(n), best
    integer :: ip(n), ir(n), j, k, kb, info3
    logical :: used(n), sorted, ties_ok
    call fill(A0)
    A1 = A0; A2 = A0
    call la_zqrpp(n, n, A1, n, ip, t1, work1, size(work1), info1)
    ! reference
    do j = 1, n
       nr(j) = dznrm2(n, A0(1,j), 1)
    enddo
    used = .false.
    do j = 1, n
       kb = 0; best = -1.d0
       do k = 1, n
          if (.not. used(k) .and. nr(k) > best) then
             best = nr(k); kb = k
          endif
       enddo
       ir(j) = kb; used(kb) = .true.
       A2(:,j) = A0(:,kb)
    enddo
    call zgeqrf(n, n, A2, n, t2, work2, size(work2), info2)
    call check('zqrpp bitwise', all(A1 == A2) .and. all(t1 == t2) .and. all(ip == ir) .and. info1 == 0 .and. info2 == 0)
    sorted = .true.
    do j = 2, n
       if (nr(ip(j)) > nr(ip(j-1))) sorted = .false.
    enddo
    call check('zqrpp sorted by norm', sorted)
    ! A0(:,ip) = Q R
    R = (0.d0, 0.d0)
    do j = 1, n
       R(1:j,j) = A1(1:j,j)
    enddo
    Q = A1
    call zungqr(n, n, n, Q, n, t1, work1, size(work1), info1)
    call check('zqrpp A(:,ipvt) = QR', maxval(abs(matmul(Q, R) - A0(:,ip))) < 1.d-13)
    ! ties: columns 2 and 4 equal, column 1 and 3 smaller: equal norms keep the input order
    A1 = A0; A1(:,4) = A1(:,2); A1(:,1) = 0.1d0*A1(:,1); A1(:,3) = 0.1d0*A1(:,3)
    call la_zqrpp(n, n, A1, n, ip, t1, work1, size(work1), info1)
    ties_ok = .false.
    do j = 1, n-1
       if (ip(j) == 2 .and. ip(j+1) == 4) ties_ok = .true.
    enddo
    call check('zqrpp ties keep order', ties_ok)
    ! NaN and Inf: info = -101, matrix not changed
    A1 = A0; A1(3,2) = transfer([9218868437227405312_8, 0_8], A1(3,2))      ! +Inf
    A2 = A1
    call la_zqrpp(n, n, A1, n, ip, t1, work1, size(work1), info1)
    call check('zqrpp Inf', info1 == -101 .and. all(A1 == A2))
    A1 = A0; A1(2,5) = transfer([9221120237041090560_8, 0_8], A1(2,5))      ! NaN
    A2 = A1
    call la_zqrpp(n, n, A1, n, ip, t1, work1, size(work1), info1)
    call check('zqrpp NaN', info1 == -101 .and. all(transfer(A1, 0_8, 2*n*n) == transfer(A2, 0_8, 2*n*n)))
    ! workspace query
    call zgeqrf(n, n, A0, n, t1, wq, -1, info3)
    call la_zqrpp(n, n, A0, n, ip, t2, work2, -1, info2)
    call check('zqrpp query', work2(1) == wq(1) .and. info2 == 0 .and. info3 == 0)
  end subroutine test_zqrpp

#if defined(ALF_LA_ONEMKL)
  !> GPU=intel: the device path (presence dispatch) of the wrappers used by the gpu_* modules, on mapped copies, against
  !> the host path of the same binary (rounding of oneMKL on the device may differ: 1e-12 relative); pivots and
  !> permutations must be equal. Skipped without an offload device.
  subroutine test_device()
    use omp_lib, only: omp_get_num_devices
    complex(dp) :: Ad(n,n), Bd(n,n), Cdv(n,n), td(n), Ah(n,n), Ch(n,n), th(n), Tr(n,n), Ccv(n,m), Cch(n,m)
    real(dp) :: Brd(n,m)
    integer :: ipd(n), iph(n), infod, infoh, i
    if (omp_get_num_devices() < 1) then
       write(*,*) 'linalg_mod test: device cases skipped (no offload device)'
       return
    endif
    Ad = A; Bd = H; Cdv = B; Ch = B; Brd = Br; Ccv = (0.d0, 0.d0)
    Tr = A
    do i = 1, n
       Tr(i,i) = Tr(i,i) + 4.d0
    enddo
    !$omp target enter data map(to: Ad, Bd, Cdv, Brd, Ccv, Tr) map(alloc: td, ipd)
    ! zgemm
    call la_zgemm('N', 'C', n, n, n, alpha, Ad, n, Bd, n, beta, Cdv, n)
    !$omp target update from(Cdv)
    call la_zgemm('N', 'C', n, n, n, alpha, A, n, H, n, beta, Ch, n)
    call check('device zgemm', rel(Cdv, Ch) < 1.d-12)
    ! zgemm_cr
    call la_zgemm_cr(n, m, n, Ad, n, Brd, n, Ccv, n)
    !$omp target update from(Ccv)
    call la_zgemm_cr(n, m, n, A, n, Br, n, Cch, n)
    call check('device zgemm_cr', rel(Ccv, Cch) < 1.d-12)
    ! ztrsm
    Ch = B; Cdv = B
    !$omp target update to(Cdv)
    call la_ztrsm('L', 'U', 'N', 'N', n, n, alpha, Tr, n, Cdv, n)
    !$omp target update from(Cdv)
    Ah = A
    do i = 1, n
       Ah(i,i) = Ah(i,i) + 4.d0
    enddo
    call la_ztrsm('L', 'U', 'N', 'N', n, n, alpha, Ah, n, Ch, n)
    call check('device ztrsm', rel(Cdv, Ch) < 1.d-12)
    ! zgeqrf, then zunmqr with its reflectors
    call la_zgeqrf(n, n, Ad, n, td, work1, size(work1), infod)
    !$omp target update from(Ad, td)
    Ah = A
    call la_zgeqrf(n, n, Ah, n, th, work2, size(work2), infoh)
    call check('device zgeqrf', rel(Ad, Ah) < 1.d-12 .and. rel1(td, th) < 1.d-12 .and. infod == 0 .and. infoh == 0)
    Cdv = B; Ch = B
    !$omp target update to(Cdv)
    call la_zunmqr('L', 'N', n, n, n, Ad, n, td, Cdv, n, work1, size(work1), infod)
    !$omp target update from(Cdv)
    call la_zunmqr('L', 'N', n, n, n, Ah, n, th, Ch, n, work2, size(work2), infoh)
    call check('device zunmqr', rel(Cdv, Ch) < 1.d-12 .and. infod == 0 .and. infoh == 0)
    ! zgetrf, zgetri
    Ad = A; Ah = A
    !$omp target update to(Ad)
    call la_zgetrf(n, n, Ad, n, ipd, infod)
    call la_zgetri(n, Ad, n, ipd, work1, size(work1), infod)
    !$omp target update from(Ad, ipd)
    call la_zgetrf(n, n, Ah, n, iph, infoh)
    call la_zgetri(n, Ah, n, iph, work2, size(work2), infoh)
    call check('device zgetrf/zgetri', all(ipd == iph) .and. rel(Ad, Ah) < 1.d-12 .and. infod == 0 .and. infoh == 0)
    ! zqrpp (ipvt valid on the host on return)
    Ad = A; Ah = A
    !$omp target update to(Ad)
    call la_zqrpp(n, n, Ad, n, ipd, td, work1, size(work1), infod)
    !$omp target update from(Ad, td)
    call la_zqrpp(n, n, Ah, n, iph, th, work2, size(work2), infoh)
    call check('device zqrpp', all(ipd == iph) .and. rel(Ad, Ah) < 1.d-12 .and. rel1(td, th) < 1.d-12 &
         & .and. infod == 0 .and. infoh == 0)
    !$omp target exit data map(delete: Ad, Bd, Cdv, Brd, Ccv, Tr, td, ipd)
  end subroutine test_device

  real(dp) function rel(X, Y)
    complex(dp), intent(in) :: X(:,:), Y(:,:)
    rel = maxval(abs(X - Y)) / maxval(abs(Y))
  end function rel

  real(dp) function rel1(x, y)
    complex(dp), intent(in) :: x(:), y(:)
    rel1 = maxval(abs(x - y)) / maxval(abs(y))
  end function rel1
#endif

  subroutine check(name, ok)
    character(len=*), intent(in) :: name
    logical, intent(in) :: ok
    if (.not. ok) then
       write(*,*) 'FAIL: ', name
       n_fail = n_fail + 1
    endif
  end subroutine check

  subroutine fill(M)
    complex(dp), intent(out) :: M(:,:)
    real(dp) :: re(size(M,1),size(M,2)), im(size(M,1),size(M,2))
    call random_number(re); call random_number(im)
    M = cmplx(re - 0.5d0, im - 0.5d0, kind=dp)
  end subroutine fill

  subroutine fillv(v)
    complex(dp), intent(out) :: v(:)
    real(dp) :: re(size(v)), im(size(v))
    call random_number(re); call random_number(im)
    v = cmplx(re - 0.5d0, im - 0.5d0, kind=dp)
  end subroutine fillv

end program test_linalg

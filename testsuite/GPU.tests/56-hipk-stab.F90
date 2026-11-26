!> @brief
!> GPU=amd (ALF_LA_HIP): kernel-level test of the HIP replacements of the gpu_stab_mod kernels.
!> @details Every kernel runs twice on the same device data, once as the OpenMP kernel of gpu_stab_mod (hk_set_enabled(.false.))
!> and once as its HIP replacement; results must agree to TOL relative (the HIP phase kernels multiply in a wave-parallel order).
!> Kernels: k_phase (cs(1), cs(2) and the S1 det cs(6), sgn +-1, cj T/F, tau from the device QR including tau = 0 entries),
!> k_det_phase (cs(4) normalised, cs(5) raw, pivot swaps, a zero pivot), k_cgr_phase (nvar 1/2, dev 0/1, bad 1), k_scale_cols,
!> k_rowdiv, k_colinv, k_tpup (default build), k_perm_rows, k_iperm_cols, k_triu, k_testscale (an Inf entry).  Orders 37 and
!> 200 (not multiples of the wave size 64).  Negative control: a perturbed tau(2) must change cs(6) beyond TOL.
!> Prints "HIP stab kernel test PASSED" or "... FAILED".
program hipk_stab_test
  use gpu_stab_mod, only: k_phase, k_det_phase, k_cgr_phase, k_scale_cols, k_rowdiv, k_colinv, k_perm_rows, k_iperm_cols, &
       &                  k_triu, k_testscale
#if !defined(STABLOG)
  use gpu_stab_mod, only: k_tpup
#endif
  use alf_hipk_mod, only: hk_set_enabled
  use linalg_mod, only: la_zgeqrf
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  implicit none
  integer, parameter :: dp = kind(0.d0)
  real(dp), parameter :: TOL = 1.d-12
  integer :: nfail, k, n
  integer, parameter :: ns(2) = [37, 200]
  complex(dp), allocatable :: A(:,:), A0(:,:), B(:,:), R(:,:), X1(:,:), X2(:,:), tau(:), D(:), D2(:)
  complex(dp) :: cs(8)
  integer, allocatable :: ipiv(:), perm(:)

  nfail = 0
  do k = 1, 2
     call run_order(ns(k))
  enddo
  if (nfail == 0) then
     write(*,'(A)') 'HIP stab kernel test PASSED'
  else
     write(*,'(A,I0,A)') 'HIP stab kernel test FAILED (', nfail, ' checks)'
  endif

contains

  subroutine rnd(X)
    complex(dp), intent(out) :: X(:,:)
    real(dp) :: re(size(X,1),size(X,2)), im(size(X,1),size(X,2))
    call random_number(re); call random_number(im)
    X = cmplx(re - 0.5d0, im - 0.5d0, dp)
  end subroutine rnd

  subroutine run_order(nn)
    integer, intent(in) :: nn
    complex(dp) :: c1(8), c2(8), c0(8), wq(1)
    complex(dp), allocatable :: wk(:)
    real(dp) :: rs(16), rs1(16)
    real(dp), allocatable :: dre(:)
    integer :: i, j, kk, info, sgn, nvar, dev, cjl
    logical :: cj
    character(len=40) :: lab

    n = nn
    if (allocated(A)) deallocate(A, A0, B, R, X1, X2, tau, D, D2, ipiv, perm)
    allocate(A(n,n), A0(n,n), B(n,n), R(n,n), X1(n,n), X2(n,n), tau(n), D(n), D2(n), ipiv(n), perm(n), dre(n))
    call rnd(A0); call rnd(B)
    do i = 1, n
       A0(i,i) = A0(i,i) + 3.d0
    enddo
    call random_number(dre)
    D = cmplx(dre + 0.5d0, 0.d0, dp)                         ! real positive scales (as the D of a UDV state)
    call rnd(R)
    D2(:) = D(:) * cmplx(1.d0, 0.3d0, dp)                    ! complex divisors for k_rowdiv
    do i = 1, n
       perm(i) = i
    enddo
    do i = n, 2, -1                                          ! deterministic shuffle
       j = 1 + mod(7*i + 3, i)
       call swap_i(perm(i), perm(j))
    enddo
    ipiv = [(i, i=1,n)]
    ipiv(2) = 5; ipiv(7) = 9; ipiv(n) = n                     ! three pivot swaps (one is the identity)
    !$omp target enter data map(alloc: A, B, R, X1, tau, D, D2, cs, rs, ipiv, perm)

    ! ---- tau and R from the device QR (own HIP QR / hipSOLVER), with a tau = 0 column ----
    A = A0
    A(2:n,1) = (0.d0, 0.d0); A(1,1) = (2.d0, 0.d0)            ! real alpha and zero x: zlarfg gives tau(1) = 0
    !$omp target update to(A)
    call la_zgeqrf(n, n, A, n, tau, wq, -1, info)
    allocate(wk(max(1,int(dble(wq(1))))))
    call la_zgeqrf(n, n, A, n, tau, wk, size(wk), info)
    !$omp target update from(tau)
    if (info /= 0 .or. tau(1) /= (0.d0,0.d0)) call fail('setup: device QR (info or tau(1) = 0 case)')
    !$omp target update to(D, D2, R, ipiv, perm, B)

    ! ---- k_phase ----
    do sgn = -1, 1, 2
       do cjl = 0, 1
          cj = (cjl == 1)
          write(lab,'(A,I0,A,I2,A,L1)') 'k_phase n=', n, ' sgn=', sgn, ' cj=', cj
          call k_on(.false.); call k_phase(n, A, n, sgn, cj, cs, tau); call get_cs(c1)
          call k_on(.true.);  call k_phase(n, A, n, sgn, cj, cs, tau); call get_cs(c2)
          call cmp_c(trim(lab)//' cs(1,2,6)', [c1(1), c1(2), c1(6)], [c2(1), c2(2), c2(6)])
       enddo
    enddo
    ! negative control: tau(2) perturbed in the HIP run must change cs(6)
    call k_phase(n, A, n, 1, .false., cs, tau); call get_cs(c0)
    tau(2) = tau(2) * cmplx(1.d0, 1.d-3, dp)
    !$omp target update to(tau)
    call k_phase(n, A, n, 1, .false., cs, tau); call get_cs(c2)
    if (abs(c2(6) - c0(6)) <= TOL) call fail('k_phase negative control (perturbed tau not seen in cs(6))')
    tau(2) = tau(2) / cmplx(1.d0, 1.d-3, dp)
    !$omp target update to(tau)

    ! ---- k_det_phase (LU-like diagonal of A0, pivot swaps; then a zero pivot) ----
    do kk = 1, 2
       A = A0
       if (kk == 2) A(3,3) = (0.d0, 0.d0)
       !$omp target update to(A)
       call k_on(.false.); call k_det_phase(n, A, n, ipiv, cs); call get_cs(c1)
       call k_on(.true.);  call k_det_phase(n, A, n, ipiv, cs); call get_cs(c2)
       write(lab,'(A,I0,A,I0)') 'k_det_phase n=', n, ' case ', kk
       call cmp_c(trim(lab)//' cs(4,5)', [c1(4), c1(5)], [c2(4), c2(5)])
       if (kk == 2 .and. abs(c2(4)) /= 0.d0) call fail('k_det_phase: zero pivot must give cs(4) = 0')
    enddo

    ! ---- k_cgr_phase (T = R factor of the QR, tau of the QR; tracked phase in cs(7)) ----
    A = A0
    !$omp target update to(A)
    call la_zgeqrf(n, n, A, n, tau, wk, size(wk), info)
    do nvar = 1, 2
       do dev = 0, 1
          cs = (0.d0, 0.d0); cs(7) = exp(cmplx(0.d0, 0.7d0, dp))
          !$omp target update to(cs)
          call k_on(.false.); call k_cgr_phase(n, A, n, tau, cs, nvar, -1, 0, exp(cmplx(0.d0, -0.4d0, dp)), dev); call get_cs(c1)
          call k_on(.true.);  call k_cgr_phase(n, A, n, tau, cs, nvar, -1, 0, exp(cmplx(0.d0, -0.4d0, dp)), dev); call get_cs(c2)
          write(lab,'(A,I0,A,I0,A,I0)') 'k_cgr_phase n=', n, ' nvar=', nvar, ' dev=', dev
          call cmp_c(trim(lab), [c1(3)], [c2(3)])
       enddo
    enddo
    call k_on(.false.); call k_cgr_phase(n, A, n, tau, cs, 1, 1, 1, (1.d0,0.d0), 0); call get_cs(c1)
    call k_on(.true.);  call k_cgr_phase(n, A, n, tau, cs, 1, 1, 1, (1.d0,0.d0), 0); call get_cs(c2)
    call cmp_c('k_cgr_phase bad=1', [c1(3)], [c2(3)])

    ! ---- element-wise kernels ----
    call mat2('k_scale_cols', 1); call mat2('k_rowdiv', 2); call mat2('k_colinv', 3); call mat2('k_triu', 4)
    call mat2('k_perm_rows', 5); call mat2('k_iperm_cols', 6)
#if !defined(STABLOG)
    call mat2('k_tpup', 7)
#endif

    ! ---- k_testscale (one Inf scale) ----
    D(n/2) = cmplx(ieee_value(1.d0, ieee_positive_inf), 0.d0, dp)
    !$omp target update to(D)
    call k_on(.false.)
    call k_testscale(n, D, rs)
    !$omp target update from(rs(1:3))
    rs1 = rs
    call k_on(.true.)
    call k_testscale(n, D, rs)
    !$omp target update from(rs(1:3))
    if (any(rs(1:3) /= rs1(1:3)) .or. nint(rs(1)) /= 1) call fail('k_testscale (Inf count, D(1), D(n))')
    write(*,'(A,I0,A,3ES12.4)') 'k_testscale n=', n, ': ', rs(1:3)

    !$omp target exit data map(delete: A, B, R, X1, tau, D, D2, cs, rs, ipiv, perm)
    call k_on(.true.)
  end subroutine run_order

    subroutine k_on(on)
      logical, intent(in) :: on
      call hk_set_enabled(on)
    end subroutine k_on

    subroutine get_cs(c)
      complex(dp), intent(out) :: c(8)
      !$omp target update from(cs)
      c = cs
    end subroutine get_cs

    ! one element-wise kernel, OpenMP then HIP, on X1 (device) started from B
    subroutine mat2(nm, kind)
      character(len=*), intent(in) :: nm
      integer, intent(in) :: kind
      integer :: pass
      do pass = 1, 2
         call k_on(pass == 2)
         X1 = B
         !$omp target update to(X1)
         select case (kind)
         case (1); call k_scale_cols(n, X1, n, D)
         case (2); call k_rowdiv(n, X1, n, D2)
         case (3); call k_colinv(n, X1, n, D)
         case (4); call k_triu(n, B, n, X1, n)
         case (5); call k_perm_rows(n, B, n, perm, X1, n)
         case (6); call k_iperm_cols(n, B, n, perm, X1, n)
#if !defined(STABLOG)
         case (7); call k_tpup(n, X1, R, D, D2)
#endif
         end select
         !$omp target update from(X1)
         if (pass == 1) X2 = X1
      enddo
      call cmp_m(nm, X1, X2, n)
    end subroutine mat2

  subroutine cmp_c(what, x, y)
    character(len=*), intent(in) :: what
    complex(dp), intent(in) :: x(:), y(:)
    real(dp) :: d
    d = maxval(abs(x - y)) / max(1.d0, maxval(abs(x)))
    write(*,'(2A,ES10.3)') what, ': max rel diff HIP vs OpenMP ', d
    if (.not. (d <= TOL)) then
       write(*,'(A,4ES12.4)') '   OpenMP: ', x(1:min(2,size(x)))
       write(*,'(A,4ES12.4)') '   HIP:    ', y(1:min(2,size(y)))
       call fail(what)
    endif
  end subroutine cmp_c

  subroutine cmp_m(what, x, y, n)
    character(len=*), intent(in) :: what
    integer, intent(in) :: n
    complex(dp), intent(in) :: x(n,n), y(n,n)
    real(dp) :: d
    d = maxval(abs(x - y)) / max(1.d0, maxval(abs(y)))
    write(*,'(A,A,I0,A,ES10.3)') what, ' n=', n, ': max rel diff HIP vs OpenMP ', d
    if (.not. (d <= TOL)) call fail(what)
  end subroutine cmp_m

  subroutine swap_i(a, b)
    integer, intent(inout) :: a, b
    integer :: t
    t = a; a = b; b = t
  end subroutine swap_i

  subroutine fail(what)
    character(len=*), intent(in) :: what
    write(*,'(2A)') 'FAIL: ', what
    nfail = nfail + 1
  end subroutine fail

end program hipk_stab_test

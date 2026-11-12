! Prototype: device-resident blocked Householder QR (LAPACK zgeqrf algorithm: zgeqr2 panel + zlarft + zlarfb) vs oneMKL zgeqrf.
!   variant A: panel = several small kernels per column step (host drives the column loop)
!   variant B: panel = ONE fused kernel (one work-group, barriers between the column steps)
!   trailing update + T: dispatched oneMKL gemm on the device (same for both variants)
! Checks the result against host MKL zgeqrf (R, V, tau).  Standalone, no ALF code.
include "mkl_blas_omp_offload.f90"
include "mkl_lapack_omp_offload.f90"
program devqr_bench
  use omp_lib
  use onemkl_blas_omp_offload_lp64, only: zgemm
  use onemkl_lapack_omp_offload_lp64, only: zgeqrf
  implicit none
  integer, parameter :: TL = 1024
  integer :: nlist(3) = (/ 576, 784, 1024 /), nblist(3) = (/ 16, 32, 64 /)
  integer :: in, ib, n, nb, rep, i, j, info, lw
  complex(8), allocatable :: A0(:,:), A(:,:), Ah(:,:), tau(:), tauh(:), wk(:), Vw(:,:), Tm(:,:), Gm(:,:), W(:,:), W2(:,:), part(:), wsh(:)
  complex(8), allocatable :: whost(:)
  real(8) :: t0, t1, tA, tB, tH, tHost, errA, errB
  integer, parameter :: NREP = 5
  character(len=16) :: ev

  call omp_set_default_device(0)
  do in = 1, 3
    n = nlist(in)
    allocate(A0(n,n), A(n,n), Ah(n,n), tau(n), tauh(n), Vw(n,128), Tm(128,128), Gm(128,128), W(128,n), W2(128,n), part(TL*2), wsh(128))
    lw = 64*n; allocate(wk(lw), whost(lw))
    call fill(A0, n)
    A = 0; tau = 0; Tm = 0; Gm = 0; W = 0; W2 = 0; Vw = 0; part = 0; wsh = 0; wk = 0
    !$omp target enter data map(to: A0) map(alloc: A, tau, Vw, Tm, Gm, W, W2, part, wsh, wk)
    ! host reference (MKL on the host, all threads of OMP/MKL_NUM_THREADS)
    Ah = A0
    t0 = omp_get_wtime()
    call zgeqrf(n, n, Ah, n, tauh, whost, lw, info)
    t1 = omp_get_wtime(); tHost = t1 - t0
    Ah = A0; call zgeqrf(n, n, Ah, n, tauh, whost, lw, info)     ! the reference result (timing above includes warm-up)
    t0 = omp_get_wtime(); Ah = A0; call zgeqrf(n, n, Ah, n, tauh, whost, lw, info); t1 = omp_get_wtime(); tHost = t1 - t0
    ! oneMKL device geqrf (hybrid)
    do rep = 0, NREP
      if (rep == 1) t0 = omp_get_wtime()
      call copy_dev(A, A0, n)
      !$omp dispatch
      call zgeqrf(n, n, A, n, tau, wk, lw, info)
    end do
    t1 = omp_get_wtime(); tH = (t1 - t0) / NREP
    write(*,'(A,I0,A,F9.2,A,F9.2,A)') 'N=', n, ' host MKL zgeqrf ', 1000*tHost, ' ms;  oneMKL device zgeqrf (dispatch) ', 1000*tH, ' ms'
    do ib = 1, 3
      nb = nblist(ib)
      call run('A', 0)
      call run('B', 1)
    end do
    !$omp target exit data map(delete: A0, A, tau, Vw, Tm, Gm, W, W2, part, wsh, wk)
    deallocate(A0, A, Ah, tau, tauh, Vw, Tm, Gm, W, W2, part, wsh, wk, whost)
  end do
contains
  subroutine run(tag, var)
    character(*) :: tag; integer :: var
    real(8) :: tt, err, errt
    do rep = 0, NREP
      if (rep == 1) t0 = omp_get_wtime()
      call copy_dev(A, A0, n)
      call qr_dev(var)
    end do
    t1 = omp_get_wtime(); tt = (t1 - t0) / NREP
    !$omp target update from(A, tau)
    err = maxval(abs(A - Ah)) / maxval(abs(Ah)); errt = maxval(abs(tau - tauh))
    write(*,'(A,I0,A,A,A,F9.2,A,ES9.2,A,ES9.2)') '   nb=', nb, ' variant ', tag, ': ', 1000*tt, ' ms   max rel |A-Ahost| ', err, '  max|tau-tauhost| ', errt
  end subroutine

  subroutine qr_dev(var)
    integer :: var
    integer :: j0, je, mrows, nc, kk
    complex(8), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0), mone = (-1.d0,0.d0)
    do j0 = 1, n, nb
      kk = min(nb, n - j0 + 1)
      je = j0 + kk - 1
      mrows = n - j0 + 1
      if (var == 0) then
        call panel_A(A, n, j0, kk, tau(j0))
      else
        call panel_B(A, n, j0, kk, tau(j0), part, wsh)
      end if
      nc = n - je
      if (nc <= 0) exit
      call build_v(A, n, Vw, j0, kk, mrows)
      call gemm_c('C', 'N', kk, kk, mrows, one, Vw, n, Vw, n, zero, Gm, 128)
      call build_t(Gm, tau(j0), Tm, kk)
      call gemm_c('C', 'N', kk, nc, mrows, one, Vw, n, A(j0,je+1), n, zero, W, 128)
      call gemm_c('C', 'N', kk, nc, kk, one, Tm, 128, W, 128, zero, W2, 128)
      call gemm_c('N', 'N', mrows, nc, kk, mone, Vw, n, W2, 128, one, A(j0,je+1), n)
    end do
  end subroutine

  subroutine gemm_c(ta, tb, m, nn, k, alpha, X, ldx, Y, ldy, beta, Z, ldz)
    character :: ta, tb
    integer :: m, nn, k, ldx, ldy, ldz
    complex(8) :: alpha, beta
    complex(8), target :: X(ldx,*), Y(ldy,*), Z(ldz,*)
    !$omp dispatch
    call zgemm(ta, tb, m, nn, k, alpha, X, ldx, Y, ldy, beta, Z, ldz)
  end subroutine

  subroutine copy_dev(D, S, m)
    integer :: m; complex(8) :: D(m,m), S(m,m); integer :: ii, jj
    !$omp target teams distribute parallel do collapse(2)
    do jj = 1, m
      do ii = 1, m
        D(ii,jj) = S(ii,jj)
      end do
    end do
  end subroutine

  subroutine build_v(Aa, lda, V, j0, kk, mrows)
    integer :: lda, j0, kk, mrows; complex(8) :: Aa(lda,*), V(lda,*); integer :: r, c
    !$omp target teams distribute parallel do collapse(2)
    do c = 1, kk
      do r = 1, mrows
        if (r < c) then
          V(r,c) = (0.d0,0.d0)
        else if (r == c) then
          V(r,c) = (1.d0,0.d0)
        else
          V(r,c) = Aa(j0 + r - 1, j0 + c - 1)
        end if
      end do
    end do
  end subroutine

  ! zlarft (forward, columnwise): T(1:i-1,i) = -tau_i * T(1:i-1,1:i-1) * G(1:i-1,i),  G = V^H V
  subroutine build_t(G, ta, T, kk)
    integer :: kk; complex(8) :: G(128,*), ta(*), T(128,*)
    complex(8) :: tmp(64), s; integer :: i, r, c
    !$omp target map(to: kk) private(tmp, s, i, r, c)
    do i = 1, kk
      do r = 1, kk
        T(r,i) = (0.d0,0.d0)
      end do
      T(i,i) = ta(i)
      do r = 1, i - 1
        tmp(r) = -ta(i) * G(r,i)
      end do
      do r = 1, i - 1
        s = (0.d0,0.d0)
        do c = r, i - 1
          s = s + T(r,c) * tmp(c)
        end do
        T(r,i) = s
      end do
    end do
    !$omp end target
  end subroutine

  ! ---------------- variant A: kernels per column step -----------------------------------------------------------------
  subroutine panel_A(Aa, lda, j0, kk, ta)
    integer :: lda, j0, kk; complex(8) :: Aa(lda,*), ta(*)
    integer :: j, je, r, c, jj
    real(8) :: ss, ar, ai, beta
    complex(8) :: alpha, inv, tj, ctj
    complex(8), allocatable :: tauh_loc(:)
    allocate(tauh_loc(kk))
    je = j0 + kk - 1
    do j = j0, je
      ss = 0.d0
      !$omp target teams distribute parallel do reduction(+:ss)
      do r = j + 1, lda
        ss = ss + real(Aa(r,j))**2 + aimag(Aa(r,j))**2
      end do
      !$omp target update from(Aa(j:j,j))
      alpha = Aa(j,j); ar = real(alpha); ai = aimag(alpha)
      if (ss == 0.d0 .and. ai == 0.d0) then
        tj = (0.d0,0.d0); inv = (1.d0,0.d0); beta = ar
      else
        beta = -sign(sqrt(ar*ar + ai*ai + ss), ar)
        tj = cmplx((beta - ar)/beta, -ai/beta, 8)
        inv = (1.d0,0.d0) / (alpha - beta)
      end if
      ctj = conjg(tj)
      tauh_loc(j-j0+1) = tj
      !$omp target teams distribute parallel do
      do r = j, lda
        if (r == j) then
          Aa(j,j) = beta
        else
          Aa(r,j) = Aa(r,j) * inv
        end if
      end do
      if (j < je) then
        !$omp target teams distribute
        do c = j + 1, je
          block
            complex(8) :: acc, vr
            acc = (0.d0,0.d0)
            !$omp parallel do reduction(+:acc) private(vr)
            do r = j, lda
              if (r == j) then; vr = (1.d0,0.d0); else; vr = Aa(r,j); end if
              acc = acc + conjg(vr) * Aa(r,c)
            end do
            wsh(c - j0 + 1) = acc
          end block
        end do
        !$omp target teams distribute parallel do collapse(2)
        do c = j + 1, je
          do r = j, lda
            block
              complex(8) :: vr
              if (r == j) then; vr = (1.d0,0.d0); else; vr = Aa(r,j); end if
              Aa(r,c) = Aa(r,c) - ctj * vr * wsh(c - j0 + 1)
            end block
          end do
        end do
      end if
    end do
    ta(1:kk) = tauh_loc(1:kk)
    !$omp target update to(ta(1:kk))
    deallocate(tauh_loc)
  end subroutine

  ! ---------------- variant B: fused panel kernel --------------------------------------------------------------------
  subroutine panel_B(Aa, lda, j0, kk, ta, pt, ws)
    integer :: lda, j0, kk; complex(8) :: Aa(lda,*), ta(*), pt(*), ws(*)
    integer :: j, je, r, c, jj, t, nq, q, nc, nqmax
    real(8) :: ss, ar, ai, beta
    complex(8) :: inv, ctj, alpha, acc, vr
    je = j0 + kk - 1
    nqmax = TL
    !$omp target teams num_teams(1) thread_limit(TL) private(ss, ar, ai, beta, inv, ctj, alpha)
    !$omp parallel num_threads(TL) private(acc, vr, c, r, t, nq, q, nc, j, je) shared(ss, inv, ctj, beta, alpha, ar, ai)
    do j = j0, j0 + kk - 1
      je = j0 + kk - 1
      !$omp single
      ss = 0.d0
      !$omp end single
      !$omp do reduction(+:ss)
      do r = j + 1, lda
        ss = ss + real(Aa(r,j))**2 + aimag(Aa(r,j))**2
      end do
      !$omp end do
      !$omp single
      alpha = Aa(j,j); ar = real(alpha); ai = aimag(alpha)
      if (ss == 0.d0 .and. ai == 0.d0) then
        ta(j-j0+1) = (0.d0,0.d0); inv = (1.d0,0.d0); beta = ar
      else
        beta = -sign(sqrt(ar*ar + ai*ai + ss), ar)
        ta(j-j0+1) = cmplx((beta - ar)/beta, -ai/beta, 8)
        inv = (1.d0,0.d0) / (alpha - beta)
      end if
      ctj = conjg(ta(j-j0+1))
      !$omp end single
      !$omp do
      do r = j + 1, lda
        Aa(r,j) = Aa(r,j) * inv
      end do
      !$omp end do
      !$omp single
      Aa(j,j) = cmplx(beta, 0.d0, 8)
      !$omp end single
      nc = je - j
      if (nc > 0) then
        nq = max(1, TL / nc)
        !$omp do
        do t = 0, nq*nc - 1
          q = mod(t, nq); c = j + 1 + t / nq
          acc = (0.d0,0.d0)
          do r = j + q, lda, nq
            if (r == j) then; vr = (1.d0,0.d0); else; vr = Aa(r,j); end if
            acc = acc + conjg(vr) * Aa(r,c)
          end do
          pt(t + 1) = acc
        end do
        !$omp end do
        !$omp do
        do c = j + 1, je
          acc = (0.d0,0.d0)
          do q = 0, nq - 1
            acc = acc + pt(q + (c - j - 1)*nq + 1)
          end do
          ws(c - j0 + 1) = acc
        end do
        !$omp end do
        !$omp do collapse(2)
        do c = j + 1, je
          do r = j, lda
            if (r == j) then; vr = (1.d0,0.d0); else; vr = Aa(r,j); end if
            Aa(r,c) = Aa(r,c) - ctj * vr * ws(c - j0 + 1)
          end do
        end do
        !$omp end do
      end if
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine

  subroutine fill(M, m_)
    integer :: m_; complex(8) :: M(m_,m_); integer :: ii, jj; real(8) :: a, b
    call srand(11)
    do jj = 1, m_
      do ii = 1, m_
        call random_number(a); call random_number(b)
        M(ii,jj) = cmplx(a - 0.5d0, b - 0.5d0, 8) * 10.d0**(-3.d0*real(jj-1,8)/real(m_-1,8))
      end do
    end do
  end subroutine
end program

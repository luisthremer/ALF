! Benchmark sections.  Each section: set up host data, map it to the device once
! (!$omp target enter data), check the device result against bench_ref, then time.
! Host<->device transfers in this file are setup/verification only (untimed), except the
! deliberate round trip of the "+rt" update variants.
module bench_sec
  use bench_la
  use bench_kern
  use bench_ref
  use bench_util
  implicit none
  external :: zgemm, zgeru
  integer :: kcur = 1
  real(8) :: opaque_one = 1d0, sink = 0d0   ! defeat dead-code elimination of the host round trip
  complex(8), parameter :: cone = (1d0, 0d0), czero = (0d0, 0d0)
contains

  ! ---------------------------------------------------------------- checkerboard B_T
  subroutine sec_cb(N)
    integer, intent(in) :: N
    integer :: L, nb, f, x, y, b
    integer, allocatable :: bi(:,:), bj(:,:)
    complex(8), allocatable :: Gd(:,:), G0(:,:), Gr(:,:)
    real(8) :: c, s, us
    L = nint(sqrt(dble(N))); nb = N / 2
    c = cosh(0.1d0); s = sinh(0.1d0)
    allocate(bi(nb,4), bj(nb,4), Gd(N,N), G0(N,N), Gr(N,N))
    do f = 1, 4
      b = 0
      do y = 0, L-1
        do x = 0, L-1
          if (f <= 2 .and. mod(x,2) == f-1) then
            b = b + 1; bi(b,f) = 1 + x + L*y; bj(b,f) = 1 + mod(x+1,L) + L*y
          else if (f >= 3 .and. mod(y,2) == f-3) then
            b = b + 1; bi(b,f) = 1 + x + L*y; bj(b,f) = 1 + x + L*mod(y+1,L)
          end if
        end do
      end do
    end do
    call rand_c(G0, N*N, 1d0); Gd = G0; Gr = G0
    !$omp target enter data map(to: Gd, bi, bj)
    call cb_apply(Gd, N, nb, bi, bj, c, s)
    !$omp target update from(Gd)       ! verification only
    call ref_cb(Gr, N, nb, bi, bj, c, s)
    call check('checkerboard B_T vs host', N, relfro(Gd, Gr, N*N), 1d-12)
    call timeit(op, us)
    call add_rec('checkerboard B_T G->B G B^-1', '8 kernels', N, us, 8d0*2d0*16d0*N*N/us*1d-3, 'GB/s')
    !$omp target exit data map(delete: Gd, bi, bj)
  contains
    subroutine op()
      call cb_apply(Gd, N, nb, bi, bj, c, s)
    end subroutine
  end subroutine

  ! ---------------------------------------------------------------- dense B_T
  subroutine sec_dense(N)
    integer, intent(in) :: N
    integer :: L
    real(8), allocatable :: E(:,:), Ei(:,:)
    complex(8), allocatable :: Ec(:,:), Eic(:,:), G0(:,:), Gd(:,:), Gr(:,:), T1(:,:), T2(:,:), T3(:,:), H(:,:)
    real(8) :: us, fl
    L = nint(sqrt(dble(N)))
    allocate(E(N,N), Ei(N,N), Ec(N,N), Eic(N,N), G0(N,N), Gd(N,N), Gr(N,N), T1(N,N), T2(N,N), T3(N,N), H(N,N))
    call ref_hop_exp(L, 0.1d0, E, Ei)
    Ec = E; Eic = Ei
    call rand_c(G0, N*N, 1d0); Gd = G0
    call zgemm('N', 'N', N, N, N, cone, Ec, N, G0, N, czero, H, N)
    call zgemm('N', 'N', N, N, N, cone, H, N, Eic, N, czero, Gr, N)
    T1 = 0; T2 = 0; T3 = 0
    !$omp target enter data map(to: E, Ei, Ec, Eic, Gd, T1, T2, T3)
    call v1()
    !$omp target update from(T2)
    call check('dense E G E^-1, ZGEMM x2 vs host', N, relfro(T2, Gr, N*N), 1d-12)
    call v2()
    !$omp target update from(T3)
    call check('dense E G E^-1, DGEMM+transpose vs host', N, relfro(T3, Gr, N*N), 1d-12)
    fl = 16d0 * dble(N)**3      ! complex-equivalent flops of two complex products
    call timeit(v1, us);  call add_rec('dense B_T G->E G E^-1', 'zgemm x2', N, us, fl/us*1d-3, 'GFLOP/s')
    call timeit(v2, us);  call add_rec('dense B_T G->E G E^-1', 'dgemm+transp x2', N, us, fl/us*1d-3, 'GFLOP/s')
    call timeit(v3, us);  call add_rec('dense B_T G->E G E^-1', 'zgemm x1', N, us, 0.5d0*fl/us*1d-3, 'GFLOP/s')
    call timeit(v4, us);  call add_rec('dense B_T G->E G E^-1', 'dgemm-R x1', N, us, 0.5d0*fl/us*1d-3, 'GFLOP/s')
    call timeit(v5, us);  call add_rec('dense B_T G->E G E^-1', 'transpose x1', N, us, 2d0*16d0*N*N/us*1d-3, 'GFLOP/s')
    !$omp target exit data map(delete: E, Ei, Ec, Eic, Gd, T1, T2, T3)
  contains
    subroutine v1()     ! two complex GEMMs with complex copies of E, E^-1: result in T2
      call la_gemm(LA_N, LA_N, N, N, N, cone, Ec, N, Gd, N, czero, T1, N)
      call la_gemm(LA_N, LA_N, N, N, N, cone, T1, N, Eic, N, czero, T2, N)
    end subroutine
    subroutine v2()     ! right: DGEMM on real view; left via E H = (H^T E)^T (E symmetric): result in T3
      call la_gemm_cr(N, N, N, Gd, N, Ei, N, T1, N)
      call k_transpose(T1, T2, N)
      call la_gemm_cr(N, N, N, T2, N, E, N, T1, N)
      call k_transpose(T1, T3, N)
    end subroutine
    subroutine v3()
      call la_gemm(LA_N, LA_N, N, N, N, cone, Ec, N, Gd, N, czero, T1, N)
    end subroutine
    subroutine v4()
      call la_gemm_cr(N, N, N, Gd, N, Ei, N, T1, N)
    end subroutine
    subroutine v5()
      call k_transpose(Gd, T1, N)
    end subroutine
  end subroutine

  ! ---------------------------------------------------------------- Op_V diagonal wrap
  subroutine sec_opv(N)
    integer, intent(in) :: N
    complex(8), allocatable :: Gd(:,:), G0(:,:), Gr(:,:)
    real(8), allocatable :: ea(:), a(:)
    real(8) :: us
    integer :: i
    allocate(Gd(N,N), G0(N,N), Gr(N,N), ea(N), a(N))
    call random_number(a); ea = exp(a - 0.5d0)
    call rand_c(G0, N*N, 1d0); Gd = G0; Gr = G0
    call ref_opv(Gr, N, ea)
    !$omp target enter data map(to: Gd, ea)
    call va()
    !$omp target update from(Gd)       ! verification only
    call check('Op_V wrap per-site launches vs host', N, relfro(Gd, Gr, N*N), 1d-12)
    Gd = G0
    !$omp target update to(Gd)         ! verification only (restore G)
    call vb()
    !$omp target update from(Gd)       ! verification only
    call check('Op_V wrap fused kernel vs host', N, relfro(Gd, Gr, N*N), 1d-12)
    call timeit(va, us); call add_rec('Op_V wrap, N sites', 'per-site x N', N, us, 2d0*16d0*N*N/us*1d-3, 'GB/s')
    call timeit(vb, us); call add_rec('Op_V wrap, N sites', 'fused', N, us, 2d0*16d0*N*N/us*1d-3, 'GB/s')
    !$omp target exit data map(delete: Gd, ea)
  contains
    subroutine va()
      integer :: i
      do i = 1, N
        call k_opv_site(Gd, N, i, ea(i))   ! ea(i): host copy is current, passed by value (scalar)
      end do
    end subroutine
    subroutine vb()
      call k_opv_all(Gd, N, ea)
    end subroutine
  end subroutine

  ! ---------------------------------------------------------------- stabilization
  subroutine sec_qr(N)
    integer, intent(in) :: N
    complex(8), allocatable :: A0(:,:), Ad(:,:), Apm(:,:), Af(:,:), Q(:,:), Cc(:,:), tau(:), work(:)
    complex(8), allocatable :: FA(:,:), QA(:,:), tauA(:), FB(:,:), QB(:,:), tauB(:)
    real(8), allocatable :: nrm(:), d3(:)
    integer, allocatable :: perm(:), jpvt(:), permB(:)
    integer(4), allocatable :: dinfo(:)
    integer :: lw, lw1, lw2, lw3, ig, j
    real(8) :: us
    allocate(A0(N,N), Ad(N,N), Apm(N,N), Af(N,N), Q(N,N), Cc(N,N), tau(N), nrm(N), d3(N))
    allocate(FA(N,N), QA(N,N), tauA(N), FB(N,N), QB(N,N), tauB(N), perm(N), permB(N), jpvt(N), dinfo(1))
    Apm = 0; Af = 0; Q = 0; Cc = 0; tau = 0; nrm = 0; perm = 0; dinfo = 0
    call ref_graded(N, 0d0, -12d0, A0)
    Ad = A0
    allocate(work(1)); work = 0
    !$omp target enter data map(to: Ad, Apm, Af, Q, Cc, tau, nrm, perm, dinfo)
    ! workspace queries need device-resident arguments, so they come after the mapping
    lw1 = la_geqrf_lwork(N, N, Ad, N)
    lw2 = la_ungqr_lwork(N, N, N, Ad, N, tau)
    lw3 = la_ormqr_lwork(LA_C, N, N, N, Ad, N, tau, Cc, N)
    lw = max(lw1, lw2, lw3, 1)
    deallocate(work); allocate(work(lw)); work = 0
    !$omp target enter data map(to: work)
    do ig = 1, 4
      select case (ig)
      case (1)
        call ref_graded(N, 0d0, -12d0, A0)
      case (2)
        call ref_graded(N, 12d0, -12d0, A0)
      case (4)     ! mildly graded (cond ~1e3): there the ticket's criterion (ii) is attainable
        call ref_graded(N, 0d0, -3d0, A0)
      case (3)     ! columns of the graded matrix scaled by random 10^(-12 u): real pivoting needed
        call ref_graded(N, 0d0, -12d0, A0)
        call random_number(d3)
        do j = 1, N
          A0(:,j) = A0(:,j) * 10d0**(-12d0*d3(j))
        end do
      end select
      Ad = A0
      !$omp target update to(Ad)         ! verification only (new input matrix)
      call qr_check(ig)
      if (ig == 1) call qr_time()
    end do
    !$omp target exit data map(delete: Ad, Apm, Af, Q, Cc, tau, nrm, perm, dinfo, work)
  contains
    subroutine pre()
      call k_colnorm(Ad, N, nrm)
      call k_sortperm(nrm, N, perm)
      call k_permcols(Ad, N, perm, Apm)
    end subroutine
    subroutine fac()
      call k_copy(Apm, Af, N*N)
      call la_geqrf(N, N, Af, N, tau, work, lw, dinfo)
    end subroutine
    subroutine gen()
      call k_copy(Af, Q, N*N)
      call la_ungqr(N, N, N, Q, N, tau, work, lw, dinfo)
    end subroutine
    subroutine app()
      call k_copy(Apm, Cc, N*N)
      call la_ormqr(LA_C, N, N, N, Af, N, tau, Cc, N, work, lw, dinfo)
    end subroutine
    subroutine tot()
      call pre(); call fac(); call gen()
    end subroutine
    subroutine hostA()
      call ref_qp3(N, A0, FA, QA, tauA, jpvt)
    end subroutine
    subroutine hostB()
      call ref_prepiv_qr(N, A0, FB, QB, tauB, permB)
    end subroutine

    subroutine qr_check(ig)
      integer, intent(in) :: ig
      complex(8), allocatable :: AP(:,:), R(:,:), QR(:,:)
      real(8) :: e, rd, rmin, rmax, ra, rb
      real(8), allocatable :: rdu(:)
      real(8) :: rtmp
      integer :: i, j, nmis
      character(len=8) :: tag
      write(tag,'(a,i1)') 'mat', ig
      allocate(AP(N,N), R(N,N), QR(N,N))
      call tot(); call app()
      !$omp target update from(perm, Af, Q, Cc)   ! verification only
      do j = 1, N
        AP(:,j) = A0(:,perm(j))
      end do
      R = 0
      do j = 1, N
        R(1:j,j) = Af(1:j,j)
      end do
      call zgemm('N', 'N', N, N, N, cone, Q, N, R, N, czero, QR, N)
      call check('QR '//tag//' (i) |AP-QR|/|A| device', N, relfro(QR, AP, N*N), 1d-12)
      call zgemm('C', 'N', N, N, N, cone, Q, N, Q, N, czero, QR, N)
      do j = 1, N
        QR(j,j) = QR(j,j) - 1d0
      end do
      call report('QR '//tag//' |Q^H Q - I|_F device', N, sqrt(sum(abs(QR)**2)))
      call check('QR '//tag//' ormqr Q^H(AP) = R', N, relfro(Cc, R, N*N), 1d-12)
      call hostB(); call hostA()
      nmis = count(perm /= permB)
      call check('QR '//tag//' perm device vs host B (#diff)', N, dble(nmis), 0.5d0)
      rd = 0d0
      do i = 1, N
        ra = abs(Af(i,i)); rb = abs(FB(i,i))
        rd = max(rd, abs(ra - rb) / rb)
      end do
      call check('QR '//tag//' (ii) |R_ii| device vs host B', N, rd, 1d-10)
      rmin = huge(1d0); rmax = 0d0
      do i = 1, N
        e = abs(FA(i,i)) / abs(FB(i,i))
        rmin = min(rmin, e); rmax = max(rmax, e)
      end do
      rtmp = 0d0
      do i = 1, N
        rtmp = max(rtmp, abs(abs(Af(i,i)) - abs(FB(i,i))) / sqrt(sum(abs(AP(:,i))**2)))
      end do
      call check('QR '//tag//' (ii-b) |dR_ii|/|col_i| dev vs B', N, rtmp, 1d-12)
      allocate(rdu(N))
      call ref_prepiv_rdiag_unblocked(N, A0, rdu)
      rtmp = 0d0
      do i = 1, N
        rtmp = max(rtmp, abs(rdu(i) - abs(FB(i,i))) / abs(FB(i,i)))
      end do
      call report('QR '//tag//' (ii) host B blocked vs unblocked', N, rtmp)
      call report('QR '//tag//' (iii) min |Rii|qp3/|Rii|prepiv', N, rmin)
      call report('QR '//tag//' (iii) max |Rii|qp3/|Rii|prepiv', N, rmax)
      call report('QR '//tag//' #columns moved by ZGEQP3 vs input', N, dble(count(jpvt /= [(i, i=1,N)])))
      call report('QR '//tag//' #columns moved by pre-pivot', N, dble(count(perm /= [(i, i=1,N)])))
    end subroutine

    subroutine qr_time()
      character(len=40) :: sec
      sec = 'stabilization QR (graded matrix)'
      call tot(); call app()      ! leaves Apm/Af/tau consistent for the partial timings
      call timeit(pre, us);   call add_rec(sec, 'dev norm+sort+perm', N, us, 0d0, '')
      call timeit(fac, us);   call add_rec(sec, 'dev geqrf', N, us, 0d0, '')
      call timeit(gen, us);   call add_rec(sec, 'dev ungqr', N, us, 0d0, '')
      call timeit(app, us);   call add_rec(sec, 'dev ormqr', N, us, 0d0, '')
      call timeit(tot, us);   call add_rec(sec, 'dev total (QR+Q)', N, us, 0d0, '')
      call timeit(hostA, us); call add_rec(sec, 'host A qp3+ungqr', N, us, 0d0, '')
      call timeit(hostB, us); call add_rec(sec, 'host B pre+geqrf+ung', N, us, 0d0, '')
    end subroutine
  end subroutine

  ! ---------------------------------------------------------------- updates G += x y^T
  subroutine sec_upd(N)
    integer, intent(in) :: N
    integer, parameter :: NU = 256
    complex(8), allocatable :: X(:,:), Y(:,:), G0(:,:), Gd(:,:), Gr(:,:)
    integer :: k, kk
    real(8) :: us
    character(len=16) :: vn
    allocate(X(N,NU), Y(NU,N), G0(N,N), Gd(N,N), Gr(N,N))
    call rand_c(G0, N*N, 1d0)
    call rand_c(X, N*NU, 0.1d0); call rand_c(Y, N*NU, 0.1d0)
    !$omp target enter data map(to: X, Y, Gd)
    ! ---- correctness: 256 sequential rank-1 updates
    Gr = G0
    do k = 1, NU
      call zgeru(N, N, cone, X(1,k), 1, Y(k,1), NU, Gr, N)
    end do
    Gd = G0
    !$omp target update to(Gd)         ! verification only (restore G)
    call r1_ker()
    !$omp target update from(Gd)       ! verification only
    call check('rank-1 x256 kernel vs host zgeru', N, relfro(Gd, Gr, N*N), 1d-12)
    Gd = G0
    !$omp target update to(Gd)         ! verification only
    call r1_nw()
    !$omp target update from(Gd)       ! verification only
    call check('rank-1 x256 kernel nowait vs host zgeru', N, relfro(Gd, Gr, N*N), 1d-12)
    Gd = G0
    !$omp target update to(Gd)         ! verification only
    call r1_blas()
    !$omp target update from(Gd)       ! verification only
    call check('rank-1 x256 la_geru vs host zgeru', N, relfro(Gd, Gr, N*N), 1d-12)
    ! ---- correctness: rank-128 as one GEMM
    Gr = G0
    call zgemm('N', 'N', N, N, 128, cone, X, N, Y, NU, cone, Gr, N)
    Gd = G0
    !$omp target update to(Gd)         ! verification only
    call la_gemm(LA_N, LA_N, N, N, 128, cone, X, N, Y, NU, cone, Gd, N)
    !$omp target update from(Gd)       ! verification only
    call check('rank-128 one ZGEMM vs host zgemm', N, relfro(Gd, Gr, N*N), 1d-12)
    ! ---- timing
    call timeit(r1_ker, us)
    call add_rec('rank-1 update [us/update]', 'kernel', N, us/NU, 2d0*16d0*N*N/(us/NU)*1d-3, 'GB/s')
    call timeit(r1_nw, us)
    call add_rec('rank-1 update [us/update]', 'kernel nowait', N, us/NU, 2d0*16d0*N*N/(us/NU)*1d-3, 'GB/s')
    call timeit(r1_blas, us)
    call add_rec('rank-1 update [us/update]', 'cublas zgeru', N, us/NU, 2d0*16d0*N*N/(us/NU)*1d-3, 'GB/s')
    call timeit(r1_ker_rt, us)
    call add_rec('rank-1 update [us/update]', 'kernel+roundtrip', N, us/NU, 2d0*16d0*N*N/(us/NU)*1d-3, 'GB/s')
    call timeit(r1_blas_rt, us)
    call add_rec('rank-1 update [us/update]', 'zgeru+roundtrip', N, us/NU, 2d0*16d0*N*N/(us/NU)*1d-3, 'GB/s')
    do kk = 0, 7
      kcur = 2**kk
      k = kcur
      call timeit(rk, us)
      write(vn,'(a,i0)') 'k=', k
      call add_rec('rank-k update [us/update]', trim(vn), N, us/k, 8d0*N*N*k/us*1d-3, 'GFLOP/s')
    end do
    !$omp target exit data map(delete: X, Y, Gd)
  contains
    subroutine r1_ker()
      integer :: k
      do k = 1, NU
        call k_ger(Gd, N, X(1,k), Y(k,1), NU)
      end do
    end subroutine
    subroutine r1_nw()
      integer :: k
      do k = 1, NU
        call k_ger_nw(Gd, N, X(1,k), Y(k,1), NU)
      end do
      !$omp taskwait
    end subroutine
    subroutine r1_blas()
      integer :: k
      do k = 1, NU
        call la_geru(N, N, cone, X(1,k), 1, Y(k,1), NU, Gd, N)
      end do
    end subroutine
    subroutine r1_ker_rt()
      integer :: k, i
      do k = 1, NU
        i = mod(k-1, N) + 1
        ! unavoidable in the real code: the accept/reject test (RNG, host) needs a ratio from G
        !$omp target update from(Gd(i,i))
        ! ... and the decision (here: the same scalar times an opaque 1.0) goes back to the device
        sink = sink + real(Gd(i,i))
        Gd(i,i) = Gd(i,i) * opaque_one
        !$omp target update to(Gd(i,i))
        call k_ger(Gd, N, X(1,k), Y(k,1), NU)
      end do
    end subroutine
    subroutine r1_blas_rt()
      integer :: k, i
      do k = 1, NU
        i = mod(k-1, N) + 1
        !$omp target update from(Gd(i,i))   ! as in r1_ker_rt
        sink = sink + real(Gd(i,i))
        Gd(i,i) = Gd(i,i) * opaque_one
        !$omp target update to(Gd(i,i))
        call la_geru(N, N, cone, X(1,k), 1, Y(k,1), NU, Gd, N)
      end do
    end subroutine
    subroutine rk()
      call la_gemm(LA_N, LA_N, N, N, kcur, cone, X, N, Y, NU, cone, Gd, N)
    end subroutine
  end subroutine

  ! ---------------------------------------------------------------- flavour batching
  subroutine sec_batch(N)
    integer, intent(in) :: N
    complex(8), allocatable :: A(:,:,:), B(:,:,:), C(:,:,:), Cr(:,:,:)
    real(8) :: us, fl, e
    integer :: f
    allocate(A(N,N,2), B(N,N,2), C(N,N,2), Cr(N,N,2))
    call rand_c(A, 2*N*N, 1d0); call rand_c(B, 2*N*N, 1d0); C = 0
    do f = 1, 2
      call zgemm('N', 'N', N, N, N, cone, A(1,1,f), N, B(1,1,f), N, czero, Cr(1,1,f), N)
    end do
    !$omp target enter data map(to: A, B, C)
    call two()
    !$omp target update from(C)         ! verification only
    call check('2x zgemm vs host', N, relfro(C, Cr, 2*N*N), 1d-12)
    C = 0
    !$omp target update to(C)           ! verification only
    call bat()
    !$omp target update from(C)         ! verification only
    call check('batched zgemm (2 flavours) vs host', N, relfro(C, Cr, 2*N*N), 1d-12)
    fl = 8d0 * dble(N)**3
    call timeit(one, us); call add_rec('flavour batching, 2 N x N ZGEMM', '1 zgemm', N, us, fl/us*1d-3, 'GFLOP/s')
    call timeit(two, us); call add_rec('flavour batching, 2 N x N ZGEMM', '2 zgemm calls', N, us, 2*fl/us*1d-3, 'GFLOP/s')
    call timeit(bat, us); call add_rec('flavour batching, 2 N x N ZGEMM', 'strided batched', N, us, 2*fl/us*1d-3, 'GFLOP/s')
    !$omp target exit data map(delete: A, B, C)
  contains
    subroutine one()
      call la_gemm(LA_N, LA_N, N, N, N, cone, A(1,1,1), N, B(1,1,1), N, czero, C(1,1,1), N)
    end subroutine
    subroutine two()
      call la_gemm(LA_N, LA_N, N, N, N, cone, A(1,1,1), N, B(1,1,1), N, czero, C(1,1,1), N)
      call la_gemm(LA_N, LA_N, N, N, N, cone, A(1,1,2), N, B(1,1,2), N, czero, C(1,1,2), N)
    end subroutine
    subroutine bat()
      call la_gemm_batched(LA_N, LA_N, N, N, N, cone, A, N, N*N, B, N, N*N, czero, C, N, N*N, 2)
    end subroutine
  end subroutine

end module bench_sec

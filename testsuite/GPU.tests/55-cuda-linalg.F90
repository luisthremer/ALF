!> @brief Device branches of linalg_mod (GPU=nvidia, linalg_cuda_mod) against the host calls on the same data.
!> @details Each wrapper with a device path is called once on host arrays (host LAPACK/BLAS) and once on copies mapped to
!> the device (cuBLAS/cuSOLVER through the presence dispatch); results agree to TOL (relative Frobenius norm). Also: the
!> pre-pivoted QR (same permutation as the host, ties keep their order, NaN -> info = -101, bitwise identical on a re-run),
!> the reflector determinant of the device QR (speedup S1: det H_i = -tau_i/conj(tau_i)) against the LU phase.
!> Argument "mixed": one operand mapped, one not; the run must stop (registered with WILL_FAIL).
module cuda_linalg_util
  implicit none
  complex(8), parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)
contains
  subroutine fill(a, seed)
    complex(8), intent(out) :: a(:,:)
    integer, intent(in) :: seed
    integer :: i, j
    do j = 1, size(a,2)
       do i = 1, size(a,1)
          ! hash-like pseudo-random entries in [-0.5, 0.5): full rank (a sum of sines would have rank <= 4)
          a(i,j) = cmplx(hrand(i, j, seed), hrand(j, i, seed + 7), kind=8)
       enddo
    enddo
  end subroutine fill
  real(8) function hrand(i, j, seed)
    integer, intent(in) :: i, j, seed
    real(8) :: x
    x = sin(12.9898d0*i + 78.233d0*j + 37.719d0*seed)*43758.5453d0
    hrand = x - floor(x) - 0.5d0
  end function hrand
  real(8) function rel(a, b)
    complex(8), intent(in) :: a(:,:), b(:,:)
    rel = sqrt(sum(abs(a - b)**2))/max(sqrt(sum(abs(b)**2)), tiny(1.d0))
  end function rel
end module cuda_linalg_util

program cuda_linalg
  use cuda_linalg_util
  use linalg_mod, only: la_zgemm, la_zgemm_cr, la_ztrsm, la_zgeqrf, la_zunmqr, la_zgetrf, la_zgetri, la_zqrpp
  use linalg_cuda_mod, only: la_cuda_sync, la_cuda_finalize
  implicit none
  real(8), parameter :: TOL = 1.d-12
  integer, parameter :: M = 37, N = 29, K = 41, LD = 45, NQ = 200
  complex(8), allocatable :: A(:,:), B(:,:), C(:,:), Ch(:,:), Ad(:,:), Bd(:,:), Cd(:,:), tau(:), taud(:), w(:)
  real(8), allocatable :: R(:,:)
  integer, allocatable :: ip(:), ipd(:)
  character(len=1), parameter :: ops(4) = ['N', 'T', 'C', 'N']
  character(len=1) :: ta, tb
  character(len=16) :: arg
  integer :: nfail, i, info, lw, t
  complex(8) :: ph_lu, ph_tr, wq(1)

  nfail = 0
  call get_command_argument(1, arg)
  if (trim(arg) == 'mixed') then
     allocate(A(N,N), B(N,N), C(N,N)); call fill(A, 1); call fill(B, 2); C = z0
     !$omp target enter data map(to: A, C)
     call la_zgemm('N', 'N', N, N, N, z1, A, N, B, N, z0, C, N)      ! B is not mapped: must terminate
     write(*,'(A)') '55-cuda-linalg mixed: NOT stopped'
     error stop 2
  endif

  ! ---- zgemm, all op pairs, lda > m
  do t = 1, 4
     ta = ops(t); tb = ops(5-t)
     allocate(A(LD,max(M,K)), B(LD,max(N,K)), C(LD,N))
     call fill(A, t); call fill(B, t+10); call fill(C, t+20)
     Ch = C
     call la_zgemm(ta, tb, M, N, K, (0.7d0,0.2d0), A, LD, B, LD, (0.3d0,-0.1d0), Ch, LD)
     !$omp target enter data map(to: A, B, C)
     call la_zgemm(ta, tb, M, N, K, (0.7d0,0.2d0), A, LD, B, LD, (0.3d0,-0.1d0), C, LD)
     !$omp target exit data map(from: C) map(delete: A, B)
     call check('zgemm '//ta//tb, rel(C(1:M,:), Ch(1:M,:)))
     deallocate(A, B, C, Ch)
  enddo

  ! ---- zgemm_cr (complex x real)
  allocate(A(LD,K), R(K,N), C(LD,N)); call fill(A, 3); C = z0
  do i = 1, N
     R(:,i) = [(sin(0.3d0*t + i), t = 1, K)]
  enddo
  Ch = C
  call la_zgemm_cr(M, N, K, A, LD, R, K, Ch, LD)
  !$omp target enter data map(to: A, R, C)
  call la_zgemm_cr(M, N, K, A, LD, R, K, C, LD)
  !$omp target exit data map(from: C) map(delete: A, R)
  call check('zgemm_cr', rel(C(1:M,:), Ch(1:M,:)))
  deallocate(A, R, C, Ch)

  ! ---- ztrsm (the variants the GPU path uses)
  call trsm('L', 'U', 'N')
  call trsm('R', 'U', 'C')
  call trsm('L', 'U', 'C')

  ! ---- zgeqrf + zunmqr (L N, L C, R N): factors and Q applied
  allocate(A(NQ,NQ), tau(NQ), B(NQ,NQ)); call fill(A, 5); call fill(B, 6)
  Ad = A; Bd = B; allocate(taud(NQ))
  call la_zgeqrf(NQ, NQ, A, NQ, tau, wq, -1, info); lw = int(dble(wq(1))); allocate(w(lw))
  call la_zgeqrf(NQ, NQ, A, NQ, tau, w, lw, info)
  !$omp target enter data map(to: Ad) map(alloc: taud)
  call la_zgeqrf(NQ, NQ, Ad, NQ, taud, w, lw, info)
  !$omp target update from(Ad, taud)
  call check('zgeqrf R', rel(Ad, A))
  call check('zgeqrf tau', rel(reshape(taud, [NQ,1]), reshape(tau, [NQ,1])))
  call unmqr('L', 'N'); call unmqr('L', 'C'); call unmqr('R', 'N')
  ! S1: phase of det(A) from the device reflectors and R, against the host LU phase
  ph_tr = z1
  do i = 1, NQ
     if (taud(i) /= z0) ph_tr = ph_tr*(-taud(i)/conjg(taud(i)))
     ph_tr = ph_tr*(Ad(i,i)/abs(Ad(i,i)))
  enddo
  call fill(A, 5); allocate(ip(NQ))
  call la_zgetrf(NQ, NQ, A, NQ, ip, info)
  ph_lu = z1
  do i = 1, NQ
     ph_lu = ph_lu*(A(i,i)/abs(A(i,i)))
     if (ip(i) /= i) ph_lu = -ph_lu
  enddo
  call check('S1 reflector det phase vs LU', abs(ph_tr - ph_lu))
  !$omp target exit data map(delete: Ad, taud)
  deallocate(A, B, Ad, Bd, tau, taud, w, ip)

  ! ---- zgetrf + zgetri (getrs on the identity on the device)
  allocate(A(NQ,NQ), ip(NQ), ipd(NQ)); call fill(A, 7)
  do i = 1, NQ
     A(i,i) = A(i,i) + 4.d0
  enddo
  Ad = A
  call la_zgetri(NQ, A, NQ, ip, wq, -1, info); lw = int(dble(wq(1))); allocate(w(lw))
  call la_zgetrf(NQ, NQ, A, NQ, ip, info)
  call la_zgetri(NQ, A, NQ, ip, w, lw, info)
  !$omp target enter data map(to: Ad) map(alloc: ipd)
  call la_zgetrf(NQ, NQ, Ad, NQ, ipd, info)
  call la_zgetri(NQ, Ad, NQ, ipd, w, lw, info)
  !$omp target exit data map(from: Ad) map(delete: ipd)
  call check('zgetrf + zgetri', rel(Ad, A))
  deallocate(A, Ad, ip, ipd, w)

  ! ---- la_zqrpp: permutation as on the host (ties: duplicated columns), factors, NaN, determinism
  allocate(A(NQ,NQ), tau(NQ), taud(NQ), ip(NQ), ipd(NQ)); call fill(A, 8)
  call ties(A)                                                        ! equal norms: the order must be kept
  Ad = A
  call la_zqrpp(NQ, NQ, A, NQ, ip, tau, wq, -1, info); lw = int(dble(wq(1))); allocate(w(lw))
  call la_zqrpp(NQ, NQ, A, NQ, ip, tau, w, lw, info)
  !$omp target enter data map(to: Ad) map(alloc: taud, ipd)
  call la_zqrpp(NQ, NQ, Ad, NQ, ipd, taud, w, lw, info)
  !$omp target update from(Ad, taud)
  call check('zqrpp info', dble(abs(info)))
  call check('zqrpp permutation (host copy)', dble(count(ip /= ipd)))
  call check('zqrpp factors', rel(Ad, A))
  B = Ad
  call fill(A, 8); call ties(A)
  Ad = A
  !$omp target update to(Ad)
  call la_zqrpp(NQ, NQ, Ad, NQ, ipd, taud, w, lw, info)
  !$omp target update from(Ad)
  call check('zqrpp bitwise re-run', merge(0.d0, 1.d0, all(Ad == B)))
  Ad(5,7) = transfer(int(z'7FF8000000000000', 8), 0.d0)
  !$omp target update to(Ad)
  call la_zqrpp(NQ, NQ, Ad, NQ, ipd, taud, w, lw, info)
  call check('zqrpp NaN -> info -101', dble(abs(info + 101)))
  !$omp target exit data map(delete: Ad, taud, ipd)

  call la_cuda_sync()
  call la_cuda_finalize()
  if (nfail > 0) then
     write(*,'(A,I0,A)') '55-cuda-linalg: FAIL (', nfail, ')'
     error stop 1
  endif
  write(*,'(A)') '55-cuda-linalg: PASS'

contains
  subroutine check(what, d)
    character(len=*), intent(in) :: what
    real(8), intent(in) :: d
    write(*,'(A,T40,ES10.3)') what, d
    if (.not. (d <= TOL)) nfail = nfail + 1
  end subroutine check

  !> columns 10 and 50 with the norm of column 3 bitwise (same squares in the same order), linearly independent
  subroutine ties(X)
    complex(8), intent(inout) :: X(:,:)
    integer :: i
    X(:,10) = conjg(X(:,3))
    X(:,50) = X(:,3)*[((-1)**i, i = 1, size(X,1))]                     ! (re/im swap would be i conj(X(:,3)): dependent)
  end subroutine ties

  subroutine trsm(side, uplo, tr)
    character(len=1), intent(in) :: side, uplo, tr
    complex(8), allocatable :: Tm(:,:), X(:,:), Xh(:,:)
    integer :: nt, j
    nt = merge(M, N, side == 'L')
    allocate(Tm(LD,nt), X(LD,N)); call fill(Tm, 11); call fill(X, 12)
    do j = 1, nt
       Tm(j,j) = Tm(j,j) + 3.d0
    enddo
    Xh = X
    call la_ztrsm(side, uplo, tr, 'N', M, N, (1.1d0,0.d0), Tm, LD, Xh, LD)
    !$omp target enter data map(to: Tm, X)
    call la_ztrsm(side, uplo, tr, 'N', M, N, (1.1d0,0.d0), Tm, LD, X, LD)
    !$omp target exit data map(from: X) map(delete: Tm)
    call check('ztrsm '//side//uplo//tr, rel(X(1:M,:), Xh(1:M,:)))
  end subroutine trsm

  subroutine unmqr(side, tr)
    character(len=1), intent(in) :: side, tr
    complex(8), allocatable :: X(:,:), Xd(:,:), wl(:)
    integer :: lwl
    X = B; Xd = Bd
    call la_zunmqr(side, tr, NQ, NQ, NQ, A, NQ, tau, X, NQ, wq, -1, info); lwl = int(dble(wq(1))); allocate(wl(lwl))
    call la_zunmqr(side, tr, NQ, NQ, NQ, A, NQ, tau, X, NQ, wl, lwl, info)
    !$omp target enter data map(to: Xd)
    call la_zunmqr(side, tr, NQ, NQ, NQ, Ad, NQ, taud, Xd, NQ, wl, lwl, info)
    !$omp target exit data map(from: Xd)
    call check('zunmqr '//side//tr, rel(Xd, X))
  end subroutine unmqr
end program cuda_linalg

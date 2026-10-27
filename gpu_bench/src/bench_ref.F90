! HOST reference implementations (plain Fortran + host BLAS/LAPACK), never offloaded.
! They define "correct" for the device results.  Same in both backends.
module bench_ref
  implicit none
  external :: zgemm, zgeru, zgeqp3, zgeqrf, zungqr, dsyev   ! host LAPACK/BLAS, plain external
contains

  subroutine ref_cb(G, N, nb, bi, bj, c, s)
    integer, intent(in) :: N, nb
    complex(8), intent(inout) :: G(N,N)
    integer, intent(in) :: bi(nb,4), bj(nb,4)
    real(8), intent(in) :: c, s
    integer :: f, b, i, j, k
    complex(8) :: gi, gj
    do f = 4, 1, -1
      do b = 1, nb
        i = bi(b,f); j = bj(b,f)
        do k = 1, N
          gi = G(i,k); gj = G(j,k)
          G(i,k) = c*gi + s*gj
          G(j,k) = s*gi + c*gj
        end do
      end do
    end do
    do f = 4, 1, -1
      do b = 1, nb
        i = bi(b,f); j = bj(b,f)
        do k = 1, N
          gi = G(k,i); gj = G(k,j)
          G(k,i) = c*gi - s*gj
          G(k,j) = -s*gi + c*gj
        end do
      end do
    end do
  end subroutine

  ! sequential per-site wrap, as ALF does
  subroutine ref_opv(G, N, ea)
    integer, intent(in) :: N
    complex(8), intent(inout) :: G(N,N)
    real(8), intent(in) :: ea(N)
    integer :: i
    do i = 1, N
      G(i,:) = G(i,:) * ea(i)
      G(:,i) = G(:,i) / ea(i)
    end do
  end subroutine

  ! E = exp(sgn*dt*T) for the periodic L x L square-lattice hopping T (t=1) via dsyev
  subroutine ref_hop_exp(L, dt, E, Einv)
    integer, intent(in) :: L
    real(8), intent(in) :: dt
    real(8), intent(out) :: E(L*L,L*L), Einv(L*L,L*L)
    real(8), allocatable :: T(:,:), w(:), work(:)
    integer :: N, x, y, i, info, j, k
    real(8) :: s1, s2
    N = L*L
    allocate(T(N,N), w(N), work(max(1,3*N)))
    T = 0d0
    do y = 0, L-1
      do x = 0, L-1
        i = 1 + x + L*y
        T(i, 1 + mod(x+1,L) + L*y) = 1d0
        T(1 + mod(x+1,L) + L*y, i) = 1d0
        T(i, 1 + x + L*mod(y+1,L)) = 1d0
        T(1 + x + L*mod(y+1,L), i) = 1d0
      end do
    end do
    call dsyev('V', 'U', N, T, N, w, work, 3*N, info)
    if (info /= 0) error stop 'dsyev failed'
    do j = 1, N
      do i = 1, N
        s1 = 0d0; s2 = 0d0
        do k = 1, N
          s1 = s1 + T(i,k)*exp(dt*w(k))*T(j,k)
          s2 = s2 + T(i,k)*exp(-dt*w(k))*T(j,k)
        end do
        E(i,j) = s1; Einv(i,j) = s2
      end do
    end do
  end subroutine

  ! ZGEQP3 + ZUNGQR on a copy: Q, R (upper triangle of F), jpvt (1-based, A(:,jpvt) = Q R)
  subroutine ref_qp3(N, A, F, Q, tau, jpvt)
    integer, intent(in) :: N
    complex(8), intent(in) :: A(N,N)
    complex(8), intent(out) :: F(N,N), Q(N,N), tau(N)
    integer, intent(out) :: jpvt(N)
    complex(8), allocatable, save :: work(:)
    real(8), allocatable, save :: rwork(:)
    integer, save :: lw = 0, nsave = 0
    complex(8) :: w1(1)
    integer :: info, lw2
    if (N /= nsave) then
      if (allocated(work)) deallocate(work, rwork)
      nsave = N
      allocate(rwork(2*N))
      call zgeqp3(N, N, F, N, jpvt, tau, w1, -1, rwork, info); lw = int(real(w1(1)))
      call zungqr(N, N, N, F, N, tau, w1, -1, info); lw2 = int(real(w1(1)))
      lw = max(lw, lw2, 2*N)
      allocate(work(lw))
    end if
    F = A
    jpvt = 0
    call zgeqp3(N, N, F, N, jpvt, tau, work, lw, rwork, info)
    Q = F
    call zungqr(N, N, N, Q, N, tau, work, lw, info)
  end subroutine

  ! stable descending column-norm order (rank counting, same rule as the device kernel)
  subroutine ref_sortperm(A, N, perm)
    integer, intent(in) :: N
    complex(8), intent(in) :: A(N,N)
    integer, intent(out) :: perm(N)
    real(8) :: nrm(N)
    integer :: j, k, r
    do j = 1, N
      nrm(j) = sum(abs(A(:,j))**2)
    end do
    do j = 1, N
      r = 1
      do k = 1, N
        if (nrm(k) > nrm(j) .or. (nrm(k) == nrm(j) .and. k < j)) r = r + 1
      end do
      perm(r) = j
    end do
  end subroutine

  ! pre-pivot + ZGEQRF + ZUNGQR on copies (host path B).  F = factored, Q = formed Q.
  subroutine ref_prepiv_qr(N, A, F, Q, tau, perm)
    integer, intent(in) :: N
    complex(8), intent(in) :: A(N,N)
    complex(8), intent(out) :: F(N,N), Q(N,N), tau(N)
    integer, intent(out) :: perm(N)
    complex(8), allocatable, save :: work(:)
    integer, save :: lw = 0, nsave = 0
    complex(8) :: w1(1)
    integer :: info, lw2, j
    if (N /= nsave) then
      if (allocated(work)) deallocate(work)
      nsave = N
      call zgeqrf(N, N, F, N, tau, w1, -1, info); lw = int(real(w1(1)))
      call zungqr(N, N, N, F, N, tau, w1, -1, info); lw2 = int(real(w1(1)))
      lw = max(lw, lw2, 2*N)
      allocate(work(lw))
    end if
    call ref_sortperm(A, N, perm)
    do j = 1, N
      F(:,j) = A(:,perm(j))
    end do
    call zgeqrf(N, N, F, N, tau, work, lw, info)
    Q = F
    call zungqr(N, N, N, Q, N, tau, work, lw, info)
  end subroutine

  ! Same pre-pivot, but UNBLOCKED ZGEQR2: only the diagonal |R_ii| is returned.  Used to show how much
  ! |R_ii| differ between two correct Householder implementations (blocked vs unblocked).
  subroutine ref_prepiv_rdiag_unblocked(N, A, rdiag)
    integer, intent(in) :: N
    complex(8), intent(in) :: A(N,N)
    real(8), intent(out) :: rdiag(N)
    complex(8), allocatable :: F(:,:), tau(:), work(:)
    integer, allocatable :: perm(:)
    integer :: j, info
    external :: zgeqr2
    allocate(F(N,N), tau(N), work(N), perm(N))
    call ref_sortperm(A, N, perm)
    do j = 1, N
      F(:,j) = A(:,perm(j))
    end do
    call zgeqr2(N, N, F, N, tau, work, info)
    do j = 1, N
      rdiag(j) = abs(F(j,j))
    end do
  end subroutine

  ! A = Q1 diag(d) R1; Q1 unitary (QR of a random matrix), R1 random upper triangular with
  ! diagonal in [1,2], d log-spaced from 10^e1 to 10^e2.
  subroutine ref_graded(N, e1, e2, A)
    integer, intent(in) :: N
    real(8), intent(in) :: e1, e2
    complex(8), intent(out) :: A(N,N)
    complex(8), allocatable :: Q1(:,:), M(:,:), tau(:), work(:)
    real(8), allocatable :: x(:,:), y(:,:)
    integer :: i, j, info, lw
    complex(8) :: one, zero
    one = 1d0; zero = 0d0
    allocate(Q1(N,N), M(N,N), tau(N), work(64*N), x(N,N), y(N,N))
    lw = 64*N
    call random_number(x); call random_number(y)
    Q1 = cmplx(x - 0.5d0, y - 0.5d0, kind=8)
    call zgeqrf(N, N, Q1, N, tau, work, lw, info)
    call zungqr(N, N, N, Q1, N, tau, work, lw, info)
    call random_number(x); call random_number(y)
    do j = 1, N
      do i = 1, N
        if (i > j) then
          M(i,j) = 0d0
        else if (i == j) then
          M(i,j) = 1d0 + x(i,j)
        else
          M(i,j) = cmplx(x(i,j) - 0.5d0, y(i,j) - 0.5d0, kind=8) / sqrt(dble(N))
        end if
      end do
    end do
    do i = 1, N
      M(i,:) = M(i,:) * 10d0**(e1 + (e2 - e1) * dble(i-1) / dble(N-1))
    end do
    call zgemm('N', 'N', N, N, N, one, Q1, N, M, N, zero, A, N)
  end subroutine

  ! random unitary matrix (Q factor of a random complex matrix)
  subroutine ref_unitary(N, Q1)
    integer, intent(in) :: N
    complex(8), intent(out) :: Q1(N,N)
    complex(8), allocatable :: tau(:), work(:)
    real(8), allocatable :: x(:,:), y(:,:)
    integer :: info, lw
    allocate(tau(N), work(64*N), x(N,N), y(N,N))
    lw = 64*N
    call random_number(x); call random_number(y)
    Q1 = cmplx(x - 0.5d0, y - 0.5d0, kind=8)
    call zgeqrf(N, N, Q1, N, tau, work, lw, info)
    call zungqr(N, N, N, Q1, N, tau, work, lw, info)
  end subroutine

  ! d(j) = 10^(e1 + (e2-e1) r(j)/(N-1)) with r a random permutation of 0..N-1 (log-spaced, random order)
  subroutine ref_randgrade(N, e1, e2, d)
    integer, intent(in) :: N
    real(8), intent(in) :: e1, e2
    real(8), intent(out) :: d(N)
    real(8) :: key(N)
    integer :: j, k, r
    call random_number(key)
    do j = 1, N
      r = 0
      do k = 1, N
        if (key(k) < key(j) .or. (key(k) == key(j) .and. k < j)) r = r + 1
      end do
      d(j) = 10d0**(e1 + (e2 - e1) * dble(r) / dble(N-1))
    end do
  end subroutine

  ! mat5: column-graded A = M diag(d), M = unitary * diag(s), s in [1,3] (cond(M) <= 3), d log-spaced
  ! 1e12 .. 1e-12 in random column order.  This is the structure of ALF's decompose (B U D).
  subroutine ref_colgraded(N, A)
    integer, intent(in) :: N
    complex(8), intent(out) :: A(N,N)
    real(8) :: d(N), s(N)
    integer :: j
    call ref_unitary(N, A)
    call random_number(s); s = 1d0 + 2d0*s
    call ref_randgrade(N, 12d0, -12d0, d)
    do j = 1, N
      A(:,j) = A(:,j) * (s(j) * d(j))
    end do
  end subroutine

  ! mat6: two-sided graded A = diag(d1) M diag(d2), d1, d2 log-spaced 1e6 .. 1e-6 in random order
  ! (entries span 1e12 .. 1e-12); structure of CGR's U_R^H U_L + D_R (V_R V_L) D_L.
  subroutine ref_twograded(N, A)
    integer, intent(in) :: N
    complex(8), intent(out) :: A(N,N)
    real(8) :: d1(N), d2(N), s(N)
    integer :: i, j
    call ref_unitary(N, A)
    call random_number(s); s = 1d0 + 2d0*s
    call ref_randgrade(N, 6d0, -6d0, d1)
    call ref_randgrade(N, 6d0, -6d0, d2)
    do j = 1, N
      do i = 1, N
        A(i,j) = d1(i) * A(i,j) * (s(j) * d2(j))
      end do
    end do
  end subroutine

end module bench_ref

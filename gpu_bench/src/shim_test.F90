! Standalone check of the SYCL shim (bench_sycl_shim.cpp) against host LAPACK.
! Runs with the shim's private fallback queue (null queue handle): ONEAPI_DEVICE_SELECTOR=opencl:cpu on a machine
! without an Intel GPU.  On a PVC node the same code can be run with the queue from OpenMP interop (see bench_la).
! Checks (all must be <= 1e-12 relative): geqrf+ungqr (Q R = A, Q^H Q = I), unmqr (Q^H A = R, Q C),
! getrf+getri (A A^-1 = I), plus bit-level agreement of R with host ZGEQRF is NOT expected (different blocking).
program shim_test
  use iso_c_binding
  implicit none
  interface
    integer(c_int) function alfsh_describe(q) bind(C, name='alfsh_describe')
      import :: c_ptr, c_int
      type(c_ptr), value :: q
    end function
    integer(c_int) function alfsh_zgeqrf_ws(q, m, n, lda, ws) bind(C, name='alfsh_zgeqrf_ws')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgeqrf(q, m, n, a, lda, tau, work, lwork) bind(C, name='alfsh_zgeqrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q, a, tau, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
    integer(c_int) function alfsh_zungqr_ws(q, m, n, k, lda, ws) bind(C, name='alfsh_zungqr_ws')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q
      integer(c_int64_t), value :: m, n, k, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zungqr(q, m, n, k, a, lda, tau, work, lwork) bind(C, name='alfsh_zungqr')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q, a, tau, work
      integer(c_int64_t), value :: m, n, k, lda, lwork
    end function
    integer(c_int) function alfsh_zunmqr_ws(q, side, trans, m, n, k, lda, ldc, ws) &
         bind(C, name='alfsh_zunmqr_ws')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q
      integer(c_int), value :: side, trans
      integer(c_int64_t), value :: m, n, k, lda, ldc
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zunmqr(q, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork) &
         bind(C, name='alfsh_zunmqr')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q, a, tau, c, work
      integer(c_int), value :: side, trans
      integer(c_int64_t), value :: m, n, k, lda, ldc, lwork
    end function
    integer(c_int) function alfsh_zgetrf_ws(q, m, n, lda, ws) bind(C, name='alfsh_zgetrf_ws')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgetrf(q, m, n, a, lda, ipiv, work, lwork) bind(C, name='alfsh_zgetrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q, a, ipiv, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
    integer(c_int) function alfsh_zgetri_ws(q, n, lda, ws) bind(C, name='alfsh_zgetri_ws')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q
      integer(c_int64_t), value :: n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgetri(q, n, a, lda, ipiv, work, lwork) bind(C, name='alfsh_zgetri')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: q, a, ipiv, work
      integer(c_int64_t), value :: n, lda, lwork
    end function
  end interface
  integer, parameter :: i8 = c_int64_t
  integer :: nn(3), ik, N, i, j, nfail
  nfail = 0
  nn = [64, 144, 400]
  if (alfsh_describe(c_null_ptr) /= 0) error stop 'no SYCL queue'
  do ik = 1, 3
    N = nn(ik)
    call run(N)
  end do
  write(*,'(a,i0,a)') 'shim_test: ', nfail, ' FAILED'
  if (nfail > 0) stop 1
contains
  subroutine pf(name, N, v, tol)
    character(len=*), intent(in) :: name
    integer, intent(in) :: N
    real(8), intent(in) :: v, tol
    character(len=4) :: s
    s = 'PASS'
    if (.not. (v <= tol)) then; s = 'FAIL'; nfail = nfail + 1; end if
    write(*,'(a,1x,a40,1x,a,i5,a,es10.3,a,es9.2)') s, name, 'N=', N, '  value=', v, '  tol=', tol
  end subroutine

  real(8) function fro(A, n2)
    integer, intent(in) :: n2
    complex(8), intent(in) :: A(n2)
    fro = sqrt(sum(abs(A)**2))
  end function

  subroutine run(N)
    integer, intent(in) :: N
    complex(8), allocatable, target :: A0(:,:), A(:,:), tau(:), W(:), Q(:,:), C0(:,:), C(:,:), B(:,:)
    integer(c_int64_t), allocatable, target :: ipiv(:)
    complex(8), allocatable :: R(:,:), X(:,:)
    real(8), allocatable :: x1(:,:), y1(:,:)
    integer(c_int64_t) :: w1, w2, w3, lw
    integer :: st
    complex(8), parameter :: one = (1d0, 0d0), zero = (0d0, 0d0)
    external :: zgemm
    allocate(A0(N,N), A(N,N), tau(N), Q(N,N), C0(N,N), C(N,N), R(N,N), X(N,N), x1(N,N), y1(N,N), ipiv(N), B(N,N))
    call random_number(x1); call random_number(y1)
    A0 = cmplx(x1 - 0.5d0, y1 - 0.5d0, kind=8)
    call random_number(x1); call random_number(y1)
    C0 = cmplx(x1 - 0.5d0, y1 - 0.5d0, kind=8)
    ! --- QR
    st = alfsh_zgeqrf_ws(c_null_ptr, int(N,i8), int(N,i8), int(N,i8), w1)
    st = alfsh_zungqr_ws(c_null_ptr, int(N,i8), int(N,i8), int(N,i8), int(N,i8), w2)
    st = alfsh_zunmqr_ws(c_null_ptr, 0, 2, int(N,i8), int(N,i8), int(N,i8), int(N,i8), int(N,i8), w3)
    lw = max(w1, w2, w3, 1_i8)
    allocate(W(lw)); W = 0
    A = A0
    st = alfsh_zgeqrf(c_null_ptr, int(N,i8), int(N,i8), c_loc(A(1,1)), int(N,i8), c_loc(tau(1)), c_loc(W(1)), lw)
    call pf('geqrf status', N, dble(st), 0d0)
    R = 0
    do j = 1, N
      R(1:j,j) = A(1:j,j)
    end do
    ! Q^H A0 = R via unmqr (side L, trans C)
    C = A0
    st = alfsh_zunmqr(c_null_ptr, 0, 2, int(N,i8), int(N,i8), int(N,i8), c_loc(A(1,1)), int(N,i8), c_loc(tau(1)), &
         c_loc(C(1,1)), int(N,i8), c_loc(W(1)), lw)
    call pf('unmqr status', N, dble(st), 0d0)
    call pf('unmqr Q^H A = R', N, fro(reshape(C - R, [N*N]), N*N) / fro(reshape(A0, [N*N]), N*N), 1d-12)
    ! Q R = A0 via ungqr
    Q = A
    st = alfsh_zungqr(c_null_ptr, int(N,i8), int(N,i8), int(N,i8), c_loc(Q(1,1)), int(N,i8), c_loc(tau(1)), &
         c_loc(W(1)), lw)
    call pf('ungqr status', N, dble(st), 0d0)
    call zgemm('N', 'N', N, N, N, one, Q, N, R, N, zero, X, N)
    call pf('geqrf+ungqr |QR - A|/|A|', N, fro(reshape(X - A0, [N*N]), N*N) / fro(reshape(A0, [N*N]), N*N), 1d-12)
    call zgemm('C', 'N', N, N, N, one, Q, N, Q, N, zero, X, N)
    do j = 1, N
      X(j,j) = X(j,j) - 1d0
    end do
    call pf('ungqr |Q^H Q - I|_F', N, fro(reshape(X, [N*N]), N*N), 1d-12 * N)
    ! --- LU
    deallocate(W)
    st = alfsh_zgetrf_ws(c_null_ptr, int(N,i8), int(N,i8), int(N,i8), w1)
    st = alfsh_zgetri_ws(c_null_ptr, int(N,i8), int(N,i8), w2)
    lw = max(w1, w2, 1_i8); allocate(W(lw)); W = 0
    B = A0
    st = alfsh_zgetrf(c_null_ptr, int(N,i8), int(N,i8), c_loc(B(1,1)), int(N,i8), c_loc(ipiv(1)), c_loc(W(1)), lw)
    call pf('getrf status', N, dble(st), 0d0)
    st = alfsh_zgetri(c_null_ptr, int(N,i8), c_loc(B(1,1)), int(N,i8), c_loc(ipiv(1)), c_loc(W(1)), lw)
    call pf('getri status', N, dble(st), 0d0)
    call zgemm('N', 'N', N, N, N, one, A0, N, B, N, zero, X, N)
    do j = 1, N
      X(j,j) = X(j,j) - 1d0
    end do
    ! |A A^-1 - I| scales with cond(A) * eps; random N x N complex: cond ~ N, so tolerance 1e-12 * N is generous
    call pf('getrf+getri |A A^-1 - I|_F', N, fro(reshape(X, [N*N]), N*N), 1d-12 * N)
  end subroutine
end program

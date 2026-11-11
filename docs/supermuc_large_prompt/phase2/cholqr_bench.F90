! Prototype benchmark: CholeskyQR2 on the device (herk/potrf/trsm via oneMKL OpenMP dispatch) vs zgeqrf+zungqr
include "mkl_blas_omp_offload.f90"
include "mkl_lapack_omp_offload.f90"
program cholqr_bench
  use omp_lib
  use onemkl_blas_omp_offload_lp64, only: zgemm, zherk, ztrsm
  use onemkl_lapack_omp_offload_lp64, only: zpotrf, zgeqrf
  implicit none
  integer, parameter :: n = 1024
  complex(8), allocatable :: X0(:,:), X(:,:), R1(:,:), R2(:,:), Gm(:,:), tau(:), wk(:)
  real(8) :: t0, t1, err, cnd
  integer :: i, j, info, rep, lw
  complex(8) :: one, zero
  real(8) :: r(2), sc
  integer :: kp
  one = (1.d0,0.d0); zero = (0.d0,0.d0)
  allocate(X0(n,n), X(n,n), R1(n,n), R2(n,n), Gm(n,n))
  call srand(7)
  do cnd = 1, 1
  end do
  do kp = 0, 5   ! conditioning of the column-normalised matrix: 10^kp ... via mixing columns
    ! X0 = Qrand * diag(s) * Qrand2 with s from 1 to 10^-(2*kp)  (cond of the scaled matrix ~ 10^(2 kp) )
    call make(X0, n, 2.d0*kp)
    ! column-normalise as ALF's pre-pivot step would see it
    do j = 1, n
      sc = sqrt(sum(abs(X0(:,j))**2)); X0(:,j) = X0(:,j)/sc
    end do
    write(*,'(A,I0,A)') '--- cond target 1e', 2*kp, ' (scaled columns)'
    !$omp target enter data map(to: X0) map(alloc: X, R1, R2, Gm)
    ! warm-up + timed CholQR2
    do rep = 0, 5
      if (rep == 1) t0 = omp_get_wtime()
      call chol2(X0, X, R1, R2, Gm, n, info)
    end do
    t1 = omp_get_wtime()
    !$omp target update from(X)
    err = orth(X, n)
    if (kp == 0) call pieces()
    write(*,'(A,F8.2,A,ES10.2,A,I0)') ' CholQR2 ms: ', 1000*(t1-t0)/5, '   |Q^H Q - 1|max: ', err, '  info=', info
    !$omp target exit data map(delete: X0, X, R1, R2, Gm)
  end do
contains
  subroutine pieces()
    integer :: r, ii, jj; real(8) :: ta, tb, tc, td; complex(8), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0)
    ta = 0; tb = 0; tc = 0; td = 0
    do r = 1, 5
      t0 = omp_get_wtime()
      !$omp dispatch
      call zherk('U', 'C', n, n, 1.d0, X0, n, 0.d0, Gm, n)
      t1 = omp_get_wtime(); ta = ta + t1 - t0
      t0 = omp_get_wtime()
      !$omp dispatch
      call zpotrf('U', n, Gm, n, info)
      t1 = omp_get_wtime(); tb = tb + t1 - t0
      t0 = omp_get_wtime()
      !$omp dispatch
      call ztrsm('R', 'U', 'N', 'N', n, n, one, Gm, n, X, n)
      t1 = omp_get_wtime(); tc = tc + t1 - t0
      t0 = omp_get_wtime()
      !$omp dispatch
      call zgemm('N','N',n,n,n,one,Gm,n,R1,n,zero,R2,n)
      t1 = omp_get_wtime(); td = td + t1 - t0
    end do
    write(*,'(A,4F8.2)') ' pieces ms (herk, potrf, trsm, gemm): ', 200*ta, 200*tb, 200*tc, 200*td
  end subroutine
  subroutine make(A, n, lg)
    complex(8) :: A(n,n); integer :: n; real(8) :: lg
    complex(8), allocatable :: U(:,:), V(:,:), tmp(:,:)
    real(8) :: s, rr(2); integer :: i, j
    allocate(U(n,n), V(n,n), tmp(n,n))
    ! random matrices -> orthogonalised by modified Gram-Schmidt on the host (n=1024: ok, once)
    call randorth(U, n); call randorth(V, n)
    do j = 1, n
      s = 10.d0**(-lg*real(j-1,8)/real(n-1,8))
      tmp(:,j) = s * V(j,:)      ! (diag(s) V^H)^T pattern; exact form is irrelevant: singular values are s
    end do
    A = matmul(U, conjg(transpose(tmp)))
  end subroutine
  subroutine randorth(U, n)
    complex(8) :: U(n,n); integer :: n; integer :: i, j, p; real(8) :: a, b
    complex(8) :: c
    do j = 1, n
      do i = 1, n
        call random_number(a); call random_number(b); U(i,j) = cmplx(a-0.5d0, b-0.5d0, 8)
      end do
    end do
    do j = 1, n
      do p = 1, j-1
        c = dot_product(U(:,p), U(:,j)); U(:,j) = U(:,j) - c*U(:,p)
      end do
      do p = 1, j-1
        c = dot_product(U(:,p), U(:,j)); U(:,j) = U(:,j) - c*U(:,p)
      end do
      U(:,j) = U(:,j)/sqrt(real(dot_product(U(:,j),U(:,j)),8))
    end do
  end subroutine
  real(8) function orth(Q, n)
    complex(8) :: Q(n,n); integer :: n
    complex(8), allocatable :: P(:,:); integer :: i
    allocate(P(n,n)); P = matmul(conjg(transpose(Q)), Q)
    do i = 1, n; P(i,i) = P(i,i) - 1.d0; end do
    orth = maxval(abs(P))
  end function
  subroutine chol2(A, Q, R1, R2, Gm, n, info)
    complex(8) :: A(n,n), Q(n,n), R1(n,n), R2(n,n), Gm(n,n); integer :: n, info
    complex(8), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0)
    integer :: i, j
    ! pass 1: G = A^H A (use gemm: herk writes one triangle, potrf reads it)
    !$omp dispatch
    call zherk('U', 'C', n, n, 1.d0, A, n, 0.d0, Gm, n)
    !$omp dispatch
    call zpotrf('U', n, Gm, n, info)           ! Gm = R1
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, n
        Q(i,j) = A(i,j)
        R1(i,j) = Gm(i,j)
      end do
    end do
    !$omp dispatch
    call ztrsm('R', 'U', 'N', 'N', n, n, one, R1, n, Q, n)     ! Q = A R1^-1
    ! pass 2
    !$omp dispatch
    call zherk('U', 'C', n, n, 1.d0, Q, n, 0.d0, Gm, n)
    !$omp dispatch
    call zpotrf('U', n, Gm, n, info)
    !$omp dispatch
    call ztrsm('R', 'U', 'N', 'N', n, n, one, Gm, n, Q, n)
    ! R = R2 R1 (triangular product) would be one more gemm; timed separately:
    !$omp dispatch
    call zgemm('N','N',n,n,n,one,Gm,n,R1,n,zero,R2,n)
  end subroutine
end program

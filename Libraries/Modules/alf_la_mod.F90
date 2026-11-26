!  Copyright (C) 2016 - 2026 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution of significant scientific interest.
!       Therefore, if you publish results produced with the ALF project, you are required to cite
!       the ALF project; see https://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!> @file alf_la_mod.F90
!> @brief Plain-array linear-algebra layer of the GPU port (ADR-002 section 1), HOST backend only.
!>
!> @details Every routine calls BLAS/LAPACK directly on contiguous arrays (the exact calls ALF makes on the CPU), so
!> that the GPU modules written against this interface run unchanged on the host (GPU=host) and, in later steps, get
!> a device backend behind the same signatures. No existing module uses it; it changes nothing in the CPU path.
!>
!> Interface rules
!>   * explicit-shape/assumed-size arrays + leading dimensions in LAPACK order, complex(8) throughout
!>     (la_gemm_cr: B real). No descriptors, no derived types. Index arrays are integer(4).
!>     LA_N, LA_T, LA_C have the values of cublasOperation_t (0, 1, 2).
!>   * Scratch memory (LAPACK work) is owned by the module and grown on demand (la_reserve_work pre-grows it),
!>     never allocated per call in steady state. The *_lwork functions return the optimal work size in complex(8)
!>     elements (LAPACK work query, info checked).
!>   * NOT thread-safe: one set of scratch buffers per process; call from one thread.
!>   * info is checked at once; any failure ends the run through Terminate_on_error. Tolerated: getrf info>0
!>     (exactly singular U; ALF ignores it).
!>   * la_sync and la_finalize exist for the device backends; on the host la_sync does nothing.
!>
!> Routine -> LAPACK/BLAS call
!>   la_gemm ZGEMM; la_gemm_cr DGEMM on the interleaved real view; la_gemm_batched ZGEMM in a loop; la_geru ZGERU;
!>   la_trsm ZTRSM; la_geqrf/_lwork ZGEQRF; la_ungqr/_lwork ZUNGQR; la_ormqr/_lwork ZUNMQR; la_getrf ZGETRF;
!>   la_getrs ZGETRS; la_qrpp column-norm pre-pivot + ZGEQRF; la_geqp3_host ZGEQP3.
!>   Not provided: la_gemm_rc (real x complex), la_getri (cold INV routines only).

#define LA_FAIL(msg) call la_fail(msg, __LINE__)

module alf_la_mod
  use iso_c_binding, only: c_int, c_loc, c_f_pointer
  use iso_fortran_env, only: error_unit
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC

  implicit none
  private
  integer, parameter, public :: LA_N = 0, LA_T = 1, LA_C = 2    ! same values as cublasOperation_t
  public :: la_init, la_finalize, la_sync, la_backend_name, la_reserve_work
  public :: la_gemm, la_gemm_cr, la_gemm_batched, la_geru, la_trsm
  public :: la_geqrf_lwork, la_geqrf, la_ungqr_lwork, la_ungqr, la_ormqr_lwork, la_ormqr
  public :: la_getrf_lwork, la_getrf, la_getrs
  public :: la_qrpp, la_geqp3_host

  logical :: inited = .false.
  complex(8), allocatable :: hws(:)       ! host LAPACK work (grown on demand)
  real(8), allocatable :: hrws(:)         ! host real work (ZGEQP3 rwork, norms)
  complex(8), allocatable :: hscr(:)      ! host scratch for the column permutation (HOST la_qrpp)

  external :: zgemm, dgemm, zgeru, ztrsm, zgeqrf, zungqr, zunmqr, zgetrf, zgetrs, zgeqp3
  real(8), external :: dznrm2

contains

  !> @brief Name of the active backend (always 'host' here).
  function la_backend_name() result(s)
    character(len=16) :: s
    s = 'host'
  end function la_backend_name

  !> Print msg and end the run through Terminate_on_error.
  subroutine la_fail(msg, line)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: line
    write(error_unit,*) 'alf_la: ', msg
    call Terminate_on_error(ERROR_GENERIC, 'alf_la_mod.F90', line)
  end subroutine la_fail


  !> @brief Initialise the layer (host: nothing to set up).
  subroutine la_init()
    inited = .true.
  end subroutine la_init

  !> @brief Release the work memory.
  subroutine la_finalize()
    if (.not. inited) return
    if (allocated(hws)) deallocate(hws)
    if (allocated(hrws)) deallocate(hrws)
    if (allocated(hscr)) deallocate(hscr)
    inited = .false.
  end subroutine la_finalize

  !> Block until all enqueued device work is finished (nothing to do on the host).
  subroutine la_sync()
  end subroutine la_sync

  !> Pre-grow the work memory (complex(8) elements).
  subroutine la_reserve_work(nelem)
    integer, intent(in) :: nelem
    call ensure_hws(nelem)
  end subroutine la_reserve_work

  subroutine ensure_hws(n)
    integer, intent(in) :: n
    if (allocated(hws)) then
      if (size(hws) >= n) return
      deallocate(hws)
    end if
    allocate(hws(max(n,1)))
  end subroutine ensure_hws

  subroutine ensure_hrws(n)
    integer, intent(in) :: n
    if (allocated(hrws)) then
      if (size(hrws) >= n) return
      deallocate(hrws)
    end if
    allocate(hrws(max(n,1)))
  end subroutine ensure_hrws

  character function opc(op)
    integer, intent(in) :: op
    select case (op)
    case (LA_N); opc = 'N'
    case (LA_T); opc = 'T'
    case default; opc = 'C'
    end select
  end function opc


  ! ---------------------------------------------------------------------------------------------
  ! BLAS
  ! ---------------------------------------------------------------------------------------------

  !> C = alpha op(A) op(B) + beta C
  subroutine la_gemm(opa, opb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(lda,*), B(ldb,*), C(ldc,*)
    if (m <= 0 .or. n <= 0) return
    call zgemm(opc(opa), opc(opb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  end subroutine la_gemm

  !> C(complex m x n) = A(complex m x k) * B(REAL k x n): one DGEMM on the interleaved real view of A and C
  !> (real 2m x k and 2m x n). lda, ldc in complex elements, ldb in real elements.
  subroutine la_gemm_cr(m, n, k, A, lda, B, ldb, C, ldc)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(8), target :: A(lda,*), C(ldc,*)
    real(8), target :: B(ldb,*)
    real(8), parameter :: one = 1.d0, zero = 0.d0
    real(8), pointer, contiguous :: Ar(:,:), Cr(:,:)
    if (m <= 0 .or. n <= 0) return
    ! real views of the complex arrays (contiguous pointers: no copy-in/out when passed on)
    call c_f_pointer(c_loc(A(1,1)), Ar, [2*lda, k])
    call c_f_pointer(c_loc(C(1,1)), Cr, [2*ldc, n])
    call dgemm('N', 'N', 2*m, n, k, one, Ar, 2*lda, B, ldb, zero, Cr, 2*ldc)
  end subroutine la_gemm_cr

  !> nb independent complex GEMMs with constant strides sA, sB, sC (in complex elements)
  subroutine la_gemm_batched(opa, opb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc, nb, sA, sB, sC
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(*), B(*), C(*)
    integer :: ib
    if (m <= 0 .or. n <= 0 .or. nb <= 0) return
    do ib = 0, nb - 1
      call zgemm(opc(opa), opc(opb), m, n, k, alpha, A(1+ib*sA), lda, B(1+ib*sB), ldb, beta, C(1+ib*sC), ldc)
    end do
  end subroutine la_gemm_batched

  !> A += alpha x y^T (no conjugation)
  subroutine la_geru(m, n, alpha, x, incx, y, incy, A, lda)
    integer, intent(in) :: m, n, incx, incy, lda
    complex(8), intent(in) :: alpha
    complex(8), target :: x(*), y(*), A(lda,*)
    if (m <= 0 .or. n <= 0) return
    call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
  end subroutine la_geru

  !> B := alpha op(A)^-1 B (side 'L') or B := alpha B op(A)^-1 (side 'R'); uplo 'U'/'L', diag 'N'/'U'.
  subroutine la_trsm(side, uplo, opa, diag, m, n, alpha, A, lda, B, ldb)
    character(len=1), intent(in) :: side, uplo, diag
    integer, intent(in) :: opa, m, n, lda, ldb
    complex(8), intent(in) :: alpha
    complex(8), target :: A(lda,*), B(ldb,*)
    if (m <= 0 .or. n <= 0) return
    call ztrsm(side, uplo, opc(opa), diag, m, n, alpha, A, lda, B, ldb)
  end subroutine la_trsm

  ! ---------------------------------------------------------------------------------------------
  ! LAPACK: QR family.  Work memory is internal; the *_lwork functions return the optimal size in
  ! complex(8) elements (for la_reserve_work).
  ! ---------------------------------------------------------------------------------------------
  !> Optimal work size of la_geqrf in complex(8) elements.
  function la_geqrf_lwork(m, n, A, lda) result(lw)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer :: lw
    complex(8) :: w(1), t1(1)
    integer :: info
    call zgeqrf(m, n, A, lda, t1, w, -1, info)
    if (info /= 0) LA_FAIL('zgeqrf work query failed')
    lw = int(dble(w(1)))
  end function la_geqrf_lwork

  !> QR factorisation A = Q R (ZGEQRF), reflectors in A and tau.
  subroutine la_geqrf(m, n, A, lda, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
    integer :: info
    if (m <= 0 .or. n <= 0) return
    lw = la_geqrf_lwork(m, n, A, lda)
    call ensure_hws(lw)
    call zgeqrf(m, n, A, lda, tau, hws, lw, info)
    if (info /= 0) LA_FAIL('zgeqrf failed')
  end subroutine la_geqrf

  !> Optimal work size of la_ungqr in complex(8) elements.
  function la_ungqr_lwork(m, n, k, A, lda, tau) result(lw)
    integer, intent(in) :: m, n, k, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
    complex(8) :: w(1)
    integer :: info
    call zungqr(m, n, k, A, lda, tau, w, -1, info)
    if (info /= 0) LA_FAIL('zungqr work query failed')
    lw = int(dble(w(1)))
  end function la_ungqr_lwork

  !> Overwrite the first n columns of the geqrf factors in A (m x n, k reflectors) by Q (m x n).
  subroutine la_ungqr(m, n, k, A, lda, tau)
    integer, intent(in) :: m, n, k, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
    integer :: info
    if (m <= 0 .or. n <= 0) return
    lw = la_ungqr_lwork(m, n, k, A, lda, tau)
    call ensure_hws(lw)
    call zungqr(m, n, k, A, lda, tau, hws, lw, info)
    if (info /= 0) LA_FAIL('zungqr failed')
  end subroutine la_ungqr

  !> Optimal work size of la_ormqr in complex(8) elements. See la_ormqr for side/trans.
  !> side 'L': C := op(Q) C ; side 'R': C := C op(Q); trans LA_N or LA_C; Q from geqrf stored in A (k reflectors).
  function la_ormqr_lwork(side, trans, m, n, k, A, lda, tau, C, ldc) result(lw)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: lw
    complex(8) :: w(1)
    integer :: info
    call zunmqr(side, opc(trans), m, n, k, A, lda, tau, C, ldc, w, -1, info)
    if (info /= 0) LA_FAIL('zunmqr work query failed')
    lw = int(dble(w(1)))
  end function la_ormqr_lwork

  !> Apply Q of la_geqrf to C (ZUNMQR); arguments as la_ormqr_lwork.
  subroutine la_ormqr(side, trans, m, n, k, A, lda, tau, C, ldc)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: lw
    integer :: info
    if (m <= 0 .or. n <= 0 .or. k <= 0) return
    lw = la_ormqr_lwork(side, trans, m, n, k, A, lda, tau, C, ldc)
    call ensure_hws(lw)
    call zunmqr(side, opc(trans), m, n, k, A, lda, tau, C, ldc, hws, lw, info)
    if (info /= 0) LA_FAIL('zunmqr failed')
  end subroutine la_ormqr

  ! ---------------------------------------------------------------------------------------------
  ! LAPACK: LU.  ipiv is integer(4), 1-based (LAPACK convention).
  ! ---------------------------------------------------------------------------------------------
  !> Work size of la_getrf: always 0.
  function la_getrf_lwork(m, n, A, lda) result(lw)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer :: lw
    lw = 0     ! LAPACK ZGETRF needs no work array (nothing to query)
  end function la_getrf_lwork

  !> LU with partial pivoting. info > 0 (exactly singular U) is tolerated, as in ALF.
  subroutine la_getrf(m, n, A, lda, ipiv)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer(c_int), target :: ipiv(*)
    integer :: info
    if (m <= 0 .or. n <= 0) return
    call zgetrf(m, n, A, lda, ipiv, info)
    if (info < 0) LA_FAIL('zgetrf: illegal argument')
  end subroutine la_getrf

  !> Solve op(A) X = B with the factors of la_getrf. trans LA_N/LA_T/LA_C.
  subroutine la_getrs(trans, n, nrhs, A, lda, ipiv, B, ldb)
    integer, intent(in) :: trans, n, nrhs, lda, ldb
    complex(8), target :: A(lda,*), B(ldb,*)
    integer(c_int), target :: ipiv(*)
    integer :: info
    if (n <= 0 .or. nrhs <= 0) return
    call zgetrs(opc(trans), n, nrhs, A, lda, ipiv, B, ldb, info)
    if (info /= 0) LA_FAIL('zgetrs failed')
  end subroutine la_getrs

  ! ---------------------------------------------------------------------------------------------
  ! Pre-pivoted QR and host QRCP (STABLOG)
  ! ---------------------------------------------------------------------------------------------

  !> Stable descending insertion sort of perm(1:n) by nrm (identical to QR_PREPIVOT variant 1 in cgr1_mod).
  subroutine sort_desc(n, nrm, perm)
    integer, intent(in) :: n
    real(8), intent(in) :: nrm(*)
    integer(c_int), intent(out) :: perm(*)
    integer :: i, j
    integer(c_int) :: kk
    do j = 1, n
      perm(j) = j
    end do
    do i = 2, n
      kk = perm(i)
      j = i - 1
      do while (j >= 1)
        if (nrm(perm(j)) >= nrm(kk)) exit
        perm(j+1) = perm(j)
        j = j - 1
      end do
      perm(j+1) = kk
    end do
  end subroutine sort_desc

  !> true if x is NaN or +-Inf (bit test on the exponent: robust under -ffast-math / finite-math-only)
  logical function not_finite(x)
    real(8), intent(in) :: x
    integer(8) :: bits
    bits = transfer(x, 0_8)
    not_finite = (iand(bits, int(z'7FF0000000000000', 8)) == int(z'7FF0000000000000', 8))
  end function not_finite

  !> QR with a column-norm pre-pivot: ipvt(1:n) = permutation sorting the columns of A (m x n, requires m >= n) by
  !> descending 2-norm (stable); then A(:,ipvt) = Q R by ZGEQRF (no further pivoting). On return A holds the
  !> geqrf factors and tau(1:n) the reflectors; ipvt has ZGEQP3 semantics (column j of A*P is
  !> column ipvt(j) of A). Same pre-pivot arithmetic as the #else branch of QR_PREPIVOT (variant 1); the factorisation
  !> is ZGEQRF, whereas variant 1 calls QDRP_decompose = ZGEQP3 with all columns fixed (t_qrpp compares the two).
  !> The separation of D from R is left to the caller (kernel).
  !> NaN/Inf in A stop the run (Terminate_on_error) after the norm transfer.
  subroutine la_qrpp(m, n, A, lda, ipvt, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer(c_int), target :: ipvt(*)
    integer :: i, j
    if (m <= 0 .or. n <= 0) return
    if (m < n) LA_FAIL('la_qrpp requires m >= n')
    ! direct finite scan of A, not relying on DZNRM2 to propagate NaN/Inf
    do j = 1, n
      do i = 1, m
        if (not_finite(dble(A(i,j))) .or. not_finite(aimag(A(i,j)))) LA_FAIL('la_qrpp: NaN or Inf in the matrix')
      end do
    end do
    call ensure_hrws(n)
    if (allocated(hscr)) then
      if (size(hscr) < m*n) deallocate(hscr)
    end if
    if (.not. allocated(hscr)) allocate(hscr(m*n))
    do j = 1, n
      hrws(j) = dznrm2(m, A(1,j), 1)
    end do
    if (any_not_finite(n, hrws)) LA_FAIL('la_qrpp: NaN or Inf in the matrix')
    call sort_desc(n, hrws, ipvt)
    do j = 1, n
      do i = 1, m
        hscr(i + (j-1)*m) = A(i, ipvt(j))
      end do
    end do
    do j = 1, n
      do i = 1, m
        A(i,j) = hscr(i + (j-1)*m)
      end do
    end do
    call la_geqrf(m, n, A, lda, tau)
  end subroutine la_qrpp

  logical function any_not_finite(n, x)
    integer, intent(in) :: n
    real(8), intent(in) :: x(*)
    integer :: i
    any_not_finite = .false.
    do i = 1, n
      if (not_finite(x(i))) any_not_finite = .true.
    end do
  end function any_not_finite

  !> ZGEQP3 with full pivoting (all columns free) on the host, for the STABLOG path. Bit-identical to ZGEQP3.
  subroutine la_geqp3_host(m, n, A, lda, ipvt, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer(c_int), target :: ipvt(*)
    integer :: lw, info
    complex(8) :: w(1)
    if (m <= 0 .or. n <= 0) return
    call ensure_hrws(2*n)
    ipvt(1:n) = 0
    call zgeqp3(m, n, A, lda, ipvt, tau, w, -1, hrws, info)
    if (info /= 0) LA_FAIL('zgeqp3 work query failed')
    lw = int(dble(w(1)))
    call ensure_hws(lw)
    call zgeqp3(m, n, A, lda, ipvt, tau, hws, lw, hrws, info)
    if (info /= 0) LA_FAIL('zgeqp3 failed')
  end subroutine la_geqp3_host

end module alf_la_mod

! Backend layer ("alf_la" preview).  Every routine takes arrays that are already DEVICE-RESIDENT
! (mapped by the caller with !$omp target enter data) in the CUDA build, and plain host arrays in
! the HOST build.  Only this module knows which backend is compiled in:
!   -DBENCH_LA_CUDA : cuBLAS / cuSOLVER through bind(C) interfaces, device pointers obtained with
!                     !$omp target data use_device_addr
!   -DBENCH_LA_HOST : host BLAS/LAPACK
!   -DBENCH_LA_ONEMKL : ifx + oneMKL (ticket 103).  BLAS (zgemm, dgemm, zgemm_batch_strided, zgeru)
!                     and, by default, the LAPACK QR family go through "!$omp dispatch" on the
!                     oneMKL OpenMP-offload Fortran interfaces (the MKL variant runs on the Intel GPU
!                     when the device is a "gen" arch, otherwise the plain host routine is called, which is
!                     the host-fallback mode used on nividia).  With -DBENCH_LA_SHIM the QR family
!                     (geqrf/ungqr/unmqr) instead calls the C++ SYCL shim bench_sycl_shim.cpp with device
!                     pointers from use_device_addr and the SYCL queue obtained from OpenMP interop.
#ifdef BENCH_LA_ONEMKL
include "mkl_blas_omp_offload.f90"
include "mkl_lapack_omp_offload.f90"
#endif
module bench_la
  use iso_c_binding

  implicit none
  private
  integer, parameter, public :: LA_N = 0, LA_T = 1, LA_C = 2    ! same values as cublasOperation_t
  public :: la_init, la_finalize, la_sync, la_backend_name
  public :: la_gemm, la_gemm_cr, la_gemm_batched, la_geru
  public :: la_geqrf_lwork, la_geqrf, la_ungqr_lwork, la_ungqr, la_ormqr_lwork, la_ormqr

#ifdef BENCH_LA_CUDA
  type(c_ptr), save :: hblas = c_null_ptr, hsolv = c_null_ptr
  ! alpha/beta are passed to cuBLAS by host pointer.  cuBLAS reads them as 16-byte-aligned
  ! cuDoubleComplex (a Fortran complex(8) PARAMETER or local is only 8-byte aligned: SIGSEGV
  ! inside libcublas, observed).  An allocated array is malloc-aligned (16 B).
  complex(8), allocatable, target, save :: cab(:)
  interface
    integer(c_int) function cublasCreate_v2(h) bind(C, name='cublasCreate_v2')
      import :: c_ptr, c_int
      type(c_ptr) :: h
    end function
    integer(c_int) function cublasDestroy_v2(h) bind(C, name='cublasDestroy_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h
    end function
    integer(c_int) function cusolverDnCreate(h) bind(C, name='cusolverDnCreate')
      import :: c_ptr, c_int
      type(c_ptr) :: h
    end function
    integer(c_int) function cusolverDnDestroy(h) bind(C, name='cusolverDnDestroy')
      import :: c_ptr, c_int
      type(c_ptr), value :: h
    end function
    integer(c_int) function cudaDeviceSynchronize() bind(C, name='cudaDeviceSynchronize')
      import :: c_int
    end function
    integer(c_int) function cublasZgemm_v2(h, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind(C, name='cublasZgemm_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, A, B, beta, C
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
    end function
    integer(c_int) function cublasDgemm_v2(h, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind(C, name='cublasDgemm_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, A, B, beta, C
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
    end function
    integer(c_int) function cublasZgemmStridedBatched(h, ta, tb, m, n, k, alpha, A, lda, sA, B, ldb, sB, &
         beta, C, ldc, sC, nb) bind(C, name='cublasZgemmStridedBatched')
      import :: c_ptr, c_int, c_long_long
      type(c_ptr), value :: h, alpha, A, B, beta, C
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc, nb
      integer(c_long_long), value :: sA, sB, sC
    end function
    integer(c_int) function cublasZgeru_v2(h, m, n, alpha, x, incx, y, incy, A, lda) &
         bind(C, name='cublasZgeru_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, x, y, A
      integer(c_int), value :: m, n, incx, incy, lda
    end function
    integer(c_int) function cusolverDnZgeqrf_bufferSize(h, m, n, A, lda, lwork) &
         bind(C, name='cusolverDnZgeqrf_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A
      integer(c_int), value :: m, n, lda
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZgeqrf(h, m, n, A, lda, tau, work, lwork, info) &
         bind(C, name='cusolverDnZgeqrf')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, tau, work, info
      integer(c_int), value :: m, n, lda, lwork
    end function
    integer(c_int) function cusolverDnZungqr_bufferSize(h, m, n, k, A, lda, tau, lwork) &
         bind(C, name='cusolverDnZungqr_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, tau
      integer(c_int), value :: m, n, k, lda
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZungqr(h, m, n, k, A, lda, tau, work, lwork, info) &
         bind(C, name='cusolverDnZungqr')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, tau, work, info
      integer(c_int), value :: m, n, k, lda, lwork
    end function
    integer(c_int) function cusolverDnZunmqr_bufferSize(h, side, trans, m, n, k, A, lda, tau, C, ldc, lwork) &
         bind(C, name='cusolverDnZunmqr_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, tau, C
      integer(c_int), value :: side, trans, m, n, k, lda, ldc
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZunmqr(h, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info) &
         bind(C, name='cusolverDnZunmqr')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, tau, C, work, info
      integer(c_int), value :: side, trans, m, n, k, lda, ldc, lwork
    end function
  end interface
#elif defined(BENCH_LA_ONEMKL)
#ifdef BENCH_LA_SHIM
  ! SYCL queue handle (sycl::queue*) obtained from OpenMP interop; null = shim creates its own queue
  ! (only meaningful in host fallback, see bench_sycl_shim.cpp)
  type(c_ptr), save :: qh = c_null_ptr
  integer(c_intptr_t), save :: iop = 0
  integer(c_int64_t), parameter :: I8 = c_int64_t
  interface
    type(c_ptr) function omp_get_interop_ptr(obj, prop, rc) bind(C, name='omp_get_interop_ptr')
      import :: c_ptr, c_intptr_t, c_int
      integer(c_intptr_t), value :: obj
      integer(c_int), value :: prop
      integer(c_int) :: rc
    end function
    integer(c_intptr_t) function omp_get_interop_int(obj, prop, rc) bind(C, name='omp_get_interop_int')
      import :: c_intptr_t, c_int
      integer(c_intptr_t), value :: obj
      integer(c_int), value :: prop
      integer(c_int) :: rc
    end function
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
  end interface
#endif
#else
  external :: zgemm, dgemm, zgeru, zgeqrf, zungqr, zunmqr
#endif

contains

  function la_backend_name() result(s)
    character(len=16) :: s
#ifdef BENCH_LA_CUDA
    s = 'cuda'
#elif defined(BENCH_LA_SHIM)
    s = 'onemkl+shim'
#elif defined(BENCH_LA_ONEMKL)
    s = 'onemkl'
#else
    s = 'host'
#endif
  end function

  subroutine chk(st, what)
    integer(c_int), intent(in) :: st
    character(len=*), intent(in) :: what
    if (st /= 0) then
      write(*,*) 'LA backend error ', st, ' in ', what
      error stop 2
    end if
  end subroutine

  subroutine la_init()
#ifdef BENCH_LA_CUDA
    allocate(cab(2))
    call chk(cublasCreate_v2(hblas), 'cublasCreate')
    call chk(cusolverDnCreate(hsolv), 'cusolverDnCreate')
#elif defined(BENCH_LA_SHIM)
    integer(c_int) :: rc, dummy
    integer(c_intptr_t) :: fr
    iop = 0
    ! OpenMP interop object for the default device; its targetsync property is the sycl::queue the
    ! OpenMP runtime uses for that device.  With no offload device iop stays 0 (omp_interop_none).
    !$omp interop init(prefer_type("sycl"), targetsync: iop)
    qh = c_null_ptr
    if (iop /= 0) then
      fr = omp_get_interop_int(iop, -1, rc)          ! omp_ipr_fr_id
      if (rc == 0 .and. fr == 4) then               ! omp_ifr_sycl
        qh = omp_get_interop_ptr(iop, -9, rc)        ! omp_ipr_targetsync
        if (rc /= 0) qh = c_null_ptr
      end if
    end if
    dummy = alfsh_describe(qh)    ! prints which queue/device the shim will use
#endif
  end subroutine

  subroutine la_finalize()
#ifdef BENCH_LA_CUDA
    call chk(cublasDestroy_v2(hblas), 'cublasDestroy')
    call chk(cusolverDnDestroy(hsolv), 'cusolverDnDestroy')
#elif defined(BENCH_LA_SHIM)
    if (iop /= 0) then
      !$omp interop destroy(iop)
    end if
#endif
  end subroutine

  ! Block until all device work is finished (needed before reading the clock).
  subroutine la_sync()
#ifdef BENCH_LA_CUDA
    call chk(cudaDeviceSynchronize(), 'cudaDeviceSynchronize')
#endif
  end subroutine

  character function opc(op)
    integer, intent(in) :: op
    select case (op)
    case (LA_N); opc = 'N'
    case (LA_T); opc = 'T'
    case default; opc = 'C'
    end select
  end function

  ! C = alpha op(A) op(B) + beta C, complex(8)
  subroutine la_gemm(opa, opb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#ifdef BENCH_LA_ONEMKL
    use onemkl_blas_omp_offload_lp64, only: zgemm
#endif
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc
    complex(8), intent(in), target :: alpha, beta
    complex(8), target :: A(lda,*), B(ldb,*), C(ldc,*)
#ifdef BENCH_LA_CUDA
    cab(1) = alpha; cab(2) = beta
    !$omp target data use_device_addr(A, B, C)
    call chk(cublasZgemm_v2(hblas, opa, opb, m, n, k, c_loc(cab(1)), c_loc(A(1,1)), lda, &
         c_loc(B(1,1)), ldb, c_loc(cab(2)), c_loc(C(1,1)), ldc), 'zgemm')
    !$omp end target data
#elif defined(BENCH_LA_ONEMKL)
    character(len=1) :: ta, tb
    ta = opc(opa); tb = opc(opb)
    !$omp dispatch
    call zgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#else
    call zgemm(opc(opa), opc(opb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#endif
  end subroutine

  ! C(complex m x n) = A(complex m x k) * B(REAL k x n): one DGEMM on the interleaved real view
  ! of A and C (real 2m x k and 2m x n).  lda, ldc in complex elements.
  subroutine la_gemm_cr(m, n, k, A, lda, B, ldb, C, ldc)
#ifdef BENCH_LA_ONEMKL
    use onemkl_blas_omp_offload_lp64, only: dgemm
#endif
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(8), target :: A(lda,*), C(ldc,*)
    real(8), target :: B(ldb,*)
    real(8), target :: one, zero
#ifdef BENCH_LA_ONEMKL
    real(8), pointer, contiguous :: Ar(:,:), Cr(:,:)
#endif
    one = 1d0; zero = 0d0
#ifdef BENCH_LA_CUDA
    !$omp target data use_device_addr(A, B, C)
    call chk(cublasDgemm_v2(hblas, LA_N, LA_N, 2*m, n, k, c_loc(one), c_loc(A(1,1)), 2*lda, &
         c_loc(B(1,1)), ldb, c_loc(zero), c_loc(C(1,1)), 2*ldc), 'dgemm')
    !$omp end target data
#elif defined(BENCH_LA_ONEMKL)
    ! real view of the complex arrays (CONTIGUOUS pointers: no copy-in/out when passed on)
    call c_f_pointer(c_loc(A(1,1)), Ar, [2*lda, k])
    call c_f_pointer(c_loc(C(1,1)), Cr, [2*ldc, n])
    !$omp dispatch
    call dgemm('N', 'N', 2*m, n, k, one, Ar, 2*lda, B, ldb, zero, Cr, 2*ldc)
#else
    call dgemm('N', 'N', 2*m, n, k, one, A, 2*lda, B, ldb, zero, C, 2*ldc)
#endif
  end subroutine

  ! nb independent complex GEMMs with constant strides (in complex elements)
  subroutine la_gemm_batched(opa, opb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
#ifdef BENCH_LA_ONEMKL
    use onemkl_blas_omp_offload_lp64, only: zgemm_batch_strided
#endif
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc, nb, sA, sB, sC
    complex(8), intent(in), target :: alpha, beta
    complex(8), target :: A(*), B(*), C(*)
#ifdef BENCH_LA_CUDA
    cab(1) = alpha; cab(2) = beta
    !$omp target data use_device_addr(A, B, C)
    call chk(cublasZgemmStridedBatched(hblas, opa, opb, m, n, k, c_loc(cab(1)), c_loc(A(1)), lda, &
         int(sA, c_long_long), c_loc(B(1)), ldb, int(sB, c_long_long), c_loc(cab(2)), c_loc(C(1)), ldc, &
         int(sC, c_long_long), nb), 'zgemmStridedBatched')
    !$omp end target data
#elif defined(BENCH_LA_ONEMKL)
    character(len=1) :: ta, tb
    ta = opc(opa); tb = opc(opb)
    !$omp dispatch
    call zgemm_batch_strided(ta, tb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
#else
    integer :: ib
    do ib = 0, nb - 1
      call zgemm(opc(opa), opc(opb), m, n, k, alpha, A(1+ib*sA), lda, B(1+ib*sB), ldb, beta, C(1+ib*sC), ldc)
    end do
#endif
  end subroutine

  ! A += alpha x y^T (no conjugation)
  subroutine la_geru(m, n, alpha, x, incx, y, incy, A, lda)
#ifdef BENCH_LA_ONEMKL
    use onemkl_blas_omp_offload_lp64, only: zgeru
#endif
    integer, intent(in) :: m, n, incx, incy, lda
    complex(8), intent(in), target :: alpha
    complex(8), target :: x(*), y(*), A(lda,*)
#ifdef BENCH_LA_CUDA
    cab(1) = alpha
    !$omp target data use_device_addr(x, y, A)
    call chk(cublasZgeru_v2(hblas, m, n, c_loc(cab(1)), c_loc(x(1)), incx, c_loc(y(1)), incy, &
         c_loc(A(1,1)), lda), 'zgeru')
    !$omp end target data
#elif defined(BENCH_LA_ONEMKL)
    !$omp dispatch
    call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
#else
    call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
#endif
  end subroutine

  ! ---- QR family.  work/tau/dinfo are caller-owned (device-resident in the CUDA build). ----
  function la_geqrf_lwork(m, n, A, lda) result(lw)
#if defined(BENCH_LA_ONEMKL) && !defined(BENCH_LA_SHIM)
    use onemkl_lapack_omp_offload_lp64, only: zgeqrf
#endif
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer :: lw
#ifdef BENCH_LA_CUDA
    integer(c_int) :: l
    !$omp target data use_device_addr(A)
    call chk(cusolverDnZgeqrf_bufferSize(hsolv, m, n, c_loc(A(1,1)), lda, l), 'geqrf_bufferSize')
    !$omp end target data
    lw = l
#elif defined(BENCH_LA_SHIM)
    integer(c_int64_t) :: ws
    call chk(alfsh_zgeqrf_ws(qh, int(m,I8), int(n,I8), int(lda,I8), ws), 'geqrf_scratchpad_size')
    lw = int(ws)
#elif defined(BENCH_LA_ONEMKL)
    ! workspace query: a plain host call (lwork = -1 only writes work(1)); no device data involved
    complex(8) :: w(1), t1(1)
    integer :: info
    call zgeqrf(m, n, A, lda, t1, w, -1, info)
    lw = int(real(w(1)))
#else
    complex(8) :: w(1)
    integer :: info
    call zgeqrf(m, n, A, lda, w, w, -1, info)
    lw = int(real(w(1)))
#endif
  end function

  subroutine la_geqrf(m, n, A, lda, tau, work, lwork, dinfo)
#if defined(BENCH_LA_ONEMKL) && !defined(BENCH_LA_SHIM)
    use onemkl_lapack_omp_offload_lp64, only: zgeqrf
#endif
    integer, intent(in) :: m, n, lda, lwork
    complex(8), target :: A(lda,*), tau(*), work(*)
    integer(4), target :: dinfo(*)
#ifdef BENCH_LA_CUDA
    !$omp target data use_device_addr(A, tau, work, dinfo)
    call chk(cusolverDnZgeqrf(hsolv, m, n, c_loc(A(1,1)), lda, c_loc(tau(1)), c_loc(work(1)), lwork, &
         c_loc(dinfo(1))), 'geqrf')
    !$omp end target data
#elif defined(BENCH_LA_SHIM)
    ! SYCL USM device pointers via use_device_addr; info comes back as the function result (no device write)
    !$omp target data use_device_addr(A, tau, work)
    call chk(alfsh_zgeqrf(qh, int(m,I8), int(n,I8), c_loc(A(1,1)), int(lda,I8), c_loc(tau(1)), &
         c_loc(work(1)), int(lwork,I8)), 'geqrf (shim)')
    !$omp end target data
    dinfo(1) = 0
#elif defined(BENCH_LA_ONEMKL)
    !$omp dispatch
    call zgeqrf(m, n, A, lda, tau, work, lwork, dinfo(1))
#else
    call zgeqrf(m, n, A, lda, tau, work, lwork, dinfo(1))
#endif
  end subroutine

  function la_ungqr_lwork(m, n, k, A, lda, tau) result(lw)
#if defined(BENCH_LA_ONEMKL) && !defined(BENCH_LA_SHIM)
    use onemkl_lapack_omp_offload_lp64, only: zungqr
#endif
    integer, intent(in) :: m, n, k, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
#ifdef BENCH_LA_CUDA
    integer(c_int) :: l
    !$omp target data use_device_addr(A, tau)
    call chk(cusolverDnZungqr_bufferSize(hsolv, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), l), &
         'ungqr_bufferSize')
    !$omp end target data
    lw = l
#elif defined(BENCH_LA_SHIM)
    integer(c_int64_t) :: ws
    call chk(alfsh_zungqr_ws(qh, int(m,I8), int(n,I8), int(k,I8), int(lda,I8), ws), 'ungqr_scratchpad_size')
    lw = int(ws)
#elif defined(BENCH_LA_ONEMKL)
    complex(8) :: w(1)
    integer :: info
    call zungqr(m, n, k, A, lda, tau, w, -1, info)     ! host workspace query
    lw = int(real(w(1)))
#else
    complex(8) :: w(1)
    integer :: info
    call zungqr(m, n, k, A, lda, tau, w, -1, info)
    lw = int(real(w(1)))
#endif
  end function

  subroutine la_ungqr(m, n, k, A, lda, tau, work, lwork, dinfo)
#if defined(BENCH_LA_ONEMKL) && !defined(BENCH_LA_SHIM)
    use onemkl_lapack_omp_offload_lp64, only: zungqr
#endif
    integer, intent(in) :: m, n, k, lda, lwork
    complex(8), target :: A(lda,*), tau(*), work(*)
    integer(4), target :: dinfo(*)
#ifdef BENCH_LA_CUDA
    !$omp target data use_device_addr(A, tau, work, dinfo)
    call chk(cusolverDnZungqr(hsolv, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), c_loc(work(1)), lwork, &
         c_loc(dinfo(1))), 'ungqr')
    !$omp end target data
#elif defined(BENCH_LA_SHIM)
    !$omp target data use_device_addr(A, tau, work)
    call chk(alfsh_zungqr(qh, int(m,I8), int(n,I8), int(k,I8), c_loc(A(1,1)), int(lda,I8), c_loc(tau(1)), &
         c_loc(work(1)), int(lwork,I8)), 'ungqr (shim)')
    !$omp end target data
    dinfo(1) = 0
#elif defined(BENCH_LA_ONEMKL)
    !$omp dispatch
    call zungqr(m, n, k, A, lda, tau, work, lwork, dinfo(1))
#else
    call zungqr(m, n, k, A, lda, tau, work, lwork, dinfo(1))
#endif
  end subroutine

  ! C := op(Q) C (side='L') with Q from geqrf stored in A.  trans: LA_N or LA_C.
  function la_ormqr_lwork(trans, m, n, k, A, lda, tau, C, ldc) result(lw)
#if defined(BENCH_LA_ONEMKL) && !defined(BENCH_LA_SHIM)
    use onemkl_lapack_omp_offload_lp64, only: zunmqr
#endif
    integer, intent(in) :: trans, m, n, k, lda, ldc
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: lw
#ifdef BENCH_LA_CUDA
    integer(c_int) :: l
    !$omp target data use_device_addr(A, tau, C)
    call chk(cusolverDnZunmqr_bufferSize(hsolv, 0, trans, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), &
         c_loc(C(1,1)), ldc, l), 'unmqr_bufferSize')
    !$omp end target data
    lw = l
#elif defined(BENCH_LA_SHIM)
    integer(c_int64_t) :: ws
    call chk(alfsh_zunmqr_ws(qh, 0, trans, int(m,I8), int(n,I8), int(k,I8), int(lda,I8), int(ldc,I8), ws), &
         'unmqr_scratchpad_size')
    lw = int(ws)
#elif defined(BENCH_LA_ONEMKL)
    complex(8) :: w(1)
    integer :: info
    character(len=1) :: tc
    tc = opc(trans)
    call zunmqr('L', tc, m, n, k, A, lda, tau, C, ldc, w, -1, info)    ! host workspace query
    lw = int(real(w(1)))
#else
    complex(8) :: w(1)
    integer :: info
    call zunmqr('L', opc(trans), m, n, k, A, lda, tau, C, ldc, w, -1, info)
    lw = int(real(w(1)))
#endif
  end function

  subroutine la_ormqr(trans, m, n, k, A, lda, tau, C, ldc, work, lwork, dinfo)
#if defined(BENCH_LA_ONEMKL) && !defined(BENCH_LA_SHIM)
    use onemkl_lapack_omp_offload_lp64, only: zunmqr
#endif
    integer, intent(in) :: trans, m, n, k, lda, ldc, lwork
    complex(8), target :: A(lda,*), tau(*), C(ldc,*), work(*)
    integer(4), target :: dinfo(*)
#ifdef BENCH_LA_CUDA
    !$omp target data use_device_addr(A, tau, C, work, dinfo)
    call chk(cusolverDnZunmqr(hsolv, 0, trans, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), &
         c_loc(C(1,1)), ldc, c_loc(work(1)), lwork, c_loc(dinfo(1))), 'unmqr')
    !$omp end target data
#elif defined(BENCH_LA_SHIM)
    !$omp target data use_device_addr(A, tau, C, work)
    call chk(alfsh_zunmqr(qh, 0, trans, int(m,I8), int(n,I8), int(k,I8), c_loc(A(1,1)), int(lda,I8), &
         c_loc(tau(1)), c_loc(C(1,1)), int(ldc,I8), c_loc(work(1)), int(lwork,I8)), 'unmqr (shim)')
    !$omp end target data
    dinfo(1) = 0
#elif defined(BENCH_LA_ONEMKL)
    character(len=1) :: tc
    tc = opc(trans)
    !$omp dispatch
    call zunmqr('L', tc, m, n, k, A, lda, tau, C, ldc, work, lwork, dinfo(1))
#else
    call zunmqr('L', opc(trans), m, n, k, A, lda, tau, C, ldc, work, lwork, dinfo(1))
#endif
  end subroutine

end module bench_la

!  Copyright (C) 2026 The ALF project
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
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> CUDA backend of linalg_mod (GPU=nvidia, -DALF_LA_CUDA): cuBLAS / cuSOLVER on device-resident arrays.
!>
!> @details
!> Compiled only for GPU=nvidia (nvfortran -mp=gpu). linalg_mod calls the cu_* routine of a wrapper first; it returns
!> done = .false. when the operands are host arrays (the host call follows unchanged) and terminates when some operands
!> are mapped on the device and others are not (a missing map). Presence is tested with omp_target_is_present on the
!> first element of each matrix operand (P2.4 step 2, decision Q-N1).
!>  - Stream: cuBLAS and cuSOLVER run on the stream of the synchronous OpenMP target regions
!>    (acc_get_cuda_stream(acc_async_sync), a runtime query of the nvfortran openacc module; gpu-port measured
!>    cudaStreamPerThread). Library calls and the caller's synchronous target regions and updates are therefore
!>    ordered without extra synchronisation; testsuite/GPU.tests/56-cuda-order checks this. Regions with nowait are
!>    not ordered (none in ALF).
!>  - Errors: library status codes are checked at once. cuSOLVER info goes to a device slot; the slots are read back
!>    (one D->H) by la_cuda_sync and when all LA_NSLOT are used, a bad one terminates the run. The info returned to
!>    the caller of a device call is 0.
!>  - Work space: the caller's LAPACK work array is not used on the device; cuSOLVER work (dws) is owned here and
!>    grown on demand. Deterministic mode (cuSOLVER) and no atomics (cuBLAS): re-runs are bitwise reproducible.
!>  - The handles are created at the first device call, on the CUDA device of omp_get_default_device().
!>  - Not thread-safe (one set of handles and slots per process).
!--------------------------------------------------------------------
module linalg_cuda_mod
  use iso_c_binding, only: c_ptr, c_int, c_size_t, c_null_ptr, c_loc, c_associated
  use iso_fortran_env, only: error_unit, int64
  use omp_lib, only: omp_target_is_present, omp_get_default_device
  use openacc, only: acc_get_cuda_stream, acc_async_sync
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  implicit none
  private
  public :: cu_zgemm, cu_zgemm_cr, cu_ztrsm, cu_zgeqrf, cu_zunmqr, cu_zgetrf, cu_zgetri
  public :: cu_qrpp_on_device, cu_colnorms, cu_permute_cols
  public :: la_cuda_sync, la_cuda_finalize

  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: LA_NSLOT = 4096                 ! device info slots
  logical, save :: inited = .false.
  type(c_ptr), save :: hblas = c_null_ptr, hsolv = c_null_ptr, hstream = c_null_ptr
  ! alpha/beta passed to cuBLAS by host pointer: 16-byte aligned (allocated) storage; a complex(8) dummy or local is only
  ! 8-byte aligned and crashed libcublas (gpu-port ticket 101)
  complex(dp), allocatable, target, save :: cab(:)
  real(dp), allocatable, target, save :: rab(:)         ! 1 and 0 for cublasDgemm (same reason)
  complex(dp), allocatable, target, save :: dws(:)      ! cuSOLVER work space (device)
  complex(dp), allocatable, target, save :: dscr(:)     ! device scratch (column permutation, identity of getri)
  real(dp), allocatable, target, save :: dnrm(:)        ! column norms of la_zqrpp (device)
  integer(c_int), allocatable, target, save :: dinfo(:) ! cuSOLVER info slots (device)
  integer, save :: islot
  character(len=16), save :: islot_what(LA_NSLOT)
  logical, save :: islot_tol(LA_NSLOT)                  ! positive info tolerated (getrf: exactly singular)

  interface
    integer(c_int) function cudaGetDevice(d) bind(C, name='cudaGetDevice')
      import :: c_int; integer(c_int) :: d
    end function
    integer(c_int) function cudaSetDevice(d) bind(C, name='cudaSetDevice')
      import :: c_int; integer(c_int), value :: d
    end function
#if defined(ALF_LA_TEST_BADSTREAM)
    integer(c_int) function cudaStreamCreate(s) bind(C, name='cudaStreamCreate')
      import :: c_ptr, c_int; type(c_ptr) :: s
    end function
#endif
    integer(c_int) function cudaStreamSynchronize(s) bind(C, name='cudaStreamSynchronize')
      import :: c_ptr, c_int; type(c_ptr), value :: s
    end function
    integer(c_int) function cublasCreate_v2(h) bind(C, name='cublasCreate_v2')
      import :: c_ptr, c_int; type(c_ptr) :: h
    end function
    integer(c_int) function cublasDestroy_v2(h) bind(C, name='cublasDestroy_v2')
      import :: c_ptr, c_int; type(c_ptr), value :: h
    end function
    integer(c_int) function cublasSetStream_v2(h, s) bind(C, name='cublasSetStream_v2')
      import :: c_ptr, c_int; type(c_ptr), value :: h, s
    end function
    integer(c_int) function cublasSetAtomicsMode(h, mode) bind(C, name='cublasSetAtomicsMode')
      import :: c_ptr, c_int; type(c_ptr), value :: h; integer(c_int), value :: mode
    end function
    integer(c_int) function cusolverDnCreate(h) bind(C, name='cusolverDnCreate')
      import :: c_ptr, c_int; type(c_ptr) :: h
    end function
    integer(c_int) function cusolverDnDestroy(h) bind(C, name='cusolverDnDestroy')
      import :: c_ptr, c_int; type(c_ptr), value :: h
    end function
    integer(c_int) function cusolverDnSetStream(h, s) bind(C, name='cusolverDnSetStream')
      import :: c_ptr, c_int; type(c_ptr), value :: h, s
    end function
    integer(c_int) function cusolverDnSetDeterministicMode(h, mode) bind(C, name='cusolverDnSetDeterministicMode')
      import :: c_ptr, c_int; type(c_ptr), value :: h; integer(c_int), value :: mode
    end function
    integer(c_int) function cublasZgemm_v2(h, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(C, name='cublasZgemm_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, a, b, beta, c
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
    end function
    integer(c_int) function cublasDgemm_v2(h, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(C, name='cublasDgemm_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, a, b, beta, c
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
    end function
    integer(c_int) function cublasZtrsm_v2(h, side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb) bind(C, name='cublasZtrsm_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, a, b
      integer(c_int), value :: side, uplo, trans, diag, m, n, lda, ldb
    end function
    integer(c_int) function cusolverDnZgeqrf_bufferSize(h, m, n, a, lda, lwork) bind(C, name='cusolverDnZgeqrf_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a
      integer(c_int), value :: m, n, lda
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZgeqrf(h, m, n, a, lda, tau, work, lwork, info) bind(C, name='cusolverDnZgeqrf')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a, tau, work, info
      integer(c_int), value :: m, n, lda, lwork
    end function
    integer(c_int) function cusolverDnZunmqr_bufferSize(h, side, trans, m, n, k, a, lda, tau, c, ldc, lwork) &
         & bind(C, name='cusolverDnZunmqr_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a, tau, c
      integer(c_int), value :: side, trans, m, n, k, lda, ldc
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZunmqr(h, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info) &
         & bind(C, name='cusolverDnZunmqr')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a, tau, c, work, info
      integer(c_int), value :: side, trans, m, n, k, lda, ldc, lwork
    end function
    integer(c_int) function cusolverDnZgetrf_bufferSize(h, m, n, a, lda, lwork) bind(C, name='cusolverDnZgetrf_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a
      integer(c_int), value :: m, n, lda
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZgetrf(h, m, n, a, lda, work, ipiv, info) bind(C, name='cusolverDnZgetrf')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a, work, ipiv, info
      integer(c_int), value :: m, n, lda
    end function
    integer(c_int) function cusolverDnZgetrs(h, trans, n, nrhs, a, lda, ipiv, b, ldb, info) bind(C, name='cusolverDnZgetrs')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, a, ipiv, b, info
      integer(c_int), value :: trans, n, nrhs, lda, ldb
    end function
  end interface

contains

  !> Stops the run with a message (backend failures and misuse)
  subroutine cu_fail(msg)
    character(len=*), intent(in) :: msg
    write(error_unit,'(2A)') 'linalg_mod CUDA backend: ', msg
    call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  end subroutine cu_fail

  subroutine chk(st, what)
    integer(c_int), intent(in) :: st
    character(len=*), intent(in) :: what
    character(len=80) :: msg
    if (st /= 0) then
      write(msg,'(A,I0,2A)') 'library status ', st, ' in ', what
      call cu_fail(trim(msg))
    endif
  end subroutine chk

  !> Creates the handles on the stream of the synchronous target regions (once)
  subroutine cu_init()
    integer(c_int) :: d
    if (inited) return
    ! the handles must live on the device of the OpenMP data (gpu_state_mod sets the default device per MPI rank)
    call chk(cudaGetDevice(d), 'cudaGetDevice')
    if (d /= omp_get_default_device()) call chk(cudaSetDevice(int(omp_get_default_device(), c_int)), 'cudaSetDevice')
    allocate(cab(2), rab(2), dinfo(LA_NSLOT))
    rab(1) = 1.d0; rab(2) = 0.d0
    dinfo = 0
    !$omp target enter data map(to: dinfo)
#if defined(ALF_LA_TEST_BADSTREAM)
    ! negative control of 56-cuda-order only: a created stream, not ordered against the target regions
    call chk(cudaStreamCreate(hstream), 'cudaStreamCreate')
#else
    hstream = transfer(acc_get_cuda_stream(acc_async_sync), c_null_ptr)
#endif
    call chk(cublasCreate_v2(hblas), 'cublasCreate')
    call chk(cusolverDnCreate(hsolv), 'cusolverDnCreate')
    call chk(cublasSetStream_v2(hblas, hstream), 'cublasSetStream')
    call chk(cusolverDnSetStream(hsolv, hstream), 'cusolverDnSetStream')
    call chk(cusolverDnSetDeterministicMode(hsolv, 1_c_int), 'cusolverDnSetDeterministicMode')
    call chk(cublasSetAtomicsMode(hblas, 0_c_int), 'cublasSetAtomicsMode')
    islot = 0
    inited = .true.
  end subroutine cu_init

  !> Waits for the device work of the stream and checks the cuSOLVER info slots (no-op before the first device call)
  subroutine la_cuda_sync()
    integer :: i
    integer(c_int), allocatable :: h(:)
    character(len=80) :: msg
    if (.not. inited) return
    call chk(cudaStreamSynchronize(hstream), 'cudaStreamSynchronize')
    if (islot == 0) return
    !$omp target update from(dinfo(1:islot))
    h = dinfo(1:islot)
    do i = 1, islot
      if (h(i) < 0 .or. (h(i) > 0 .and. .not. islot_tol(i))) then
        write(msg,'(3A,I0)') 'cuSOLVER ', trim(islot_what(i)), ' failed, info = ', h(i)
        call cu_fail(trim(msg))
      endif
    enddo
    islot = 0
  end subroutine la_cuda_sync

  !> Last sync and slot check, destroys the handles and frees the device scratch
  subroutine la_cuda_finalize()
    if (.not. inited) return
    call la_cuda_sync()
    call chk(cublasDestroy_v2(hblas), 'cublasDestroy')
    call chk(cusolverDnDestroy(hsolv), 'cusolverDnDestroy')
    !$omp target exit data map(delete: dinfo)
    deallocate(dinfo, cab, rab)
    call free_c(dws)
    call free_c(dscr)
    if (allocated(dnrm)) then
      !$omp target exit data map(delete: dnrm)
      deallocate(dnrm)
    endif
    inited = .false.
  end subroutine la_cuda_finalize

  subroutine free_c(x)
    complex(dp), allocatable, intent(inout) :: x(:)
    if (.not. allocated(x)) return
    !$omp target exit data map(delete: x)
    deallocate(x)
  end subroutine free_c

  !> Device array x of at least n elements (old contents are lost)
  subroutine grow_c(x, n)
    complex(dp), allocatable, intent(inout) :: x(:)
    integer, intent(in) :: n
    if (allocated(x)) then
      if (size(x) >= n) return
      call la_cuda_sync()                                ! the old buffer may still be in use
      call free_c(x)
    endif
    allocate(x(max(n,1)))
    !$omp target enter data map(alloc: x)
  end subroutine grow_c

  !> Next info slot (reads all slots back first when they are used up)
  integer function next_slot(what, tol)
    character(len=*), intent(in) :: what
    logical, intent(in) :: tol
    if (islot >= LA_NSLOT) call la_cuda_sync()
    islot = islot + 1
    islot_what(islot) = what
    islot_tol(islot) = tol
    next_slot = islot
  end function next_slot

  !> .true. if all operands are on the device, .false. if none is; terminates on a mix. Null pointers are skipped.
  logical function cu_on_device(rname, p)
    character(len=*), intent(in) :: rname
    type(c_ptr), intent(in) :: p(:)
    integer :: i, np, nh
    np = 0; nh = 0
    do i = 1, size(p)
      if (.not. c_associated(p(i))) cycle
      if (omp_target_is_present(p(i), omp_get_default_device()) /= 0) then
        np = np + 1
      else
        nh = nh + 1
      endif
    enddo
    if (np > 0 .and. nh > 0) call cu_fail(rname//': some operands are mapped on the device, others are not (missing map)')
    cu_on_device = (np > 0)
    if (cu_on_device) call cu_init()
  end function cu_on_device

  !> Device test of the operands of la_zqrpp (a, tau, ipvt), as cu_on_device
  logical function cu_qrpp_on_device(a, tau, ipvt)
    complex(dp), target, intent(in) :: a(*), tau(*)
    integer(c_int), target, intent(in) :: ipvt(*)
    cu_qrpp_on_device = cu_on_device('la_zqrpp', [c_loc(a(1)), c_loc(tau(1)), c_loc(ipvt(1))])
  end function cu_qrpp_on_device

  integer(c_int) function op(c)
    character(len=1), intent(in) :: c
    select case (c)
    case ('N', 'n'); op = 0
    case ('T', 't'); op = 1
    case ('C', 'c'); op = 2
    case default; op = -1; call cu_fail('unknown op '//c)
    end select
  end function op

  integer(c_int) function sidec(c)
    character(len=1), intent(in) :: c
    sidec = merge(1_c_int, 0_c_int, c == 'R' .or. c == 'r')      ! CUBLAS_SIDE_LEFT = 0, RIGHT = 1
  end function sidec

  !> C = alpha op(A) op(B) + beta C on the device (cublasZgemm)
  subroutine cu_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, done)
    character(len=1), intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(dp), intent(in) :: alpha, beta
    complex(dp), target, intent(in) :: a(lda,*), b(ldb,*)
    complex(dp), target :: c(ldc,*)
    logical, intent(out) :: done
    done = .false.
    if (m <= 0 .or. n <= 0) return
    if (.not. cu_on_device('la_zgemm', [c_loc(a(1,1)), c_loc(b(1,1)), c_loc(c(1,1))])) return
    cab(1) = alpha; cab(2) = beta
    !$omp target data use_device_addr(a, b, c)
    call chk(cublasZgemm_v2(hblas, op(transa), op(transb), m, n, k, c_loc(cab(1)), c_loc(a(1,1)), lda, &
         & c_loc(b(1,1)), ldb, c_loc(cab(2)), c_loc(c(1,1)), ldc), 'zgemm')
    !$omp end target data
    done = .true.
  end subroutine cu_zgemm

  !> C(complex) = A(complex) B(real) as one cublasDgemm on the interleaved real view (as la_zgemm_cr)
  subroutine cu_zgemm_cr(m, n, k, a, lda, b, ldb, c, ldc, done)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(dp), target, intent(in) :: a(lda,*)
    complex(dp), target :: c(ldc,*)
    real(dp), target, intent(in) :: b(ldb,*)
    logical, intent(out) :: done
    done = .false.
    if (m <= 0 .or. n <= 0) return
    if (.not. cu_on_device('la_zgemm_cr', [c_loc(a(1,1)), c_loc(b(1,1)), c_loc(c(1,1))])) return
    !$omp target data use_device_addr(a, b, c)
    call chk(cublasDgemm_v2(hblas, 0_c_int, 0_c_int, 2*m, n, k, c_loc(rab(1)), c_loc(a(1,1)), 2*lda, &
         & c_loc(b(1,1)), ldb, c_loc(rab(2)), c_loc(c(1,1)), 2*ldc), 'dgemm')
    !$omp end target data
    done = .true.
  end subroutine cu_zgemm_cr

  !> Triangular solve on the device (cublasZtrsm)
  subroutine cu_ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, done)
    character(len=1), intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(dp), intent(in) :: alpha
    complex(dp), target, intent(in) :: a(lda,*)
    complex(dp), target :: b(ldb,*)
    logical, intent(out) :: done
    integer(c_int) :: cu, cd
    done = .false.
    if (m <= 0 .or. n <= 0) return
    if (.not. cu_on_device('la_ztrsm', [c_loc(a(1,1)), c_loc(b(1,1))])) return
    cu = merge(1_c_int, 0_c_int, uplo == 'U' .or. uplo == 'u')    ! CUBLAS_FILL_MODE_LOWER = 0, UPPER = 1
    cd = merge(1_c_int, 0_c_int, diag == 'U' .or. diag == 'u')    ! CUBLAS_DIAG_NON_UNIT = 0, UNIT = 1
    cab(1) = alpha
    !$omp target data use_device_addr(a, b)
    call chk(cublasZtrsm_v2(hblas, sidec(side), cu, op(transa), cd, m, n, c_loc(cab(1)), c_loc(a(1,1)), lda, &
         & c_loc(b(1,1)), ldb), 'ztrsm')
    !$omp end target data
    done = .true.
  end subroutine cu_ztrsm

  !> QR factorisation on the device (cusolverDnZgeqrf; reflectors in the LAPACK convention)
  subroutine cu_zgeqrf(m, n, a, lda, tau, done)
    integer, intent(in) :: m, n, lda
    complex(dp), target :: a(lda,*), tau(*)
    logical, intent(out) :: done
    integer(c_int) :: lw
    integer :: s
    done = .false.
    if (m <= 0 .or. n <= 0) return
    if (.not. cu_on_device('la_zgeqrf', [c_loc(a(1,1)), c_loc(tau(1))])) return
    !$omp target data use_device_addr(a)
    call chk(cusolverDnZgeqrf_bufferSize(hsolv, m, n, c_loc(a(1,1)), lda, lw), 'zgeqrf_bufferSize')
    !$omp end target data
    call grow_c(dws, int(lw))
    s = next_slot('zgeqrf', .false.)
    !$omp target data use_device_addr(a, tau, dws, dinfo)
    call chk(cusolverDnZgeqrf(hsolv, m, n, c_loc(a(1,1)), lda, c_loc(tau(1)), c_loc(dws(1)), lw, c_loc(dinfo(s))), 'zgeqrf')
    !$omp end target data
    done = .true.
  end subroutine cu_zgeqrf

  !> C = op(Q) C or C op(Q) on the device (cusolverDnZunmqr)
  subroutine cu_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, done)
    character(len=1), intent(in) :: side, trans
    integer, intent(in) :: m, n, k, lda, ldc
    complex(dp), target, intent(in) :: a(lda,*), tau(*)
    complex(dp), target :: c(ldc,*)
    logical, intent(out) :: done
    integer(c_int) :: lw
    integer :: s
    done = .false.
    if (m <= 0 .or. n <= 0 .or. k <= 0) return
    if (.not. cu_on_device('la_zunmqr', [c_loc(a(1,1)), c_loc(tau(1)), c_loc(c(1,1))])) return
    !$omp target data use_device_addr(a, tau, c)
    call chk(cusolverDnZunmqr_bufferSize(hsolv, sidec(side), op(trans), m, n, k, c_loc(a(1,1)), lda, c_loc(tau(1)), &
         & c_loc(c(1,1)), ldc, lw), 'zunmqr_bufferSize')
    !$omp end target data
    call grow_c(dws, int(lw))
    s = next_slot('zunmqr', .false.)
    !$omp target data use_device_addr(a, tau, c, dws, dinfo)
    call chk(cusolverDnZunmqr(hsolv, sidec(side), op(trans), m, n, k, c_loc(a(1,1)), lda, c_loc(tau(1)), &
         & c_loc(c(1,1)), ldc, c_loc(dws(1)), lw, c_loc(dinfo(s))), 'zunmqr')
    !$omp end target data
    done = .true.
  end subroutine cu_zunmqr

  !> LU factorisation on the device (cusolverDnZgetrf); info > 0 (exactly singular) is tolerated as by the callers
  subroutine cu_zgetrf(m, n, a, lda, ipiv, done)
    integer, intent(in) :: m, n, lda
    complex(dp), target :: a(lda,*)
    integer(c_int), target :: ipiv(*)
    logical, intent(out) :: done
    integer(c_int) :: lw
    integer :: s
    done = .false.
    if (m <= 0 .or. n <= 0) return
    if (.not. cu_on_device('la_zgetrf', [c_loc(a(1,1)), c_loc(ipiv(1))])) return
    !$omp target data use_device_addr(a)
    call chk(cusolverDnZgetrf_bufferSize(hsolv, m, n, c_loc(a(1,1)), lda, lw), 'zgetrf_bufferSize')
    !$omp end target data
    call grow_c(dws, int(lw))
    s = next_slot('zgetrf', .true.)
    !$omp target data use_device_addr(a, dws, ipiv, dinfo)
    call chk(cusolverDnZgetrf(hsolv, m, n, c_loc(a(1,1)), lda, c_loc(dws(1)), c_loc(ipiv(1)), c_loc(dinfo(s))), 'zgetrf')
    !$omp end target data
    done = .true.
  end subroutine cu_zgetrf

  !> A^-1 from the LU factors on the device. cuSOLVER has no getri: getrs on the identity in dscr, then copied into A.
  !> Rounding differs from ZGETRI in the last bits.
  subroutine cu_zgetri(n, a, lda, ipiv, done)
    integer, intent(in) :: n, lda
    complex(dp), target :: a(lda,n)
    integer(c_int), target, intent(in) :: ipiv(*)
    logical, intent(out) :: done
    integer :: i, j, s
    done = .false.
    if (n <= 0) return
    if (.not. cu_on_device('la_zgetri', [c_loc(a(1,1)), c_loc(ipiv(1))])) return
    call grow_c(dscr, n*n)
    call ident(n, dscr)
    s = next_slot('zgetrs', .false.)
    !$omp target data use_device_addr(a, ipiv, dscr, dinfo)
    call chk(cusolverDnZgetrs(hsolv, 0_c_int, n, n, c_loc(a(1,1)), lda, c_loc(ipiv(1)), c_loc(dscr(1)), n, &
         & c_loc(dinfo(s))), 'zgetrs')
    !$omp end target data
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, n
        a(i,j) = dscr(i + (j-1)*n)
      enddo
    enddo
    done = .true.
  end subroutine cu_zgetri

  subroutine ident(n, x)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: x(n*n)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, n
        if (i == j) then                         ! no merge on complex: nvfortran drops such a kernel (CDSELECT)
          x(i + (j-1)*n) = (1.d0, 0.d0)
        else
          x(i + (j-1)*n) = (0.d0, 0.d0)
        endif
      enddo
    enddo
  end subroutine ident

  !> nrm(1:n) (host) = 2-norms of the columns of the device matrix A(1:m,1:n); NaN if a column has a NaN or Inf entry.
  !> One thread per column with a serial sum: fixed summation order, so pivots do not flip from run to run.
  subroutine cu_colnorms(m, n, a, lda, nrm)
    integer, intent(in) :: m, n, lda
    complex(dp), intent(in) :: a(lda,n)
    real(dp), intent(out) :: nrm(n)
    integer :: i, j, nf
    real(dp) :: mx, s, v
    integer(int64) :: b1, b2
    integer(int64), parameter :: mask = int(z'7FF0000000000000', int64)
    if (allocated(dnrm)) then
      if (size(dnrm) < n) then
        !$omp target exit data map(delete: dnrm)
        deallocate(dnrm)
      endif
    endif
    if (.not. allocated(dnrm)) then
      allocate(dnrm(n))
      !$omp target enter data map(alloc: dnrm)
    endif
    !$omp target teams distribute parallel do private(i, mx, s, v, nf, b1, b2)
    do j = 1, n
      mx = 0.d0
      nf = 0
      do i = 1, m
        b1 = transfer(dble(a(i,j)), 0_int64)
        b2 = transfer(aimag(a(i,j)), 0_int64)
        if (iand(b1, mask) == mask .or. iand(b2, mask) == mask) nf = nf + 1
        mx = max(mx, abs(dble(a(i,j))), abs(aimag(a(i,j))))
      enddo
      if (nf > 0) then
        dnrm(j) = transfer(int(z'7FF8000000000000', int64), 0.d0)
      else if (mx > 0.d0) then
        s = 0.d0
        do i = 1, m
          v = dble(a(i,j))/mx
          s = s + v*v
          v = aimag(a(i,j))/mx
          s = s + v*v
        enddo
        dnrm(j) = mx*sqrt(s)
      else
        dnrm(j) = 0.d0
      endif
    enddo
    ! D->H: n column norms for the host sort of la_zqrpp
    !$omp target update from(dnrm(1:n))
    nrm(1:n) = dnrm(1:n)
  end subroutine cu_colnorms

  !> Device A(:,j) <- A(:,ipvt(j)), j = 1..n; ipvt is uploaded to the device copy first (valid on host and device after)
  subroutine cu_permute_cols(m, n, a, lda, ipvt)
    integer, intent(in) :: m, n, lda
    complex(dp), intent(inout) :: a(lda,n)
    integer(c_int), intent(inout) :: ipvt(n)
    integer :: i, j
    ! H->D: the n-entry permutation
    !$omp target update to(ipvt(1:n))
    call grow_c(dscr, m*n)
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, m
        dscr(i + (j-1)*m) = a(i, ipvt(j))
      enddo
    enddo
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, m
        a(i,j) = dscr(i + (j-1)*m)
      enddo
    enddo
  end subroutine cu_permute_cols
end module linalg_cuda_mod

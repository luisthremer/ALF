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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your changes as modified.

include "mkl_blas_omp_offload.f90"

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> oneMKL device backend of linalg_mod (GPU=intel, -DALF_LA_ONEMKL; P2.4 step 2).
!>
!> @details
!> Every om_<name> has the arguments of la_<name> and returns .true. when it did the work on the device, .false. when
!> the caller must do the plain host call. The decision (presence dispatch, P2.4 step 2 rule 1): the matrix operands
!> are checked with omp_target_is_present (first and last element touched): all present -> device, none -> host,
!> a mix -> the run stops (a host array handed to the device would be read as garbage). Without an offload device
!> everything stays on the host.
!>  - BLAS (zgemm, the real-view dgemm of la_zgemm_cr, ztrsm): oneMKL OpenMP "dispatch" variants on the mapped arrays;
!>    synchronous (no nowait: ifx task creation costs ~26 us per op, measured on gpu-port), so they are ordered with
!>    the synchronous target regions of the gpu_* modules.
!>  - LAPACK (zgeqrf, zunmqr, zgetrf, zgetri): oneMKL USM LAPACK through the SYCL shim linalg_sycl_shim.cpp on the
!>    sycl::queue of the OpenMP runtime (OpenMP interop, prefer_type("sycl"), set up at the first device call). The
!>    shim waits for every call, so info is returned to the caller at once. lwork = -1 queries stay on the host (rule 1):
!>    the device calls ignore the caller's work array and use the module scratch dws (oneMKL scratchpad size).
!>  - zqrpp: column norms on the device (fixed summation order), host NaN test and the stable insertion sort of the
!>    host branch, ipvt H->D (valid on host and device afterwards), device column permutation, device zgeqrf.
!>  - zgeqp3 has no device path (the callers bracket it with target update).
!> Ported from fork/gpu-port alf_la_om_mod (Libraries/Modules/alf_la_mod.F90 :125-978) without its info slots and test
!> hooks.
!--------------------------------------------------------------------
module linalg_onemkl_mod
  use iso_c_binding, only: c_ptr, c_loc, c_int, c_int64_t, c_intptr_t, c_char, c_null_char, c_null_ptr, c_associated, &
       & c_f_pointer
  use omp_lib, only: omp_get_default_device, omp_get_num_devices, omp_target_is_present, omp_interop_kind, &
       & omp_interop_none, omp_ifr_sycl
  use onemkl_blas_omp_offload_lp64, only: zgemm, dgemm, ztrsm
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: om_zgemm, om_zgemm_cr, om_ztrsm, om_zgeqrf, om_zunmqr, om_zgetrf, om_zgetri, om_zqrpp, om_finalize, om_get_queue

  integer, save :: ndev = -1            ! number of offload devices, queried once
  logical, save :: inited = .false.     ! interop queue handed to the shim
  integer(omp_interop_kind), save :: iop = omp_interop_none
  type(c_ptr), save :: qsave = c_null_ptr  ! the sycl::queue* handed to the shim (also used by gpu_lowrank's block path)
  complex(kind=kind(0.d0)), allocatable, target, save :: dws(:)    ! device scratchpad of the shim calls
  integer(c_int64_t), allocatable, target, save :: dpiv(:)         ! device int64 pivots (oneMKL USM API)
  complex(kind=kind(0.d0)), allocatable, target, save :: dscr(:)   ! column permutation scratch of om_zqrpp
  real(kind=kind(0.d0)), allocatable, target, save :: dnrm(:)      ! column norms of om_zqrpp

  ! omp_get_interop_int/_ptr are missing from omp_lib of ifx 2025.x; the runtime exports them (OpenMP 5.1 values)
  integer(c_int), parameter :: omp_ipr_fr_id = -1, omp_ipr_targetsync = -9
  interface
    integer(c_intptr_t) function omp_get_interop_int(interop, property_id, ret_code) bind(C, name='omp_get_interop_int')
      import :: c_intptr_t, c_int
      integer(c_intptr_t), value :: interop
      integer(c_int), value :: property_id
      integer(c_int), intent(out) :: ret_code
    end function
    type(c_ptr) function omp_get_interop_ptr(interop, property_id, ret_code) bind(C, name='omp_get_interop_ptr')
      import :: c_intptr_t, c_int, c_ptr
      integer(c_intptr_t), value :: interop
      integer(c_int), value :: property_id
      integer(c_int), intent(out) :: ret_code
    end function
    integer(c_int) function alfsh_init(q) bind(C, name='alfsh_init')
      import :: c_ptr, c_int
      type(c_ptr), value :: q
    end function
    subroutine alfsh_finalize() bind(C, name='alfsh_finalize')
    end subroutine
    integer(c_int) function alfsh_device_info(buf, n, is_l0) bind(C, name='alfsh_device_info')
      import :: c_char, c_int
      character(kind=c_char) :: buf(*)
      integer(c_int), value :: n
      integer(c_int) :: is_l0
    end function
    integer(c_int) function alfsh_ptr_kind(p) bind(C, name='alfsh_ptr_kind')
      import :: c_ptr, c_int
      type(c_ptr), value :: p
    end function
    integer(c_int) function alfsh_wait() bind(C, name='alfsh_wait')
      import :: c_int
    end function
    integer(c_int) function alfsh_zgeqrf_ws(m, n, lda, ws) bind(C, name='alfsh_zgeqrf_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgeqrf(m, n, a, lda, tau, work, lwork) bind(C, name='alfsh_zgeqrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
    integer(c_int) function alfsh_zunmqr_ws(side, trans, m, n, k, lda, ldc, ws) bind(C, name='alfsh_zunmqr_ws')
      import :: c_int, c_int64_t
      integer(c_int), value :: side, trans
      integer(c_int64_t), value :: m, n, k, lda, ldc
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork) &
         bind(C, name='alfsh_zunmqr')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, c, work
      integer(c_int), value :: side, trans
      integer(c_int64_t), value :: m, n, k, lda, ldc, lwork
    end function
    integer(c_int) function alfsh_zgetrf_ws(m, n, lda, ws) bind(C, name='alfsh_zgetrf_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgetrf(m, n, a, lda, ipiv, work, lwork) bind(C, name='alfsh_zgetrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, ipiv, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
    integer(c_int) function alfsh_zgetri_ws(n, lda, ws) bind(C, name='alfsh_zgetri_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgetri(n, a, lda, ipiv, work, lwork) bind(C, name='alfsh_zgetri')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, ipiv, work
      integer(c_int64_t), value :: n, lda, lwork
    end function
  end interface

contains

  !> number of elements from the first to the last one touched by a rows x cols block with leading dimension ld
  integer function ext2(rows, cols, ld)
    integer, intent(in) :: rows, cols, ld
    ext2 = 0
    if (rows > 0 .and. cols > 0) ext2 = (cols-1)*ld + rows
  end function ext2

  !> first and last element (p1, pn) of a complex operand of nel elements (sequence association: no bounds problem)
  subroutine zends(a, nel, p1, pn)
    complex(kind=kind(0.d0)), target, intent(in) :: a(*)
    integer, intent(in) :: nel
    type(c_ptr), intent(out) :: p1, pn
    p1 = c_loc(a(1)); pn = c_loc(a(nel))
  end subroutine zends

  !> same for a real operand
  subroutine dends(a, nel, p1, pn)
    real(kind=kind(0.d0)), target, intent(in) :: a(*)
    integer, intent(in) :: nel
    type(c_ptr), intent(out) :: p1, pn
    p1 = c_loc(a(1)); pn = c_loc(a(nel))
  end subroutine dends

  !> same for an integer operand
  subroutine iends(a, nel, p1, pn)
    integer, target, intent(in) :: a(*)
    integer, intent(in) :: nel
    type(c_ptr), intent(out) :: p1, pn
    p1 = c_loc(a(1)); pn = c_loc(a(nel))
  end subroutine iends

  subroutine om_fail(msg, line)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: line
    write(error_unit,*) 'linalg_mod (oneMKL): ', msg
    call Terminate_on_error(ERROR_GENERIC, 'linalg_onemkl_mod.F90', line)
  end subroutine om_fail

  !> Hands the sycl::queue of the OpenMP runtime (default device) to the shim; checks that it is a Level Zero GPU and
  !> that device memory of the OpenMP default device is device USM of the queue's own device (catches a tile mismatch).
  subroutine om_init()
    type(c_ptr) :: qh
    integer(c_intptr_t) :: fr
    integer(c_int) :: rc, is_l0
    integer :: dev, i
    character(kind=c_char) :: nm(128)
    character(len=128) :: name
    complex(kind=kind(0.d0)), allocatable, target :: probe(:)
    if (inited) return
    dev = omp_get_default_device()
    iop = omp_interop_none
    !$omp interop init(prefer_type("sycl"), targetsync: iop) device(dev)
    if (iop == omp_interop_none) call om_fail('OpenMP interop object not available', __LINE__)
    fr = omp_get_interop_int(iop, omp_ipr_fr_id, rc)
    qh = c_null_ptr
    if (rc == 0 .and. fr == omp_ifr_sycl) qh = omp_get_interop_ptr(iop, omp_ipr_targetsync, rc)
    if (.not. c_associated(qh)) call om_fail('no sycl::queue from OpenMP interop', __LINE__)
    if (alfsh_init(qh) /= 0) call om_fail('shim initialisation failed', __LINE__)
    qsave = qh
    nm = c_null_char
    is_l0 = 0
    if (alfsh_device_info(nm, 128_c_int, is_l0) /= 0) call om_fail('cannot query the SYCL device', __LINE__)
    name = ''
    do i = 1, 127
      if (nm(i) == c_null_char) exit
      name(i:i) = nm(i)
    end do
    write(*,'(2A)') ' linalg_mod: oneMKL device: ', trim(name)
    if (is_l0 /= 1) call om_fail('the SYCL queue is not on a Level Zero GPU', __LINE__)
    allocate(probe(1))
    probe = (0.d0, 0.d0)
    !$omp target enter data map(to: probe)
    !$omp target data use_device_addr(probe)
    rc = alfsh_ptr_kind(c_loc(probe(1)))
    !$omp end target data
    !$omp target exit data map(delete: probe)
    if (rc /= 2) call om_fail('OpenMP default device and SYCL queue device differ', __LINE__)
    inited = .true.
  end subroutine om_init

  !> The sycl::queue* of the OpenMP runtime (default device) that every shim call uses; initialises the backend if needed.
  !> Used by the asynchronous block path of gpu_lowrank_mod (libalfsb, linalg_sycl_blk.cpp), which must use the same in-order queue.
  type(c_ptr) function om_get_queue()
    call om_init()
    om_get_queue = qsave
  end function om_get_queue

  !> End of run: frees the device scratch and the interop object (nothing if the backend was never used)
  subroutine om_finalize()
    if (.not. inited) return
    call om_sync()
    if (allocated(dws)) then
      !$omp target exit data map(delete: dws)
      deallocate(dws)
    end if
    if (allocated(dpiv)) then
      !$omp target exit data map(delete: dpiv)
      deallocate(dpiv)
    end if
    if (allocated(dscr)) then
      !$omp target exit data map(delete: dscr)
      deallocate(dscr)
    end if
    if (allocated(dnrm)) then
      !$omp target exit data map(delete: dnrm)
      deallocate(dnrm)
    end if
    call alfsh_finalize()
    !$omp interop destroy(iop)
    iop = omp_interop_none
    qsave = c_null_ptr
    inited = .false.
  end subroutine om_finalize

  subroutine om_sync()
    if (alfsh_wait() /= 0) call om_fail('queue wait failed', __LINE__)
  end subroutine om_sync

  !> device scratch of at least n elements (grown on demand; the queue is idle before a free)
  subroutine ensure_dws(n)
    integer(c_int64_t), intent(in) :: n
    if (allocated(dws)) then
      if (size(dws, kind=c_int64_t) >= n) return
      call om_sync()
      !$omp target exit data map(delete: dws)
      deallocate(dws)
    end if
    allocate(dws(max(n, 1_c_int64_t)))
    !$omp target enter data map(alloc: dws)
  end subroutine ensure_dws

  subroutine ensure_dpiv(n)
    integer, intent(in) :: n
    if (allocated(dpiv)) then
      if (size(dpiv) >= n) return
      call om_sync()
      !$omp target exit data map(delete: dpiv)
      deallocate(dpiv)
    end if
    allocate(dpiv(max(n, 1)))
    !$omp target enter data map(alloc: dpiv)
  end subroutine ensure_dpiv

  subroutine ensure_dscr(n)
    integer, intent(in) :: n
    if (allocated(dscr)) then
      if (size(dscr) >= n) return
      call om_sync()
      !$omp target exit data map(delete: dscr)
      deallocate(dscr)
    end if
    allocate(dscr(max(n, 1)))
    !$omp target enter data map(alloc: dscr)
  end subroutine ensure_dscr

  subroutine ensure_dnrm(n)
    integer, intent(in) :: n
    if (allocated(dnrm)) then
      if (size(dnrm) >= n) return
      call om_sync()
      !$omp target exit data map(delete: dnrm)
      deallocate(dnrm)
    end if
    allocate(dnrm(max(n, 1)))
    !$omp target enter data map(alloc: dnrm)
  end subroutine ensure_dnrm

  !> scratchpad size of a shim query; a failed query stops the run
  subroutine ws_check(rc, what, line)
    integer(c_int), intent(in) :: rc
    character(len=*), intent(in) :: what
    integer, intent(in) :: line
    if (rc /= 0) call om_fail(what//': scratchpad size query failed', line)
  end subroutine ws_check

  !> Presence dispatch: p1(i), pn(i) = first and last element of operand i. .true. if all are on the default device,
  !> .false. if none is (or there is no offload device); a mix stops the run.
  logical function on_device(rname, p1, pn)
    character(len=*), intent(in) :: rname
    type(c_ptr), intent(in) :: p1(:), pn(:)
    integer(c_int) :: dev
    integer :: i, np
    if (ndev < 0) ndev = omp_get_num_devices()
    on_device = .false.
    if (ndev < 1 .or. size(p1) == 0) return
    dev = int(omp_get_default_device(), c_int)
    np = 0
    do i = 1, size(p1)
      if (omp_target_is_present(p1(i), dev) /= 0) np = np + 1
      if (omp_target_is_present(pn(i), dev) /= 0) np = np + 1
    end do
    if (np /= 0 .and. np /= 2*size(p1)) then
      write(error_unit,*) 'linalg_mod (oneMKL): ', rname, ': some operands are on the device, some on the host'
      call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
    end if
    on_device = (np > 0)
  end function on_device

  !> zgemm on the device if a, b, c are there (see the module header)
  logical function om_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    character(len=1), intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(kind=kind(0.d0)), intent(in) :: alpha, beta
    complex(kind=kind(0.d0)), target, intent(in) :: a(lda,*), b(ldb,*)
    complex(kind=kind(0.d0)), target, intent(inout) :: c(ldc,*)
    integer :: na, nb, nc
    type(c_ptr) :: p1(3), pn(3)
    om_zgemm = .false.
    if (m <= 0 .or. n <= 0) return
    if (transa == 'N' .or. transa == 'n') then
      na = ext2(m, k, lda)
    else
      na = ext2(k, m, lda)
    end if
    if (transb == 'N' .or. transb == 'n') then
      nb = ext2(k, n, ldb)
    else
      nb = ext2(n, k, ldb)
    end if
    nc = ext2(m, n, ldc)
    call zends(c, nc, p1(3), pn(3))
    if (na > 0) then
      call zends(a, na, p1(1), pn(1))
      call zends(b, nb, p1(2), pn(2))
      if (.not. on_device('la_zgemm', p1, pn)) return
    else
      if (.not. on_device('la_zgemm', p1(3:3), pn(3:3))) return
    end if
    !$omp dispatch
    call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    om_zgemm = .true.
  end function om_zgemm

  !> C(complex m x n) = A(complex m x k) * B(real k x n) on the device (one dgemm on the real views, as la_zgemm_cr)
  logical function om_zgemm_cr(m, n, k, a, lda, b, ldb, c, ldc)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(kind=kind(0.d0)), target, intent(in) :: a(lda,*)
    real(kind=kind(0.d0)), target, intent(in) :: b(ldb,*)
    complex(kind=kind(0.d0)), target, intent(inout) :: c(ldc,*)
    real(kind=kind(0.d0)), pointer, contiguous :: ar(:,:), cr(:,:)
    real(kind=kind(0.d0)), parameter :: one = 1.d0, zero = 0.d0
    type(c_ptr) :: p1(3), pn(3)
    om_zgemm_cr = .false.
    if (m <= 0 .or. n <= 0 .or. k <= 0) return
    call zends(a, ext2(m, k, lda), p1(1), pn(1))
    call dends(b, ext2(k, n, ldb), p1(2), pn(2))
    call zends(c, ext2(m, n, ldc), p1(3), pn(3))
    if (.not. on_device('la_zgemm_cr', p1, pn)) return
    call c_f_pointer(c_loc(a(1,1)), ar, [2*lda, k])
    call c_f_pointer(c_loc(c(1,1)), cr, [2*ldc, n])
    !$omp dispatch
    call dgemm('N', 'N', 2*m, n, k, one, ar, 2*lda, b, ldb, zero, cr, 2*ldc)
    om_zgemm_cr = .true.
  end function om_zgemm_cr

  !> ztrsm on the device if a and b are there
  logical function om_ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    character(len=1), intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(kind=kind(0.d0)), intent(in) :: alpha
    complex(kind=kind(0.d0)), target, intent(in) :: a(lda,*)
    complex(kind=kind(0.d0)), target, intent(inout) :: b(ldb,*)
    integer :: na
    type(c_ptr) :: p1(2), pn(2)
    om_ztrsm = .false.
    if (m <= 0 .or. n <= 0) return
    if (side == 'L' .or. side == 'l') then
      na = ext2(m, m, lda)
    else
      na = ext2(n, n, lda)
    end if
    call zends(a, na, p1(1), pn(1))
    call zends(b, ext2(m, n, ldb), p1(2), pn(2))
    if (.not. on_device('la_ztrsm', p1, pn)) return
    !$omp dispatch
    call ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    om_ztrsm = .true.
  end function om_ztrsm

  !> zgeqrf on the device if a and tau are there (lwork = -1: host query, see the module header)
  logical function om_zgeqrf(m, n, a, lda, tau, lwork, info)
    integer, intent(in) :: m, n, lda, lwork
    complex(kind=kind(0.d0)), target, intent(inout) :: a(lda,*), tau(*)
    integer, intent(out) :: info
    type(c_ptr) :: p1(2), pn(2)
    om_zgeqrf = .false.
    if (lwork == -1 .or. m <= 0 .or. n <= 0) return
    call zends(a, ext2(m, n, lda), p1(1), pn(1))
    call zends(tau, min(m, n), p1(2), pn(2))
    if (.not. on_device('la_zgeqrf', p1, pn)) return
    call dev_zgeqrf(m, n, a, lda, tau, info)
    om_zgeqrf = .true.
  end function om_zgeqrf

  !> device zgeqrf (a, tau present); info from oneMKL
  subroutine dev_zgeqrf(m, n, a, lda, tau, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target, intent(inout) :: a(lda,*), tau(*)
    integer, intent(out) :: info
    integer(c_int64_t) :: ws
    call om_init()
    ws = 0
    call ws_check(alfsh_zgeqrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws), 'zgeqrf', __LINE__)
    call ensure_dws(ws)
    !$omp target data use_device_addr(a, tau, dws)
    info = alfsh_zgeqrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), c_loc(tau(1)), &
         & c_loc(dws(1)), size(dws, kind=c_int64_t))
    !$omp end target data
  end subroutine dev_zgeqrf

  !> zunmqr on the device if a, tau, c are there
  logical function om_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, lwork, info)
    character(len=1), intent(in) :: side, trans
    integer, intent(in) :: m, n, k, lda, ldc, lwork
    complex(kind=kind(0.d0)), target, intent(in) :: a(lda,*), tau(*)
    complex(kind=kind(0.d0)), target, intent(inout) :: c(ldc,*)
    integer, intent(out) :: info
    type(c_ptr) :: p1(3), pn(3)
    integer(c_int64_t) :: ws
    integer(c_int) :: cs, ct
    om_zunmqr = .false.
    if (lwork == -1 .or. m <= 0 .or. n <= 0 .or. k <= 0) return
    cs = 0
    if (side == 'R' .or. side == 'r') cs = 1
    ct = 0
    if (trans == 'T' .or. trans == 't') ct = 1
    if (trans == 'C' .or. trans == 'c') ct = 2
    if (cs == 0) then
      call zends(a, ext2(m, k, lda), p1(1), pn(1))
    else
      call zends(a, ext2(n, k, lda), p1(1), pn(1))
    end if
    call zends(tau, k, p1(2), pn(2))
    call zends(c, ext2(m, n, ldc), p1(3), pn(3))
    if (.not. on_device('la_zunmqr', p1, pn)) return
    call om_init()
    ws = 0
    call ws_check(alfsh_zunmqr_ws(cs, ct, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), int(lda,c_int64_t), &
         & int(ldc,c_int64_t), ws), 'zunmqr', __LINE__)
    call ensure_dws(ws)
    !$omp target data use_device_addr(a, tau, c, dws)
    info = alfsh_zunmqr(cs, ct, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), &
         & c_loc(tau(1)), c_loc(c(1,1)), int(ldc,c_int64_t), c_loc(dws(1)), size(dws, kind=c_int64_t))
    !$omp end target data
    om_zunmqr = .true.
  end function om_zunmqr

  !> zgetrf on the device if a and ipiv are there; ipiv (int32, device) is filled from the int64 pivots of oneMKL
  logical function om_zgetrf(m, n, a, lda, ipiv, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target, intent(inout) :: a(lda,*)
    integer, target, intent(inout) :: ipiv(*)
    integer, intent(out) :: info
    type(c_ptr) :: p1(2), pn(2)
    integer(c_int64_t) :: ws
    om_zgetrf = .false.
    if (m <= 0 .or. n <= 0) return
    call zends(a, ext2(m, n, lda), p1(1), pn(1))
    call iends(ipiv, min(m, n), p1(2), pn(2))
    if (.not. on_device('la_zgetrf', p1, pn)) return
    call om_init()
    ws = 0
    call ws_check(alfsh_zgetrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws), 'zgetrf', __LINE__)
    call ensure_dws(ws)
    call ensure_dpiv(min(m, n))
    !$omp target data use_device_addr(a, dws, dpiv)
    info = alfsh_zgetrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), c_loc(dpiv(1)), &
         & c_loc(dws(1)), size(dws, kind=c_int64_t))
    !$omp end target data
    call piv_from64(min(m, n), ipiv)
    om_zgetrf = .true.
  end function om_zgetrf

  !> zgetri on the device if a and ipiv are there (lwork = -1: host query)
  logical function om_zgetri(n, a, lda, ipiv, lwork, info)
    integer, intent(in) :: n, lda, lwork
    complex(kind=kind(0.d0)), target, intent(inout) :: a(lda,*)
    integer, target, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    type(c_ptr) :: p1(2), pn(2)
    integer(c_int64_t) :: ws
    om_zgetri = .false.
    if (lwork == -1 .or. n <= 0) return
    call zends(a, ext2(n, n, lda), p1(1), pn(1))
    call iends(ipiv, n, p1(2), pn(2))
    if (.not. on_device('la_zgetri', p1, pn)) return
    call om_init()
    ws = 0
    call ws_check(alfsh_zgetri_ws(int(n,c_int64_t), int(lda,c_int64_t), ws), 'zgetri', __LINE__)
    call ensure_dws(ws)
    call ensure_dpiv(n)
    call piv_to64(n, ipiv)
    !$omp target data use_device_addr(a, dws, dpiv)
    info = alfsh_zgetri(int(n,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), c_loc(dpiv(1)), c_loc(dws(1)), &
         & size(dws, kind=c_int64_t))
    !$omp end target data
    om_zgetri = .true.
  end function om_zgetri

  !> ipiv(1:n) (int32, device) = dpiv(1:n) (int64, device)
  subroutine piv_from64(n, ipiv)
    integer, intent(in) :: n
    integer, intent(inout) :: ipiv(n)
    integer :: i
    !$omp target teams distribute parallel do map(tofrom: ipiv(1:n))
    do i = 1, n
      ipiv(i) = int(dpiv(i))
    end do
  end subroutine piv_from64

  !> dpiv(1:n) (int64, device) = ipiv(1:n) (int32, device)
  subroutine piv_to64(n, ipiv)
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    integer :: i
    !$omp target teams distribute parallel do map(to: ipiv(1:n))
    do i = 1, n
      dpiv(i) = int(ipiv(i), c_int64_t)
    end do
  end subroutine piv_to64

  !> la_zqrpp on the device if a, tau, ipvt are there (m >= n > 0 and lwork /= -1 checked by the caller). Same contract
  !> as the host branch: info = -101 and A unchanged if A has a NaN or Inf entry, else the info of zgeqrf; ipvt valid on
  !> the host and on the device.
  logical function om_zqrpp(m, n, a, lda, ipvt, tau, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target, intent(inout) :: a(lda,*), tau(*)
    integer, target, intent(inout) :: ipvt(*)
    integer, intent(out) :: info
    type(c_ptr) :: p1(3), pn(3)
    integer :: i, j, kk
    om_zqrpp = .false.
    call zends(a, ext2(m, n, lda), p1(1), pn(1))
    call zends(tau, n, p1(2), pn(2))
    call iends(ipvt, n, p1(3), pn(3))
    if (.not. on_device('la_zqrpp', p1, pn)) return
    om_zqrpp = .true.
    call ensure_dnrm(n)
    call ensure_dscr(m*n)
    call k_colnorms(m, n, a, lda, dnrm)
    ! D->H: n column norms for the host sort
    !$omp target update from(dnrm(1:n))
    do j = 1, n
      if (not_finite(dnrm(j))) then
        info = -101
        return
      end if
    end do
    ! stable descending insertion sort of the column norms (as the host branch of la_zqrpp)
    do j = 1, n
      ipvt(j) = j
    end do
    do i = 2, n
      kk = ipvt(i)
      j = i - 1
      do while (j >= 1)
        if (dnrm(ipvt(j)) >= dnrm(kk)) exit
        ipvt(j+1) = ipvt(j)
        j = j - 1
      end do
      ipvt(j+1) = kk
    end do
    ! H->D: the n-entry permutation
    call upd_to_i(n, ipvt)
    call k_permute_cols(m, n, a, lda, ipvt, dscr)
    call dev_zgeqrf(m, n, a, lda, tau, info)
  end function om_zqrpp

  logical function not_finite(x)
    real(kind=kind(0.d0)), intent(in) :: x
    integer(kind=8), parameter :: expmask = int(z'7FF0000000000000', 8)
    not_finite = (iand(transfer(x, 0_8), expmask) == expmask)
  end function not_finite

  subroutine upd_to_i(n, ipvt)
    integer, intent(in) :: n
    integer, intent(in) :: ipvt(n)
    !$omp target update to(ipvt(1:n))
  end subroutine upd_to_i

  !> nrm(j) = 2-norm of column j of the device matrix A (scaled by the largest modulus; NaN if an entry is NaN or Inf)
  subroutine k_colnorms(m, n, a, lda, nrm)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), intent(in) :: a(lda,n)
    real(kind=kind(0.d0)), intent(inout) :: nrm(n)
    integer :: i, j, nf
    real(kind=kind(0.d0)) :: mx, s, v
    integer(kind=8) :: b1, b2
    integer(kind=8), parameter :: expmask = int(z'7FF0000000000000', 8)
    !$omp target teams distribute parallel do private(mx, s, v, nf, b1, b2, i) map(tofrom: nrm(1:n))
    do j = 1, n
      mx = 0.d0
      nf = 0
      do i = 1, m
        b1 = transfer(dble(a(i,j)), 0_8)
        b2 = transfer(aimag(a(i,j)), 0_8)
        if (iand(b1, expmask) == expmask .or. iand(b2, expmask) == expmask) nf = nf + 1
        mx = max(mx, abs(dble(a(i,j))), abs(aimag(a(i,j))))
      end do
      s = 0.d0
      if (nf > 0) then
        nrm(j) = transfer(int(z'7FF8000000000000', 8), 0.d0)
      else if (mx > 0.d0) then
        do i = 1, m
          v = dble(a(i,j)) / mx
          s = s + v*v
          v = aimag(a(i,j)) / mx
          s = s + v*v
        end do
        nrm(j) = mx * sqrt(s)
      else
        nrm(j) = 0.d0
      end if
    end do
  end subroutine k_colnorms

  !> A(:,j) = A(:,ipvt(j)) on the device, through scr (m*n)
  subroutine k_permute_cols(m, n, a, lda, ipvt, scr)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), intent(inout) :: a(lda,n)
    integer, intent(in) :: ipvt(n)
    complex(kind=kind(0.d0)), intent(inout) :: scr(m*n)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2) map(tofrom: scr(1:m*n))
    do j = 1, n
      do i = 1, m
        scr(i + (j-1)*m) = a(i, ipvt(j))
      end do
    end do
    !$omp target teams distribute parallel do collapse(2) map(tofrom: scr(1:m*n))
    do j = 1, n
      do i = 1, m
        a(i,j) = scr(i + (j-1)*m)
      end do
    end do
  end subroutine k_permute_cols

end module linalg_onemkl_mod

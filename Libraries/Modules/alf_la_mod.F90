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
!> @brief Plain-array linear-algebra layer for the GPU port (ADR-002 section 1, ticket 202).
!>
!> Backends, chosen by a preprocessor macro:
!>   (none) / ALF_LA_HOST : the exact BLAS/LAPACK calls ALF makes today (default; GPU=none).
!>   ALF_LA_CUDA          : cuBLAS / cuSOLVER through bind(C) interfaces on DEVICE-RESIDENT arrays.
!>                          Needs -mp=gpu (nvfortran); the caller maps the arrays (!$omp target enter data) and
!>                          alf_la gets the device addresses itself (!$omp target data use_device_addr).
!>   ALF_LA_ONEMKL        : not implemented yet (ticket 302). Compiles (GPU=intel builds), but la_init and every
!>                          HOST code path stop the run (LA_GUARD), so host code never touches device data.
!>
!> Interface rules
!>   * explicit-shape/assumed-size arrays + leading dimensions in LAPACK order, complex(8) throughout
!>     (la_gemm_cr: B real). No descriptors, no derived types. Index arrays are integer(4).
!>     LA_N, LA_T, LA_C have the values of cublasOperation_t (0, 1, 2).
!>   * Scratch memory (LAPACK work, cuSOLVER info) is owned by the module: host work and device work are
!>     grown on demand (la_reserve_work can pre-grow), never allocated per call in steady state. The *_lwork
!>     functions return the optimal work size in complex(8) elements (HOST: LAPACK query, info checked;
!>     CUDA: cuSOLVER bufferSize, a host call that needs no completion of earlier device work).
!>   * NOT thread-safe: one set of handles, scratch buffers and info slots per process; call from one thread.
!>
!> Streams and ordering (CUDA)
!>   * Both library handles are set (la_init) to the stream on which nvfortran's OpenMP offload runtime runs
!>     SYNCHRONOUS target regions, obtained with acc_get_cuda_stream(acc_async_sync) from the nvfortran `openacc` module
!>     (a runtime query only; no OpenACC directives, the code is OpenMP target). Measured on nividia, nvfortran 26.9,
!>     A30, driver 615.71: that call returns the CUDA handle 0x2 = cudaStreamPerThread, so library calls and the caller's
!>     synchronous target regions are in-order on one stream (ADR-002's intent). What did NOT work, each shown by
!>     t_order (RAW/WAR/WAW errors of order 1): a stream from cudaStreamCreate; the legacy NULL stream (passed RAW by
!>     accident only); and the stream acc_get_cuda_stream(acc_async_noval) returns, which is a different,
!>     NON-BLOCKING stream (cudaStreamGetFlags = 1; cudaStreamQuery on it reports "not ready" while a `nowait`
!>     region runs) used by `nowait` regions. Hence regions with `nowait` (and OpenMP async/depend) are NOT ordered
!>     against alf_la calls: la_sync first. `!$omp interop` is not implemented in nvfortran 26.9 (W-0287), so the
!>     interop route is unavailable; the OpenACC runtime query is the substitute.
!>   * Single GPU per process: the OpenACC "current device" queried above must equal the OpenMP default device (it does
!>     with one visible GPU; bind ranks to GPUs with CUDA_VISIBLE_DEVICES, not with omp_set_default_device).
!>     -DALF_LA_DEBUG prints the stream handle once in la_init (0x2 expected).
!>   * -DALF_LA_TEST_BADSTREAM (negative control, tests only): la_init uses a cudaStreamCreate stream instead; the
!>     cudaStreamCreate/Destroy interfaces are kept for this purpose. t_order must FAIL then.
!>   * Safe fallback: compile with -DALF_LA_SYNC_ALWAYS and every library call is bracketed by la_sync.
!>   * Canary: testsuite/alf_la.tests/t_order.F90 (RAW, WAR and WAW hazards between alf_la calls and caller regions,
!>     with no host-blocking read in between). Re-run it for every new compiler/driver version.
!>   * Library calls are asynchronous; la_sync is the only intended blocking call (plus a rare amortised one every
!>     LA_NSLOT solver calls, see "info", and the documented transfers of la_qrpp / la_geqp3_host). alf_la's own
!>     helper kernels and transfers synchronise the stream first. To read results D->H a caller must la_sync first
!>     if it wants failures reported before the data are used.
!>   * info: HOST checks it immediately. CUDA: every solver call writes into its own slot of a device info
!>     array; slots are read back (one D->H copy) by la_sync and when the slots are used up; any nonzero
!>     slot ends the run via Terminate_on_error. Slots are not reset on the device (cuSOLVER overwrites its slot on
!>     every call); islot restarts at 0 after each read-back. Tolerated: getrf info>0 (exactly singular U; ALF ignores it).
!>   * cuBLAS alpha/beta are passed as host pointers to 16-byte-aligned module storage (a complex(8)
!>     parameter/local is only 8-byte aligned and crashed libcublas, ticket 101).
!>
!> Routine -> ALF call sites that the GPU path will use (file:line as of 2026-10-27; verified with grep).
!> Hot-path default (CGR, ...) / STABLOG share these.
!>   la_gemm        ZGEMM  Prog/cgr1_mod.F90:223,586,634 (+108,118,139,140 STAB1 variants); Prog/wrapul_mod.F90:93,99;
!>                         Prog/upgrade_mod.F90:276,280 (op N/T for the submatrix blocks); Prog/Global_mod.F90:873,929;
!>                         Prog/WaveFunction_mod.F90:93; Libraries/Modules/Mat_subroutines_mod.F90:334,338,497,501
!>   la_gemm_cr     DGEMM on the interleaved view: dense B_T (G*E, real E) - ticket 101; no ALF call site yet
!>   la_gemm_batched cublasZgemmStridedBatched: the two flavours of the dense/checkerboard wrap (ticket 101)
!>   la_geru        ZGERU  Prog/upgrade_mod.F90:269 (rank-1 update of G per accepted flip)
!>   la_geqrf/_lwork ZGEQRF Libraries/Modules/mymats_mod.F90:963,973,1076,1085 (UDV); inside ZGEQP3 for Prog/QDRP_decompose_mod.F90:78,82
!>   la_ungqr/_lwork ZUNGQR Prog/udv_state_mod.F90:576; Libraries/Modules/mymats_mod.F90:990,1100
!>   la_ormqr/_lwork ZUNMQR Prog/cgr1_mod.F90:556 ('L','C'), 602 ('R','N'); Prog/cgr2_2_mod.F90:190 ('L','N')
!>   la_trsm        ZTRSM  Prog/cgr1_mod.F90:571 ('L','U','N','N'), 619 ('R','U','C','N'); Prog/cgr2_2_mod.F90:182 ('L','U','C','N')
!>   la_getrf/_lwork ZGETRF Prog/cgr1_mod.F90:124,141,687; Prog/cgr2_2_mod.F90:120; Prog/Global_mod.F90:876
!>   la_getrs       ZGETRS Prog/cgr1_mod.F90:164,697; Prog/cgr2_2_mod.F90:121
!>   la_qrpp        column-norm pre-pivot + ZGEQRF = variant 1 of QR_PREPIVOT (Prog/cgr1_mod.F90:312ff,
!>                  Prog/cgr2_2_mod.F90:405ff and 610ff), instead of ZGEQP3 in Prog/QDRP_decompose_mod.F90:78,82
!>   la_geqp3_host  ZGEQP3 on the host (STABLOG path, ticket 108): QDRP_decompose_mod.F90:78,82
!> Not in alf_la (small element-wise steps: kernels of tickets 203-206): ZLAPMR (cgr1_mod.F90:574), ZLAPMT (:622),
!>   ZSCAL of RHS rows/columns by 1/D or exp(-L) (cgr1_mod.F90:526-540,607-614), the separation of D from R
!>   (QDRP_decompose_mod.F90: |R_ii|, row scaling), the det/phase loops (reflector determinants, permutation sign),
!>   construction of RHS = U^dagger / copy and scaling loops, MMULT-style scalings with D, DZNRM2 of columns on the host
!>   (replaced by the device norm kernel inside la_qrpp).
!>   Not provided: la_gemm_rc (real x complex): no DGEMM view exists for it without transposes (ticket 101 used
!>   transpose kernels); no ALF call site. la_getri: its only call sites are the cold INV routines (mymats_mod.F90:575,623,677).

#if defined(ALF_LA_ONEMKL) && defined(ALF_LA_CUDA)
#error "alf_la_mod: define at most one of ALF_LA_ONEMKL, ALF_LA_CUDA"
#endif
! ALF_LA_ONEMKL (GPU=intel): backend not implemented yet (ticket 302). The module compiles with the HOST code paths so
! that GPU=intel builds; la_init stops the run, so no device data is ever processed by host code.
#if defined(ALF_LA_CUDA) && defined(ALF_LA_HOST)
#error "alf_la_mod: define at most one of ALF_LA_HOST, ALF_LA_CUDA"
#endif

#define LA_FAIL(msg) call la_fail(msg, __LINE__)
! LA_GUARD heads every HOST code path of a public routine: under ALF_LA_ONEMKL (GPU=intel, backend not yet
! implemented) it stops the run, so host BLAS/LAPACK can never touch device-resident data; otherwise it is empty.
#if defined(ALF_LA_ONEMKL)
#define LA_GUARD LA_FAIL('ONEMKL backend not implemented yet (ticket 302)')
#else
#define LA_GUARD
#endif

#if defined(ALF_LA_ONEMKL)
! ---------------------------------------------------------------------------------------------------------------
! ALF_LA_ONEMKL (GPU=intel, ticket 302a): the oneMKL backend, a module of its own in this file so that the HOST and CUDA
! code below keeps its text and its line numbers (every ONEMKL-only insertion is followed by a #line directive).
!   BLAS     : oneMKL OpenMP "!$omp dispatch" variants on the device-resident (mapped) arrays.
!   LAPACK   : geqrf, ungqr, unmqr, getrf, getrs through the C++ SYCL shim (alf_la_sycl_shim.cpp), on the sycl::queue of the
!              OpenMP offload runtime (OpenMP interop, prefer_type("sycl")); the shim waits for every call, so the
!              caller's synchronous target regions, the dispatched BLAS calls and the shim calls are ordered.
!   zgeqp3   : no device implementation in oneMKL: host ZGEQP3 with D->H / H->D copies (la_geqp3_host).
!   info     : every solver call records its status in a slot; slots are checked by la_sync and when used up
!              (same contract as the CUDA backend).
!   Test hooks (environment, read here only): ALF_LA_TEST_BREAK_ORDER, ALF_LA_TEST_POISON_INFO, ALF_LA_TEST_WRONG_LD.
! ---------------------------------------------------------------------------------------------------------------
include "mkl_blas_omp_offload.f90"
module alf_la_om_mod
  use iso_c_binding
  use iso_fortran_env, only: error_unit
  use omp_lib, only: omp_interop_kind, omp_get_interop_ptr, omp_get_interop_int, omp_ipr_fr_id, omp_ipr_targetsync, &
       omp_interop_none, omp_get_default_device, omp_get_num_devices, omp_target_is_present, omp_ifr_sycl
  use onemkl_blas_omp_offload_lp64, only: zgemm, dgemm, zgemm_batch_strided, zgeru, ztrsm
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  implicit none
  private
  integer, parameter :: OM_NSLOT = 4096                      ! = LA_NSLOT
  public :: om_init, om_finalize, om_sync, om_reserve, om_poison_last
  public :: om_gemm, om_gemm_cr, om_gemm_batched, om_geru, om_trsm
  public :: om_geqrf_lwork, om_geqrf, om_ungqr_lwork, om_ungqr, om_ormqr_lwork, om_ormqr
  public :: om_getrf_lwork, om_getrf, om_getrs, om_qrpp, om_geqp3_host

#define OM_FAIL(msg) call om_fail(msg, __LINE__)
  logical, save :: oinited = .false.
  integer(omp_interop_kind), save :: iop = omp_interop_none
  logical, save :: hook_break = .false., hook_poison = .false., hook_wrongld = .false.
  complex(8), allocatable, target, save :: dws(:)              ! device scratch of the shim solvers
  integer, save :: dws_n = 0
  ! dpiv: device int64 pivots (oneMKL USM API); module-owned scratch, mapped by ensure_dpiv, so it needs no guard
  integer(c_int64_t), allocatable, target, save :: dpiv(:)
  integer, save :: dpiv_n = 0
  complex(8), allocatable, target, save :: dscr(:)             ! column-permutation scratch of om_qrpp
  integer, save :: dscr_n = 0
  real(8), allocatable, target, save :: dnrm(:)                ! column norms of om_qrpp
  integer, save :: dnrm_n = 0
  complex(8), allocatable, save :: hws(:)                      ! host work of ZGEQP3
  real(8), allocatable, save :: hrws(:)
  integer, save :: oinfo(OM_NSLOT), islot = 0
  character(len=16), save :: islot_what(OM_NSLOT)
  logical, save :: islot_tol(OM_NSLOT)

  interface
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
    subroutine alfsh_release() bind(C, name='alfsh_release')   ! opens the BREAK_ORDER gate (error exit path)
    end subroutine
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
    integer(c_int) function alfsh_zungqr_ws(m, n, k, lda, ws) bind(C, name='alfsh_zungqr_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, k, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zungqr(m, n, k, a, lda, tau, work, lwork) bind(C, name='alfsh_zungqr')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, work
      integer(c_int64_t), value :: m, n, k, lda, lwork
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
    integer(c_int) function alfsh_zgetrs_ws(trans, n, nrhs, lda, ldb, ws) bind(C, name='alfsh_zgetrs_ws')
      import :: c_int, c_int64_t
      integer(c_int), value :: trans
      integer(c_int64_t), value :: n, nrhs, lda, ldb
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, work, lwork) bind(C, name='alfsh_zgetrs')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, ipiv, b, work
      integer(c_int), value :: trans
      integer(c_int64_t), value :: n, nrhs, lda, ldb, lwork
    end function
    integer(c_int) function alfsh_zgemm_unordered(ta, tb, m, n, k, ar, ai, a, lda, b, ldb, br, bi, c, ldc) &
         bind(C, name='alfsh_zgemm_unordered')
      import :: c_ptr, c_int, c_int64_t, c_double
      integer(c_int), value :: ta, tb
      integer(c_int64_t), value :: m, n, k, lda, ldb, ldc
      real(c_double), value :: ar, ai, br, bi
      type(c_ptr), value :: a, b, c
    end function
  end interface
  external :: zgeqp3
contains

  subroutine om_fail(msg, line)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: line
    write(error_unit,*) 'alf_la: ', msg
    ! the test-hook gate of the unordered GEMMs must not outlive an error exit (a closed gate would hang the runtime
    ! shutdown); the shim also releases it in alfsh_finalize and from an atexit handler
    call alfsh_release()
    call Terminate_on_error(ERROR_GENERIC, 'alf_la_mod.F90[onemkl]', line)
  end subroutine om_fail

  ! ---- host-side presence guard ----------------------------------------------------------------------------------
  ! Every array handed to a dispatched oneMKL call or to a shim call must be device-resident. An unmapped array makes
  ! the runtime pass a HOST address to the device (hard-reset the Arc workstation in T0), so the guard checks the first
  ! and the last element with omp_target_is_present (default device) on the host, before anything is queued.
  ! nel <= 0: nothing is read or written, nothing to check.
  subroutine need_present(p1, pn, rname, aname)
    type(c_ptr), intent(in) :: p1, pn
    character(len=*), intent(in) :: rname, aname
    character(len=160) :: msg
    integer(c_int) :: dev, r1, rn
    dev = int(omp_get_default_device(), c_int)
    r1 = omp_target_is_present(p1, dev)
    rn = omp_target_is_present(pn, dev)
    if (r1 == 0 .or. rn == 0) then
      write(msg,'(5A)') trim(rname), ': argument ', trim(aname), ' is not present on the device, nothing was queued'
      OM_FAIL(trim(msg))
    end if
  end subroutine need_present

  subroutine need_i(a, nel, rname, aname)
    integer(c_int), target :: a(*)
    integer(8), intent(in) :: nel
    character(len=*), intent(in) :: rname, aname
    if (nel <= 0) return
    call need_present(c_loc(a(1)), c_loc(a(nel)), rname, aname)
  end subroutine need_i

  subroutine need_z(a, nel, rname, aname)
    complex(8), target :: a(*)
    integer(8), intent(in) :: nel
    character(len=*), intent(in) :: rname, aname
    if (nel <= 0) return
    call need_present(c_loc(a(1)), c_loc(a(nel)), rname, aname)
  end subroutine need_z

  subroutine need_d(a, nel, rname, aname)
    real(8), target :: a(*)
    integer(8), intent(in) :: nel
    character(len=*), intent(in) :: rname, aname
    if (nel <= 0) return
    call need_present(c_loc(a(1)), c_loc(a(nel)), rname, aname)
  end subroutine need_d

  !> number of elements from the first to the last one touched by a rows x cols block with leading dimension ld
  integer(8) function ext2(rows, cols, ld)
    integer, intent(in) :: rows, cols, ld
    ext2 = 0
    if (rows > 0 .and. cols > 0) ext2 = int(cols-1,8)*int(ld,8) + int(rows,8)
  end function ext2

  !> same for a strided vector of n elements
  integer(8) function extv(n, inc)
    integer, intent(in) :: n, inc
    extv = 0
    if (n > 0) extv = int(n-1,8)*int(abs(inc),8) + 1
  end function extv

  logical function hook_armed(name)
    character(len=*), intent(in) :: name
    character(len=8) :: v
    integer :: st
    v = ''
    call get_environment_variable(name, v, status=st)
    hook_armed = (st == 0 .and. trim(v) == '1')
    if (hook_armed) write(*,'(2A)') ' alf_la test hook armed: ', name
  end function hook_armed

  subroutine om_init()
    type(c_ptr) :: qh
    integer(c_intptr_t) :: fr
    integer(c_int) :: rc, is_l0
    integer :: dev
    character(kind=c_char) :: nm(128)
    character(len=128) :: name
    integer :: i
    complex(8), allocatable, target :: probe(:)
    if (oinited) return
    if (omp_get_num_devices() < 1) OM_FAIL('ONEMKL backend: no OpenMP offload device')
    dev = omp_get_default_device()
    iop = omp_interop_none
    ! the sycl::queue the OpenMP runtime uses for the default device (see the header)
    !$omp interop init(prefer_type("sycl"), targetsync: iop) device(dev)
    if (iop == omp_interop_none) OM_FAIL('ONEMKL backend: OpenMP interop object not available')
    fr = omp_get_interop_int(iop, omp_ipr_fr_id, rc)
    qh = c_null_ptr
    if (rc == 0 .and. fr == omp_ifr_sycl) qh = omp_get_interop_ptr(iop, omp_ipr_targetsync, rc)
    if (.not. c_associated(qh)) OM_FAIL('ONEMKL backend: no sycl::queue from OpenMP interop')
    if (alfsh_init(qh) /= 0) OM_FAIL('ONEMKL backend: shim initialisation failed')
    nm = c_null_char
    is_l0 = 0
    rc = alfsh_device_info(nm, 128_c_int, is_l0)
    if (rc /= 0) OM_FAIL('ONEMKL backend: cannot query the SYCL device')
    name = ''
    do i = 1, 127
      if (nm(i) == c_null_char) exit
      name(i:i) = nm(i)
    end do
    write(*,'(2A)') ' alf_la: ONEMKL device: ', trim(name)
    if (is_l0 /= 1) OM_FAIL('ONEMKL backend: the SYCL queue is not on a Level Zero device')
    ! the OpenMP default device and the queue's device must be the same: a device pointer of the OpenMP device must be
    ! device USM of the queue's device itself (not only of its context: another device or tile of one context is refused)
    allocate(probe(1)); probe = (0.d0, 0.d0)
    !$omp target enter data map(to: probe)
    !$omp target data use_device_addr(probe)
    rc = alfsh_ptr_kind(c_loc(probe(1)))
    !$omp end target data
    !$omp target exit data map(delete: probe)
    if (rc /= 2) OM_FAIL('ONEMKL backend: OpenMP default device and SYCL queue device differ')
    hook_break = hook_armed('ALF_LA_TEST_BREAK_ORDER')
    hook_poison = hook_armed('ALF_LA_TEST_POISON_INFO')
    hook_wrongld = hook_armed('ALF_LA_TEST_WRONG_LD')
    islot = 0
    oinited = .true.
  end subroutine om_init

  subroutine om_finalize()
    if (.not. oinited) return
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
    dws_n = 0; dpiv_n = 0; dscr_n = 0; dnrm_n = 0
    if (allocated(hws)) deallocate(hws)
    if (allocated(hrws)) deallocate(hrws)
    call alfsh_finalize()
    !$omp interop destroy(iop)
    iop = omp_interop_none
    oinited = .false.
  end subroutine om_finalize

  subroutine om_sync()
    if (alfsh_wait() /= 0) OM_FAIL('ONEMKL backend: queue wait failed')
    call check_slots()
  end subroutine om_sync

  subroutine om_reserve(nelem)
    integer, intent(in) :: nelem
    call ensure_dws(nelem)
  end subroutine om_reserve

  ! ---- scratch ---------------------------------------------------------------------------------------------------
  subroutine ensure_dws(n)
    integer, intent(in) :: n
    if (n <= dws_n) return
    call om_sync()
    if (allocated(dws)) then
      !$omp target exit data map(delete: dws)
      deallocate(dws)
    end if
    allocate(dws(max(n,1)))
    !$omp target enter data map(alloc: dws)
    dws_n = max(n,1)
  end subroutine ensure_dws

  subroutine ensure_dpiv(n)
    integer, intent(in) :: n
    if (n <= dpiv_n) return
    call om_sync()
    if (allocated(dpiv)) then
      !$omp target exit data map(delete: dpiv)
      deallocate(dpiv)
    end if
    allocate(dpiv(max(n,1)))
    !$omp target enter data map(alloc: dpiv)
    dpiv_n = max(n,1)
  end subroutine ensure_dpiv

  subroutine ensure_dscr(n)
    integer, intent(in) :: n
    if (n <= dscr_n) return
    call om_sync()
    if (allocated(dscr)) then
      !$omp target exit data map(delete: dscr)
      deallocate(dscr)
    end if
    allocate(dscr(max(n,1)))
    !$omp target enter data map(alloc: dscr)
    dscr_n = max(n,1)
  end subroutine ensure_dscr

  subroutine ensure_dnrm(n)
    integer, intent(in) :: n
    if (n <= dnrm_n) return
    call om_sync()
    if (allocated(dnrm)) then
      !$omp target exit data map(delete: dnrm)
      deallocate(dnrm)
    end if
    allocate(dnrm(max(n,1)))
    !$omp target enter data map(alloc: dnrm)
    dnrm_n = max(n,1)
  end subroutine ensure_dnrm

  ! ---- info slots (host side: the shim returns the status of every call) --------------------------------------------
  integer function next_slot(what, tol)
    character(len=*), intent(in) :: what
    logical, intent(in) :: tol
    if (islot >= OM_NSLOT) call check_slots()
    islot = islot + 1
    islot_what(islot) = what
    islot_tol(islot) = tol
    oinfo(islot) = 0
    next_slot = islot
  end function next_slot

  subroutine check_slots()
    integer :: i
    character(len=80) :: msg
    if (hook_poison) then      ! test hook: the slots are never read back, failures stay unnoticed
      islot = 0
      return
    end if
    do i = 1, islot
      if (oinfo(i) < 0 .or. (oinfo(i) > 0 .and. .not. islot_tol(i))) then
        write(msg,'(3A,I0)') 'oneMKL ', trim(islot_what(i)), ' failed, info = ', oinfo(i)
        OM_FAIL(trim(msg))
      end if
    end do
    islot = 0
  end subroutine check_slots

  subroutine om_poison_last()      ! test hook of t_info: make the most recent slot look failed
    if (islot < 1) return
    oinfo(islot) = -3
  end subroutine om_poison_last

  ! ---- BLAS (dispatch) ---------------------------------------------------------------------------------------------
  character function opc(op)
    integer, intent(in) :: op
    select case (op)
    case (0); opc = 'N'
    case (1); opc = 'T'
    case default; opc = 'C'
    end select
  end function opc

  subroutine guard_gemm(rname, opa, opb, m, n, k, A, lda, B, ldb, C, ldc)
    character(len=*), intent(in) :: rname
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc
    complex(8), target :: A(*), B(*), C(*)
    if (opa == 0) then
      call need_z(A, ext2(m, k, lda), rname, 'A')
    else
      call need_z(A, ext2(k, m, lda), rname, 'A')
    end if
    if (opb == 0) then
      call need_z(B, ext2(k, n, ldb), rname, 'B')
    else
      call need_z(B, ext2(n, k, ldb), rname, 'B')
    end if
    call need_z(C, ext2(m, n, ldc), rname, 'C')
  end subroutine guard_gemm

  subroutine om_gemm(opa, opb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(lda,*), B(ldb,*), C(ldc,*)
    character(len=1) :: ta, tb
    integer :: lda_eff, rc
    ta = opc(opa); tb = opc(opb)
    call guard_gemm('la_gemm', opa, opb, m, n, k, A, lda, B, ldb, C, ldc)
    lda_eff = lda
    if (hook_wrongld .and. opa == 0 .and. lda > m) lda_eff = m     ! test hook: a valid but wrong leading dimension
    if (hook_break) then
      ! test hook: the GEMM goes to a SECOND queue and nobody waits for it (la_sync / la_finalize do)
      !$omp target data use_device_addr(A, B, C)
      rc = alfsh_zgemm_unordered(opa, opb, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), dble(alpha), aimag(alpha), &
           c_loc(A(1,1)), int(lda_eff,c_int64_t), c_loc(B(1,1)), int(ldb,c_int64_t), dble(beta), aimag(beta), &
           c_loc(C(1,1)), int(ldc,c_int64_t))
      !$omp end target data
      if (rc /= 0) OM_FAIL('zgemm (unordered test hook) failed')
      return
    end if
    !$omp dispatch
    call zgemm(ta, tb, m, n, k, alpha, A, lda_eff, B, ldb, beta, C, ldc)
  end subroutine om_gemm

  subroutine om_gemm_cr(m, n, k, A, lda, B, ldb, C, ldc)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(8), target :: A(lda,*), C(ldc,*)
    real(8), target :: B(ldb,*)
    real(8), target :: one, zero
    real(8), pointer, contiguous :: Ar(:,:), Cr(:,:)
    one = 1.d0; zero = 0.d0
    call need_z(A, ext2(m, k, lda), 'la_gemm_cr', 'A')
    call need_d(B, ext2(k, n, ldb), 'la_gemm_cr', 'B')
    call need_z(C, ext2(m, n, ldc), 'la_gemm_cr', 'C')
    ! real views of the complex arrays (contiguous pointers: no copy when passed on)
    call c_f_pointer(c_loc(A(1,1)), Ar, [2*lda, k])
    call c_f_pointer(c_loc(C(1,1)), Cr, [2*ldc, n])
    !$omp dispatch
    call dgemm('N', 'N', 2*m, n, k, one, Ar, 2*lda, B, ldb, zero, Cr, 2*ldc)
  end subroutine om_gemm_cr

  subroutine om_gemm_batched(opa, opb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc, nb, sA, sB, sC
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(*), B(*), C(*)
    character(len=1) :: ta, tb
    ta = opc(opa); tb = opc(opb)
    if (nb > 0) then
      if (opa == 0) then
        call need_z(A, int(nb-1,8)*int(sA,8) + ext2(m, k, lda), 'la_gemm_batched', 'A')
      else
        call need_z(A, int(nb-1,8)*int(sA,8) + ext2(k, m, lda), 'la_gemm_batched', 'A')
      end if
      if (opb == 0) then
        call need_z(B, int(nb-1,8)*int(sB,8) + ext2(k, n, ldb), 'la_gemm_batched', 'B')
      else
        call need_z(B, int(nb-1,8)*int(sB,8) + ext2(n, k, ldb), 'la_gemm_batched', 'B')
      end if
      call need_z(C, int(nb-1,8)*int(sC,8) + ext2(m, n, ldc), 'la_gemm_batched', 'C')
    end if
    !$omp dispatch
    call zgemm_batch_strided(ta, tb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
  end subroutine om_gemm_batched

  subroutine om_geru(m, n, alpha, x, incx, y, incy, A, lda)
    integer, intent(in) :: m, n, incx, incy, lda
    complex(8), intent(in) :: alpha
    complex(8), target :: x(*), y(*), A(lda,*)
    call need_z(x, extv(m, incx), 'la_geru', 'x')
    call need_z(y, extv(n, incy), 'la_geru', 'y')
    call need_z(A, ext2(m, n, lda), 'la_geru', 'A')
    !$omp dispatch
    call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
  end subroutine om_geru

  subroutine om_trsm(side, uplo, opa, diag, m, n, alpha, A, lda, B, ldb)
    character(len=1), intent(in) :: side, uplo, diag
    integer, intent(in) :: opa, m, n, lda, ldb
    complex(8), intent(in) :: alpha
    complex(8), target :: A(lda,*), B(ldb,*)
    character(len=1) :: ta
    ta = opc(opa)
    if (side == 'L' .or. side == 'l') then
      call need_z(A, ext2(m, m, lda), 'la_trsm', 'A')
    else
      call need_z(A, ext2(n, n, lda), 'la_trsm', 'A')
    end if
    call need_z(B, ext2(m, n, ldb), 'la_trsm', 'B')
    !$omp dispatch
    call ztrsm(side, uplo, ta, diag, m, n, alpha, A, lda, B, ldb)
  end subroutine om_trsm

  ! ---- QR family (shim) --------------------------------------------------------------------------------------------
  ! A failed work-size query is recorded like a failed solver call (it is reported by la_sync).
  integer function ws_result(rc, ws, what)
    integer(c_int), intent(in) :: rc
    integer(c_int64_t), intent(in) :: ws
    character(len=*), intent(in) :: what
    integer :: s
    if (rc /= 0) then
      s = next_slot(what, .false.)
      oinfo(s) = -1
      ws_result = 1
    else
      ws_result = int(max(ws, 1_c_int64_t))
    end if
  end function ws_result

  function om_geqrf_lwork(m, n, lda) result(lw)
    integer, intent(in) :: m, n, lda
    integer :: lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    ws = 0
    rc = alfsh_zgeqrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws)
    lw = ws_result(rc, ws, 'geqrf_ws')
  end function om_geqrf_lwork

  subroutine om_geqrf(m, n, A, lda, tau, lw)
    integer, intent(in) :: m, n, lda, lw
    complex(8), target :: A(lda,*), tau(*)
    integer :: s
    call ensure_dws(lw)
    call need_z(A, ext2(m, n, lda), 'la_geqrf', 'A')
    call need_z(tau, int(min(m,n),8), 'la_geqrf', 'tau')
    call need_z(dws, int(lw,8), 'la_geqrf', 'work')
    s = next_slot('geqrf', .false.)
    !$omp target data use_device_addr(A, tau, dws)
    oinfo(s) = alfsh_zgeqrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), c_loc(tau(1)), &
         c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
  end subroutine om_geqrf

  function om_ungqr_lwork(m, n, k, lda) result(lw)
    integer, intent(in) :: m, n, k, lda
    integer :: lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    ws = 0
    rc = alfsh_zungqr_ws(int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), int(lda,c_int64_t), ws)
    lw = ws_result(rc, ws, 'ungqr_ws')
  end function om_ungqr_lwork

  subroutine om_ungqr(m, n, k, A, lda, tau, lw)
    integer, intent(in) :: m, n, k, lda, lw
    complex(8), target :: A(lda,*), tau(*)
    integer :: s
    call ensure_dws(lw)
    call need_z(A, ext2(m, n, lda), 'la_ungqr', 'A')
    call need_z(tau, int(k,8), 'la_ungqr', 'tau')
    call need_z(dws, int(lw,8), 'la_ungqr', 'work')
    s = next_slot('ungqr', .false.)
    !$omp target data use_device_addr(A, tau, dws)
    oinfo(s) = alfsh_zungqr(int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), &
         c_loc(tau(1)), c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
  end subroutine om_ungqr

  function om_ormqr_lwork(side, trans, m, n, k, lda, ldc) result(lw)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc
    integer :: lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc, cs
    cs = 0; if (side == 'R' .or. side == 'r') cs = 1
    ws = 0
    rc = alfsh_zunmqr_ws(cs, int(trans,c_int), int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), int(lda,c_int64_t), &
         int(ldc,c_int64_t), ws)
    lw = ws_result(rc, ws, 'unmqr_ws')
  end function om_ormqr_lwork

  subroutine om_ormqr(side, trans, m, n, k, A, lda, tau, C, ldc, lw)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc, lw
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: s
    integer(c_int) :: cs
    cs = 0; if (side == 'R' .or. side == 'r') cs = 1
    call ensure_dws(lw)
    if (cs == 0) then
      call need_z(A, ext2(m, k, lda), 'la_ormqr', 'A')
    else
      call need_z(A, ext2(n, k, lda), 'la_ormqr', 'A')
    end if
    call need_z(tau, int(k,8), 'la_ormqr', 'tau')
    call need_z(C, ext2(m, n, ldc), 'la_ormqr', 'C')
    call need_z(dws, int(lw,8), 'la_ormqr', 'work')
    s = next_slot('unmqr', .false.)
    !$omp target data use_device_addr(A, tau, C, dws)
    oinfo(s) = alfsh_zunmqr(cs, int(trans,c_int), int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), c_loc(A(1,1)), &
         int(lda,c_int64_t), c_loc(tau(1)), c_loc(C(1,1)), int(ldc,c_int64_t), c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
  end subroutine om_ormqr

  ! ---- LU (shim; ipiv: the caller's int32 device array <-> int64 device scratch) -----------------------------------
  function om_getrf_lwork(m, n, lda) result(lw)
    integer, intent(in) :: m, n, lda
    integer :: lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    ws = 0
    rc = alfsh_zgetrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws)
    lw = ws_result(rc, ws, 'getrf_ws')
  end function om_getrf_lwork

  subroutine om_getrf(m, n, A, lda, ipiv)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer(c_int), target :: ipiv(*)
    integer :: s, lw, i, mn
    mn = min(m, n)
    ! ipiv is checked before anything runs: the int32 <-> int64 conversion regions below would map an absent array
    ! implicitly (a silent temporary mapping)
    call need_i(ipiv, int(mn,8), 'la_getrf', 'ipiv')
    lw = om_getrf_lwork(m, n, lda)
    call ensure_dws(lw)
    call ensure_dpiv(mn)
    call need_z(A, ext2(m, n, lda), 'la_getrf', 'A')
    call need_z(dws, int(lw,8), 'la_getrf', 'work')
    s = next_slot('getrf', .true.)
    !$omp target data use_device_addr(A, dws, dpiv)
    oinfo(s) = alfsh_zgetrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), c_loc(dpiv(1)), &
         c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
    !$omp target teams distribute parallel do map(tofrom: ipiv(1:mn))
    do i = 1, mn
      ipiv(i) = int(dpiv(i), c_int)
    end do
  end subroutine om_getrf

  subroutine om_getrs(trans, n, nrhs, A, lda, ipiv, B, ldb)
    integer, intent(in) :: trans, n, nrhs, lda, ldb
    complex(8), target :: A(lda,*), B(ldb,*)
    integer(c_int), target :: ipiv(*)
    integer :: s, lw, i
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    call need_i(ipiv, int(n,8), 'la_getrs', 'ipiv')      ! before any target region (see om_getrf)
    ws = 0
    rc = alfsh_zgetrs_ws(int(trans,c_int), int(n,c_int64_t), int(nrhs,c_int64_t), int(lda,c_int64_t), int(ldb,c_int64_t), ws)
    lw = ws_result(rc, ws, 'getrs_ws')
    call ensure_dws(lw)
    call ensure_dpiv(n)
    !$omp target teams distribute parallel do map(tofrom: ipiv(1:n))
    do i = 1, n
      dpiv(i) = int(ipiv(i), c_int64_t)
    end do
    call need_z(A, ext2(n, n, lda), 'la_getrs', 'A')
    call need_z(B, ext2(n, nrhs, ldb), 'la_getrs', 'B')
    call need_z(dws, int(lw,8), 'la_getrs', 'work')
    s = next_slot('getrs', .false.)
    !$omp target data use_device_addr(A, B, dws, dpiv)
    oinfo(s) = alfsh_zgetrs(int(trans,c_int), int(n,c_int64_t), int(nrhs,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), &
         c_loc(dpiv(1)), c_loc(B(1,1)), int(ldb,c_int64_t), c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
  end subroutine om_getrs

  ! ---- pre-pivoted QR ----------------------------------------------------------------------------------------------
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

  logical function not_finite(x)
    real(8), intent(in) :: x
    integer(8) :: bits
    bits = transfer(x, 0_8)
    not_finite = (iand(bits, int(z'7FF0000000000000', 8)) == int(z'7FF0000000000000', 8))
  end function not_finite

  !> nrm(j) = 2-norm of column j (scaled; NaN if an entry is NaN/Inf), all device-resident
  subroutine k_colnorms(m, n, A, lda, nrm)
    integer, intent(in) :: m, n, lda
    complex(8), intent(in) :: A(lda,*)
    real(8), intent(inout) :: nrm(*)
    integer :: i, j, nf
    real(8) :: mx, s, v
    integer(8) :: b1, b2, mask
    mask = int(z'7FF0000000000000', 8)
    !$omp target teams distribute parallel do private(mx, s, v, nf, b1, b2, i) map(tofrom: nrm(1:n))
    do j = 1, n
      mx = 0.d0
      nf = 0
      do i = 1, m
        b1 = transfer(dble(A(i,j)), 0_8)
        b2 = transfer(aimag(A(i,j)), 0_8)
        if (iand(b1, mask) == mask .or. iand(b2, mask) == mask) nf = nf + 1
        mx = max(mx, abs(dble(A(i,j))), abs(aimag(A(i,j))))
      end do
      s = 0.d0
      if (nf > 0) then
        nrm(j) = transfer(int(z'7FF8000000000000', 8), 0.d0)
      else if (mx > 0.d0) then
        do i = 1, m
          v = dble(A(i,j)) / mx
          s = s + v*v
          v = aimag(A(i,j)) / mx
          s = s + v*v
        end do
        nrm(j) = mx * sqrt(s)
      else
        nrm(j) = 0.d0
      end if
    end do
  end subroutine k_colnorms

  ! ipvt, A and scr are device-resident: om_qrpp checks ipvt and A on the host before calling this (scr is module scratch)
  subroutine k_permute_cols(m, n, A, lda, ipvt, scr)
    integer, intent(in) :: m, n, lda
    complex(8), intent(inout) :: A(lda,*)
    integer(c_int), intent(in) :: ipvt(*)
    complex(8), intent(inout) :: scr(*)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2) map(tofrom: scr(1:m*n))
    do j = 1, n
      do i = 1, m
        scr(i + (j-1)*m) = A(i, ipvt(j))
      end do
    end do
    !$omp target teams distribute parallel do collapse(2) map(tofrom: scr(1:m*n))
    do j = 1, n
      do i = 1, m
        A(i,j) = scr(i + (j-1)*m)
      end do
    end do
  end subroutine k_permute_cols

  subroutine om_qrpp(m, n, A, lda, ipvt, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer(c_int), target :: ipvt(*)
    integer :: i
    logical :: bad
    ! presence of every array argument first: the target regions below (colnorms, permutation) and the H->D update of
    ! ipvt would otherwise map an absent array implicitly
    call need_z(A, ext2(m, n, lda), 'la_qrpp', 'A')
    call need_z(tau, int(min(m,n),8), 'la_qrpp', 'tau')
    call need_i(ipvt, int(n,8), 'la_qrpp', 'ipvt')
    call ensure_dnrm(n)
    call ensure_dscr(m*n)
    call k_colnorms(m, n, A, lda, dnrm)
    ! D->H: n column norms for the host sort (ticket 202: la_qrpp)
    !$omp target update from(dnrm(1:n))
    bad = .false.
    do i = 1, n
      if (not_finite(dnrm(i))) bad = .true.
    end do
    if (bad) OM_FAIL('la_qrpp: NaN or Inf in the matrix')
    call sort_desc(n, dnrm, ipvt)
    ! H->D: the n-entry permutation (ticket 202: la_qrpp)
    call upd_to_i(ipvt, n)
    call k_permute_cols(m, n, A, lda, ipvt, dscr)
  end subroutine om_qrpp

  subroutine upd_from_z(a, nel)
    integer, intent(in) :: nel
    complex(8) :: a(*)
    !$omp target update from(a(1:nel))
  end subroutine upd_from_z
  subroutine upd_to_z(a, nel)
    integer, intent(in) :: nel
    complex(8) :: a(*)
    !$omp target update to(a(1:nel))
  end subroutine upd_to_z
  subroutine upd_to_i(a, nel)
    integer, intent(in) :: nel
    integer(c_int) :: a(*)
    !$omp target update to(a(1:nel))
  end subroutine upd_to_i

  !> ZGEQP3 on the host with D->H / H->D copies (no QRCP in oneMKL on the device)
  subroutine om_geqp3_host(m, n, A, lda, ipvt, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer(c_int), target :: ipvt(*)
    integer :: lw, info, nel
    complex(8) :: w(1)
    nel = (n-1)*lda + m
    ! the D->H / H->D updates below would not fail on an absent array: check before touching anything
    call need_z(A, int(nel,8), 'la_geqp3_host', 'A')
    call need_z(tau, int(min(m,n),8), 'la_geqp3_host', 'tau')
    call need_i(ipvt, int(n,8), 'la_geqp3_host', 'ipvt')
    call om_sync()
    ! D->H: ZGEQP3 has no device implementation (STABLOG pivoting, ticket 108)
    call upd_from_z(A, nel)
    ! host work arrays grow only when needed (as the HOST and CUDA paths do)
    if (allocated(hrws)) then
      if (size(hrws) < 2*n) deallocate(hrws)
    end if
    if (.not. allocated(hrws)) allocate(hrws(max(2*n,1)))
    ipvt(1:n) = 0
    call zgeqp3(m, n, A, lda, ipvt, tau, w, -1, hrws, info)
    if (info /= 0) OM_FAIL('zgeqp3 work query failed')
    lw = int(dble(w(1)))
    if (allocated(hws)) then
      if (size(hws) < lw) deallocate(hws)
    end if
    if (.not. allocated(hws)) allocate(hws(max(lw,1)))
    call zgeqp3(m, n, A, lda, ipvt, tau, hws, lw, hrws, info)
    if (info /= 0) OM_FAIL('zgeqp3 failed')
    ! H->D: factors, reflectors and permutation back to the device (STABLOG pivoting, ticket 108)
    call upd_to_z(A, nel)
    call upd_to_z(tau, min(m,n))
    call upd_to_i(ipvt, n)
  end subroutine om_geqp3_host

end module alf_la_om_mod
#endif
#line 120
module alf_la_mod
  use iso_c_binding
  use iso_fortran_env, only: error_unit
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
#if defined(ALF_LA_CUDA)
  use openacc, only: acc_get_cuda_stream, acc_async_sync   ! runtime query of the OpenMP stream only
#endif

#if defined(ALF_LA_ONEMKL)
  use alf_la_om_mod
#endif
#line 128
  implicit none
  private
  integer, parameter, public :: LA_N = 0, LA_T = 1, LA_C = 2    ! same values as cublasOperation_t
  integer, parameter, public :: LA_NSLOT = 4096                 ! device info slots (CUDA)
  public :: la_init, la_finalize, la_sync, la_backend_name, la_reserve_work
  public :: la_gemm, la_gemm_cr, la_gemm_batched, la_geru, la_trsm
  public :: la_geqrf_lwork, la_geqrf, la_ungqr_lwork, la_ungqr, la_ormqr_lwork, la_ormqr
  public :: la_getrf_lwork, la_getrf, la_getrs
  public :: la_qrpp, la_geqp3_host
  public :: la_debug_poison_last_info                        ! test hook (t_info)

  logical, save :: inited = .false.
  complex(8), allocatable, save :: hws(:)       ! host LAPACK work (grown on demand)
  real(8), allocatable, save :: hrws(:)         ! host real work (ZGEQP3 rwork, norms)
  complex(8), allocatable, save :: hscr(:)      ! host scratch for the column permutation (HOST la_qrpp)

#if defined(ALF_LA_CUDA)
  type(c_ptr), save :: hblas = c_null_ptr, hsolv = c_null_ptr, hstream = c_null_ptr
  logical, save :: lib_pending = .false.        ! library work enqueued on hstream and not yet waited for
  ! alpha/beta passed to cuBLAS by host pointer; 16-byte aligned (malloc) storage, see header.
  complex(8), allocatable, target, save :: cab(:)
  ! device scratch (all mapped with target enter data in the module, grown on demand)
  complex(8), allocatable, target, save :: dws(:)    ! cuSOLVER workspace
  integer, save :: dws_n = 0
  complex(8), allocatable, target, save :: dscr(:)   ! column-permutation scratch of la_qrpp
  integer, save :: dscr_n = 0
  real(8), allocatable, target, save :: dnrm(:)      ! column norms of la_qrpp
  integer, save :: dnrm_n = 0
  integer(c_int), allocatable, target, save :: dinfo(:)   ! LA_NSLOT cuSOLVER info slots (device)
  integer, save :: islot = 0
  character(len=16), save :: islot_what(LA_NSLOT)
  logical, save :: islot_tol(LA_NSLOT)          ! positive info tolerated (getrf)
  interface
    integer(c_int) function cublasCreate_v2(h) bind(C, name='cublasCreate_v2')
      import :: c_ptr, c_int; type(c_ptr) :: h
    end function
    integer(c_int) function cublasDestroy_v2(h) bind(C, name='cublasDestroy_v2')
      import :: c_ptr, c_int; type(c_ptr), value :: h
    end function
    integer(c_int) function cublasSetStream_v2(h, s) bind(C, name='cublasSetStream_v2')
      import :: c_ptr, c_int; type(c_ptr), value :: h, s
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
      import :: c_ptr, c_int; type(c_ptr), value :: h; integer(c_int), value :: mode   ! mode 1 = CUSOLVER_DETERMINISTIC_RESULTS
    end function
    integer(c_int) function cublasSetAtomicsMode(h, mode) bind(C, name='cublasSetAtomicsMode')
      import :: c_ptr, c_int; type(c_ptr), value :: h; integer(c_int), value :: mode   ! mode 0 = CUBLAS_ATOMICS_NOT_ALLOWED
    end function
    integer(c_int) function cudaStreamCreate(s) bind(C, name='cudaStreamCreate')
      import :: c_ptr, c_int; type(c_ptr) :: s
    end function
    integer(c_int) function cudaStreamDestroy(s) bind(C, name='cudaStreamDestroy')
      import :: c_ptr, c_int; type(c_ptr), value :: s
    end function
    integer(c_int) function cudaStreamSynchronize(s) bind(C, name='cudaStreamSynchronize')
      import :: c_ptr, c_int
      type(c_ptr), value :: s
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
    integer(c_int) function cublasZtrsm_v2(h, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb) &
         bind(C, name='cublasZtrsm_v2')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, alpha, A, B
      integer(c_int), value :: side, uplo, trans, diag, m, n, lda, ldb
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
    integer(c_int) function cusolverDnZgetrf_bufferSize(h, m, n, A, lda, lwork) &
         bind(C, name='cusolverDnZgetrf_bufferSize')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A
      integer(c_int), value :: m, n, lda
      integer(c_int) :: lwork
    end function
    integer(c_int) function cusolverDnZgetrf(h, m, n, A, lda, work, ipiv, info) &
         bind(C, name='cusolverDnZgetrf')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, work, ipiv, info
      integer(c_int), value :: m, n, lda
    end function
    integer(c_int) function cusolverDnZgetrs(h, trans, n, nrhs, A, lda, ipiv, B, ldb, info) &
         bind(C, name='cusolverDnZgetrs')
      import :: c_ptr, c_int
      type(c_ptr), value :: h, A, ipiv, B, info
      integer(c_int), value :: trans, n, nrhs, lda, ldb
    end function
  end interface
#else
  external :: zgemm, dgemm, zgeru, ztrsm, zgeqrf, zungqr, zunmqr, zgetrf, zgetrs, zgeqp3
  real(8), external :: dznrm2
#endif

contains

  function la_backend_name() result(s)
    character(len=16) :: s
#if defined(ALF_LA_ONEMKL)
    s = 'onemkl'
#else
#line 296
#if defined(ALF_LA_CUDA)
    s = 'cuda'
#else
    s = 'host'
#endif
#endif
#line 301
  end function la_backend_name

  subroutine la_fail(msg, line)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: line
    write(error_unit,*) 'alf_la: ', msg
    call Terminate_on_error(ERROR_GENERIC, 'alf_la_mod.F90', line)
  end subroutine la_fail

#if defined(ALF_LA_CUDA)
  subroutine chk(st, what)
    integer(c_int), intent(in) :: st
    character(len=*), intent(in) :: what
    character(len=80) :: msg
    if (st /= 0) then
      write(msg,'(A,I0,2A)') 'library status ', st, ' in ', what
      LA_FAIL(trim(msg))
    end if
  end subroutine chk
#endif

  subroutine la_init()
#if defined(ALF_LA_CUDA)
    if (inited) return
    allocate(cab(2))
    allocate(dinfo(LA_NSLOT))
    dinfo = 0
    !$omp target enter data map(to: dinfo)
#if defined(ALF_LA_TEST_BADSTREAM)
    ! NEGATIVE CONTROL for t_order only (never for production): an ordinary created stream, NOT ordered against the
    ! OpenMP target regions; t_order must fail with this (see testsuite/alf_la.tests/t_order.F90)
    call chk(cudaStreamCreate(hstream), 'cudaStreamCreate')
#else
    hstream = transfer(acc_get_cuda_stream(acc_async_sync), c_null_ptr)    ! stream of synchronous OpenMP regions, see header
#endif
#if defined(ALF_LA_DEBUG)
    write(error_unit,'(A,Z0)') 'alf_la: library stream handle 0x', transfer(hstream, 0_c_intptr_t)
#endif
    call chk(cublasCreate_v2(hblas), 'cublasCreate')
    call chk(cusolverDnCreate(hsolv), 'cusolverDnCreate')
    call chk(cublasSetStream_v2(hblas, hstream), 'cublasSetStream')
    call chk(cusolverDnSetStream(hsolv, hstream), 'cusolverDnSetStream')
    call chk(cusolverDnSetDeterministicMode(hsolv, 1_c_int), 'cusolverDnSetDeterministicMode')   ! ticket 212: bitwise reproducible
    call chk(cublasSetAtomicsMode(hblas, 0_c_int), 'cublasSetAtomicsMode')                       ! library results (no atomics)
    islot = 0
    lib_pending = .false.
#endif
#if defined(ALF_LA_ONEMKL)
    ! ONEMKL backend (ticket 302a): interop queue, oneMKL dispatch + SYCL shim, test hooks (module alf_la_om_mod)
    ! a failed initialisation (no device, queue/device mismatch) ends the run in om_init
    call om_init()
#endif
    inited = .true.
  end subroutine la_init

  subroutine la_finalize()
    if (.not. inited) return
#if defined(ALF_LA_ONEMKL)
    call om_finalize()
#endif
#line 356
#if defined(ALF_LA_CUDA)
    call la_sync()
    call chk(cublasDestroy_v2(hblas), 'cublasDestroy')
    call chk(cusolverDnDestroy(hsolv), 'cusolverDnDestroy')
#if defined(ALF_LA_TEST_BADSTREAM)
    call chk(cudaStreamDestroy(hstream), 'cudaStreamDestroy')
#endif
    !$omp target exit data map(delete: dinfo)
    deallocate(dinfo, cab)
    if (allocated(dws)) then
      !$omp target exit data map(delete: dws)
      deallocate(dws)
    end if
    if (allocated(dscr)) then
      !$omp target exit data map(delete: dscr)
      deallocate(dscr)
    end if
    if (allocated(dnrm)) then
      !$omp target exit data map(delete: dnrm)
      deallocate(dnrm)
    end if
    dws_n = 0; dscr_n = 0; dnrm_n = 0
#endif
    if (allocated(hws)) deallocate(hws)
    if (allocated(hrws)) deallocate(hrws)
    if (allocated(hscr)) deallocate(hscr)
    inited = .false.
  end subroutine la_finalize

  !> Block until all enqueued device work is finished and check the solver info slots.
  subroutine la_sync()
#if defined(ALF_LA_ONEMKL)
    call om_sync()
#endif
#line 387
#if defined(ALF_LA_CUDA)
    call chk(cudaStreamSynchronize(hstream), 'cudaStreamSynchronize')
    lib_pending = .false.
    call check_info_slots()
#endif
  end subroutine la_sync

  !> Pre-grow the work memory (complex(8) elements): host and, with CUDA, device.
  subroutine la_reserve_work(nelem)
    integer, intent(in) :: nelem
#if defined(ALF_LA_ONEMKL)
    call om_reserve(nelem)      ! the host workspace of the HOST/CUDA paths is unused under ONEMKL
#else
#line 397
    call ensure_hws(nelem)
#endif
#line 398
#if defined(ALF_LA_CUDA)
    call ensure_dws(nelem)
#endif
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

#if defined(ALF_LA_CUDA)
  ! ---------------------------------------------------------------------------------------------
  ! CUDA-only helpers
  ! ---------------------------------------------------------------------------------------------
  subroutine pre_lib()
#if defined(ALF_LA_SYNC_ALWAYS)
    call la_sync()
#endif
  end subroutine pre_lib

  subroutine post_lib()
    lib_pending = .true.
#if defined(ALF_LA_SYNC_ALWAYS)
    call la_sync()
#endif
  end subroutine post_lib

  !> Wait for enqueued library work before alf_la's own synchronous target regions / transfers.
  subroutine sync_pend()
    if (lib_pending) then
      call chk(cudaStreamSynchronize(hstream), 'cudaStreamSynchronize')
      lib_pending = .false.
    end if
  end subroutine sync_pend

  subroutine ensure_dws(n)
    integer, intent(in) :: n
    if (n <= dws_n) return
    call la_sync()     ! the old buffer may still be in use
    if (allocated(dws)) then
      !$omp target exit data map(delete: dws)
      deallocate(dws)
    end if
    allocate(dws(max(n,1)))
    !$omp target enter data map(alloc: dws)
    dws_n = max(n,1)
  end subroutine ensure_dws

  subroutine ensure_dscr(n)
    integer, intent(in) :: n
    if (n <= dscr_n) return
    call la_sync()
    if (allocated(dscr)) then
      !$omp target exit data map(delete: dscr)
      deallocate(dscr)
    end if
    allocate(dscr(max(n,1)))
    !$omp target enter data map(alloc: dscr)
    dscr_n = max(n,1)
  end subroutine ensure_dscr

  subroutine ensure_dnrm(n)
    integer, intent(in) :: n
    if (n <= dnrm_n) return
    call la_sync()
    if (allocated(dnrm)) then
      !$omp target exit data map(delete: dnrm)
      deallocate(dnrm)
    end if
    allocate(dnrm(max(n,1)))
    !$omp target enter data map(alloc: dnrm)
    dnrm_n = max(n,1)
  end subroutine ensure_dnrm

  !> Reserve the next info slot (and read all slots back first when they are used up).
  integer function next_slot(what, tol)
    character(len=*), intent(in) :: what
    logical, intent(in) :: tol
    if (islot >= LA_NSLOT) then
      call chk(cudaStreamSynchronize(hstream), 'cudaStreamSynchronize')
      lib_pending = .false.
      call check_info_slots()
    end if
    islot = islot + 1
    islot_what(islot) = what
    islot_tol(islot) = tol
    next_slot = islot
  end function next_slot

  subroutine check_info_slots()
    integer, allocatable :: h(:)
    integer :: i
    character(len=80) :: msg
    if (islot == 0) return
    allocate(h(islot))
    ! D->H: read the cuSOLVER info slots to detect failures (amortised: la_sync or every LA_NSLOT calls)
    !$omp target update from(dinfo(1:islot))
    h = dinfo(1:islot)
    do i = 1, islot
      if (h(i) < 0 .or. (h(i) > 0 .and. .not. islot_tol(i))) then
        write(msg,'(3A,I0)') 'cuSOLVER ', trim(islot_what(i)), ' failed, info = ', h(i)
        LA_FAIL(trim(msg))
      end if
    end do
    islot = 0
    deallocate(h)
  end subroutine check_info_slots

  ! D<->H helpers on flat 1-D views (explicit element counts: exact extent of a sub-block with lda > m)
  subroutine upd_from_z(a, nel)
    integer, intent(in) :: nel
    complex(8) :: a(*)
    !$omp target update from(a(1:nel))
  end subroutine upd_from_z
  subroutine upd_to_z(a, nel)
    integer, intent(in) :: nel
    complex(8) :: a(*)
    !$omp target update to(a(1:nel))
  end subroutine upd_to_z
  subroutine upd_from_d(a, nel)
    integer, intent(in) :: nel
    real(8) :: a(*)
    !$omp target update from(a(1:nel))
  end subroutine upd_from_d
  subroutine upd_to_i(a, nel)
    integer, intent(in) :: nel
    integer(c_int) :: a(*)
    !$omp target update to(a(1:nel))
  end subroutine upd_to_i

  !> nrm(j) = 2-norm of column j of the m x n matrix A (scaled to avoid overflow; a NaN/Inf entry gives NaN,
  !> which la_qrpp detects on the host). All device-resident.
  subroutine k_colnorms(m, n, A, lda, nrm)
    integer, intent(in) :: m, n, lda
    complex(8), intent(in) :: A(lda,*)
    real(8), intent(inout) :: nrm(*)
    integer :: i, j, nf
    real(8) :: mx, s, v
    integer(8) :: b1, b2, mask
    mask = int(z'7FF0000000000000', 8)
    ! ticket 212: one thread per column with serial loops, so the floating-point sum has a fixed order (an OpenMP
    ! device reduction combines in an unspecified order and made near-tied pivots flip from run to run)
    !$omp target teams distribute parallel do private(i, mx, s, v, nf, b1, b2)
    do j = 1, n
      mx = 0.d0
      nf = 0
      do i = 1, m
        b1 = transfer(dble(A(i,j)), 0_8)
        b2 = transfer(aimag(A(i,j)), 0_8)
        if (iand(b1, mask) == mask .or. iand(b2, mask) == mask) nf = nf + 1    ! NaN or Inf entry
        mx = max(mx, abs(dble(A(i,j))), abs(aimag(A(i,j))))
      end do
      s = 0.d0
      if (nf > 0) then
        nrm(j) = transfer(int(z'7FF8000000000000', 8), 0.d0)                  ! quiet NaN: detected on the host
      else if (mx > 0.d0) then
        do i = 1, m
          v = dble(A(i,j)) / mx
          s = s + v*v
          v = aimag(A(i,j)) / mx
          s = s + v*v
        end do
        nrm(j) = mx * sqrt(s)
      else
        nrm(j) = 0.d0
      end if
    end do
  end subroutine k_colnorms

  !> A(:,j) <- A(:,ipvt(j)), j = 1..n, through the scratch array scr(m*n).
  subroutine k_permute_cols(m, n, A, lda, ipvt, scr)
    integer, intent(in) :: m, n, lda
    complex(8), intent(inout) :: A(lda,*)
    integer(c_int), intent(in) :: ipvt(*)
    complex(8), intent(inout) :: scr(*)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, m
        scr(i + (j-1)*m) = A(i, ipvt(j))
      end do
    end do
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, m
        A(i,j) = scr(i + (j-1)*m)
      end do
    end do
  end subroutine k_permute_cols
#endif

  ! ---------------------------------------------------------------------------------------------
  ! BLAS
  ! ---------------------------------------------------------------------------------------------

  !> C = alpha op(A) op(B) + beta C
  subroutine la_gemm(opa, opb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(lda,*), B(ldb,*), C(ldc,*)
    if (m <= 0 .or. n <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_gemm(opa, opb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#else
#line 621
#if defined(ALF_LA_CUDA)
    cab(1) = alpha; cab(2) = beta
    call pre_lib()
    !$omp target data use_device_addr(A, B, C)
    call chk(cublasZgemm_v2(hblas, opa, opb, m, n, k, c_loc(cab(1)), c_loc(A(1,1)), lda, &
         c_loc(B(1,1)), ldb, c_loc(cab(2)), c_loc(C(1,1)), ldc), 'zgemm')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call zgemm(opc(opa), opc(opb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#endif
#endif
#line 633
  end subroutine la_gemm

  !> C(complex m x n) = A(complex m x k) * B(REAL k x n): one DGEMM on the interleaved real view of A and C
  !> (real 2m x k and 2m x n). lda, ldc in complex elements, ldb in real elements.
  subroutine la_gemm_cr(m, n, k, A, lda, B, ldb, C, ldc)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(8), target :: A(lda,*), C(ldc,*)
    real(8), target :: B(ldb,*)
    real(8), target, save :: one = 1.d0, zero = 0.d0
#if !defined(ALF_LA_CUDA)
    real(8), pointer, contiguous :: Ar(:,:), Cr(:,:)
#endif
    if (m <= 0 .or. n <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_gemm_cr(m, n, k, A, lda, B, ldb, C, ldc)
#else
#line 646
#if defined(ALF_LA_CUDA)
    call pre_lib()
    !$omp target data use_device_addr(A, B, C)
    call chk(cublasDgemm_v2(hblas, LA_N, LA_N, 2*m, n, k, c_loc(one), c_loc(A(1,1)), 2*lda, &
         c_loc(B(1,1)), ldb, c_loc(zero), c_loc(C(1,1)), 2*ldc), 'dgemm')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    ! real views of the complex arrays (contiguous pointers: no copy-in/out when passed on)
    call c_f_pointer(c_loc(A(1,1)), Ar, [2*lda, k])
    call c_f_pointer(c_loc(C(1,1)), Cr, [2*ldc, n])
    call dgemm('N', 'N', 2*m, n, k, one, Ar, 2*lda, B, ldb, zero, Cr, 2*ldc)
#endif
#endif
#line 660
  end subroutine la_gemm_cr

  !> nb independent complex GEMMs with constant strides sA, sB, sC (in complex elements)
  subroutine la_gemm_batched(opa, opb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc, nb, sA, sB, sC
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(*), B(*), C(*)
#if !defined(ALF_LA_CUDA)
    integer :: ib
#endif
    if (m <= 0 .or. n <= 0 .or. nb <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_gemm_batched(opa, opb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
#else
#line 671
#if defined(ALF_LA_CUDA)
    cab(1) = alpha; cab(2) = beta
    call pre_lib()
    !$omp target data use_device_addr(A, B, C)
    call chk(cublasZgemmStridedBatched(hblas, opa, opb, m, n, k, c_loc(cab(1)), c_loc(A(1)), lda, &
         int(sA, c_long_long), c_loc(B(1)), ldb, int(sB, c_long_long), c_loc(cab(2)), c_loc(C(1)), ldc, &
         int(sC, c_long_long), nb), 'zgemmStridedBatched')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    do ib = 0, nb - 1
      call zgemm(opc(opa), opc(opb), m, n, k, alpha, A(1+ib*sA), lda, B(1+ib*sB), ldb, beta, C(1+ib*sC), ldc)
    end do
#endif
#endif
#line 686
  end subroutine la_gemm_batched

  !> A += alpha x y^T (no conjugation)
  subroutine la_geru(m, n, alpha, x, incx, y, incy, A, lda)
    integer, intent(in) :: m, n, incx, incy, lda
    complex(8), intent(in) :: alpha
    complex(8), target :: x(*), y(*), A(lda,*)
    if (m <= 0 .or. n <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_geru(m, n, alpha, x, incx, y, incy, A, lda)
#else
#line 694
#if defined(ALF_LA_CUDA)
    cab(1) = alpha
    call pre_lib()
    !$omp target data use_device_addr(x, y, A)
    call chk(cublasZgeru_v2(hblas, m, n, c_loc(cab(1)), c_loc(x(1)), incx, c_loc(y(1)), incy, &
         c_loc(A(1,1)), lda), 'zgeru')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call zgeru(m, n, alpha, x, incx, y, incy, A, lda)
#endif
#endif
#line 706
  end subroutine la_geru

  !> B := alpha op(A)^-1 B (side 'L') or B := alpha B op(A)^-1 (side 'R'); uplo 'U'/'L', diag 'N'/'U'.
  subroutine la_trsm(side, uplo, opa, diag, m, n, alpha, A, lda, B, ldb)
    character(len=1), intent(in) :: side, uplo, diag
    integer, intent(in) :: opa, m, n, lda, ldb
    complex(8), intent(in) :: alpha
    complex(8), target :: A(lda,*), B(ldb,*)
#if defined(ALF_LA_CUDA)
    integer(c_int) :: cs, cu, cd
#endif
    if (m <= 0 .or. n <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_trsm(side, uplo, opa, diag, m, n, alpha, A, lda, B, ldb)
#else
#line 718
#if defined(ALF_LA_CUDA)
    cs = 0; if (side == 'R' .or. side == 'r') cs = 1               ! CUBLAS_SIDE_LEFT/RIGHT
    cu = 0; if (uplo == 'U' .or. uplo == 'u') cu = 1               ! CUBLAS_FILL_MODE_LOWER/UPPER
    cd = 0; if (diag == 'U' .or. diag == 'u') cd = 1               ! CUBLAS_DIAG_NON_UNIT/UNIT
    cab(1) = alpha
    call pre_lib()
    !$omp target data use_device_addr(A, B)
    call chk(cublasZtrsm_v2(hblas, cs, cu, opa, cd, m, n, c_loc(cab(1)), c_loc(A(1,1)), lda, &
         c_loc(B(1,1)), ldb), 'ztrsm')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call ztrsm(side, uplo, opc(opa), diag, m, n, alpha, A, lda, B, ldb)
#endif
#endif
#line 733
  end subroutine la_trsm

  ! ---------------------------------------------------------------------------------------------
  ! LAPACK: QR family.  Work memory is internal; the *_lwork functions return the optimal size in
  ! complex(8) elements (for la_reserve_work) and are always a host-side query.
  ! ---------------------------------------------------------------------------------------------
  function la_geqrf_lwork(m, n, A, lda) result(lw)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer :: lw
#if defined(ALF_LA_ONEMKL)
    lw = om_geqrf_lwork(m, n, lda)
#else
#line 743
#if defined(ALF_LA_CUDA)
    integer(c_int) :: l
    call pre_lib()
    !$omp target data use_device_addr(A)
    call chk(cusolverDnZgeqrf_bufferSize(hsolv, m, n, c_loc(A(1,1)), lda, l), 'geqrf_bufferSize')
    !$omp end target data
    lw = l
#else
    complex(8) :: w(1), t1(1)
    integer :: info
    LA_GUARD
    call zgeqrf(m, n, A, lda, t1, w, -1, info)
    if (info /= 0) LA_FAIL('zgeqrf work query failed')
    lw = int(dble(w(1)))
#endif
#endif
#line 758
  end function la_geqrf_lwork

  subroutine la_geqrf(m, n, A, lda, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
#if defined(ALF_LA_CUDA)
    integer :: s
#else
    integer :: info
#endif
    if (m <= 0 .or. n <= 0) return
    lw = la_geqrf_lwork(m, n, A, lda)
#if defined(ALF_LA_ONEMKL)
    call om_geqrf(m, n, A, lda, tau, lw)
#else
#line 771
#if defined(ALF_LA_CUDA)
    call ensure_dws(max(lw,1))
    s = next_slot('geqrf', .false.)
    call pre_lib()
    !$omp target data use_device_addr(A, tau, dws, dinfo)
    call chk(cusolverDnZgeqrf(hsolv, m, n, c_loc(A(1,1)), lda, c_loc(tau(1)), c_loc(dws(1)), lw, &
         c_loc(dinfo(s))), 'geqrf')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call ensure_hws(lw)
    call zgeqrf(m, n, A, lda, tau, hws, lw, info)
    if (info /= 0) LA_FAIL('zgeqrf failed')
#endif
#endif
#line 786
  end subroutine la_geqrf

  function la_ungqr_lwork(m, n, k, A, lda, tau) result(lw)
    integer, intent(in) :: m, n, k, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
#if defined(ALF_LA_ONEMKL)
    lw = om_ungqr_lwork(m, n, k, lda)
#else
#line 792
#if defined(ALF_LA_CUDA)
    integer(c_int) :: l
    call pre_lib()
    !$omp target data use_device_addr(A, tau)
    call chk(cusolverDnZungqr_bufferSize(hsolv, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), l), &
         'ungqr_bufferSize')
    !$omp end target data
    lw = l
#else
    complex(8) :: w(1)
    integer :: info
    LA_GUARD
    call zungqr(m, n, k, A, lda, tau, w, -1, info)
    if (info /= 0) LA_FAIL('zungqr work query failed')
    lw = int(dble(w(1)))
#endif
#endif
#line 808
  end function la_ungqr_lwork

  !> Overwrite the first n columns of the geqrf factors in A (m x n, k reflectors) by Q (m x n).
  subroutine la_ungqr(m, n, k, A, lda, tau)
    integer, intent(in) :: m, n, k, lda
    complex(8), target :: A(lda,*), tau(*)
    integer :: lw
#if defined(ALF_LA_CUDA)
    integer :: s
#else
    integer :: info
#endif
    if (m <= 0 .or. n <= 0) return
    lw = la_ungqr_lwork(m, n, k, A, lda, tau)
#if defined(ALF_LA_ONEMKL)
    call om_ungqr(m, n, k, A, lda, tau, lw)
#else
#line 822
#if defined(ALF_LA_CUDA)
    call ensure_dws(max(lw,1))
    s = next_slot('ungqr', .false.)
    call pre_lib()
    !$omp target data use_device_addr(A, tau, dws, dinfo)
    call chk(cusolverDnZungqr(hsolv, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), c_loc(dws(1)), lw, &
         c_loc(dinfo(s))), 'ungqr')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call ensure_hws(lw)
    call zungqr(m, n, k, A, lda, tau, hws, lw, info)
    if (info /= 0) LA_FAIL('zungqr failed')
#endif
#endif
#line 837
  end subroutine la_ungqr

  !> side 'L': C := op(Q) C ; side 'R': C := C op(Q); trans LA_N or LA_C; Q from geqrf stored in A (k reflectors).
  function la_ormqr_lwork(side, trans, m, n, k, A, lda, tau, C, ldc) result(lw)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: lw
#if defined(ALF_LA_ONEMKL)
    lw = om_ormqr_lwork(side, trans, m, n, k, lda, ldc)
#else
#line 845
#if defined(ALF_LA_CUDA)
    integer(c_int) :: l, cs
    cs = 0; if (side == 'R' .or. side == 'r') cs = 1
    call pre_lib()
    !$omp target data use_device_addr(A, tau, C)
    call chk(cusolverDnZunmqr_bufferSize(hsolv, cs, trans, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), &
         c_loc(C(1,1)), ldc, l), 'unmqr_bufferSize')
    !$omp end target data
    lw = l
#else
    complex(8) :: w(1)
    integer :: info
    LA_GUARD
    call zunmqr(side, opc(trans), m, n, k, A, lda, tau, C, ldc, w, -1, info)
    if (info /= 0) LA_FAIL('zunmqr work query failed')
    lw = int(dble(w(1)))
#endif
#endif
#line 862
  end function la_ormqr_lwork

  subroutine la_ormqr(side, trans, m, n, k, A, lda, tau, C, ldc)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: lw
#if defined(ALF_LA_CUDA)
    integer :: s
    integer(c_int) :: cs
#else
    integer :: info
#endif
    if (m <= 0 .or. n <= 0 .or. k <= 0) return
    lw = la_ormqr_lwork(side, trans, m, n, k, A, lda, tau, C, ldc)
#if defined(ALF_LA_ONEMKL)
    call om_ormqr(side, trans, m, n, k, A, lda, tau, C, ldc, lw)
#else
#line 877
#if defined(ALF_LA_CUDA)
    cs = 0; if (side == 'R' .or. side == 'r') cs = 1
    call ensure_dws(max(lw,1))
    s = next_slot('unmqr', .false.)
    call pre_lib()
    !$omp target data use_device_addr(A, tau, C, dws, dinfo)
    call chk(cusolverDnZunmqr(hsolv, cs, trans, m, n, k, c_loc(A(1,1)), lda, c_loc(tau(1)), &
         c_loc(C(1,1)), ldc, c_loc(dws(1)), lw, c_loc(dinfo(s))), 'unmqr')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call ensure_hws(lw)
    call zunmqr(side, opc(trans), m, n, k, A, lda, tau, C, ldc, hws, lw, info)
    if (info /= 0) LA_FAIL('zunmqr failed')
#endif
#endif
#line 893
  end subroutine la_ormqr

  ! ---------------------------------------------------------------------------------------------
  ! LAPACK: LU.  ipiv is integer(4), 1-based (LAPACK convention; cuSOLVER uses the same), device-resident with CUDA.
  ! ---------------------------------------------------------------------------------------------
  function la_getrf_lwork(m, n, A, lda) result(lw)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer :: lw
#if defined(ALF_LA_ONEMKL)
    lw = om_getrf_lwork(m, n, lda)
#else
#line 902
#if defined(ALF_LA_CUDA)
    integer(c_int) :: l
    call pre_lib()
    !$omp target data use_device_addr(A)
    call chk(cusolverDnZgetrf_bufferSize(hsolv, m, n, c_loc(A(1,1)), lda, l), 'getrf_bufferSize')
    !$omp end target data
    lw = l
#else
    LA_GUARD
    lw = 0     ! LAPACK ZGETRF needs no work array (nothing to query)
#endif
#endif
#line 913
  end function la_getrf_lwork

  !> LU with partial pivoting. info > 0 (exactly singular U) is tolerated, as in ALF.
  subroutine la_getrf(m, n, A, lda, ipiv)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer(c_int), target :: ipiv(*)
#if defined(ALF_LA_CUDA)
    integer :: lw, s
#else
    integer :: info
#endif
    if (m <= 0 .or. n <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_getrf(m, n, A, lda, ipiv)
#else
#line 926
#if defined(ALF_LA_CUDA)
    lw = la_getrf_lwork(m, n, A, lda)
    call ensure_dws(max(lw,1))
    s = next_slot('getrf', .true.)
    call pre_lib()
    !$omp target data use_device_addr(A, ipiv, dws, dinfo)
    call chk(cusolverDnZgetrf(hsolv, m, n, c_loc(A(1,1)), lda, c_loc(dws(1)), c_loc(ipiv(1)), &
         c_loc(dinfo(s))), 'getrf')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call zgetrf(m, n, A, lda, ipiv, info)
    if (info < 0) LA_FAIL('zgetrf: illegal argument')
#endif
#endif
#line 941
  end subroutine la_getrf

  !> Solve op(A) X = B with the factors of la_getrf. trans LA_N/LA_T/LA_C.
  subroutine la_getrs(trans, n, nrhs, A, lda, ipiv, B, ldb)
    integer, intent(in) :: trans, n, nrhs, lda, ldb
    complex(8), target :: A(lda,*), B(ldb,*)
    integer(c_int), target :: ipiv(*)
#if defined(ALF_LA_CUDA)
    integer :: s
#else
    integer :: info
#endif
    if (n <= 0 .or. nrhs <= 0) return
#if defined(ALF_LA_ONEMKL)
    call om_getrs(trans, n, nrhs, A, lda, ipiv, B, ldb)
#else
#line 954
#if defined(ALF_LA_CUDA)
    s = next_slot('getrs', .false.)
    call pre_lib()
    !$omp target data use_device_addr(A, ipiv, B, dinfo)
    call chk(cusolverDnZgetrs(hsolv, trans, n, nrhs, c_loc(A(1,1)), lda, c_loc(ipiv(1)), &
         c_loc(B(1,1)), ldb, c_loc(dinfo(s))), 'getrs')
    !$omp end target data
    call post_lib()
#else
    LA_GUARD
    call zgetrs(opc(trans), n, nrhs, A, lda, ipiv, B, ldb, info)
    if (info /= 0) LA_FAIL('zgetrs failed')
#endif
#endif
#line 967
  end subroutine la_getrs

  ! ---------------------------------------------------------------------------------------------
  ! Pre-pivoted QR (ticket 104 / 108 variant 1) and host QRCP (STABLOG)
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
  !> Norms: DZNRM2 (HOST) / a scaled device kernel (CUDA); they can differ in the last bits, so columns whose norms
  !> tie to rounding may be ordered differently on the two backends - an equally valid pivot order.
  !> NaN/Inf in A stop the run (Terminate_on_error) after the norm transfer.
  !> CUDA: transfers: n norms D->H and the n-entry permutation H->D (sort on the host, the only unavoidable host step);
  !> ipvt, tau and A are device-resident, ipvt also valid on the host.
  subroutine la_qrpp(m, n, A, lda, ipvt, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer(c_int), target :: ipvt(*)
#if !defined(ALF_LA_CUDA)
    integer :: i, j
#endif
    if (m <= 0 .or. n <= 0) return
    if (m < n) LA_FAIL('la_qrpp requires m >= n')
#if defined(ALF_LA_ONEMKL)
    call om_qrpp(m, n, A, lda, ipvt, tau)
#else
#line 1023
#if defined(ALF_LA_CUDA)
    call ensure_dnrm(n)
    call ensure_dscr(m*n)
    call sync_pend()
    call k_colnorms(m, n, A, lda, dnrm)
    ! D->H: n column norms for the host sort (ticket 202: la_qrpp)
    call upd_from_d(dnrm, n)
    if (any_not_finite(n, dnrm)) LA_FAIL('la_qrpp: NaN or Inf in the matrix')
    call sort_desc(n, dnrm, ipvt)
    ! H->D: the n-entry permutation (ticket 202: la_qrpp)
    call upd_to_i(ipvt, n)
    call k_permute_cols(m, n, A, lda, ipvt, dscr)
#else
    LA_GUARD
    ! direct finite scan of A (as the CUDA kernel does), not relying on DZNRM2 to propagate NaN/Inf
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
#endif
#endif
#line 1064
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

  !> Test hook: make the info slot of the most recent solver call look failed (CUDA), or fail at once (HOST).
  !> Used by t_info to show that a failure is detected without la_sync (slot recycling) and at la_sync.
  subroutine la_debug_poison_last_info()
#if defined(ALF_LA_ONEMKL)
    call om_poison_last()
#else
#line 1080
#if defined(ALF_LA_CUDA)
    if (islot < 1) return
    call sync_pend()
    dinfo(islot) = -3
    !$omp target update to(dinfo(islot:islot))
#else
    LA_FAIL('geqrf failed (poisoned info, test hook)')
#endif
#endif
#line 1088
  end subroutine la_debug_poison_last_info

  !> ZGEQP3 with full pivoting (all columns free) on the HOST, for STABLOG's CGR / CGR2_2 (ticket 108).
  !> CUDA: A is copied D->H, factored on the host, and A, tau, ipvt are copied H->D (documented transfer:
  !> no device QRCP in cuSOLVER; ticket 202/108). ipvt, tau also valid on the host. Bit-identical to ZGEQP3.
  subroutine la_geqp3_host(m, n, A, lda, ipvt, tau)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*), tau(*)
    integer(c_int), target :: ipvt(*)
    integer :: lw, info, nel
    complex(8) :: w(1)
    if (m <= 0 .or. n <= 0) return
    nel = (n-1)*lda + m
#if defined(ALF_LA_ONEMKL)
    call om_geqp3_host(m, n, A, lda, ipvt, tau)
#else
#line 1101
#if defined(ALF_LA_CUDA)
    call sync_pend()
    ! D->H: ZGEQP3 has no device implementation (STABLOG pivoting, ticket 108)
    call upd_from_z(A, nel)
#else
    LA_GUARD
#endif
    call ensure_hrws(2*n)
    ipvt(1:n) = 0
    call zgeqp3(m, n, A, lda, ipvt, tau, w, -1, hrws, info)
    if (info /= 0) LA_FAIL('zgeqp3 work query failed')
    lw = int(dble(w(1)))
    call ensure_hws(lw)
    call zgeqp3(m, n, A, lda, ipvt, tau, hws, lw, hrws, info)
    if (info /= 0) LA_FAIL('zgeqp3 failed')
#if defined(ALF_LA_CUDA)
    ! H->D: factors, reflectors and permutation back to the device (STABLOG pivoting, ticket 108)
    call upd_to_z(A, nel)
    call upd_to_z(tau, min(m,n))
    call upd_to_i(ipvt, n)
#endif
#endif
#line 1122
  end subroutine la_geqp3_host

end module alf_la_mod

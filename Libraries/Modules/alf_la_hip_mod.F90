! ALF_LA_HIP (GPU=amd): HIP backend of alf_la for AMD GPUs.
!   BLAS   : hipBLAS (C shim alf_la_hip_shim.cpp) on the device-resident (mapped) arrays, device pointers via use_device_addr.
!   LAPACK : hipSOLVER geqrf, ungqr, unmqr, getrf, getrs (int32 pivots, used directly on the caller's device array).
!   zgeqp3 : host ZGEQP3 with D->H / H->D copies (la_geqp3_host).
!   Ordering: every shim call ends with hipDeviceSynchronize, because the OpenMP target regions (libgomp GCN plugin) run on
!   queues of their own; solver info is read back inside the call and checked by la_sync (same contract as the other backends).
!   Derived from alf_la_om_mod (oneMKL backend, alf_la_mod.F90); same om_* interface.
#if defined(ALF_LA_HIP)
module alf_la_hip_mod
  use iso_c_binding
  use iso_fortran_env, only: error_unit
  use omp_lib, only: omp_get_default_device, omp_get_num_devices, omp_target_is_present
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
  logical, save :: hook_break = .false., hook_poison = .false., hook_wrongld = .false.
  complex(8), allocatable, target, save :: dws(:)              ! device scratch of the shim solvers
  integer, save :: dws_n = 0
  ! dpiv: device int64 pivots (oneMKL USM API); module-owned scratch, mapped by ensure_dpiv, so it needs no guard
  integer(c_int), allocatable, target, save :: dpiv(:)
  integer, save :: dpiv_n = 0
  complex(8), allocatable, target, save :: dscr(:)             ! column-permutation scratch of om_qrpp
  integer, save :: dscr_n = 0
  real(8), allocatable, target, save :: dnrm(:)                ! column norms of om_qrpp
  integer, save :: dnrm_n = 0
  complex(8), allocatable, save :: hws(:)                      ! host work of ZGEQP3
  real(8), allocatable, save :: hrws(:)
  integer, save :: oinfo(OM_NSLOT), islot = 0
  logical, save :: use_hq = .true.                         ! blocked HIP Householder QR / Q-apply (alf_hip_qr.hip.cpp); ALF_HIP_QR=0 -> hipSOLVER
  character(len=16), save :: islot_what(OM_NSLOT)
  logical, save :: islot_tol(OM_NSLOT)

  interface
    integer(c_int) function alfhp_init(dev, name, nname) bind(C, name='alfhp_init')
      import :: c_char, c_int
      integer(c_int), value :: dev, nname
      character(kind=c_char) :: name(*)
    end function
    subroutine alfhp_finalize() bind(C, name='alfhp_finalize')
    end subroutine
    integer(c_int) function alfhp_wait() bind(C, name='alfhp_wait')
      import :: c_int
    end function
    integer(c_int) function alfhp_zgemm(ta, tb, m, n, k, ar, ai, a, lda, b, ldb, br, bi, c, ldc) bind(C, name='alfhp_zgemm')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
      real(c_double), value :: ar, ai, br, bi
      type(c_ptr), value :: a, b, c
    end function
    integer(c_int) function alfhp_dgemm(m, n, k, a, lda, b, ldb, c, ldc) bind(C, name='alfhp_dgemm')
      import :: c_ptr, c_int
      integer(c_int), value :: m, n, k, lda, ldb, ldc
      type(c_ptr), value :: a, b, c
    end function
    integer(c_int) function alfhp_zgemm_batched(ta, tb, m, n, k, ar, ai, a, lda, sa, b, ldb, sb, br, bi, c, ldc, sc, nb) &
         bind(C, name='alfhp_zgemm_batched')
      import :: c_ptr, c_int, c_int64_t, c_double
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc, nb
      integer(c_int64_t), value :: sa, sb, sc
      real(c_double), value :: ar, ai, br, bi
      type(c_ptr), value :: a, b, c
    end function
    integer(c_int) function alfhp_zgeru(m, n, ar, ai, x, incx, y, incy, a, lda) bind(C, name='alfhp_zgeru')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: m, n, incx, incy, lda
      real(c_double), value :: ar, ai
      type(c_ptr), value :: x, y, a
    end function
    integer(c_int) function alfhp_ztrsm(side, uplo, op, diag, m, n, ar, ai, a, lda, b, ldb) bind(C, name='alfhp_ztrsm')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: side, uplo, diag, op, m, n, lda, ldb
      real(c_double), value :: ar, ai
      type(c_ptr), value :: a, b
    end function
    integer(c_int) function alfhp_zgeqrf_ws(m, n, lda, ws) bind(C, name='alfhp_zgeqrf_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfhp_zgeqrf(m, n, a, lda, tau, work, lwork) bind(C, name='alfhp_zgeqrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
    integer(c_int) function alfhp_zungqr_ws(m, n, k, lda, ws) bind(C, name='alfhp_zungqr_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, k, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfhp_zungqr(m, n, k, a, lda, tau, work, lwork) bind(C, name='alfhp_zungqr')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, work
      integer(c_int64_t), value :: m, n, k, lda, lwork
    end function
    integer(c_int) function alfhp_zunmqr_ws(side, trans, m, n, k, lda, ldc, ws) bind(C, name='alfhp_zunmqr_ws')
      import :: c_int, c_int64_t
      integer(c_int), value :: side, trans
      integer(c_int64_t), value :: m, n, k, lda, ldc
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfhp_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork) bind(C, name='alfhp_zunmqr')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, c, work
      integer(c_int), value :: side, trans
      integer(c_int64_t), value :: m, n, k, lda, ldc, lwork
    end function
    integer(c_int) function alfhp_zgetrf_ws(m, n, lda, ws) bind(C, name='alfhp_zgetrf_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfhp_zgetrf(m, n, a, lda, ipiv, work, lwork) bind(C, name='alfhp_zgetrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, ipiv, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
    integer(c_int) function alfhq_geqrf(m, n, a, lda, tau) bind(C, name='alfhq_geqrf')
      import :: c_ptr, c_int
      type(c_ptr), value :: a, tau
      integer(c_int), value :: m, n, lda
    end function
    integer(c_int) function alfhq_unmqr(side, trans, m, n, k, a, lda, tau, c, ldc) bind(C, name='alfhq_unmqr')
      import :: c_ptr, c_int
      type(c_ptr), value :: a, tau, c
      integer(c_int), value :: side, trans, m, n, k, lda, ldc
    end function
    integer(c_int) function alfhp_zgetrs_ws(trans, n, nrhs, lda, ldb, ws) bind(C, name='alfhp_zgetrs_ws')
      import :: c_int, c_int64_t
      integer(c_int), value :: trans
      integer(c_int64_t), value :: n, nrhs, lda, ldb
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfhp_zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, work, lwork) bind(C, name='alfhp_zgetrs')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, ipiv, b, work
      integer(c_int), value :: trans
      integer(c_int64_t), value :: n, nrhs, lda, ldb, lwork
    end function
  end interface
  external :: zgeqp3
contains

  subroutine om_fail(msg, line)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: line
    write(error_unit,*) 'alf_la: ', msg
    ! the test-hook gate of the unordered GEMMs must not outlive an error exit (a closed gate would hang the runtime
    ! shutdown); the shim also releases it in alfhp_finalize and from an atexit handler
    call Terminate_on_error(ERROR_GENERIC, 'alf_la_mod.F90[hip]', line)
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
    integer :: dev, rc, i
    character(kind=c_char) :: nm(128)
    character(len=128) :: name
    if (oinited) return
    if (omp_get_num_devices() < 1) OM_FAIL('HIP backend: no OpenMP offload device')
    dev = omp_get_default_device()
    nm = c_null_char
    rc = alfhp_init(int(dev, c_int), nm, 128_c_int)
    if (rc /= 0) then
      write(error_unit,*) 'alf_la: alfhp_init failed, code ', rc, ' device ', dev, ' of ', omp_get_num_devices()
      OM_FAIL('HIP backend: initialisation failed')
    end if
    name = ''
    do i = 1, 127
      if (nm(i) == c_null_char) exit
      name(i:i) = nm(i)
    end do
    write(*,'(2A)') ' alf_la: HIP device: ', trim(name)
    islot = 0
    block
      character(len=8) :: ev
      integer :: evl, evs
      call get_environment_variable('ALF_HIP_QR', ev, evl, evs)
      if (evs == 0 .and. evl >= 1) use_hq = (ev(1:1) /= '0')
    end block
    oinited = .true.
  end subroutine om_init

  subroutine om_finalize()
    if (.not. oinited) return
    call om_sync()
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
    dws_n = 0; dpiv_n = 0; dscr_n = 0; dnrm_n = 0
    if (allocated(hws)) deallocate(hws)
    if (allocated(hrws)) deallocate(hrws)
    call alfhp_finalize()
    oinited = .false.
  end subroutine om_finalize

  subroutine om_sync()
    if (alfhp_wait() /= 0) OM_FAIL('HIP backend: device synchronisation failed')
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
        write(msg,'(3A,I0)') 'hipSOLVER ', trim(islot_what(i)), ' failed, info = ', oinfo(i)
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
    integer :: rc
    call guard_gemm('la_gemm', opa, opb, m, n, k, A, lda, B, ldb, C, ldc)
    !$omp target data use_device_addr(A, B, C)
    rc = alfhp_zgemm(int(opa,c_int), int(opb,c_int), int(m,c_int), int(n,c_int), int(k,c_int), dble(alpha), aimag(alpha), &
         c_loc(A(1,1)), int(lda,c_int), c_loc(B(1,1)), int(ldb,c_int), dble(beta), aimag(beta), c_loc(C(1,1)), int(ldc,c_int))
    !$omp end target data
    if (rc /= 0) OM_FAIL('hipblasZgemm failed')
  end subroutine om_gemm

  subroutine om_gemm_cr(m, n, k, A, lda, B, ldb, C, ldc)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(8), target :: A(lda,*), C(ldc,*)
    real(8), target :: B(ldb,*)
    integer :: rc
    call need_z(A, ext2(m, k, lda), 'la_gemm_cr', 'A')
    call need_d(B, ext2(k, n, ldb), 'la_gemm_cr', 'B')
    call need_z(C, ext2(m, n, ldc), 'la_gemm_cr', 'C')
    ! real views of the complex arrays: C = A_real_view(2m x k) * B(k x n)
    !$omp target data use_device_addr(A, B, C)
    rc = alfhp_dgemm(int(2*m,c_int), int(n,c_int), int(k,c_int), c_loc(A(1,1)), int(2*lda,c_int), c_loc(B(1,1)), int(ldb,c_int), &
         c_loc(C(1,1)), int(2*ldc,c_int))
    !$omp end target data
    if (rc /= 0) OM_FAIL('hipblasDgemm failed')
  end subroutine om_gemm_cr

  subroutine om_gemm_batched(opa, opb, m, n, k, alpha, A, lda, sA, B, ldb, sB, beta, C, ldc, sC, nb)
    integer, intent(in) :: opa, opb, m, n, k, lda, ldb, ldc, nb, sA, sB, sC
    complex(8), intent(in) :: alpha, beta
    complex(8), target :: A(*), B(*), C(*)
    integer :: rc
    if (nb <= 0) return
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
    !$omp target data use_device_addr(A, B, C)
    rc = alfhp_zgemm_batched(int(opa,c_int), int(opb,c_int), int(m,c_int), int(n,c_int), int(k,c_int), dble(alpha), aimag(alpha), &
         c_loc(A(1)), int(lda,c_int), int(sA,c_int64_t), c_loc(B(1)), int(ldb,c_int), int(sB,c_int64_t), dble(beta), aimag(beta), &
         c_loc(C(1)), int(ldc,c_int), int(sC,c_int64_t), int(nb,c_int))
    !$omp end target data
    if (rc /= 0) OM_FAIL('hipblasZgemmStridedBatched failed')
  end subroutine om_gemm_batched

  subroutine om_geru(m, n, alpha, x, incx, y, incy, A, lda)
    integer, intent(in) :: m, n, incx, incy, lda
    complex(8), intent(in) :: alpha
    complex(8), target :: x(*), y(*), A(lda,*)
    integer :: rc
    call need_z(x, extv(m, incx), 'la_geru', 'x')
    call need_z(y, extv(n, incy), 'la_geru', 'y')
    call need_z(A, ext2(m, n, lda), 'la_geru', 'A')
    !$omp target data use_device_addr(x, y, A)
    rc = alfhp_zgeru(int(m,c_int), int(n,c_int), dble(alpha), aimag(alpha), c_loc(x(1)), int(incx,c_int), c_loc(y(1)), &
         int(incy,c_int), c_loc(A(1,1)), int(lda,c_int))
    !$omp end target data
    if (rc /= 0) OM_FAIL('hipblasZgeru failed')
  end subroutine om_geru

  subroutine om_trsm(side, uplo, opa, diag, m, n, alpha, A, lda, B, ldb)
    character(len=1), intent(in) :: side, uplo, diag
    integer, intent(in) :: opa, m, n, lda, ldb
    complex(8), intent(in) :: alpha
    complex(8), target :: A(lda,*), B(ldb,*)
    integer :: rc
    integer(c_int) :: cs, cu, cd
    if (side == 'L' .or. side == 'l') then
      call need_z(A, ext2(m, m, lda), 'la_trsm', 'A')
    else
      call need_z(A, ext2(n, n, lda), 'la_trsm', 'A')
    end if
    call need_z(B, ext2(m, n, ldb), 'la_trsm', 'B')
    !$omp target data use_device_addr(A, B)
    cs = 0; if (side == 'R' .or. side == 'r') cs = 1          ! 0 left, 1 right
    cu = 0; if (uplo == 'U' .or. uplo == 'u') cu = 1          ! 0 lower, 1 upper
    cd = 0; if (diag == 'U' .or. diag == 'u') cd = 1          ! 0 non-unit, 1 unit
    rc = alfhp_ztrsm(cs, cu, int(opa,c_int), cd, int(m,c_int), int(n,c_int), dble(alpha), aimag(alpha), c_loc(A(1,1)), &
         int(lda,c_int), c_loc(B(1,1)), int(ldb,c_int))
    !$omp end target data
    if (rc /= 0) OM_FAIL('hipblasZtrsm failed')
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
    rc = alfhp_zgeqrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws)
    lw = ws_result(rc, ws, 'geqrf_ws')
  end function om_geqrf_lwork

  subroutine om_geqrf(m, n, A, lda, tau, lw)
    integer, intent(in) :: m, n, lda, lw
    complex(8), target :: A(lda,*), tau(*)
    integer :: s, rc
    call ensure_dws(lw)
    call need_z(A, ext2(m, n, lda), 'la_geqrf', 'A')
    call need_z(tau, int(min(m,n),8), 'la_geqrf', 'tau')
    call need_z(dws, int(lw,8), 'la_geqrf', 'work')
    s = next_slot('geqrf', .false.)
    rc = 3
    if (use_hq) then
      !$omp target data use_device_addr(A, tau)
      rc = alfhq_geqrf(int(m,c_int), int(n,c_int), c_loc(A(1,1)), int(lda,c_int), c_loc(tau(1)))
      !$omp end target data
    end if
    if (rc == 3) then       ! size not supported by the HIP QR (m > 2304) or disabled: hipSOLVER
      !$omp target data use_device_addr(A, tau, dws)
      oinfo(s) = alfhp_zgeqrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), c_loc(tau(1)), &
           c_loc(dws(1)), int(lw,c_int64_t))
      !$omp end target data
    else
      oinfo(s) = rc
    end if
  end subroutine om_geqrf

  function om_ungqr_lwork(m, n, k, lda) result(lw)
    integer, intent(in) :: m, n, k, lda
    integer :: lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    ws = 0
    rc = alfhp_zungqr_ws(int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), int(lda,c_int64_t), ws)
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
    oinfo(s) = alfhp_zungqr(int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), &
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
    rc = alfhp_zunmqr_ws(cs, int(trans,c_int), int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), int(lda,c_int64_t), &
         int(ldc,c_int64_t), ws)
    lw = ws_result(rc, ws, 'unmqr_ws')
  end function om_ormqr_lwork

  subroutine om_ormqr(side, trans, m, n, k, A, lda, tau, C, ldc, lw)
    character(len=1), intent(in) :: side
    integer, intent(in) :: trans, m, n, k, lda, ldc, lw
    complex(8), target :: A(lda,*), tau(*), C(ldc,*)
    integer :: s
    integer(c_int) :: cs, rc
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
    rc = 3
    if (use_hq .and. trans /= 1) then      ! trans = 1 (plain transpose) is not a unitary operation: hipSOLVER rejects it
      !$omp target data use_device_addr(A, tau, C)
      rc = alfhq_unmqr(cs, int(trans,c_int), int(m,c_int), int(n,c_int), int(k,c_int), c_loc(A(1,1)), int(lda,c_int), &
           c_loc(tau(1)), c_loc(C(1,1)), int(ldc,c_int))
      !$omp end target data
    end if
    if (rc == 3) then
      !$omp target data use_device_addr(A, tau, C, dws)
      oinfo(s) = alfhp_zunmqr(cs, int(trans,c_int), int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), c_loc(A(1,1)), &
           int(lda,c_int64_t), c_loc(tau(1)), c_loc(C(1,1)), int(ldc,c_int64_t), c_loc(dws(1)), int(lw,c_int64_t))
      !$omp end target data
    else
      oinfo(s) = rc
    end if
  end subroutine om_ormqr

  ! ---- LU (shim; ipiv: the caller's int32 device array <-> int64 device scratch) -----------------------------------
  function om_getrf_lwork(m, n, lda) result(lw)
    integer, intent(in) :: m, n, lda
    integer :: lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    ws = 0
    rc = alfhp_zgetrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws)
    lw = ws_result(rc, ws, 'getrf_ws')
  end function om_getrf_lwork

  subroutine om_getrf(m, n, A, lda, ipiv)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(lda,*)
    integer(c_int), target :: ipiv(*)
    integer :: s, lw, mn
    mn = min(m, n)
    call need_i(ipiv, int(mn,8), 'la_getrf', 'ipiv')
    lw = om_getrf_lwork(m, n, lda)
    call ensure_dws(lw)
    call need_z(A, ext2(m, n, lda), 'la_getrf', 'A')
    call need_z(dws, int(lw,8), 'la_getrf', 'work')
    s = next_slot('getrf', .true.)
    !$omp target data use_device_addr(A, dws, ipiv)
    oinfo(s) = alfhp_zgetrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), c_loc(ipiv(1)), &
         c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
  end subroutine om_getrf

  subroutine om_getrs(trans, n, nrhs, A, lda, ipiv, B, ldb)
    integer, intent(in) :: trans, n, nrhs, lda, ldb
    complex(8), target :: A(lda,*), B(ldb,*)
    integer(c_int), target :: ipiv(*)
    integer :: s, lw
    integer(c_int64_t) :: ws
    integer(c_int) :: rc
    call need_i(ipiv, int(n,8), 'la_getrs', 'ipiv')
    ws = 0
    rc = alfhp_zgetrs_ws(int(trans,c_int), int(n,c_int64_t), int(nrhs,c_int64_t), int(lda,c_int64_t), int(ldb,c_int64_t), ws)
    lw = ws_result(rc, ws, 'getrs_ws')
    call ensure_dws(lw)
    call need_z(A, ext2(n, n, lda), 'la_getrs', 'A')
    call need_z(B, ext2(n, nrhs, ldb), 'la_getrs', 'B')
    call need_z(dws, int(lw,8), 'la_getrs', 'work')
    s = next_slot('getrs', .false.)
    !$omp target data use_device_addr(A, B, dws, ipiv)
    oinfo(s) = alfhp_zgetrs(int(trans,c_int), int(n,c_int64_t), int(nrhs,c_int64_t), c_loc(A(1,1)), int(lda,c_int64_t), &
         c_loc(ipiv(1)), c_loc(B(1,1)), int(ldb,c_int64_t), c_loc(dws(1)), int(lw,c_int64_t))
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
    complex(8), intent(in) :: A(lda,n)
    real(8), intent(inout) :: nrm(n)
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
    complex(8), intent(inout) :: A(lda,n)
    integer(c_int), intent(in) :: ipvt(n)
    complex(8), intent(inout) :: scr(m*n)
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

end module alf_la_hip_mod
#endif

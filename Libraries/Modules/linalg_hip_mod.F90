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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> GPU=amd (ALF_LA_HIP): the HIP device path of linalg_mod (hipBLAS, hipSOLVER, own HIP QR) on device-resident arrays.
!
!> @details
!> linalg_mod calls lh_dev first (presence dispatch, P2.4 decision Q1): all matrix operands present on the default
!> device -> the lh_z* routine here; none present -> the host BLAS/LAPACK call; mixed -> the run stops.  Device addresses
!> come from use_device_addr.  Ported from gpu-port alf_la_hip_mod (om_*), with these changes: info is returned to the
!> caller (hipSOLVER calls block and read their info back), lh_zqrpp includes the ZGEQRF (the contract of la_zqrpp), and
!> lh_zgetri is ZGETRS on the identity (hipSOLVER has no getri).
!>  - Ordering: the OpenMP target regions (libgomp GCN plugin) run on their own queues.  hipSOLVER calls synchronise before
!>    they return.  hipBLAS calls and the own HIP kernels/QR only mark the work as pending when the patched plugin has the
!>    presync hook (deferred synchronisation, ALF_HIP_ASYNC=0 switches it off); the plugin then synchronises before every
!>    OpenMP launch and device copy (alf_hip_kernels.hip.cpp).
!>  - QR: own blocked Householder QR / Q-apply for m <= 2304 (alf_hip_qr.hip.cpp; LAPACK reflector convention, so the
!>    tracked det(U) of gpu_stab_mod holds); larger m, trans = 'T', or ALF_HIP_QR=0: hipSOLVER.
!>  - Device scratch (solver work, norms) is owned here and grows on demand; the caller's work/lwork are not used.
!--------------------------------------------------------------------
#if defined(ALF_LA_HIP)
module linalg_hip_mod
  use iso_c_binding, only: c_int, c_int64_t, c_double, c_char, c_null_char, c_ptr, c_loc
  use iso_fortran_env, only: error_unit
  use omp_lib, only: omp_get_default_device, omp_get_num_devices, omp_target_is_present
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  use alf_hipk_mod, only: hk_copyz, hk_ident, hk_colnorms, hk_permute_cols
  implicit none
  private
  public :: lh_pz, lh_pd, lh_pi, lh_dev
  public :: lh_zgemm, lh_zgemm_cr, lh_ztrsm, lh_zgeqrf, lh_zunmqr, lh_zqrpp, lh_zgetrf, lh_zgetri

  logical, save :: inited = .false.
  logical, save :: use_hq = .true.                             ! own HIP QR (ALF_HIP_QR=0 -> hipSOLVER)
  complex(kind=kind(0.d0)), allocatable, target, save :: dws(:)   ! device work of the hipSOLVER calls
  complex(kind=kind(0.d0)), allocatable, target, save :: dscr(:)  ! device n x n scratch of lh_zgetri
  real(kind=kind(0.d0)), allocatable, target, save :: dnrm(:)     ! column norms of lh_zqrpp

  interface
    integer(c_int) function alfhp_init(dev, name, nname) bind(C, name='alfhp_init')
      import :: c_char, c_int
      integer(c_int), value :: dev, nname
      character(kind=c_char) :: name(*)
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
  end interface

contains

  subroutine lh_fail(msg, line)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: line
    write(error_unit,'(2A)') ' linalg (HIP): ', msg
    call Terminate_on_error(ERROR_GENERIC, __FILE__, line)
  end subroutine lh_fail

  !> .true. if a(1) is present on the default device
  logical function lh_pz(a)
    complex(kind=kind(0.d0)), intent(in), target :: a(*)
    lh_pz = omp_target_is_present(c_loc(a(1)), omp_get_default_device()) /= 0
  end function lh_pz

  !> .true. if a(1) is present on the default device (real array)
  logical function lh_pd(a)
    real(kind=kind(0.d0)), intent(in), target :: a(*)
    lh_pd = omp_target_is_present(c_loc(a(1)), omp_get_default_device()) /= 0
  end function lh_pd

  !> .true. if a(1) is present on the default device (integer array)
  logical function lh_pi(a)
    integer, intent(in), target :: a(*)
    lh_pi = omp_target_is_present(c_loc(a(1)), omp_get_default_device()) /= 0
  end function lh_pi

  !> Presence dispatch: .true. if every operand is on the device, .false. if none is; a mix stops the run (host BLAS on
  !> device data or the reverse would give silently wrong results).
  logical function lh_dev(p, rname)
    logical, intent(in) :: p(:)
    character(len=*), intent(in) :: rname
    lh_dev = all(p)
    if (lh_dev .or. .not. any(p)) return
    call lh_fail(rname//': some matrix operands are on the device and some are not', __LINE__)
  end function lh_dev

  subroutine lh_init()
    integer :: dev, rc, i, evl, evs
    character(kind=c_char) :: nm(128)
    character(len=128) :: name
    character(len=8) :: ev
    if (inited) return
    if (omp_get_num_devices() < 1) call lh_fail('no OpenMP offload device', __LINE__)
    dev = omp_get_default_device()
    nm = c_null_char
    rc = alfhp_init(int(dev, c_int), nm, 128_c_int)
    if (rc /= 0) then
       write(error_unit,'(A,I0,A,I0)') ' linalg (HIP): alfhp_init failed, code ', rc, ', device ', dev
       call lh_fail('initialisation failed', __LINE__)
    endif
    name = ''
    do i = 1, 127
       if (nm(i) == c_null_char) exit
       name(i:i) = nm(i)
    enddo
    call get_environment_variable('ALF_HIP_QR', ev, evl, evs)
    if (evs == 0 .and. evl >= 1) use_hq = (ev(1:1) /= '0')
    if (use_hq) then
       write(*,'(3A)') ' linalg (HIP): device ', trim(name), '; QR: own HIP QR for m <= 2304, else hipSOLVER'
    else
       write(*,'(3A)') ' linalg (HIP): device ', trim(name), '; QR: hipSOLVER (ALF_HIP_QR=0)'
    endif
    inited = .true.
  end subroutine lh_init

  subroutine ensure_z(a, n)
    complex(kind=kind(0.d0)), allocatable, target, intent(inout) :: a(:)
    integer, intent(in) :: n
    if (allocated(a)) then
       if (size(a) >= n) return
       !$omp target exit data map(delete: a)
       deallocate(a)
    endif
    allocate(a(max(n,1)))
    !$omp target enter data map(alloc: a)
  end subroutine ensure_z

  !> API failure codes of the shim (>= 998) stop the run; a solver info is returned
  integer function api(rc, what)
    integer(c_int), intent(in) :: rc
    character(len=*), intent(in) :: what
    character(len=80) :: msg
    if (rc >= 998) then
       write(msg,'(2A,I0)') what, ' failed (HIP API), code ', rc
       call lh_fail(trim(msg), __LINE__)
    endif
    api = int(rc)
  end function api

  integer(c_int) function opn(op)
    character(len=1), intent(in) :: op
    select case (op)
    case ('N', 'n'); opn = 0
    case ('T', 't'); opn = 1
    case default;    opn = 2
    end select
  end function opn

  integer function ws_size(rc, ws, what)
    integer(c_int), intent(in) :: rc
    integer(c_int64_t), intent(in) :: ws
    character(len=*), intent(in) :: what
    if (rc /= 0) call lh_fail(what//': work-size query failed', __LINE__)
    ws_size = int(max(ws, 1_c_int64_t))
  end function ws_size

  !> C = alpha op(A) op(B) + beta C (hipBLAS)
  subroutine lh_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    character(len=1), intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(kind=kind(0.d0)), intent(in) :: alpha, beta
    complex(kind=kind(0.d0)), target :: a(lda,*), b(ldb,*), c(ldc,*)
    integer(c_int) :: rc
    call lh_init()
    !$omp target data use_device_addr(a, b, c)
    rc = alfhp_zgemm(opn(transa), opn(transb), int(m,c_int), int(n,c_int), int(k,c_int), dble(alpha), aimag(alpha), &
         & c_loc(a(1,1)), int(lda,c_int), c_loc(b(1,1)), int(ldb,c_int), dble(beta), aimag(beta), c_loc(c(1,1)), int(ldc,c_int))
    !$omp end target data
    if (rc /= 0) call lh_fail('hipblasZgemm failed', __LINE__)
  end subroutine lh_zgemm

  !> C(complex m x n) = A(complex m x k) * B(real k x n): hipBLAS DGEMM on the interleaved real view
  subroutine lh_zgemm_cr(m, n, k, a, lda, b, ldb, c, ldc)
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(kind=kind(0.d0)), target :: a(lda,*), c(ldc,*)
    real(kind=kind(0.d0)), target :: b(ldb,*)
    integer(c_int) :: rc
    if (m <= 0 .or. n <= 0) return
    call lh_init()
    !$omp target data use_device_addr(a, b, c)
    rc = alfhp_dgemm(int(2*m,c_int), int(n,c_int), int(k,c_int), c_loc(a(1,1)), int(2*lda,c_int), c_loc(b(1,1)), &
         & int(ldb,c_int), c_loc(c(1,1)), int(2*ldc,c_int))
    !$omp end target data
    if (rc /= 0) call lh_fail('hipblasDgemm failed', __LINE__)
  end subroutine lh_zgemm_cr

  !> Triangular solve with several right-hand sides (hipBLAS)
  subroutine lh_ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    character(len=1), intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(kind=kind(0.d0)), intent(in) :: alpha
    complex(kind=kind(0.d0)), target :: a(lda,*), b(ldb,*)
    integer(c_int) :: rc, cs, cu, cdg
    call lh_init()
    cs = merge(1, 0, side == 'R' .or. side == 'r')
    cu = merge(1, 0, uplo == 'U' .or. uplo == 'u')
    cdg = merge(1, 0, diag == 'U' .or. diag == 'u')
    !$omp target data use_device_addr(a, b)
    rc = alfhp_ztrsm(cs, cu, opn(transa), cdg, int(m,c_int), int(n,c_int), dble(alpha), aimag(alpha), c_loc(a(1,1)), &
         & int(lda,c_int), c_loc(b(1,1)), int(ldb,c_int))
    !$omp end target data
    if (rc /= 0) call lh_fail('hipblasZtrsm failed', __LINE__)
  end subroutine lh_ztrsm

  !> QR factorisation (ZGEQRF semantics): own HIP QR, or hipSOLVER
  subroutine lh_zgeqrf(m, n, a, lda, tau, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target :: a(lda,*), tau(*)
    integer, intent(out) :: info
    integer(c_int) :: rc
    integer(c_int64_t) :: ws
    integer :: lw
    call lh_init()
    rc = 3
    if (use_hq) then
       !$omp target data use_device_addr(a, tau)
       rc = alfhq_geqrf(int(m,c_int), int(n,c_int), c_loc(a(1,1)), int(lda,c_int), c_loc(tau(1)))
       !$omp end target data
    endif
    if (rc == 3) then           ! size not supported by the HIP QR (m > 2304) or switched off: hipSOLVER
       ws = 0
       lw = ws_size(alfhp_zgeqrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws), ws, 'zgeqrf')
       call ensure_z(dws, lw)
       !$omp target data use_device_addr(a, tau, dws)
       rc = alfhp_zgeqrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), c_loc(tau(1)), &
            & c_loc(dws(1)), int(lw,c_int64_t))
       !$omp end target data
    endif
    info = api(rc, 'zgeqrf')
  end subroutine lh_zgeqrf

  !> C = op(Q) C or C op(Q) (ZUNMQR semantics): own HIP Q-apply, or hipSOLVER
  subroutine lh_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, info)
    character(len=1), intent(in) :: side, trans
    integer, intent(in) :: m, n, k, lda, ldc
    complex(kind=kind(0.d0)), target :: a(lda,*), tau(*), c(ldc,*)
    integer, intent(out) :: info
    integer(c_int) :: rc, cs, ct
    integer(c_int64_t) :: ws
    integer :: lw
    call lh_init()
    cs = merge(1, 0, side == 'R' .or. side == 'r')
    ct = opn(trans)
    rc = 3
    if (use_hq .and. ct /= 1) then      ! trans = 'T' is not unitary: hipSOLVER only
       !$omp target data use_device_addr(a, tau, c)
       rc = alfhq_unmqr(cs, ct, int(m,c_int), int(n,c_int), int(k,c_int), c_loc(a(1,1)), int(lda,c_int), c_loc(tau(1)), &
            & c_loc(c(1,1)), int(ldc,c_int))
       !$omp end target data
    endif
    if (rc == 3) then
       ws = 0
       lw = ws_size(alfhp_zunmqr_ws(cs, ct, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), int(lda,c_int64_t), &
            & int(ldc,c_int64_t), ws), ws, 'zunmqr')
       call ensure_z(dws, lw)
       !$omp target data use_device_addr(a, tau, c, dws)
       rc = alfhp_zunmqr(cs, ct, int(m,c_int64_t), int(n,c_int64_t), int(k,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), &
            & c_loc(tau(1)), c_loc(c(1,1)), int(ldc,c_int64_t), c_loc(dws(1)), int(lw,c_int64_t))
       !$omp end target data
    endif
    info = api(rc, 'zunmqr')
  end subroutine lh_zunmqr

  !> Column-norm pre-pivot + ZGEQRF (la_zqrpp contract): norms on the device, n norms D->H, the stable descending insertion
  !> sort of the host la_zqrpp, ipvt on the host and H->D, columns permuted on the device, device QR.  info = -101 if A has a
  !> NaN or Inf.
  subroutine lh_zqrpp(m, n, a, lda, ipvt, tau, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target :: a(lda,*), tau(*)
    integer, target, intent(inout) :: ipvt(*)
    integer, intent(out) :: info
    integer :: i, j, kk
    integer(kind=8) :: bits
    info = 0
    if (m <= 0 .or. n <= 0) return
    if (m < n) then
       info = -1
       return
    endif
    call lh_init()
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
    call hk_colnorms(m, n, a, lda, dnrm)
    ! D->H: n column norms for the host sort
    !$omp target update from(dnrm(1:n))
    do j = 1, n
       bits = transfer(dnrm(j), 0_8)
       if (iand(bits, int(z'7FF0000000000000', 8)) == int(z'7FF0000000000000', 8)) then
          info = -101
          return
       endif
    enddo
    do j = 1, n
       ipvt(j) = j
    enddo
    do i = 2, n
       kk = ipvt(i)
       j = i - 1
       do while (j >= 1)
          if (dnrm(ipvt(j)) >= dnrm(kk)) exit
          ipvt(j+1) = ipvt(j)
          j = j - 1
       enddo
       ipvt(j+1) = kk
    enddo
    ! H->D: the n-entry permutation
    !$omp target update to(ipvt(1:n))
    call hk_permute_cols(m, n, a, lda, ipvt)
    call lh_zgeqrf(m, n, a, lda, tau, info)
  end subroutine lh_zqrpp

  !> LU factorisation with partial pivoting (hipSOLVER, device ipiv); info > 0 (exactly singular) is returned as on the host
  subroutine lh_zgetrf(m, n, a, lda, ipiv, info)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target :: a(lda,*)
    integer, target :: ipiv(*)
    integer, intent(out) :: info
    integer(c_int64_t) :: ws
    integer :: lw
    integer(c_int) :: rc
    call lh_init()
    ws = 0
    lw = ws_size(alfhp_zgetrf_ws(int(m,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), ws), ws, 'zgetrf')
    call ensure_z(dws, lw)
    !$omp target data use_device_addr(a, ipiv, dws)
    rc = alfhp_zgetrf(int(m,c_int64_t), int(n,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), c_loc(ipiv(1)), &
         & c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
    info = api(rc, 'zgetrf')
  end subroutine lh_zgetrf

  !> Inverse from the LU factors: ZGETRS on the identity in a device scratch, copied back to A
  subroutine lh_zgetri(n, a, lda, ipiv, info)
    integer, intent(in) :: n, lda
    complex(kind=kind(0.d0)), target :: a(lda,*)
    integer, target :: ipiv(*)
    integer, intent(out) :: info
    integer(c_int64_t) :: ws
    integer :: lw
    integer(c_int) :: rc
    info = 0
    if (n <= 0) return
    ! shortcut: A must be stored contiguously (lda = n; the only caller, gpu_tau_mod, does); a 2-D copy kernel lifts this
    if (lda /= n) call lh_fail('la_zgetri on the device needs lda = n', __LINE__)
    call lh_init()
    call ensure_z(dscr, n*n)
    call hk_ident(n, n, dscr)
    ws = 0
    lw = ws_size(alfhp_zgetrs_ws(0_c_int, int(n,c_int64_t), int(n,c_int64_t), int(lda,c_int64_t), int(n,c_int64_t), ws), &
         & ws, 'zgetrs')
    call ensure_z(dws, lw)
    !$omp target data use_device_addr(a, ipiv, dscr, dws)
    rc = alfhp_zgetrs(0_c_int, int(n,c_int64_t), int(n,c_int64_t), c_loc(a(1,1)), int(lda,c_int64_t), c_loc(ipiv(1)), &
         & c_loc(dscr(1)), int(n,c_int64_t), c_loc(dws(1)), int(lw,c_int64_t))
    !$omp end target data
    info = api(rc, 'zgetrs')
    if (info /= 0) return
    call hk_copyz(n*n, dscr, a)
  end subroutine lh_zgetri

end module linalg_hip_mod
#endif

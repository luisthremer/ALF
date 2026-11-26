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
!> Ported from fork/gpu-port alf_la_om_mod (Libraries/Modules/alf_la_mod.F90 :125-978) without its info slots and test
!> hooks.
!--------------------------------------------------------------------
module linalg_onemkl_mod
  use iso_c_binding, only: c_ptr, c_loc, c_int, c_f_pointer
  use omp_lib, only: omp_get_default_device, omp_get_num_devices, omp_target_is_present
  use onemkl_blas_omp_offload_lp64, only: zgemm, dgemm, ztrsm
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: om_zgemm, om_zgemm_cr, om_ztrsm

  integer, save :: ndev = -1            ! number of offload devices, queried once

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

end module linalg_onemkl_mod

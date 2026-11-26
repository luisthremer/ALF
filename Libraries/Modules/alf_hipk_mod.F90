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
!> GPU=amd (ALF_LA_HIP): HIP implementations of small device kernels (alf_hip_kernels.hip.cpp) behind plain subroutines.
!
!> @details
!> A libgomp-GCN OpenMP target region with a parallel loop costs ~150 us to launch, a HIP kernel ~15 us.  The arguments are
!> the host addresses of device-resident (mapped) arrays; the device addresses come from use_device_addr.  Every routine
!> reproduces the OpenMP kernel it replaces (same arguments and results); a failed launch ends the run.  Compiled only for
!> GPU=amd.
!--------------------------------------------------------------------
#if defined(ALF_LA_HIP)
module alf_hipk_mod
  use iso_c_binding, only: c_int, c_long, c_ptr, c_loc
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  implicit none
  private
  public :: hk_copyz, hk_ident, hk_colnorms, hk_permute_cols

  interface
    integer(c_int) function alfhk_copyz(n, a, b) bind(C, name='alfhk_copyz')
      import :: c_ptr, c_int, c_long
      integer(c_long), value :: n
      type(c_ptr), value :: a, b
    end function
    integer(c_int) function alfhk_ident(n, ld, U) bind(C, name='alfhk_ident')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U
    end function
    integer(c_int) function alfhk_colnorms(m, n, A, lda, nrm) bind(C, name='alfhk_colnorms')
      import :: c_ptr, c_int
      integer(c_int), value :: m, n, lda
      type(c_ptr), value :: A, nrm
    end function
    integer(c_int) function alfhk_permute_cols(m, n, A, lda, ipvt) bind(C, name='alfhk_permute_cols')
      import :: c_ptr, c_int
      integer(c_int), value :: m, n, lda
      type(c_ptr), value :: A, ipvt
    end function
  end interface

contains

  subroutine chk(rc, what)
    integer(c_int), intent(in) :: rc
    character(len=*), intent(in) :: what
    if (rc /= 0) then
       write(*,*) 'alf_hipk: ', what, ' failed, rc = ', rc
       call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
    endif
  end subroutine chk

  !> b(1:n) = a(1:n)
  subroutine hk_copyz(n, a, b)
    integer, intent(in) :: n
    complex(kind=kind(0.d0)), target :: a(*), b(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(a, b)
    rc = alfhk_copyz(int(n,c_long), c_loc(a(1)), c_loc(b(1)))
    !$omp end target data
    call chk(rc, 'copyz')
  end subroutine hk_copyz

  !> U(1:n,1:n) = identity (leading dimension ld)
  subroutine hk_ident(n, ld, U)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: U(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U)
    rc = alfhk_ident(int(n,c_int), int(ld,c_int), c_loc(U(1)))
    !$omp end target data
    call chk(rc, 'ident')
  end subroutine hk_ident

  !> nrm(j) = 2-norm of column j of A(m,n) (scaled); NaN if the column holds a NaN or Inf
  subroutine hk_colnorms(m, n, A, lda, nrm)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target :: A(*)
    real(kind=kind(0.d0)), target :: nrm(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, nrm)
    rc = alfhk_colnorms(int(m,c_int), int(n,c_int), c_loc(A(1)), int(lda,c_int), c_loc(nrm(1)))
    !$omp end target data
    call chk(rc, 'colnorms')
  end subroutine hk_colnorms

  !> A(:,j) <- A(:,ipvt(j)), j = 1..n (device ipvt)
  subroutine hk_permute_cols(m, n, A, lda, ipvt)
    integer, intent(in) :: m, n, lda
    complex(kind=kind(0.d0)), target :: A(*)
    integer, target :: ipvt(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, ipvt)
    rc = alfhk_permute_cols(int(m,c_int), int(n,c_int), c_loc(A(1)), int(lda,c_int), c_loc(ipvt(1)))
    !$omp end target data
    call chk(rc, 'permute_cols')
  end subroutine hk_permute_cols

end module alf_hipk_mod
#endif

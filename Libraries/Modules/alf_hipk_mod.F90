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
  use iso_c_binding, only: c_int, c_long, c_double, c_ptr, c_loc
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  implicit none
  private
  public :: hk_copyz, hk_ident, hk_colnorms, hk_permute_cols
  public :: hk_enabled, hk_set_enabled
  public :: hk_setd, hk_setl, hk_scale_cols, hk_diagabs, hk_divrows, hk_phase, hk_scale_row1, hk_scale_col1, hk_triu, &
       &    hk_perm, hk_perm_vec, hk_ct, hk_rowdiv, hk_colinv, hk_tpup, hk_det_phase, hk_cgr_phase, hk_testscale, hk_prec
  public :: hk_scatter_rows, hk_vslice_r, hk_blk_gather, hk_blk_scale, hk_sub_flush_gather
  public :: hk_corr, hk_scalars

  ! HIP replacements on (default) or off (ALF_HIP_KERNELS=0: the OpenMP kernels run, for A/B checks); -1 = not read yet
  integer, save :: hk_state = -1

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
    integer(c_int) function alfhk_setd(n, D) bind(C, name='alfhk_setd')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: D
    end function
    integer(c_int) function alfhk_setl(n, L) bind(C, name='alfhk_setl')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: L
    end function
    integer(c_int) function alfhk_scale_cols(n, ld, U, D) bind(C, name='alfhk_scale_cols')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, D
    end function
    integer(c_int) function alfhk_diagabs(n, ld, U, D) bind(C, name='alfhk_diagabs')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, D
    end function
    integer(c_int) function alfhk_divrows(n, ld, U, D) bind(C, name='alfhk_divrows')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, D
    end function
    integer(c_int) function alfhk_phase(n, ld, U, sgn, cj, tau, cs) bind(C, name='alfhk_phase')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld, sgn, cj
      type(c_ptr), value :: U, tau, cs
    end function
    integer(c_int) function alfhk_scale_row1(n, ld, U, cs) bind(C, name='alfhk_scale_row1')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, cs
    end function
    integer(c_int) function alfhk_scale_col1(n, ld, U, cs) bind(C, name='alfhk_scale_col1')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, cs
    end function
    integer(c_int) function alfhk_triu(n, ld, U, ldr, R) bind(C, name='alfhk_triu')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld, ldr
      type(c_ptr), value :: U, R
    end function
    integer(c_int) function alfhk_perm(mode, n, X, ldx, perm, W, ldw) bind(C, name='alfhk_perm')
      import :: c_ptr, c_int
      integer(c_int), value :: mode, n, ldx, ldw
      type(c_ptr), value :: X, perm, W
    end function
    integer(c_int) function alfhk_perm_vec(n, L, perm, W) bind(C, name='alfhk_perm_vec')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: L, perm, W
    end function
    integer(c_int) function alfhk_ct(n, A, lda, B, ldb) bind(C, name='alfhk_ct')
      import :: c_ptr, c_int
      integer(c_int), value :: n, lda, ldb
      type(c_ptr), value :: A, B
    end function
    integer(c_int) function alfhk_rowdiv(n, lda, A, D) bind(C, name='alfhk_rowdiv')
      import :: c_ptr, c_int
      integer(c_int), value :: n, lda
      type(c_ptr), value :: A, D
    end function
    integer(c_int) function alfhk_colinv(n, lda, A, D) bind(C, name='alfhk_colinv')
      import :: c_ptr, c_int
      integer(c_int), value :: n, lda
      type(c_ptr), value :: A, D
    end function
    integer(c_int) function alfhk_tpup(n, T, R, Dr, Dl) bind(C, name='alfhk_tpup')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: T, R, Dr, Dl
    end function
    integer(c_int) function alfhk_det_phase(n, ld, A, ipiv, cs) bind(C, name='alfhk_det_phase')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: A, ipiv, cs
    end function
    integer(c_int) function alfhk_cgr_phase(n, ld, T, tau, cs, nvar, sgn, bad, p0r, p0i, dev) bind(C, name='alfhk_cgr_phase')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: n, ld, nvar, sgn, bad, dev
      real(c_double), value :: p0r, p0i
      type(c_ptr), value :: T, tau, cs
    end function
    integer(c_int) function alfhk_testscale(n, D, r, bad) bind(C, name='alfhk_testscale')
      import :: c_ptr, c_int
      integer(c_int), value :: n, bad
      type(c_ptr), value :: D, r
    end function
    integer(c_int) function alfhk_prec(n, A, B, r) bind(C, name='alfhk_prec')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: A, B, r
    end function
    integer(c_int) function alfhk_scatter_rows(nd, ncol, P, M, ld, W, ldw) bind(C, name='alfhk_scatter_rows')
      import :: c_ptr, c_int
      integer(c_int), value :: nd, ncol, ld, ldw
      type(c_ptr), value :: P, M, W
    end function
    integer(c_int) function alfhk_vslice_r(nd, ncol, ld, nf, t, cj, desc, ldptr, nop, ldn, vr_ptr, vr_ent, V_meta, V_type, &
         & nsig, V_exp, M) bind(C, name='alfhk_vslice_r')
      import :: c_ptr, c_int
      integer(c_int), value :: nd, ncol, ld, nf, t, cj, desc, ldptr, nop, ldn
      type(c_ptr), value :: vr_ptr, vr_ent, V_meta, V_type, nsig, V_exp, M
    end function
    integer(c_int) function alfhk_blk_gather(nfl, nfa, nfb, ld, b, nd, G, bsk, gbk) bind(C, name='alfhk_blk_gather')
      import :: c_ptr, c_int
      integer(c_int), value :: nfl, nfa, nfb, ld, b, nd
      type(c_ptr), value :: G, bsk, gbk
    end function
    integer(c_int) function alfhk_blk_scale(nfl, nfa, nfb, ld, b, nd, ldc, G, bsk, cbk) bind(C, name='alfhk_blk_scale')
      import :: c_ptr, c_int
      integer(c_int), value :: nfl, nfa, nfb, ld, b, nd, ldc
      type(c_ptr), value :: G, bsk, cbk
    end function
    integer(c_int) function alfhk_sub_flush_gather(ld, nrow, k, nf, ldslot, G, S, C, Qt, slot) bind(C, name='alfhk_sub_flush_gather')
      import :: c_ptr, c_int
      integer(c_int), value :: ld, nrow, k, nf, ldslot
      type(c_ptr), value :: G, S, C, Qt, slot
    end function
    integer(c_int) function alfhk_corr(kind, k, n, norb, nfl, ns, sq3, wr, wi, ld, g, jc, inv, acc, acc0) bind(C, name='alfhk_corr')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: kind, k, n, norb, nfl, ld
      real(c_double), value :: ns, sq3, wr, wi
      type(c_ptr), value :: g, jc, inv, acc, acc0
    end function
    integer(c_int) function alfhk_scalars(nk, hi, hj, hf, hc, n, norb, nfl, dec, ld, g, inv, uo, hres) bind(C, name='alfhk_scalars')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: nk, n, norb, nfl, dec, ld
      type(c_ptr), value :: hi, hj, hf, hc, g, inv, uo
      real(c_double) :: hres(6)
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

  !> .true. if the HIP replacements are on (default; ALF_HIP_KERNELS=0 switches them off, read once)
  logical function hk_enabled()
    character(len=8) :: v
    integer :: l, st
    if (hk_state < 0) then
       hk_state = 1
       call get_environment_variable('ALF_HIP_KERNELS', v, l, st)
       if (st == 0 .and. l >= 1) then
          if (v(1:1) == '0') hk_state = 0
       endif
       if (hk_state == 0) write(*,'(A)') ' alf_hipk: HIP kernels OFF (ALF_HIP_KERNELS=0), OpenMP kernels run'
    endif
    hk_enabled = (hk_state == 1)
  end function hk_enabled

  !> switch the HIP replacements on or off (tests: the OpenMP and the HIP kernel on the same data)
  subroutine hk_set_enabled(on)
    logical, intent(in) :: on
    hk_state = merge(1, 0, on)
  end subroutine hk_set_enabled

  subroutine hk_setd(n, D)
    integer, intent(in) :: n
    complex(kind=kind(0.d0)), target :: D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(D)
    rc = alfhk_setd(int(n,c_int), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'setd')
  end subroutine hk_setd

  subroutine hk_setl(n, L)
    integer, intent(in) :: n
    real(kind=kind(0.d0)), target :: L(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(L)
    rc = alfhk_setl(int(n,c_int), c_loc(L(1)))
    !$omp end target data
    call chk(rc, 'setl')
  end subroutine hk_setl

  subroutine hk_scale_cols(n, ld, U, D)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: U(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, D)
    rc = alfhk_scale_cols(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'scale_cols')
  end subroutine hk_scale_cols

  subroutine hk_diagabs(n, ld, U, D)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: U(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, D)
    rc = alfhk_diagabs(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'diagabs')
  end subroutine hk_diagabs

  subroutine hk_divrows(n, ld, U, D)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: U(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, D)
    rc = alfhk_divrows(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'divrows')
  end subroutine hk_divrows

  subroutine hk_phase(n, ld, U, sgn, cj, tau, cs)
    integer, intent(in) :: n, ld, sgn
    logical, intent(in) :: cj
    complex(kind=kind(0.d0)), target :: U(*), tau(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, tau, cs)
    rc = alfhk_phase(int(n,c_int), int(ld,c_int), c_loc(U(1)), int(sgn,c_int), int(merge(1,0,cj),c_int), c_loc(tau(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'phase')
  end subroutine hk_phase

  subroutine hk_scale_row1(n, ld, U, cs)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: U(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, cs)
    rc = alfhk_scale_row1(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'scale_row1')
  end subroutine hk_scale_row1

  subroutine hk_scale_col1(n, ld, U, cs)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: U(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, cs)
    rc = alfhk_scale_col1(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'scale_col1')
  end subroutine hk_scale_col1

  subroutine hk_triu(n, ld, U, ldr, R)
    integer, intent(in) :: n, ld, ldr
    complex(kind=kind(0.d0)), target :: U(*), R(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, R)
    rc = alfhk_triu(int(n,c_int), int(ld,c_int), c_loc(U(1)), int(ldr,c_int), c_loc(R(1)))
    !$omp end target data
    call chk(rc, 'triu')
  end subroutine hk_triu

  subroutine hk_perm(mode, n, X, ldx, perm, W, ldw)
    integer, intent(in) :: mode, n, ldx, ldw
    complex(kind=kind(0.d0)), target :: X(*), W(*)
    integer, target :: perm(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(X, perm, W)
    rc = alfhk_perm(int(mode,c_int), int(n,c_int), c_loc(X(1)), int(ldx,c_int), c_loc(perm(1)), c_loc(W(1)), int(ldw,c_int))
    !$omp end target data
    call chk(rc, 'perm')
  end subroutine hk_perm

  subroutine hk_perm_vec(n, L, perm, W)
    integer, intent(in) :: n
    real(kind=kind(0.d0)), target :: L(*), W(*)
    integer, target :: perm(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(L, perm, W)
    rc = alfhk_perm_vec(int(n,c_int), c_loc(L(1)), c_loc(perm(1)), c_loc(W(1)))
    !$omp end target data
    call chk(rc, 'perm_vec')
  end subroutine hk_perm_vec

  subroutine hk_ct(n, A, lda, B, ldb)
    integer, intent(in) :: n, lda, ldb
    complex(kind=kind(0.d0)), target :: A(*), B(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, B)
    rc = alfhk_ct(int(n,c_int), c_loc(A(1)), int(lda,c_int), c_loc(B(1)), int(ldb,c_int))
    !$omp end target data
    call chk(rc, 'ct')
  end subroutine hk_ct

  subroutine hk_rowdiv(n, lda, A, D)
    integer, intent(in) :: n, lda
    complex(kind=kind(0.d0)), target :: A(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, D)
    rc = alfhk_rowdiv(int(n,c_int), int(lda,c_int), c_loc(A(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'rowdiv')
  end subroutine hk_rowdiv

  subroutine hk_colinv(n, lda, A, D)
    integer, intent(in) :: n, lda
    complex(kind=kind(0.d0)), target :: A(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, D)
    rc = alfhk_colinv(int(n,c_int), int(lda,c_int), c_loc(A(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'colinv')
  end subroutine hk_colinv

  subroutine hk_tpup(n, T, R, Dr, Dl)
    integer, intent(in) :: n
    complex(kind=kind(0.d0)), target :: T(*), R(*), Dr(*), Dl(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(T, R, Dr, Dl)
    rc = alfhk_tpup(int(n,c_int), c_loc(T(1)), c_loc(R(1)), c_loc(Dr(1)), c_loc(Dl(1)))
    !$omp end target data
    call chk(rc, 'tpup')
  end subroutine hk_tpup

  subroutine hk_det_phase(n, ld, A, ipiv, cs)
    integer, intent(in) :: n, ld
    complex(kind=kind(0.d0)), target :: A(*), cs(*)
    integer, target :: ipiv(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, ipiv, cs)
    rc = alfhk_det_phase(int(n,c_int), int(ld,c_int), c_loc(A(1)), c_loc(ipiv(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'det_phase')
  end subroutine hk_det_phase

  subroutine hk_cgr_phase(n, ld, T, tau, cs, nvar, sgn, bad, ph0, dev)
    integer, intent(in) :: n, ld, nvar, sgn, bad, dev
    complex(kind=kind(0.d0)), intent(in) :: ph0
    complex(kind=kind(0.d0)), target :: T(*), tau(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(T, tau, cs)
    rc = alfhk_cgr_phase(int(n,c_int), int(ld,c_int), c_loc(T(1)), c_loc(tau(1)), c_loc(cs(1)), int(nvar,c_int), int(sgn,c_int), &
         & int(bad,c_int), real(ph0,c_double), aimag(ph0), int(dev,c_int))
    !$omp end target data
    call chk(rc, 'cgr_phase')
  end subroutine hk_cgr_phase

  subroutine hk_testscale(n, D, r, bad)
    integer, intent(in) :: n, bad
    complex(kind=kind(0.d0)), target :: D(*)
    real(kind=kind(0.d0)), target :: r(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(D, r)
    rc = alfhk_testscale(int(n,c_int), c_loc(D(1)), c_loc(r(1)), int(bad,c_int))
    !$omp end target data
    call chk(rc, 'testscale')
  end subroutine hk_testscale

  subroutine hk_prec(n, A, B, r)
    integer, intent(in) :: n
    complex(kind=kind(0.d0)), target :: A(*), B(*)
    real(kind=kind(0.d0)), target :: r(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, B, r)
    rc = alfhk_prec(int(n,c_int), c_loc(A(1)), c_loc(B(1)), c_loc(r(1)))
    !$omp end target data
    call chk(rc, 'prec')
  end subroutine hk_prec

  !> M(P(l), j) = W(l, j)  (k_scatter_rows of gpu_wrap_mod)
  subroutine hk_scatter_rows(nd, ncol, P, M, ld, W, ldw)
    integer, intent(in) :: nd, ncol, ld, ldw
    integer, target :: P(*)
    complex(kind=kind(0.d0)), target :: M(*), W(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(P, M, W)
    rc = alfhk_scatter_rows(int(nd,c_int), int(ncol,c_int), c_loc(P(1)), c_loc(M(1)), int(ld,c_int), c_loc(W(1)), int(ldw,c_int))
    !$omp end target data
    call chk(rc, 'scatter_rows')
  end subroutine hk_scatter_rows

  !> the Op_V slice of k_vslice_R (gpu_wrap_mod); ldptr, nop, ldn: leading extents of vr_ptr, V_meta (2nd), nsigma_d
  subroutine hk_vslice_r(nd, ncol, ld, nf, t, cj, desc, ldptr, nop, ldn, vr_ptr, vr_ent, V_meta, V_type, nsig, V_exp, M)
    integer, intent(in) :: nd, ncol, ld, nf, t, ldptr, nop, ldn
    logical, intent(in) :: cj, desc
    integer, target :: vr_ptr(*), vr_ent(*), V_meta(*), V_type(*)
    complex(kind=kind(0.d0)), target :: nsig(*), V_exp(*), M(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(vr_ptr, vr_ent, V_meta, V_type, nsig, V_exp, M)
    rc = alfhk_vslice_r(int(nd,c_int), int(ncol,c_int), int(ld,c_int), int(nf,c_int), int(t,c_int), int(merge(1,0,cj),c_int), &
         & int(merge(1,0,desc),c_int), int(ldptr,c_int), int(nop,c_int), int(ldn,c_int), c_loc(vr_ptr(1)), c_loc(vr_ent(1)), &
         & c_loc(V_meta(1)), c_loc(V_type(1)), c_loc(nsig(1)), c_loc(V_exp(1)), c_loc(M(1)))
    !$omp end target data
    call chk(rc, 'vslice_r')
  end subroutine hk_vslice_r

  !> k_blk_gather of gpu_lowrank_mod
  subroutine hk_blk_gather(nfl, nfa, nfb, ld, b, nd, G, bsk, gbk)
    integer, intent(in) :: nfl, nfa, nfb, ld, b, nd
    complex(kind=kind(0.d0)), target :: G(*), gbk(*)
    integer, target :: bsk(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(G, bsk, gbk)
    rc = alfhk_blk_gather(int(nfl,c_int), int(nfa,c_int), int(nfb,c_int), int(ld,c_int), int(b,c_int), int(nd,c_int), &
         & c_loc(G(1)), c_loc(bsk(1)), c_loc(gbk(1)))
    !$omp end target data
    call chk(rc, 'blk_gather')
  end subroutine hk_blk_gather

  !> k_blk_scale of gpu_lowrank_mod
  subroutine hk_blk_scale(nfl, nfa, nfb, ld, b, nd, ldc, G, bsk, cbk)
    integer, intent(in) :: nfl, nfa, nfb, ld, b, nd, ldc
    complex(kind=kind(0.d0)), target :: G(*), cbk(*)
    integer, target :: bsk(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(G, bsk, cbk)
    rc = alfhk_blk_scale(int(nfl,c_int), int(nfa,c_int), int(nfb,c_int), int(ld,c_int), int(b,c_int), int(nd,c_int), &
         & int(ldc,c_int), c_loc(G(1)), c_loc(bsk(1)), c_loc(cbk(1)))
    !$omp end target data
    call chk(rc, 'blk_scale')
  end subroutine hk_blk_scale

  !> k_sub_flush_gather of gpu_lowrank_mod
  subroutine hk_sub_flush_gather(G, ld, nrow, S, k, nf, C, Qt, slot, ldslot)
    integer, intent(in) :: ld, nrow, k, nf, ldslot
    complex(kind=kind(0.d0)), target :: G(*), C(*), Qt(*)
    integer, target :: S(*), slot(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(G, S, C, Qt, slot)
    rc = alfhk_sub_flush_gather(int(ld,c_int), int(nrow,c_int), int(k,c_int), int(nf,c_int), int(ldslot,c_int), c_loc(G(1)), &
         & c_loc(S(1)), c_loc(C(1)), c_loc(Qt(1)), c_loc(slot(1)))
    !$omp end target data
    call chk(rc, 'sub_flush_gather')
  end subroutine hk_sub_flush_gather

  !> equal-time correlation of kind kind_k (k_corr of gpu_obs_mod) into acc(:,1,:,:,k..) and acc0(:,k..)
  subroutine hk_corr(kind_k, k, n, norb, nfl, ns, sq3, wr, wi, ld, g, jc, inv, acc, acc0)
    integer, intent(in) :: kind_k, k, n, norb, nfl, ld
    real(kind=kind(0.d0)), intent(in) :: ns, sq3, wr, wi
    complex(kind=kind(0.d0)), target :: g(*), acc(*), acc0(*)
    integer, target :: jc(*), inv(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(g, jc, inv, acc, acc0)
    rc = alfhk_corr(int(kind_k,c_int), int(k,c_int), int(n,c_int), int(norb,c_int), int(nfl,c_int), ns, sq3, wr, wi, &
         & int(ld,c_int), c_loc(g(1)), c_loc(jc(1)), c_loc(inv(1)), c_loc(acc(1)), c_loc(acc0(1)))
    !$omp end target data
    call chk(rc, 'corr')
  end subroutine hk_corr

  !> hres(1:6) = re/im of Zkin (before the N_SUN factor), Zpot, Zrho (before N_SUN) of measure_on (gpu_obs_mod); 48 B D->H
  subroutine hk_scalars(nk, hi, hj, hf, hc, n, norb, nfl, dec, ld, g, inv, uo, hres)
    integer, intent(in) :: nk, n, norb, nfl, dec, ld
    integer, target :: hi(*), hj(*), hf(*), inv(*)
    complex(kind=kind(0.d0)), target :: hc(*), g(*)
    real(kind=kind(0.d0)), target :: uo(*)
    real(kind=kind(0.d0)), intent(out) :: hres(6)
    integer(c_int) :: rc
    !$omp target data use_device_addr(hi, hj, hf, hc, g, inv, uo)
    rc = alfhk_scalars(int(nk,c_int), c_loc(hi(1)), c_loc(hj(1)), c_loc(hf(1)), c_loc(hc(1)), int(n,c_int), int(norb,c_int), &
         & int(nfl,c_int), int(dec,c_int), int(ld,c_int), c_loc(g(1)), c_loc(inv(1)), c_loc(uo(1)), hres)
    !$omp end target data
    call chk(rc, 'scalars')
  end subroutine hk_scalars

end module alf_hipk_mod
#endif

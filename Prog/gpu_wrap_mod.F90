!  Copyright (C) 2016 - 2025 The ALF project
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
!> @brief
!> GPU port, ticket 204 (ADR-002): device versions of the hopping exponential B_T (Hop_mod_mmthr, _m1, mmthl, mmthlc,
!> mmthl_m1) and of the Op_V operators of one time slice (Op_mmultL/R, Op_Wrapup/Op_Wrapdo), acting on DEVICE-RESIDENT
!> complex matrices, built on the flat arrays of gpu_state_mod and on alf_la.  Compiled only with -DALF_GPU.
!> Nothing in the sweep calls these yet (tickets 205-207 will).
!>
!> Conventions
!>   * every matrix argument is an explicit-shape/assumed-size device array M(ld,*) that the caller has mapped; the
!>     routines never allocate on the device (scratch = gpu_state scr_c(Ndim,Ndim,2)) and use synchronous target
!>     regions only (ADR-002 amended); alf_la calls are ordered with them on the same stream.
!>   * left routines (mmthr, mmthr_m1, mmthlc): M is Ndim x ncol, rows are sites.  Right routines (mmthl, mmthl_m1):
!>     M is nrow x Ndim (nrow <= Ndim for the dense/gather scratch), columns are sites.
!>   * the field of operator n at slice t is read on the device from nsigma_d(n,t) (the caller keeps it in sync).
!>
!> Hopping kernels (see gpu_state_mod header for groups)
!>   checkerboard (all elements nd <= 8), default (gpu_cb_fused = .false.): one kernel launch per group (launches per B_T
!>     = number of groups, 4 for the square-lattice checkerboard), threads over (element, tile of gpu_tl columns / rows);
!>     the element data (site list, matrix) are read once per tile.  gpu_cb_fused = .true.: ONE launch per B_T, teams over
!>     tiles of columns (rows), the groups in the order of the CPU routine inside the team with a team barrier between
!>     them (columns/rows of M are independent, so this is the sequential result).  Measured slower on the A30; kept for
!>     comparison.
!>   element arithmetic: nd <= 8 hand-written sums in the order of ZDSLSYMM/ZSLHEMM (l = 1..nd); real elements use real
!>     arithmetic (what the promoted product is); nd = 1 scales by the real value.
!>   dense (element nd > 8): gather (if P is not the identity) + GEMM into scratch + scatter.  Complex: ZGEMM
!>     (la_gemm).  Real, right: la_gemm_cr (DGEMM on the interleaved view).  Real, left: gpu_dense_left_mode = 1:
!>     complex copy of A + la_gemm (4x the flops); = 2: transpose M, la_gemm_cr with A^T, transpose back.
!>
!> Op_V kernels: tables V_exp of gpu_state_mod, index = field value from nsigma_d; per-operator routines (Op_mmultR/L,
!> Op_Wrapup/Wrapdo semantics, the L/R skip of tiny |g| honoured, the wrap ones ignoring it, as ALF does) and fused
!> routines for a whole slice that apply, for each row (column) of the matrix, the factors of all operators in ALF's
!> operator order (list vr_ptr/vr_ent), i.e. the same sequence of scalar products per element as ALF.
!>
!> Self-test: ALF_GPU_SELFTEST=2 (gpu_wrap_hook, called from main after gpu_state_setup); timing: ALF_GPU_TIMING=1;
!> ALF_GPU_EXIT=1 ends the run after the hook; ALF_GPU_FAULT=hop|hopc|vexp|nan|order injects a fault into the device
!> copy (negative controls, the self-test must fail).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_wrap_mod

  Use omp_lib
  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod
  Use Natural_Constants, only: Eps_machine
  Use Operator_mod
  Use Fields_mod
  Use Hamiltonian_main
  Use QMC_runtime_var
  Use Hop_mod
  Use alf_la_mod
  Use gpu_state_mod

  Implicit none

  Private
  Public :: gpu_mmthr, gpu_mmthr_m1, gpu_mmthlc, gpu_mmthl, gpu_mmthl_m1
  Public :: gpu_Vop_R, gpu_Vop_L, gpu_Vop_wrap, gpu_Vslice_R, gpu_Vslice_L, gpu_Vwrap_all
  Public :: gpu_wrapur_slices, gpu_wrapul_slices, gpu_wrapgr_up, gpu_wrapgr_do
  Public :: gpu_wrap_hook, gpu_wrap_selftest, gpu_wrap_timing
  Public :: gpu_cb_fused, gpu_dense_left_mode, gpu_ct, gpu_rbt, gpu_tl, gpu_tr, gpu_vb, gpu_n_launch, gpu_n_lib

  Logical, Save :: gpu_cb_fused = .false.     ! .true.: one launch per B_T (slower, see k_cb_left); default: one launch per group
  Integer, Save :: gpu_dense_left_mode = 0     ! 0 auto (v1 for nd <= 400, else v2), 1 complex copy + ZGEMM, 2 transposes + DGEMM (measured, ticket 204)
  Integer(int64), Save :: gpu_n_launch = 0, gpu_n_lib = 0    ! counters: target-region launches / alf_la calls issued here
  Logical, Save :: la_up = .false.
  Logical, Save :: fault_order = .false.
  Integer, Save :: gpu_tl = 1, gpu_tr = 2   ! columns (left) / rows (right) per thread of the per-group kernels (tuned on the A30)
  Integer, Save :: gpu_vb = 8                 ! columns/rows per thread of the fused Op_V slice kernels
  Integer, Save :: gpu_ct = 4, gpu_rbt = 8   ! columns (left) / rows (right) per team of the fused kernels (tuned, see gpu_wrap_timing)
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

#if defined(ALF_LA_ONEMKL)
  Logical, Parameter :: la_available = .false.
  ! measured on the Arc B580 (ifx 2025.3, spir64 JIT): the single-launch kernel (team barrier between groups inside one
  ! launch) gives wrong results (rel diff 0.2); the per-group kernels are exact.  Fused mode is therefore NVIDIA-only.
  Logical, Parameter :: fused_ok = .false.
#else
  Logical, Parameter :: la_available = .true.
  Logical, Parameter :: fused_ok = .true.
#endif

Contains

!====================================================================
!  Hopping
!====================================================================
  Subroutine ensure_la()
    Implicit none
    if (.not. la_up) then
       call la_init()
       la_up = .true.
    endif
  end subroutine ensure_la

!--------------------------------------------------------------------
  ! 1 mmthr: left, mat, elements descending;  2 mmthr_m1: left, invmat, ascending;  3 mmthl: right, mat, ascending;
  ! 4 mmthlc: left, mat, ascending;  5 mmthl_m1: right, invmat, descending
  Subroutine gpu_mmthr(M, ld, ncol, nf)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    call hop_left(M, ld, ncol, nf, 1, .true.)
  end subroutine gpu_mmthr

  Subroutine gpu_mmthr_m1(M, ld, ncol, nf)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    call hop_left(M, ld, ncol, nf, 2, .false.)
  end subroutine gpu_mmthr_m1

  Subroutine gpu_mmthlc(M, ld, ncol, nf)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    call hop_left(M, ld, ncol, nf, 1, .false.)
  end subroutine gpu_mmthlc

  Subroutine gpu_mmthl(M, ld, nrow, nf)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    call hop_right(M, ld, nrow, nf, 1, .false.)
  end subroutine gpu_mmthl

  Subroutine gpu_mmthl_m1(M, ld, nrow, nf)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    call hop_right(M, ld, nrow, nf, 2, .true.)
  end subroutine gpu_mmthl_m1

!--------------------------------------------------------------------
  Subroutine hop_left(M, ld, ncol, nf, k, desc)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, k
    Logical, intent(in) :: desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: ga, gb, g, first, last, step, e
    Logical :: dsc

    ga = hop_gstart(nf); gb = hop_gstart(nf+1) - 1
    if (gb < ga .or. ncol <= 0) return
    dsc = desc
    if (fault_order) dsc = .not. desc
    if (dsc) then
       first = gb; last = ga; step = -1
    else
       first = ga; last = gb; step = 1
    endif
    if (.not. hop_big(nf)) then
       if (gpu_cb_fused .and. fused_ok) then
          call k_cb_left(M, ld, ncol, first, last, step, k, gpu_ct)
          gpu_n_launch = gpu_n_launch + 1
       else
          do g = first, last, step
             call k_cb_left_g(M, ld, ncol, g, k, gpu_tl)
             gpu_n_launch = gpu_n_launch + 1
          enddo
       endif
    else
       do g = first, last, step
          e = hop_gptr(g)
          if (hop_meta(2,e) > 8) then
             call dense_left(M, ld, ncol, e, k)
          else
             call k_cb_left_g(M, ld, ncol, g, k, gpu_tl)
             gpu_n_launch = gpu_n_launch + 1
          endif
       enddo
    endif
  end subroutine hop_left

  Subroutine hop_right(M, ld, nrow, nf, k, desc)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf, k
    Logical, intent(in) :: desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: ga, gb, g, first, last, step, e
    Logical :: dsc

    ga = hop_gstart(nf); gb = hop_gstart(nf+1) - 1
    if (gb < ga .or. nrow <= 0) return
    dsc = desc
    if (fault_order) dsc = .not. desc
    if (dsc) then
       first = gb; last = ga; step = -1
    else
       first = ga; last = gb; step = 1
    endif
    if (.not. hop_big(nf)) then
       if (gpu_cb_fused .and. fused_ok) then
          call k_cb_right(M, ld, nrow, first, last, step, k, gpu_rbt)
          gpu_n_launch = gpu_n_launch + 1
       else
          do g = first, last, step
             call k_cb_right_g(M, ld, nrow, g, k, gpu_tr)
             gpu_n_launch = gpu_n_launch + 1
          enddo
       endif
    else
       do g = first, last, step
          e = hop_gptr(g)
          if (hop_meta(2,e) > 8) then
             call dense_right(M, ld, nrow, e, k)
          else
             call k_cb_right_g(M, ld, nrow, g, k, gpu_tr)
             gpu_n_launch = gpu_n_launch + 1
          endif
       enddo
    endif
  end subroutine hop_right

!--------------------------------------------------------------------
!> Element bodies (device routines): apply element e (matrix k) to the columns j1..j2 of M (left: M(P,j) = A M(P,j)) or
!> to the rows i1..i2 of M (right: M(i,P) = M(i,P) A).  The element data are read once per call and reused for all
!> columns/rows of the range.  Arrays are passed explicitly so that the routines only see device pointers.
  Subroutine el_left_d(M, ld, j1, j2, e, k, hmeta, hP, hR, hC)
    Implicit none
    !$omp declare target
    Integer, intent(in) :: ld, j1, j2, e, k
    Integer, intent(in) :: hmeta(4,*), hP(*)
    Real    (Kind=Kind(0.d0)), intent(in) :: hR(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: hC(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: nd, po, mo, i, l, r, r1, r2, j
    Real    (Kind=Kind(0.d0)) :: a, zr, zi, a11, a12, a21, a22
    Complex (Kind=Kind(0.d0)) :: w(8), z, c11, c12, c21, c22, w1, w2

    nd = hmeta(2,e); po = hmeta(3,e); mo = hmeta(4,e) + (k-1)*nd*nd
    if (hmeta(1,e) == 1) then
       if (nd == 1) then
          r = hP(po+1)
          a = hR(mo+1)
          do j = j1, j2
             M(r,j) = cmplx(a*real(M(r,j)), a*aimag(M(r,j)), kind(0.d0))
          enddo
       else if (nd == 2) then
          r1 = hP(po+1); r2 = hP(po+2)
          a11 = hR(mo+1); a21 = hR(mo+2); a12 = hR(mo+3); a22 = hR(mo+4)
          do j = j1, j2
             w1 = M(r1,j); w2 = M(r2,j)
             zr = a11*real(w1) + a12*real(w2); zi = a11*aimag(w1) + a12*aimag(w2)
             M(r1,j) = cmplx(zr, zi, kind(0.d0))
             zr = a21*real(w1) + a22*real(w2); zi = a21*aimag(w1) + a22*aimag(w2)
             M(r2,j) = cmplx(zr, zi, kind(0.d0))
          enddo
       else
          do j = j1, j2
             do l = 1, nd
                w(l) = M(hP(po+l), j)
             enddo
             do i = 1, nd
                a = hR(mo + i)
                zr = a*real(w(1)); zi = a*aimag(w(1))
                do l = 2, nd
                   a = hR(mo + i + (l-1)*nd)
                   zr = zr + a*real(w(l)); zi = zi + a*aimag(w(l))
                enddo
                M(hP(po+i), j) = cmplx(zr, zi, kind(0.d0))
             enddo
          enddo
       endif
    else
       if (nd == 1) then
          r = hP(po+1)
          c11 = hC(mo+1)
          do j = j1, j2
             M(r,j) = c11 * M(r,j)
          enddo
       else if (nd == 2) then
          r1 = hP(po+1); r2 = hP(po+2)
          c11 = hC(mo+1); c21 = hC(mo+2); c12 = hC(mo+3); c22 = hC(mo+4)
          do j = j1, j2
             w1 = M(r1,j); w2 = M(r2,j)
             M(r1,j) = c11*w1 + c12*w2
             M(r2,j) = c21*w1 + c22*w2
          enddo
       else
          do j = j1, j2
             do l = 1, nd
                w(l) = M(hP(po+l), j)
             enddo
             do i = 1, nd
                z = hC(mo + i) * w(1)
                do l = 2, nd
                   z = z + hC(mo + i + (l-1)*nd) * w(l)
                enddo
                M(hP(po+i), j) = z
             enddo
          enddo
       endif
    endif
  end subroutine el_left_d

  Subroutine el_right_d(M, ld, i1, i2, istep, e, k, hmeta, hP, hR, hC)
    Implicit none
    !$omp declare target
    Integer, intent(in) :: ld, i1, i2, istep, e, k
    Integer, intent(in) :: hmeta(4,*), hP(*)
    Real    (Kind=Kind(0.d0)), intent(in) :: hR(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: hC(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: nd, po, mo, j, l, c, c1, c2, i
    Real    (Kind=Kind(0.d0)) :: a, zr, zi, a11, a12, a21, a22
    Complex (Kind=Kind(0.d0)) :: w(8), z, c11, c12, c21, c22, w1, w2

    nd = hmeta(2,e); po = hmeta(3,e); mo = hmeta(4,e) + (k-1)*nd*nd
    if (hmeta(1,e) == 1) then
       if (nd == 1) then
          c = hP(po+1)
          a = hR(mo+1)
          do i = i1, i2, istep
             M(i,c) = cmplx(a*real(M(i,c)), a*aimag(M(i,c)), kind(0.d0))
          enddo
       else if (nd == 2) then
          c1 = hP(po+1); c2 = hP(po+2)
          a11 = hR(mo+1); a21 = hR(mo+2); a12 = hR(mo+3); a22 = hR(mo+4)
          do i = i1, i2, istep
             w1 = M(i,c1); w2 = M(i,c2)
             zr = a11*real(w1) + a21*real(w2); zi = a11*aimag(w1) + a21*aimag(w2)
             M(i,c1) = cmplx(zr, zi, kind(0.d0))
             zr = a12*real(w1) + a22*real(w2); zi = a12*aimag(w1) + a22*aimag(w2)
             M(i,c2) = cmplx(zr, zi, kind(0.d0))
          enddo
       else
          do i = i1, i2, istep
             do l = 1, nd
                w(l) = M(i, hP(po+l))
             enddo
             do j = 1, nd
                a = hR(mo + 1 + (j-1)*nd)
                zr = a*real(w(1)); zi = a*aimag(w(1))
                do l = 2, nd
                   a = hR(mo + l + (j-1)*nd)
                   zr = zr + a*real(w(l)); zi = zi + a*aimag(w(l))
                enddo
                M(i, hP(po+j)) = cmplx(zr, zi, kind(0.d0))
             enddo
          enddo
       endif
    else
       if (nd == 1) then
          c = hP(po+1)
          c11 = hC(mo+1)
          do i = i1, i2, istep
             M(i,c) = M(i,c) * c11
          enddo
       else if (nd == 2) then
          c1 = hP(po+1); c2 = hP(po+2)
          c11 = hC(mo+1); c21 = hC(mo+2); c12 = hC(mo+3); c22 = hC(mo+4)
          do i = i1, i2, istep
             w1 = M(i,c1); w2 = M(i,c2)
             M(i,c1) = w1*c11 + w2*c21
             M(i,c2) = w1*c12 + w2*c22
          enddo
       else
          do i = i1, i2, istep
             do l = 1, nd
                w(l) = M(i, hP(po+l))
             enddo
             do j = 1, nd
                z = w(1) * hC(mo + 1 + (j-1)*nd)
                do l = 2, nd
                   z = z + w(l) * hC(mo + l + (j-1)*nd)
                enddo
                M(i, hP(po+j)) = z
             enddo
          enddo
       endif
    endif
  end subroutine el_right_d

!--------------------------------------------------------------------
!> FUSED, left (gpu_cb_fused = .true.): groups ga..gb (step gstep) for all columns in ONE launch; teams over tiles of CT
!> columns, inside a team the groups run in sequence (the parallel loop of a group ends in a team barrier).  Measured
!> slower than one launch per group on the A30 (see BENCHMARKS / ticket 204), kept for comparison.
  Subroutine k_cb_left(M, ld, ncol, ga, gb, gstep, k, ct)
    Implicit none
    Integer, intent(in) :: ld, ncol, ga, gb, gstep, k, ct
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: it, ntile, j0, jj, j, g, e, e1, e2

    ntile = (ncol + ct - 1) / ct
    !$omp target teams distribute private(j0, g, e1, e2)
    do it = 1, ntile
       j0 = (it-1)*ct
       do g = ga, gb, gstep
          e1 = hop_gptr(g); e2 = hop_gptr(g+1) - 1
          !$omp parallel do collapse(2) private(j)
          do jj = 1, ct
             do e = e1, e2
                j = j0 + jj
                if (j <= ncol) call el_left_d(M, ld, j, j, e, k, hop_meta, hop_P, hop_Rval, hop_Cval)
             enddo
          enddo
       enddo
    enddo
  end subroutine k_cb_left

!> FUSED, right: teams over tiles of RBT rows.
  Subroutine k_cb_right(M, ld, nrow, ga, gb, gstep, k, rbt)
    Implicit none
    Integer, intent(in) :: ld, nrow, ga, gb, gstep, k, rbt
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: it, ntile, i0, ii, i, g, e, e1, e2

    ntile = (nrow + rbt - 1) / rbt
    !$omp target teams distribute private(i0, g, e1, e2)
    do it = 1, ntile
       i0 = (it-1)*rbt
       do g = ga, gb, gstep
          e1 = hop_gptr(g); e2 = hop_gptr(g+1) - 1
          !$omp parallel do collapse(2) private(i)
          do e = e1, e2
             do ii = 1, rbt
                i = i0 + ii
                if (i <= nrow) call el_right_d(M, ld, i, i, 1, e, k, hop_meta, hop_P, hop_Rval, hop_Cval)
             enddo
          enddo
       enddo
    enddo
  end subroutine k_cb_right

!> PER GROUP (one launch per group), flat: threads over (column tile, element) / (element, row tile); a thread handles
!> tl columns (left: consecutive) or rows (right: strided by the tile count), the element data are read once per thread.
  Subroutine k_cb_left_g(M, ld, ncol, g, k, tl)
    Implicit none
    Integer, intent(in) :: ld, ncol, g, k, tl
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: jt, ntile, e, e1, e2, j1, j2
    e1 = hop_gptr(g); e2 = hop_gptr(g+1) - 1
    ntile = (ncol + tl - 1) / tl
    !$omp target teams distribute parallel do collapse(2) private(j1, j2)
    do jt = 1, ntile
       do e = e1, e2
          j1 = (jt-1)*tl + 1
          j2 = min(ncol, jt*tl)
          call el_left_d(M, ld, j1, j2, e, k, hop_meta, hop_P, hop_Rval, hop_Cval)
       enddo
    enddo
  end subroutine k_cb_left_g

  Subroutine k_cb_right_g(M, ld, nrow, g, k, tl)
    Implicit none
    Integer, intent(in) :: ld, nrow, g, k, tl
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: it, ntile, e, e1, e2
    e1 = hop_gptr(g); e2 = hop_gptr(g+1) - 1
    ntile = (nrow + tl - 1) / tl
    ! thread (e, it) handles the rows it, it+ntile, it+2*ntile, ...: neighbouring threads touch neighbouring rows
    !$omp target teams distribute parallel do collapse(2)
    do e = e1, e2
       do it = 1, ntile
          call el_right_d(M, ld, it, nrow, ntile, e, k, hop_meta, hop_P, hop_Rval, hop_Cval)
       enddo
    enddo
  end subroutine k_cb_right_g

!--------------------------------------------------------------------
!  Dense elements
!--------------------------------------------------------------------
  Subroutine dense_left(M, ld, ncol, e, k)
    Implicit none
    Integer, intent(in) :: ld, ncol, e, k
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: nd, po, mo, xo, pid, mode
    nd = hop_meta(2,e); po = hop_meta(3,e); mo = hop_meta(4,e) + (k-1)*nd*nd
    xo = hop_aux(1,e); pid = hop_aux(2,e)
    if (nd > Ndim .or. ncol > Ndim) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    if (.not. la_available) then
       write(error_unit,*) 'gpu_wrap: dense hopping needs alf_la (ticket 302 for GPU=intel)'
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    endif
    call ensure_la()
    mode = gpu_dense_left_mode
    if (mode == 0) mode = merge(1, 2, nd <= 400)
    if (hop_meta(1,e) == 1 .and. mode == 2) then
       ! (A M)^T = M^T A^T: transpose-gather, DGEMM on the interleaved view with A^T, transpose-scatter
       call k_tgather(nd, ncol, hop_P(po+1), M, ld, scr_c(1,1,1))               ! scr1(j,l) = M(P(l),j), ld = ncol
       call la_gemm_cr(ncol, nd, nd, scr_c(1,1,1), ncol, hop_At(xo + (k-1)*nd*nd + 1), nd, scr_c(1,1,2), ncol)
       call k_tscatter(nd, ncol, hop_P(po+1), M, ld, scr_c(1,1,2))              ! M(P(i),j) = scr2(j,i)
       gpu_n_launch = gpu_n_launch + 2; gpu_n_lib = gpu_n_lib + 1
    else
       if (pid == 1) then
          call gemm_left(e, k, nd, ncol, mo, xo, M, ld, scr_c(1,1,2))
       else
          call k_gather_rows(nd, ncol, hop_P(po+1), M, ld, scr_c(1,1,1), Ndim)
          call gemm_left(e, k, nd, ncol, mo, xo, scr_c(1,1,1), Ndim, scr_c(1,1,2))
          gpu_n_launch = gpu_n_launch + 1
       endif
       call k_scatter_rows(nd, ncol, hop_P(po+1), M, ld, scr_c(1,1,2), Ndim)
       gpu_n_launch = gpu_n_launch + 1; gpu_n_lib = gpu_n_lib + 1
    endif
  end subroutine dense_left

  ! C(nd x ncol, ld Ndim) = A * B(nd x ncol, ld ldb) with A the element matrix as complex (complex element or the
  ! complex copy of the real one)
  Subroutine gemm_left(e, k, nd, ncol, mo, xo, B, ldb, C)
    Implicit none
    Integer, intent(in) :: e, k, nd, ncol, mo, xo, ldb
    Complex (Kind=Kind(0.d0)), intent(inout) :: B(ldb,*), C(Ndim,*)
    if (hop_meta(1,e) == 2) then
       call la_gemm(LA_N, LA_N, nd, ncol, nd, z1, hop_Cval(mo+1), nd, B, ldb, z0, C, Ndim)
    else
       call la_gemm(LA_N, LA_N, nd, ncol, nd, z1, hop_Ac(xo + (k-1)*nd*nd + 1), nd, B, ldb, z0, C, Ndim)
    endif
  end subroutine gemm_left

  Subroutine dense_right(M, ld, nrow, e, k)
    Implicit none
    Integer, intent(in) :: ld, nrow, e, k
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: nd, po, mo, pid
    nd = hop_meta(2,e); po = hop_meta(3,e); mo = hop_meta(4,e) + (k-1)*nd*nd
    pid = hop_aux(2,e)
    if (nd > Ndim .or. nrow > Ndim) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    if (.not. la_available) then
       write(error_unit,*) 'gpu_wrap: dense hopping needs alf_la (ticket 302 for GPU=intel)'
       call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    endif
    call ensure_la()
    if (pid == 1) then
       call gemm_right(e, nd, nrow, mo, M, ld, scr_c(1,1,2))
    else
       call k_gather_cols(nrow, nd, hop_P(po+1), M, ld, scr_c(1,1,1))          ! scr1(i,l) = M(i,P(l)), ld = nrow
       call gemm_right(e, nd, nrow, mo, scr_c(1,1,1), nrow, scr_c(1,1,2))
       gpu_n_launch = gpu_n_launch + 1
    endif
    call k_scatter_cols(nrow, nd, hop_P(po+1), M, ld, scr_c(1,1,2))             ! M(i,P(j)) = scr2(i,j), ld = nrow
    gpu_n_launch = gpu_n_launch + 1; gpu_n_lib = gpu_n_lib + 1
  end subroutine dense_right

  ! C(nrow x nd, ld nrow) = B(nrow x nd, ld ldb) * A
  Subroutine gemm_right(e, nd, nrow, mo, B, ldb, C)
    Implicit none
    Integer, intent(in) :: e, nd, nrow, mo, ldb
    Complex (Kind=Kind(0.d0)), intent(inout) :: B(ldb,*), C(nrow,*)
    if (hop_meta(1,e) == 2) then
       call la_gemm(LA_N, LA_N, nrow, nd, nd, z1, B, ldb, hop_Cval(mo+1), nd, z0, C, nrow)
    else
       call la_gemm_cr(nrow, nd, nd, B, ldb, hop_Rval(mo+1), nd, C, nrow)
    endif
  end subroutine gemm_right

  Subroutine k_gather_rows(nd, ncol, P, M, ld, W, ldw)
    Implicit none
    Integer, intent(in) :: nd, ncol, ld, ldw
    Integer, intent(in) :: P(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,*)
    Integer :: j, l
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, ncol
       do l = 1, nd
          W(l,j) = M(P(l), j)
       enddo
    enddo
  end subroutine k_gather_rows

  Subroutine k_scatter_rows(nd, ncol, P, M, ld, W, ldw)
    Implicit none
    Integer, intent(in) :: nd, ncol, ld, ldw
    Integer, intent(in) :: P(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(ldw,*)
    Integer :: j, l
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, ncol
       do l = 1, nd
          M(P(l), j) = W(l,j)
       enddo
    enddo
  end subroutine k_scatter_rows

  ! W(j,l) = M(P(l), j)   (W is ncol x nd, ld ncol)
  Subroutine k_tgather(nd, ncol, P, M, ld, W)
    Implicit none
    Integer, intent(in) :: nd, ncol, ld
    Integer, intent(in) :: P(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ncol,*)
    Integer :: j, l
    !$omp target teams distribute parallel do collapse(2)
    do l = 1, nd
       do j = 1, ncol
          W(j,l) = M(P(l), j)
       enddo
    enddo
  end subroutine k_tgather

  ! M(P(l), j) = W(j,l)
  Subroutine k_tscatter(nd, ncol, P, M, ld, W)
    Implicit none
    Integer, intent(in) :: nd, ncol, ld
    Integer, intent(in) :: P(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(ncol,*)
    Integer :: j, l
    !$omp target teams distribute parallel do collapse(2)
    do l = 1, nd
       do j = 1, ncol
          M(P(l), j) = W(j,l)
       enddo
    enddo
  end subroutine k_tscatter

  ! W(i,l) = M(i, P(l))   (W is nrow x nd, ld nrow)
  Subroutine k_gather_cols(nrow, nd, P, M, ld, W)
    Implicit none
    Integer, intent(in) :: nrow, nd, ld
    Integer, intent(in) :: P(*)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(nrow,*)
    Integer :: i, l
    !$omp target teams distribute parallel do collapse(2)
    do l = 1, nd
       do i = 1, nrow
          W(i,l) = M(i, P(l))
       enddo
    enddo
  end subroutine k_gather_cols

  Subroutine k_scatter_cols(nrow, nd, P, M, ld, W)
    Implicit none
    Integer, intent(in) :: nrow, nd, ld
    Integer, intent(in) :: P(*)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(nrow,*)
    Integer :: i, l
    !$omp target teams distribute parallel do collapse(2)
    do l = 1, nd
       do i = 1, nrow
          M(i, P(l)) = W(i,l)
       enddo
    enddo
  end subroutine k_scatter_cols

!====================================================================
!  Op_V
!====================================================================
!> Op_mmultR(Mat, Op_V(n,nf), nsigma(n,t), cop, t): rows P(i) of M (ld x ncol) scaled by E_exp(i,sp).
  Subroutine gpu_Vop_R(M, ld, ncol, n, nf, t, cop)
    Implicit none
    Integer, intent(in) :: ld, ncol, n, nf, t
    Character (len=1), intent(in) :: cop
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Logical :: cj
    if (V_meta(4,n,nf) == 1 .or. ncol <= 0) return        ! Op_mmultR returns early for |g| < Eps_machine
    cj = (cop == 'c' .or. cop == 'C')
    call k_vop_R(M, ld, ncol, n, nf, t, cj, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vop_R

  Subroutine k_vop_R(M, ld, ncol, n, nf, t, cj, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, ncol, n, nf, t, N1, po, eo, nt
    Logical, intent(in) :: cj
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: i, j, sp, r
    Complex (Kind=Kind(0.d0)) :: z
    !$omp target teams distribute parallel do collapse(2) private(sp, r, z)
    do j = 1, ncol
       do i = 1, N1
          sp = nint(real(nsigma_d(n,t)))
          z = V_exp(eo + i + (sp+nt)*N1)
          if (cj) z = conjg(z)
          r = V_P(po+i)
          M(r,j) = z * M(r,j)
       enddo
    enddo
  end subroutine k_vop_R

!> Op_mmultL(Mat, Op_V(n,nf), nsigma(n,t), cop, t, sign): columns P(i) of M (nrow x ld) scaled by E_exp(i, sign*sp).
  Subroutine gpu_Vop_L(M, ld, nrow, n, nf, t, cop, sign)
    Implicit none
    Integer, intent(in) :: ld, nrow, n, nf, t, sign
    Character (len=1), intent(in) :: cop
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Logical :: cj
    if (V_meta(4,n,nf) == 1 .or. nrow <= 0) return        ! Op_mmultL returns early for |g| < Eps_machine
    cj = (cop == 'c' .or. cop == 'C')
    call k_vop_L(M, ld, nrow, n, nf, t, cj, sign, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vop_L

  Subroutine k_vop_L(M, ld, nrow, n, nf, t, cj, sign, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nrow, n, nf, t, sign, N1, po, eo, nt
    Logical, intent(in) :: cj
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: i, c, sp, cc
    Complex (Kind=Kind(0.d0)) :: z
    !$omp target teams distribute parallel do collapse(2) private(sp, cc, z)
    do c = 1, N1
       do i = 1, nrow
          sp = sign*nint(real(nsigma_d(n,t)))
          z = V_exp(eo + c + (sp+nt)*N1)
          if (cj) z = conjg(z)
          cc = V_P(po+c)
          M(i,cc) = z * M(i,cc)
       enddo
    enddo
  end subroutine k_vop_L

!> Op_Wrapup (dir = +1) / Op_Wrapdo (dir = -1) for operator n at slice t on G (ld x Ndim): rows by E_exp(i, dir*sp), then
!> columns by E_exp(i, -dir*sp).  N_Type 2 applies nothing for a diagonal operator.  No skip check (as ALF).
  Subroutine gpu_Vop_wrap(G, ld, n, nf, t, dir, ntype)
    Implicit none
    Integer, intent(in) :: ld, n, nf, t, dir, ntype
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    if (ntype /= 1) return
    call k_vop_wrap_rows(G, ld, Ndim, n, t, dir, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    call k_vop_wrap_cols(G, ld, Ndim, n, t, dir, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    gpu_n_launch = gpu_n_launch + 2
  end subroutine gpu_Vop_wrap

  Subroutine k_vop_wrap_rows(G, ld, nd, n, t, dir, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nd, n, t, dir, N1, po, eo, nt
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Integer :: i, j, sp, r
    Complex (Kind=Kind(0.d0)) :: z
    !$omp target teams distribute parallel do collapse(2) private(sp, r, z)
    do j = 1, nd
       do i = 1, N1
          sp = nint(real(nsigma_d(n,t)))
          z = V_exp(eo + i + (dir*sp+nt)*N1)
          r = V_P(po+i)
          G(r,j) = z * G(r,j)
       enddo
    enddo
  end subroutine k_vop_wrap_rows

  Subroutine k_vop_wrap_cols(G, ld, nd, n, t, dir, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nd, n, t, dir, N1, po, eo, nt
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Integer :: i, c, sp, cc
    Complex (Kind=Kind(0.d0)) :: z
    !$omp target teams distribute parallel do collapse(2) private(sp, cc, z)
    do c = 1, N1
       do i = 1, nd
          sp = nint(real(nsigma_d(n,t)))
          z = V_exp(eo + c + (-dir*sp+nt)*N1)
          cc = V_P(po+c)
          G(i,cc) = z * G(i,cc)
       enddo
    enddo
  end subroutine k_vop_wrap_cols

!--------------------------------------------------------------------
!> Whole slice, rows: for every operator n (ascending, or descending if desc) Op_mmultR(M, Op_V(n,nf), nsigma(n,t),
!> cop, t), in one kernel.  Per row r the factors of the operators containing r are applied to each element in that
!> order (skipped operators excluded): the same product sequence as the loop of calls.
  Subroutine gpu_Vslice_R(M, ld, ncol, nf, t, cop, desc)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, t
    Character (len=1), intent(in) :: cop
    Logical, intent(in) :: desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (ncol <= 0) return
    call k_vslice_R(M, ld, ncol, Ndim, nf, t, (cop == 'c' .or. cop == 'C'), desc, gpu_vb)
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vslice_R

  Subroutine k_vslice_R(M, ld, ncol, nd, nf, t, cj, desc, vb)
    Implicit none
    Integer, intent(in) :: ld, ncol, nd, nf, t, vb
    Logical, intent(in) :: cj, desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: jb, nblk, j, j1, j2, r, pa, pb, p, q, n, i, sp, nact
    Complex (Kind=Kind(0.d0)) :: x, z, zq(4)

    nblk = (ncol + vb - 1) / vb
    ! thread (r, column block): the factors of row r are looked up once per thread when there are at most 4 of them
    ! (always in the models with one operator per site); more: per element, the same products in the same order
    !$omp target teams distribute parallel do collapse(2) private(j1, j2, j, pa, pb, p, q, n, i, sp, nact, x, z, zq)
    do jb = 1, nblk
       do r = 1, nd
          pa = vr_ptr(r,nf); pb = vr_ptr(r+1,nf) - 1
          j1 = (jb-1)*vb + 1; j2 = min(ncol, jb*vb)
          if (pb >= pa) then
             if (pb - pa < 4) then
                nact = 0
                do q = 0, pb - pa
                   if (desc) then
                      p = pb - q
                   else
                      p = pa + q
                   endif
                   n = vr_ent(1,p); i = vr_ent(2,p)
                   if (V_meta(4,n,nf) /= 1) then
                      sp = nint(real(nsigma_d(n,t)))
                      z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                      if (cj) z = conjg(z)
                      nact = nact + 1
                      zq(nact) = z
                   endif
                enddo
                do j = j1, j2
                   x = M(r,j)
                   do q = 1, nact
                      x = zq(q) * x
                   enddo
                   M(r,j) = x
                enddo
             else
                do j = j1, j2
                   x = M(r,j)
                   do q = 0, pb - pa
                      if (desc) then
                         p = pb - q
                      else
                         p = pa + q
                      endif
                      n = vr_ent(1,p); i = vr_ent(2,p)
                      if (V_meta(4,n,nf) /= 1) then
                         sp = nint(real(nsigma_d(n,t)))
                         z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                         if (cj) z = conjg(z)
                         x = z * x
                      endif
                   enddo
                   M(r,j) = x
                enddo
             endif
          endif
       enddo
    enddo
  end subroutine k_vslice_R

!> Whole slice, columns: Op_mmultL(M, Op_V(n,nf), nsigma(n,t), cop, t, sign) for all n in order.
  Subroutine gpu_Vslice_L(M, ld, nrow, nf, t, cop, sign, desc)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf, t, sign
    Character (len=1), intent(in) :: cop
    Logical, intent(in) :: desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (nrow <= 0) return
    call k_vslice_L(M, ld, nrow, Ndim, nf, t, (cop == 'c' .or. cop == 'C'), sign, desc, gpu_vb)
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vslice_L

  Subroutine k_vslice_L(M, ld, nrow, nd, nf, t, cj, sign, desc, vb)
    Implicit none
    Integer, intent(in) :: ld, nrow, nd, nf, t, sign, vb
    Logical, intent(in) :: cj, desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: ib, nblk, c, i1, i2, i0, pa, pb, p, q, n, i, sp, nact
    Complex (Kind=Kind(0.d0)) :: x, z, zq(4)

    nblk = (nrow + vb - 1) / vb
    !$omp target teams distribute parallel do collapse(2) private(i1, i2, i0, pa, pb, p, q, n, i, sp, nact, x, z, zq)
    do ib = 1, nblk
       do c = 1, nd
          pa = vr_ptr(c,nf); pb = vr_ptr(c+1,nf) - 1
          i1 = (ib-1)*vb + 1; i2 = min(nrow, ib*vb)
          if (pb >= pa) then
             if (pb - pa < 4) then
                nact = 0
                do q = 0, pb - pa
                   if (desc) then
                      p = pb - q
                   else
                      p = pa + q
                   endif
                   n = vr_ent(1,p); i = vr_ent(2,p)
                   if (V_meta(4,n,nf) /= 1) then
                      sp = sign*nint(real(nsigma_d(n,t)))
                      z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                      if (cj) z = conjg(z)
                      nact = nact + 1
                      zq(nact) = z
                   endif
                enddo
                do i0 = i1, i2
                   x = M(i0,c)
                   do q = 1, nact
                      x = zq(q) * x
                   enddo
                   M(i0,c) = x
                enddo
             else
                do i0 = i1, i2
                   x = M(i0,c)
                   do q = 0, pb - pa
                      if (desc) then
                         p = pb - q
                      else
                         p = pa + q
                      endif
                      n = vr_ent(1,p); i = vr_ent(2,p)
                      if (V_meta(4,n,nf) /= 1) then
                         sp = sign*nint(real(nsigma_d(n,t)))
                         z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                         if (cj) z = conjg(z)
                         x = z * x
                      endif
                   enddo
                   M(i0,c) = x
                enddo
             endif
          endif
       enddo
    enddo
  end subroutine k_vslice_L

!> Whole slice wrap, G (ld x Ndim): the loop over n of Op_Wrapup (dir = +1, n ascending) or Op_Wrapdo (dir = -1, n
!> descending) with N_Type = 1, in one kernel.  Element (r,c): the row factors of the operators containing r and the
!> column factors of those containing c, merged in operator order (row before column for the same operator).
  Subroutine gpu_Vwrap_all(G, ld, nf, t, dir)
    Implicit none
    Integer, intent(in) :: ld, nf, t, dir
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    call k_vwrap_all(G, ld, Ndim, nf, t, dir)
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vwrap_all

  Subroutine k_vwrap_all(G, ld, nd, nf, t, dir)
    Implicit none
    Integer, intent(in) :: ld, nd, nf, t, dir
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Integer, parameter :: big = 1073741824
    Integer :: r, c, pr, pre, pc, pce, ir, ic, nr, nc, n, i, sp, cnt
    Complex (Kind=Kind(0.d0)) :: x, z
    !$omp target teams distribute parallel do collapse(2) private(pr, pre, pc, pce, ir, ic, nr, nc, n, i, sp, cnt, x, z)
    do c = 1, nd
       do r = 1, nd
          pr = vr_ptr(r,nf); pre = vr_ptr(r+1,nf) - 1
          pc = vr_ptr(c,nf); pce = vr_ptr(c+1,nf) - 1
          if (pre >= pr .or. pce >= pc) then
             x = G(r,c)
             if (dir > 0) then
                ir = pr; ic = pc
                do cnt = 1, (pre - pr + 1) + (pce - pc + 1)
                   nr = big; nc = big
                   if (ir <= pre) nr = vr_ent(1,ir)
                   if (ic <= pce) nc = vr_ent(1,ic)
                   if (nr <= nc) then
                      n = nr; i = vr_ent(2,ir); ir = ir + 1
                      sp = nint(real(nsigma_d(n,t)))
                      z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                   else
                      n = nc; i = vr_ent(2,ic); ic = ic + 1
                      sp = nint(real(nsigma_d(n,t)))
                      z = V_exp(V_meta(3,n,nf) + i + (-sp+V_type(n))*V_meta(1,n,nf))
                   endif
                   x = z * x
                enddo
             else
                ir = pre; ic = pce
                do cnt = 1, (pre - pr + 1) + (pce - pc + 1)
                   nr = -1; nc = -1
                   if (ir >= pr) nr = vr_ent(1,ir)
                   if (ic >= pc) nc = vr_ent(1,ic)
                   if (nr >= nc) then
                      n = nr; i = vr_ent(2,ir); ir = ir - 1
                      sp = nint(real(nsigma_d(n,t)))
                      z = V_exp(V_meta(3,n,nf) + i + (-sp+V_type(n))*V_meta(1,n,nf))
                   else
                      n = nc; i = vr_ent(2,ic); ic = ic - 1
                      sp = nint(real(nsigma_d(n,t)))
                      z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                   endif
                   x = z * x
                enddo
             endif
             G(r,c) = x
          endif
       enddo
    enddo
  end subroutine k_vwrap_all

!====================================================================
!  Composite slice routines (what wrapur / wrapul / WRAPGRUP / WRAPGRDO do without the decomposition / updates)
!====================================================================
!> WRAPUR propagation: for t = ntau+1..ntau1:  U = B_T U, then Op_mmultR 'n' for n = 1..N_op.   U is Ndim x ncol.
  Subroutine gpu_wrapur_slices(U, ld, ncol, nf, ntau, ntau1)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, ntau, ntau1
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Integer :: t
    do t = ntau + 1, ntau1
       call gpu_mmthr(U, ld, ncol, nf)
       call gpu_Vslice_R(U, ld, ncol, nf, t, 'n', .false.)
    enddo
  end subroutine gpu_wrapur_slices

!> WRAPUL propagation: for t = ntau1..ntau+1 (descending): Op_mmultR 'c' for n = N_op..1, then U = mmthlc(U).
  Subroutine gpu_wrapul_slices(U, ld, ncol, nf, ntau, ntau1)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, ntau, ntau1
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(ld,*)
    Integer :: t
    do t = ntau1, ntau + 1, -1
       call gpu_Vslice_R(U, ld, ncol, nf, t, 'c', .true.)
       call gpu_mmthlc(U, ld, ncol, nf)
    enddo
  end subroutine gpu_wrapul_slices

!> WRAPGRUP without updates: G = B(ntau+1) G B(ntau+1)^-1  (mmthr, mmthl_m1, then Op_Wrapup N_Type 1 for n = 1..N_op)
  Subroutine gpu_wrapgr_up(G, ld, nf, ntau)
    Implicit none
    Integer, intent(in) :: ld, nf, ntau
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    call gpu_mmthr(G, ld, Ndim, nf)
    call gpu_mmthl_m1(G, ld, Ndim, nf)
    call gpu_Vwrap_all(G, ld, nf, ntau + 1, 1)
  end subroutine gpu_wrapgr_up

!> WRAPGRDO without updates: G = B(ntau)^-1 G B(ntau)  (Op_Wrapdo N_Type 1 for n = N_op..1, then mmthl, mmthr_m1)
  Subroutine gpu_wrapgr_do(G, ld, nf, ntau)
    Implicit none
    Integer, intent(in) :: ld, nf, ntau
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    call gpu_Vwrap_all(G, ld, nf, ntau, -1)
    call gpu_mmthl(G, ld, Ndim, nf)
    call gpu_mmthr_m1(G, ld, Ndim, nf)
  end subroutine gpu_wrapgr_do

!====================================================================
!  Hook, self-test, timing
!====================================================================
  Logical function env_is(name, val)
    Implicit none
    Character (len=*), intent(in) :: name, val
    Character (len=64) :: e
    Integer :: st
    e = ''
    call get_environment_variable(name, e, status=st)
    env_is = (st == 0 .and. trim(e) == val)
  end function env_is

  Subroutine env_int(name, v)
    Implicit none
    Character (len=*), intent(in) :: name
    Integer, intent(inout) :: v
    Character (len=32) :: e
    Integer :: st, iv, ios
    e = ''
    call get_environment_variable(name, e, status=st)
    if (st == 0 .and. len_trim(e) > 0) then
       read(e,*,iostat=ios) iv
       if (ios == 0 .and. iv >= 1) v = iv
    endif
  end subroutine env_int

!> Called from main after gpu_state_setup (guarded by ALF_GPU and gpu_active).
  Subroutine gpu_wrap_hook()
    Implicit none
    Character (len=64) :: f
    Integer :: st
    Logical :: did
    did = .false.
    call env_int('ALF_GPU_CT', gpu_ct)         ! tuning of the fused hopping kernels (timing study)
    call env_int('ALF_GPU_RBT', gpu_rbt)
    call env_int('ALF_GPU_TL', gpu_tl)
    call env_int('ALF_GPU_TR', gpu_tr)
    call env_int('ALF_GPU_VB', gpu_vb)
    f = ''
    call get_environment_variable('ALF_GPU_FAULT', f, status=st)
    if (st == 0 .and. len_trim(f) > 0 .and. env_is('ALF_GPU_SELFTEST', '2')) call inject_fault(trim(f))
    if (env_is('ALF_GPU_SELFTEST', '2')) then
       call gpu_wrap_selftest()
       did = .true.
    endif
    if (env_is('ALF_GPU_TIMING', '1')) then
       call gpu_wrap_timing()
       did = .true.
    endif
    if (la_up) then
       call la_finalize()
       la_up = .false.
    endif
    if (did .and. env_is('ALF_GPU_EXIT', '1')) then
       call gpu_state_finalize()
       write(output_unit,'(A)') ' GPU wrap hook: ALF_GPU_EXIT=1, stopping after the test'
       stop
    endif
  end subroutine gpu_wrap_hook

  ! negative controls: change the DEVICE copy only (the CPU reference stays correct)
  Subroutine inject_fault(f)
    Implicit none
    Character (len=*), intent(in) :: f
    Real (Kind=Kind(0.d0)) :: nan
    nan = transfer(9221120237041090560_int64, 0.d0)
    select case (f)
    case ('hop')
       hop_Rval(1) = hop_Rval(1)*(1.d0 + 1.d-9)
       !$omp target update to(hop_Rval)
    case ('hopc')
       hop_Cval(1) = hop_Cval(1) + cmplx(1.d-9, 0.d0, kind(0.d0))
       !$omp target update to(hop_Cval)
    case ('vexp')      ! 1e-6: a 1e-9 change of one table entry is below the 1e-10 Frobenius tolerance (one row of N)
       V_exp(1) = V_exp(1)*(1.d0 + 1.d-6)
       !$omp target update to(V_exp)
    case ('nan')
       hop_Rval(1) = nan
       !$omp target update to(hop_Rval)
    case ('order')
       fault_order = .true.
    case default
       write(error_unit,'(A,A)') ' ALF_GPU_FAULT: unknown value ', f
    end select
    call info_line('GPU wrap self-test: fault injected into the device copy: '//f)
  end subroutine inject_fault

!--------------------------------------------------------------------
  ! relative Frobenius difference |A-B|_F / |B|_F; any non-finite entry gives huge(1.d0) (no tolerance can pass)
  Real (Kind=Kind(0.d0)) function frob_rel(A, B)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:), B(:,:)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: sd, sb, dr, di
    sd = 0.d0; sb = 0.d0
    do j = 1, size(A,2)
       do i = 1, size(A,1)
          if (.not. (finite_c(A(i,j)) .and. finite_c(B(i,j)))) then
             frob_rel = huge(1.d0); return
          endif
          dr = real(A(i,j)) - real(B(i,j)); di = aimag(A(i,j)) - aimag(B(i,j))
          sd = sd + dr*dr + di*di
          sb = sb + real(B(i,j))**2 + aimag(B(i,j))**2
       enddo
    enddo
    frob_rel = sqrt(sd) / max(sqrt(sb), tiny(1.d0))
  end function frob_rel

  ! NaN-safe maximum: a NaN or Inf measure propagates as huge
  Subroutine acc_max(m, d)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: m
    Real (Kind=Kind(0.d0)), intent(in) :: d
    if (.not. is_finite(d)) then
       m = huge(1.d0)
    else if (d > m) then
       m = d
    endif
  end subroutine acc_max

  Subroutine say(line)
    Implicit none
    Character (len=*), intent(in) :: line
    call info_line(line)
    write(output_unit,'(1X,A)') trim(line)
  end subroutine say

  Integer function test_ncolp()
    Implicit none
    test_ncolp = max(1, Ndim/2 + 1)      ! stand-in for N_part (finite-T runs have N_part = Ndim, projector runs less)
  end function test_ncolp

!--------------------------------------------------------------------
!> Apply hopping routine op (1 mmthr, 2 mmthr_m1, 3 mmthl, 4 mmthlc, 5 mmthl_m1) on the host (CPU routine of Hop_mod).
  Subroutine cpu_hop(op, nf, M)
    Implicit none
    Integer, intent(in) :: op, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    select case (op)
    case (1); call Hop_mod_mmthr   (M, nf, 1)
    case (2); call Hop_mod_mmthr_m1(M, nf, 1)
    case (3); call Hop_mod_mmthl   (M, nf, 1)
    case (4); call Hop_mod_mmthlc  (M, nf, 1)
    case (5); call Hop_mod_mmthl_m1(M, nf, 1)
    end select
  end subroutine cpu_hop

  Subroutine dev_hop(op, nf, D, ld, nn)
    Implicit none
    Integer, intent(in) :: op, nf, ld, nn
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(ld,*)
    select case (op)
    case (1); call gpu_mmthr   (D, ld, nn, nf)
    case (2); call gpu_mmthr_m1(D, ld, nn, nf)
    case (3); call gpu_mmthl   (D, ld, nn, nf)
    case (4); call gpu_mmthlc  (D, ld, nn, nf)
    case (5); call gpu_mmthl_m1(D, ld, nn, nf)
    end select
  end subroutine dev_hop

  Subroutine gpu_wrap_selftest()
    Implicit none
    Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-10
    Integer(int64) :: s
    Integer :: nf, op, ic, ish, nrow, ncl, nbad, ncase, nsp, nops, nsyn
    Real (Kind=Kind(0.d0)) :: dh(5,2), dv(8), dsl(4), dwa(2), dwo(2), dskip, dsyn
    Complex (Kind=Kind(0.d0)), allocatable :: A0(:,:), A1(:,:), A2(:,:), Dt(:,:)
    Character (len=512) :: line
    Character (len=10), parameter :: opname(5) = (/ 'mmthr     ', 'mmthr_m1  ', 'mmthl     ', 'mmthlc    ', 'mmthl_m1  ' /)
    Integer :: nw, ncolp

    s = 20251029_int64
    ncolp = test_ncolp()
    nw = get_Nwrap()
    call say('GPU wrap self-test (ticket 204) start')
    gpu_n_launch = 0; gpu_n_lib = 0

    ! ---- inputs finite ----
    nbad = 0
    nbad = n_nonfinite_r(size(hop_Rval), hop_Rval) + n_nonfinite_c(size(hop_Cval), hop_Cval) &
         & + n_nonfinite_c(size(V_exp), V_exp)
    if (nbad /= 0) call selftest_fail('non-finite value in the mirror (wrap self-test)')

    ! ---- structure report ----
    call report_structure()

    ! ---- hopping: 5 routines, two code-path combinations, two shapes ----
    ! combination 1: cb fused (1 launch per B_T) + dense left variant 1 (complex copy, ZGEMM);
    ! combination 2: cb one launch per group + dense left variant 2 (transposes, DGEMM)
    dh = 0.d0; ncase = 0
    do nf = 1, N_FL
       do op = 1, 5
          do ish = 1, 2
             if (ish == 1) then
                nrow = Ndim; ncl = Ndim
             else
                nrow = ncolp; ncl = ncolp
             endif
             ! left routines (1,2,4): matrix Ndim x ncl, right (3,5): nrow x Ndim
             if (op == 3 .or. op == 5) then
                allocate(A0(nrow,Ndim), A1(nrow,Ndim), A2(nrow,Ndim), Dt(nrow,Ndim))
             else
                allocate(A0(Ndim,ncl), A1(Ndim,ncl), A2(Ndim,ncl), Dt(Ndim,ncl))
             endif
             call rand_mat(A0, s)
             A1 = A0
             call cpu_hop(op, nf, A1)
             !$omp target enter data map(alloc: Dt)
             do ic = 1, 2
                if (ic == 1) then
                   gpu_cb_fused = fused_ok;  gpu_dense_left_mode = 1
                else
                   gpu_cb_fused = .false.; gpu_dense_left_mode = 0
                endif
                if (hop_big(nf) .and. .not. la_available) cycle
                Dt = A0
                !$omp target update to(Dt)       ! self-test: upload the random input
                if (op == 3 .or. op == 5) then
                   call dev_hop(op, nf, Dt, size(Dt,1), size(Dt,1))      ! right: nrow rows
                else
                   call dev_hop(op, nf, Dt, size(Dt,1), size(Dt,2))      ! left: ncol columns
                endif
                !$omp target update from(Dt)     ! self-test: download the result for the comparison
                A2 = Dt
                call acc_max(dh(op,ic), frob_rel(A2, A1))
                ncase = ncase + 1
             enddo
             !$omp target exit data map(delete: Dt)
             deallocate(A0, A1, A2, Dt)
          enddo
       enddo
    enddo
    gpu_cb_fused = .false.; gpu_dense_left_mode = 0
    do op = 1, 5
       write(line,'(A,A,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU wrap self-test hopping ', opname(op), &
            & ': max rel Frobenius diff (cb fused / dense-left v1) ', dh(op,1), ', (cb per-group / dense-left v2) ', dh(op,2), &
            & ' (tol ', tol
       call say(line)
       if (.not. (dh(op,1) <= tol .and. dh(op,2) <= tol)) call selftest_fail('hopping '//trim(opname(op)))
    enddo
    if (.not. la_available) call say('GPU wrap self-test: dense hopping paths SKIPPED (alf_la backend for this build not available)')

    ! ---- Op_V: per operator ----
    call test_vop(s, dv, dwa, dwo, dskip, nsp, nops)
    write(line,'(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU wrap self-test Op_V per operator (', nops, &
         & ' operators of ', N_FL*size(Op_V,1), ' tested, ', nsp, ' operator/field cases): Op_mmultR n ', dv(1), ', c ', dv(2), &
         & '; Op_mmultL n ', dv(3), ', c ', dv(4), ' (tol ', tol
    call say(line)
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU wrap self-test Op_V per operator: Op_Wrapup ', dwa(1), &
         & ', Op_Wrapdo ', dwo(1), ', N_Type 2 (no-op) ', max(dwa(2),dwo(2)), '; tiny-g copy (L/R skip honoured, wrap ignores) ', &
         & dskip, ' (tol ', tol
    call say(line)
    if (.not. (dv(1) <= tol .and. dv(2) <= tol .and. dv(3) <= tol .and. dv(4) <= tol .and. dwa(1) <= tol .and. dwo(1) <= tol &
         & .and. dwa(2) <= tol .and. dwo(2) <= tol .and. dskip <= tol)) call selftest_fail('Op_V per operator')

    ! ---- fused slice kernels ----
    call test_slice(s, dsl, dwa(1), dwo(1))
    write(line,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU wrap self-test fused slice: Op_mmultR rows (4 variants) ', &
         & dsl(1), ', Op_mmultL columns (8 variants) ', dsl(2), ', tiny-g skip ', dsl(3), '; slice wraps up/down ', &
         & max(dwa(1),dwo(1)), ' (tol ', tol
    call say(line)
    if (.not. (dsl(1) <= tol .and. dsl(2) <= tol .and. dsl(3) <= tol .and. dwa(1) <= tol .and. dwo(1) <= tol)) &
         & call selftest_fail('fused Op_V slice kernels')

    ! ---- synthetic operators with overlapping P (the Hubbard operators are single-site and disjoint) ----
    call test_synthetic_V(s, dsyn, nsyn)
    if (nsyn > 0) then
       write(line,'(A,I0,A,ES10.3,A,ES8.1)') 'GPU wrap self-test synthetic overlapping operators (', nsyn, &
            & ' operators, shared sites, mixed field types, one skipped; fused R/L/wrap kernels vs sequential per-operator products) max rel diff ', dsyn, ' (tol ', tol
       call say(line)
       if (.not. (dsyn <= tol)) call selftest_fail('fused kernels with overlapping operators')
    else
       call say('GPU wrap self-test synthetic overlapping operators: SKIPPED (too few operators/table entries in this model)')
    endif

    ! ---- sequences over Nwrap slices ----
    call test_sequences(s, tol, nw, ncolp)

    write(line,'(A,I0,A,I0)') 'GPU wrap self-test device launches issued ', gpu_n_launch, ', alf_la calls ', gpu_n_lib
    call say(line)
    call say('GPU wrap self-test PASSED')
  end subroutine gpu_wrap_selftest

!--------------------------------------------------------------------
  Integer function n_nonfinite_r(n, x)
    Implicit none
    Integer, intent(in) :: n
    Real (Kind=Kind(0.d0)), intent(in) :: x(n)
    Integer :: i
    n_nonfinite_r = 0
    do i = 1, n
       if (.not. is_finite(x(i))) n_nonfinite_r = n_nonfinite_r + 1
    enddo
  end function n_nonfinite_r

  Integer function n_nonfinite_c(n, x)
    Implicit none
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(in) :: x(n)
    Integer :: i
    n_nonfinite_c = 0
    do i = 1, n
       if (.not. finite_c(x(i))) n_nonfinite_c = n_nonfinite_c + 1
    enddo
  end function n_nonfinite_c

!--------------------------------------------------------------------
!> Groups, launches per B_T and nd histogram per flavour (and the measured launch count of one mmthr call).
  Subroutine report_structure()
    Implicit none
    Integer :: nf, e, nel, ng, nsmall, nbig, ndmax
    Integer(int64) :: l0
    Complex (Kind=Kind(0.d0)), allocatable :: Dt(:,:)
    Character (len=512) :: line

    do nf = 1, N_FL
       nel = hop_start(nf+1) - hop_start(nf)
       ng = hop_gstart(nf+1) - hop_gstart(nf)
       nsmall = 0; nbig = 0; ndmax = 0
       do e = hop_start(nf), hop_start(nf+1)-1
          if (hop_meta(2,e) > 8) then
             nbig = nbig + 1
          else
             nsmall = nsmall + 1
          endif
          if (hop_meta(2,e) > ndmax) ndmax = hop_meta(2,e)
       enddo
       write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'GPU wrap structure flavour ', nf, ': elements ', nel, ' (small nd<=8: ', &
            & nsmall, ', dense nd>8: ', nbig, '), max nd ', ndmax, ', groups ', ng, ' (checkerboard groups = launches per B_T if per-group)'
       call say(line)
    enddo
    if (.not. hop_big(1)) then
       allocate(Dt(Ndim,Ndim)); Dt = z1
       !$omp target enter data map(to: Dt)
       l0 = gpu_n_launch; gpu_cb_fused = .true.;  call gpu_mmthr(Dt, Ndim, Ndim, 1)
       write(line,'(A,I0,A)') 'GPU wrap structure: launches per B_T (mmthr), fused mode: ', gpu_n_launch - l0, ' (1 team per column, groups in sequence)'
       call say(line)
       l0 = gpu_n_launch; gpu_cb_fused = .false.; call gpu_mmthr(Dt, Ndim, Ndim, 1)
       write(line,'(A,I0,A)') 'GPU wrap structure: launches per B_T (mmthr), per-group mode: ', gpu_n_launch - l0, ' (= number of groups)'
       call say(line)
       gpu_cb_fused = .false.
       !$omp target exit data map(delete: Dt)
       deallocate(Dt)
    endif
  end subroutine report_structure

!--------------------------------------------------------------------
!> Per-operator Op_V routines vs Op_mmultR/L, Op_Wrapup/do.  The field of operator n at slice t0 is set on the device
!> to each allowed value (nsigma_d(n,t0), restored afterwards).
  Subroutine test_vop(s, dv, dwa, dwo, dskip, ncases, nops)
    Implicit none
    Integer(int64), intent(inout) :: s
    Real (Kind=Kind(0.d0)), intent(out) :: dv(8), dwa(2), dwo(2), dskip
    Integer, intent(out) :: ncases, nops
    Integer :: nf, n, sp, nt, i, nstride, t0
    Complex (Kind=Kind(0.d0)) :: zsp, old
    Complex (Kind=Kind(0.d0)), allocatable :: A0(:,:), A1(:,:), A2(:,:), Dt(:,:)
    Character (len=1) :: cop
    Type (Operator) :: Otmp
    Real (Kind=Kind(0.d0)) :: d

    t0 = 1
    dv = 0.d0; dwa = 0.d0; dwo = 0.d0; dskip = 0.d0; ncases = 0; nops = 0
    nstride = max(1, size(Op_V,1)/48)
    allocate(A0(Ndim,Ndim), A1(Ndim,Ndim), A2(Ndim,Ndim), Dt(Ndim,Ndim))
    !$omp target enter data map(alloc: Dt)
    do nf = 1, N_FL
       do n = 1, size(Op_V,1)
          if (mod(n-1, nstride) /= 0 .and. n /= size(Op_V,1)) cycle
          if (nf == 1) nops = nops + 1
          nt = V_type(n)
          old = nsigma_d(n,t0)
          do sp = -nt, nt
             if (sp == 0) cycle
             zsp = cmplx(real(sp,kind(0.d0)), 0.d0, kind(0.d0))
             nsigma_d(n,t0) = zsp
             !$omp target update to(nsigma_d(n:n,t0:t0))   ! self-test: set the field value of this case on the device
             ncases = ncases + 1
             do i = 1, 2
                if (i == 1) cop = 'n'
                if (i == 2) cop = 'c'
                ! rows
                call rand_mat(A0, s); A1 = A0
                call Op_mmultR(A1, Op_V(n,nf), zsp, cop, t0)
                call dev_run_R(A0, Dt, n, nf, t0, cop, A2)
                call acc_max(dv(i), frob_rel(A2, A1))
                ! columns, both signs
                call rand_mat(A0, s); A1 = A0
                call Op_mmultL(A1, Op_V(n,nf), zsp, cop, t0, 1)
                call dev_run_L(A0, Dt, n, nf, t0, cop, 1, A2)
                d = frob_rel(A2, A1)
                A1 = A0
                call Op_mmultL(A1, Op_V(n,nf), zsp, cop, t0, -1)
                call dev_run_L(A0, Dt, n, nf, t0, cop, -1, A2)
                call acc_max(d, frob_rel(A2, A1))
                call acc_max(dv(2+i), d)
             enddo
             ! wraps, N_Type 1 and 2
             call rand_mat(A0, s); A1 = A0
             call Op_Wrapup(A1, Op_V(n,nf), zsp, Ndim, 1, t0)
             call dev_run_W(A0, Dt, n, nf, t0, 1, 1, A2); call acc_max(dwa(1), frob_rel(A2, A1))
             A1 = A0
             call Op_Wrapup(A1, Op_V(n,nf), zsp, Ndim, 2, t0)
             call dev_run_W(A0, Dt, n, nf, t0, 1, 2, A2); call acc_max(dwa(2), frob_rel(A2, A1))
             A1 = A0
             call Op_Wrapdo(A1, Op_V(n,nf), zsp, Ndim, 1, t0)
             call dev_run_W(A0, Dt, n, nf, t0, -1, 1, A2); call acc_max(dwo(1), frob_rel(A2, A1))
             A1 = A0
             call Op_Wrapdo(A1, Op_V(n,nf), zsp, Ndim, 2, t0)
             call dev_run_W(A0, Dt, n, nf, t0, -1, 2, A2); call acc_max(dwo(2), frob_rel(A2, A1))
          enddo
          nsigma_d(n,t0) = old
          !$omp target update to(nsigma_d(n:n,t0:t0))
       enddo
    enddo

    ! tiny-g copy of operator (1,1): L/R skip honoured (device flag set), wraps ignore it
    Otmp = Op_V(1,1)
    Otmp%g = cmplx(1.d-20, 0.d0, kind(0.d0))
    nt = V_type(1)
    old = nsigma_d(1,t0)
    V_meta(4,1,1) = 1
    do sp = -nt, nt
       if (sp == 0) cycle
       zsp = cmplx(real(sp,kind(0.d0)), 0.d0, kind(0.d0))
       nsigma_d(1,t0) = zsp
       !$omp target update to(nsigma_d(1:1,t0:t0))
       call rand_mat(A0, s)
       A1 = A0; call Op_mmultR(A1, Otmp, zsp, 'n', t0)
       call dev_run_R(A0, Dt, 1, 1, t0, 'n', A2); call acc_max(dskip, frob_rel(A2, A1))
       if (frob_rel(A1, A0) /= 0.d0) call selftest_fail('Op_mmultR did not return early for tiny g')
       A1 = A0; call Op_mmultL(A1, Otmp, zsp, 'n', t0, 1)
       call dev_run_L(A0, Dt, 1, 1, t0, 'n', 1, A2); call acc_max(dskip, frob_rel(A2, A1))
       A1 = A0; call Op_Wrapup(A1, Otmp, zsp, Ndim, 1, t0)
       call dev_run_W(A0, Dt, 1, 1, t0, 1, 1, A2); call acc_max(dskip, frob_rel(A2, A1))
       if (frob_rel(A1, A0) == 0.d0) call selftest_fail('Op_Wrapup skipped for tiny g')
    enddo
    V_meta(4,1,1) = v_skip_of(Op_V(1,1))
    nsigma_d(1,t0) = old
    !$omp target update to(nsigma_d(1:1,t0:t0))
    !$omp target exit data map(delete: Dt)
    deallocate(A0, A1, A2, Dt)
  end subroutine test_vop

  Integer function v_skip_of(Op)
    Implicit none
    Type (Operator), intent(in) :: Op
    v_skip_of = 0
    if ( abs(Op%g) < Eps_machine ) v_skip_of = 1               ! same test as Op_mmultL/R
  end function v_skip_of

  ! The three helpers work on Ndim x Ndim matrices (explicit shape: the device array is addressed, never copied)
  Subroutine dev_run_R(A0, Dt, n, nf, t, cop, Aout)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A0(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dt(Ndim,Ndim), Aout(Ndim,Ndim)
    Integer, intent(in) :: n, nf, t
    Character (len=1), intent(in) :: cop
    Dt = A0
    !$omp target update to(Dt)
    call gpu_Vop_R(Dt, Ndim, Ndim, n, nf, t, cop)
    !$omp target update from(Dt)
    Aout = Dt
  end subroutine dev_run_R

  Subroutine dev_run_L(A0, Dt, n, nf, t, cop, sign, Aout)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A0(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dt(Ndim,Ndim), Aout(Ndim,Ndim)
    Integer, intent(in) :: n, nf, t, sign
    Character (len=1), intent(in) :: cop
    Dt = A0
    !$omp target update to(Dt)
    call gpu_Vop_L(Dt, Ndim, Ndim, n, nf, t, cop, sign)
    !$omp target update from(Dt)
    Aout = Dt
  end subroutine dev_run_L

  Subroutine dev_run_W(A0, Dt, n, nf, t, dir, ntype, Aout)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: A0(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: Dt(Ndim,Ndim), Aout(Ndim,Ndim)
    Integer, intent(in) :: n, nf, t, dir, ntype
    Dt = A0
    !$omp target update to(Dt)
    call gpu_Vop_wrap(Dt, Ndim, n, nf, t, dir, ntype)
    !$omp target update from(Dt)
    Aout = Dt
  end subroutine dev_run_W

!--------------------------------------------------------------------
!> Fused slice kernels vs the CPU loops over all operators, for three slices, all flavours, square and non-square shapes.
  Subroutine test_slice(s, dsl, dwu, dwd)
    Implicit none
    Integer(int64), intent(inout) :: s
    Real (Kind=Kind(0.d0)), intent(out) :: dsl(4), dwu, dwd
    Integer :: nf, it, t, n, ish, iv, sg, nrow, ncl, ncolp
    Complex (Kind=Kind(0.d0)), allocatable :: A0(:,:), A1(:,:), A2(:,:), Dt(:,:)
    Character (len=1) :: cop
    Logical :: desc
    Integer :: ts(3)
    Type (Operator) :: Otmp

    ncolp = test_ncolp()
    dsl = 0.d0; dwu = 0.d0; dwd = 0.d0
    ts = (/ 1, (Ltrot+1)/2, Ltrot /)
    do nf = 1, N_FL
       do it = 1, 3
          t = ts(it)
          do ish = 1, 2
             if (ish == 1) then
                nrow = Ndim; ncl = Ndim
             else
                nrow = ncolp; ncl = ncolp
             endif
             if (nf == 1 .and. it == 1 .and. ish == 1) call say('GPU wrap self-test progress: fused slice rows')
             ! rows (matrix Ndim x ncl)
             allocate(A0(Ndim,ncl), A1(Ndim,ncl), A2(Ndim,ncl), Dt(Ndim,ncl))
             !$omp target enter data map(alloc: Dt)
             do iv = 1, 4
                if (iv == 1 .or. iv == 3) then; cop = 'n'; else; cop = 'c'; endif
                desc = (iv >= 3)
                call rand_mat(A0, s); A1 = A0
                if (.not. desc) then
                   do n = 1, size(Op_V,1);  call Op_mmultR(A1, Op_V(n,nf), nsigma%f(n,t), cop, t); enddo
                else
                   do n = size(Op_V,1), 1, -1; call Op_mmultR(A1, Op_V(n,nf), nsigma%f(n,t), cop, t); enddo
                endif
                Dt = A0
                !$omp target update to(Dt)
                call gpu_Vslice_R(Dt, Ndim, ncl, nf, t, cop, desc)
                !$omp target update from(Dt)
                A2 = Dt
                call acc_max(dsl(1), frob_rel(A2, A1))
             enddo
             !$omp target exit data map(delete: Dt)
             deallocate(A0, A1, A2, Dt)
             if (nf == 1 .and. it == 1 .and. ish == 1) call say('GPU wrap self-test progress: fused slice columns')
             ! columns (matrix nrow x Ndim)
             allocate(A0(nrow,Ndim), A1(nrow,Ndim), A2(nrow,Ndim), Dt(nrow,Ndim))
             !$omp target enter data map(alloc: Dt)
             do iv = 1, 8
                if (mod(iv,2) == 1) then; cop = 'n'; else; cop = 'c'; endif
                if (iv <= 4) then; sg = 1; else; sg = -1; endif
                desc = (mod((iv-1)/2, 2) == 1)
                call rand_mat(A0, s); A1 = A0
                if (.not. desc) then
                   do n = 1, size(Op_V,1);  call Op_mmultL(A1, Op_V(n,nf), nsigma%f(n,t), cop, t, sg); enddo
                else
                   do n = size(Op_V,1), 1, -1; call Op_mmultL(A1, Op_V(n,nf), nsigma%f(n,t), cop, t, sg); enddo
                endif
                Dt = A0
                !$omp target update to(Dt)
                call gpu_Vslice_L(Dt, nrow, nrow, nf, t, cop, sg, desc)
                !$omp target update from(Dt)
                A2 = Dt
                call acc_max(dsl(2), frob_rel(A2, A1))
             enddo
             !$omp target exit data map(delete: Dt)
             deallocate(A0, A1, A2, Dt)
          enddo
          if (nf == 1 .and. it == 1) call say('GPU wrap self-test progress: fused slice wraps')
          ! slice wraps (G is Ndim x Ndim)
          allocate(A0(Ndim,Ndim), A1(Ndim,Ndim), A2(Ndim,Ndim), Dt(Ndim,Ndim))
          !$omp target enter data map(alloc: Dt)
          call rand_mat(A0, s); A1 = A0
          do n = 1, size(Op_V,1);  call Op_Wrapup(A1, Op_V(n,nf), nsigma%f(n,t), Ndim, 1, t); enddo
          Dt = A0
          !$omp target update to(Dt)
          call gpu_Vwrap_all(Dt, Ndim, nf, t, 1)
          !$omp target update from(Dt)
          A2 = Dt
          call acc_max(dwu, frob_rel(A2, A1))
          A1 = A0
          do n = size(Op_V,1), 1, -1;  call Op_Wrapdo(A1, Op_V(n,nf), nsigma%f(n,t), Ndim, 1, t); enddo
          Dt = A0
          !$omp target update to(Dt)
          call gpu_Vwrap_all(Dt, Ndim, nf, t, -1)
          !$omp target update from(Dt)
          A2 = Dt
          call acc_max(dwd, frob_rel(A2, A1))
          !$omp target exit data map(delete: Dt)
          deallocate(A0, A1, A2, Dt)
       enddo
    enddo

    ! tiny-g in the fused kernels: operator 1 of flavour 1 skipped by R/L (device flag), still applied by the wrap
    Otmp = Op_V(1,1)
    Otmp%g = cmplx(1.d-20, 0.d0, kind(0.d0))
    allocate(A0(Ndim,Ndim), A1(Ndim,Ndim), A2(Ndim,Ndim), Dt(Ndim,Ndim))
    !$omp target enter data map(alloc: Dt)
    V_meta(4,1,1) = 1
    !$omp target update to(V_meta)
    t = ts(1)
    call rand_mat(A0, s); A1 = A0
    call Op_mmultR(A1, Otmp, nsigma%f(1,t), 'n', t)
    do n = 2, size(Op_V,1);  call Op_mmultR(A1, Op_V(n,1), nsigma%f(n,t), 'n', t); enddo
    Dt = A0
    !$omp target update to(Dt)
    call gpu_Vslice_R(Dt, Ndim, Ndim, 1, t, 'n', .false.)
    !$omp target update from(Dt)
    A2 = Dt
    call acc_max(dsl(3), frob_rel(A2, A1))
    A1 = A0
    call Op_mmultL(A1, Otmp, nsigma%f(1,t), 'n', t, 1)
    do n = 2, size(Op_V,1);  call Op_mmultL(A1, Op_V(n,1), nsigma%f(n,t), 'n', t, 1); enddo
    Dt = A0
    !$omp target update to(Dt)
    call gpu_Vslice_L(Dt, Ndim, Ndim, 1, t, 'n', 1, .false.)
    !$omp target update from(Dt)
    A2 = Dt
    call acc_max(dsl(3), frob_rel(A2, A1))
    V_meta(4,1,1) = v_skip_of(Op_V(1,1))
    !$omp target update to(V_meta)
    !$omp target exit data map(delete: Dt)
    deallocate(A0, A1, A2, Dt)
  end subroutine test_slice

!--------------------------------------------------------------------
!> Synthetic operator set with overlapping site lists, mixed field types and one skipped operator, written over the
!> gpu_state tables (host and device) for the duration of the test and restored afterwards.  The fused kernels
!> (rows, columns, whole-slice wraps) must reproduce the sequential per-operator products bit for bit, as must the
!> chained per-operator device routines.  Flavour 1 only; flavour 2 gets empty lists.
  Subroutine test_synthetic_V(s, dmax, nsyn)
    Implicit none
    Integer(int64), intent(inout) :: s
    Real (Kind=Kind(0.d0)), intent(out) :: dmax
    Integer, intent(out) :: nsyn
    Integer, parameter :: ns = 6
    Integer, parameter :: n1s(ns) = (/ 3, 2, 4, 1, 3, 2 /), nts(ns) = (/ 1, 2, 1, 2, 1, 2 /)
    Integer :: pool, n, i, j, poff, eoff, p, r, nf2, t, iv, sg, dirv, k, sp, q
    Integer, allocatable :: Pcopy(:), cnt(:)
    Integer, allocatable :: Vm_s(:,:,:), Vt_s(:), VP_s(:), vrp_s(:,:), vre_s(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: Ve_s(:), ns_s(:,:)
    Complex (Kind=Kind(0.d0)), allocatable :: A0(:,:), A1(:,:), A2(:,:), A3(:,:), Dt(:,:)
    Character (len=1) :: cop
    Real (Kind=Kind(0.d0)) :: x, y

    dmax = 0.d0; nsyn = 0
    if (size(Op_V,1) < ns .or. Ndim < 4) return
    poff = sum(n1s)
    eoff = 0
    do n = 1, ns
       eoff = eoff + n1s(n)*(2*nts(n)+1)
    enddo
    if (size(V_P) < poff .or. size(V_exp) < eoff .or. size(vr_ent,2) < poff .or. size(V_type) < ns) return
    nsyn = ns
    t = 1
    ! save
    allocate(Vm_s(size(V_meta,1),size(V_meta,2),size(V_meta,3)), Vt_s(size(V_type)), VP_s(size(V_P)), Ve_s(size(V_exp)))
    allocate(vrp_s(size(vr_ptr,1),size(vr_ptr,2)), vre_s(size(vr_ent,1),size(vr_ent,2)), ns_s(size(nsigma_d,1),size(nsigma_d,2)))
    Vm_s = V_meta; Vt_s = V_type; VP_s = V_P; Ve_s = V_exp; vrp_s = vr_ptr; vre_s = vr_ent; ns_s = nsigma_d
    ! build
    pool = min(Ndim, 8)
    poff = 0; eoff = 0
    do n = 1, ns
       V_type(n) = nts(n)
       V_meta(1,n,1) = n1s(n); V_meta(2,n,1) = poff; V_meta(3,n,1) = eoff
       V_meta(4,n,1) = 0
       if (n == 3) V_meta(4,n,1) = 1                      ! skipped by Op_mmultL/R (tiny coupling), not by the wraps
       ! n1s distinct sites out of the pool: partial Fisher-Yates
       allocate(Pcopy(pool))
       do i = 1, pool; Pcopy(i) = i; enddo
       do i = 1, n1s(n)
          j = i + int((rnd(s) + 0.5d0)*real(pool - i + 1, kind(0.d0)))
          if (j > pool) j = pool
          k = Pcopy(i); Pcopy(i) = Pcopy(j); Pcopy(j) = k
          V_P(poff+i) = Pcopy(i)
       enddo
       deallocate(Pcopy)
       do i = 1, n1s(n)*(2*nts(n)+1)
          x = 0.6d0 + abs(rnd(s)) + 0.3d0*rnd(s)
          y = 0.4d0*rnd(s)
          V_exp(eoff+i) = cmplx(x, y, kind(0.d0))
       enddo
       poff = poff + n1s(n); eoff = eoff + n1s(n)*(2*nts(n)+1)
    enddo
    ! per-site lists, flavour 1; flavour 2 empty (other flavours as well)
    allocate(cnt(Ndim))
    cnt = 0
    do n = 1, ns
       do i = 1, n1s(n)
          r = V_P(V_meta(2,n,1)+i); cnt(r) = cnt(r) + 1
       enddo
    enddo
    p = 1
    do r = 1, Ndim
       vr_ptr(r,1) = p; p = p + cnt(r)
    enddo
    vr_ptr(Ndim+1,1) = p
    do nf2 = 2, size(vr_ptr,2)
       vr_ptr(:,nf2) = p
    enddo
    do r = 1, Ndim
       cnt(r) = vr_ptr(r,1)
    enddo
    do n = 1, ns
       do i = 1, n1s(n)
          r = V_P(V_meta(2,n,1)+i)
          vr_ent(1,cnt(r)) = n; vr_ent(2,cnt(r)) = i; cnt(r) = cnt(r) + 1
       enddo
    enddo
    deallocate(cnt)
    do n = 1, ns
       k = 1 + int((rnd(s) + 0.5d0)*real(2*nts(n), kind(0.d0)))      ! 1..2nt
       sp = k - nts(n) - 1; if (sp >= 0) sp = sp + 1                  ! -nt..-1, 1..nt
       nsigma_d(n,t) = cmplx(real(sp,kind(0.d0)), 0.d0, kind(0.d0))
    enddo
    ! device copies (self-test only: synthetic data, restored below)
    !$omp target update to(V_meta, V_type, V_P, V_exp, vr_ptr, vr_ent, nsigma_d)

    allocate(A0(Ndim,Ndim), A1(Ndim,Ndim), A2(Ndim,Ndim), A3(Ndim,Ndim), Dt(Ndim,Ndim))
    !$omp target enter data map(alloc: Dt)
    ! rows / columns, all variants
    do iv = 0, 7
       if (mod(iv,2) == 0) then; cop = 'n'; else; cop = 'c'; endif
       dirv = merge(1, 0, mod(iv/2,2) == 1)                           ! 1 = descending
       sg = merge(-1, 1, iv >= 4)
       ! rows (sign irrelevant), use iv < 4 for R and all for L
       if (iv < 4) then
          call rand_mat(A0, s); A1 = A0
          call ref_slice('R', A1, cop, 1, dirv == 1, t)
          Dt = A0
          !$omp target update to(Dt)
          call gpu_Vslice_R(Dt, Ndim, Ndim, 1, t, cop, dirv == 1)
          !$omp target update from(Dt)
          A2 = Dt
          call acc_max(dmax, frob_rel(A2, A1))
          ! chained per-operator routines
          A3 = A0
          Dt = A0
          !$omp target update to(Dt)
          do q = 1, ns
             n = merge(ns+1-q, q, dirv == 1)
             call gpu_Vop_R(Dt, Ndim, Ndim, n, 1, t, cop)
          enddo
          !$omp target update from(Dt)
          A3 = Dt
          call acc_max(dmax, frob_rel(A3, A1))
       endif
       call rand_mat(A0, s); A1 = A0
       call ref_slice('L', A1, cop, sg, dirv == 1, t)
       Dt = A0
       !$omp target update to(Dt)
       call gpu_Vslice_L(Dt, Ndim, Ndim, 1, t, cop, sg, dirv == 1)
       !$omp target update from(Dt)
       A2 = Dt
       call acc_max(dmax, frob_rel(A2, A1))
       Dt = A0
       !$omp target update to(Dt)
       do q = 1, ns
          n = merge(ns+1-q, q, dirv == 1)
          call gpu_Vop_L(Dt, Ndim, Ndim, n, 1, t, cop, sg)
       enddo
       !$omp target update from(Dt)
       A3 = Dt
       call acc_max(dmax, frob_rel(A3, A1))
    enddo
    ! slice wraps, both directions; chained per-operator wraps
    do iv = 1, 2
       dirv = merge(1, -1, iv == 1)
       call rand_mat(A0, s); A1 = A0
       call ref_wrap(A1, dirv, t)
       Dt = A0
       !$omp target update to(Dt)
       call gpu_Vwrap_all(Dt, Ndim, 1, t, dirv)
       !$omp target update from(Dt)
       A2 = Dt
       call acc_max(dmax, frob_rel(A2, A1))
       Dt = A0
       !$omp target update to(Dt)
       do q = 1, ns
          n = merge(q, ns+1-q, dirv == 1)
          call gpu_Vop_wrap(Dt, Ndim, n, 1, t, dirv, 1)
       enddo
       !$omp target update from(Dt)
       A3 = Dt
       call acc_max(dmax, frob_rel(A3, A1))
    enddo
    !$omp target exit data map(delete: Dt)
    deallocate(A0, A1, A2, A3, Dt)
    ! restore host and device
    V_meta = Vm_s; V_type = Vt_s; V_P = VP_s; V_exp = Ve_s; vr_ptr = vrp_s; vr_ent = vre_s; nsigma_d = ns_s
    !$omp target update to(V_meta, V_type, V_P, V_exp, vr_ptr, vr_ent, nsigma_d)
    deallocate(Vm_s, Vt_s, VP_s, Ve_s, vrp_s, vre_s, ns_s)
  end subroutine test_synthetic_V

  ! Sequential reference from the flat tables (flavour 1, operators 1..6 in order or reversed), host, as ZSCAL loops
  Subroutine ref_slice(side, M, cop, sg, desc, t)
    Implicit none
    Character (len=1), intent(in) :: side, cop
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer, intent(in) :: sg, t
    Logical, intent(in) :: desc
    Integer :: q, n, i, N1, po, eo, nt, sp
    Complex (Kind=Kind(0.d0)) :: z
    do q = 1, 6
       n = q
       if (desc) n = 7 - q
       if (V_meta(4,n,1) == 1) cycle
       N1 = V_meta(1,n,1); po = V_meta(2,n,1); eo = V_meta(3,n,1); nt = V_type(n)
       sp = sg*nint(real(nsigma_d(n,t)))
       do i = 1, N1
          z = V_exp(eo + i + (sp+nt)*N1)
          if (cop == 'c') z = conjg(z)
          if (side == 'R') then
             M(V_P(po+i),:) = z*M(V_P(po+i),:)
          else
             M(:,V_P(po+i)) = z*M(:,V_P(po+i))
          endif
       enddo
    enddo
  end subroutine ref_slice

  Subroutine ref_wrap(M, dir, t)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(:,:)
    Integer, intent(in) :: dir, t
    Integer :: q, n, i, N1, po, eo, nt, sp
    do q = 1, 6
       n = q
       if (dir < 0) n = 7 - q
       N1 = V_meta(1,n,1); po = V_meta(2,n,1); eo = V_meta(3,n,1); nt = V_type(n)
       sp = nint(real(nsigma_d(n,t)))
       do i = 1, N1
          M(V_P(po+i),:) = V_exp(eo + i + (dir*sp+nt)*N1) * M(V_P(po+i),:)
       enddo
       do i = 1, N1
          M(:,V_P(po+i)) = V_exp(eo + i + (-dir*sp+nt)*N1) * M(:,V_P(po+i))
       enddo
    enddo
  end subroutine ref_wrap

!--------------------------------------------------------------------
!> wrapur / wrapul / WRAPGR sequences over Nwrap slices (and a longer one over all Ltrot slices), device vs CPU, error
!> after every slice (growth), round trip up/down.
  Subroutine test_sequences(s, tol, nw, ncolp)
    Implicit none
    Integer(int64), intent(inout) :: s
    Real (Kind=Kind(0.d0)), intent(in) :: tol
    Integer, intent(in) :: nw, ncolp
    Integer :: nf, nsl, t, ic
    Real (Kind=Kind(0.d0)) :: dur(0:Ltrot), dul(0:Ltrot), dgu(0:Ltrot), dgd(0:Ltrot), dmax(5), dlong(4), dround, dret
    Real (Kind=Kind(0.d0)) :: drt(0:2*Ltrot), dsens(0:2*Ltrot)
    Complex (Kind=Kind(0.d0)), allocatable :: U0(:,:), U1(:,:), U2(:,:), Dt(:,:), G0(:,:), G1(:,:), Gp(:,:)
    Character (len=1024) :: line

    nsl = min(max(nw,1), Ltrot)
    dmax = 0.d0; dlong = 0.d0; dround = 0.d0
    do nf = 1, N_FL
       if (hop_big(nf) .and. .not. la_available) cycle
       do ic = 1, 2
          if (ic == 1) then
             gpu_cb_fused = fused_ok;  gpu_dense_left_mode = 1
          else
             gpu_cb_fused = .false.; gpu_dense_left_mode = 0
          endif
          dur = 0.d0; dul = 0.d0; dgu = 0.d0; dgd = 0.d0; drt = 0.d0; dsens = 0.d0
          ! ---- wrapur: U (Ndim x ncolp), slices 1..Ltrot one at a time, compare after each
          allocate(U0(Ndim,ncolp), U1(Ndim,ncolp), U2(Ndim,ncolp), Dt(Ndim,ncolp))
          call rand_mat(U0, s)
          U1 = U0
          !$omp target enter data map(alloc: Dt)
          Dt = U0
          !$omp target update to(Dt)
          do t = 1, Ltrot
             call cpu_wrapur(U1, nf, t-1, t)
             call gpu_wrapur_slices(Dt, Ndim, ncolp, nf, t-1, t)
             !$omp target update from(Dt)    ! self-test: compare after every slice
             U2 = Dt
             dur(t) = frob_rel(U2, U1)
             if (.not. is_finite(dur(t))) dur(t) = huge(1.d0)
          enddo
          !$omp target exit data map(delete: Dt)
          deallocate(U0, U1, U2, Dt)
          ! ---- wrapul: U (Ndim x ncolp), slices Ltrot..1 descending
          allocate(U0(Ndim,ncolp), U1(Ndim,ncolp), U2(Ndim,ncolp), Dt(Ndim,ncolp))
          call rand_mat(U0, s)
          U1 = U0
          !$omp target enter data map(alloc: Dt)
          Dt = U0
          !$omp target update to(Dt)
          do t = Ltrot, 1, -1
             call cpu_wrapul(U1, nf, t-1, t)
             call gpu_wrapul_slices(Dt, Ndim, ncolp, nf, t-1, t)
             !$omp target update from(Dt)
             U2 = Dt
             dul(Ltrot-t+1) = frob_rel(U2, U1)
          enddo
          !$omp target exit data map(delete: Dt)
          deallocate(U0, U1, U2, Dt)
          ! ---- WRAPGRUP chain (slices 1..Ltrot) and WRAPGRDO chain (Ltrot..1), each from a fresh G (Ndim x Ndim)
          allocate(G0(Ndim,Ndim), G1(Ndim,Ndim), U2(Ndim,Ndim), Dt(Ndim,Ndim))
          !$omp target enter data map(alloc: Dt)
          call rand_mat(G0, s)
          G1 = G0; Dt = G0
          !$omp target update to(Dt)
          do t = 1, Ltrot
             call cpu_wrapgr_up(G1, nf, t-1)
             call gpu_wrapgr_up(Dt, Ndim, nf, t-1)
             !$omp target update from(Dt)
             U2 = Dt
             dgu(t) = frob_rel(U2, G1)
          enddo
          call rand_mat(G0, s)
          G1 = G0; Dt = G0
          !$omp target update to(Dt)
          do t = Ltrot, 1, -1
             call cpu_wrapgr_do(G1, nf, t)
             call gpu_wrapgr_do(Dt, Ndim, nf, t)
             !$omp target update from(Dt)
             U2 = Dt
             dgd(Ltrot-t+1) = frob_rel(U2, G1)
          enddo
          ! ---- round trip: Nwrap up wraps then the Nwrap down wraps back (CPU, device, and a CPU run whose start matrix
          ! is perturbed by 1e-16 relative: the intrinsic sensitivity of the chain to rounding)
          call rand_mat(G0, s)
          G1 = G0; Dt = G0
          !$omp target update to(Dt)
          allocate(Gp(Ndim,Ndim))
          do t = 1, Ndim
             Gp(:,t) = G0(:,t) * (1.d0 + 1.d-16*(1.d0 + rnd(s)))
          enddo
          drt = 0.d0; dsens = 0.d0
          do t = 1, nsl
             call cpu_wrapgr_up(G1, nf, t-1); call cpu_wrapgr_up(Gp, nf, t-1)
             call gpu_wrapgr_up(Dt, Ndim, nf, t-1)
             !$omp target update from(Dt)
             U2 = Dt
             drt(t) = frob_rel(U2, G1); dsens(t) = frob_rel(Gp, G1)
          enddo
          do t = nsl, 1, -1
             call cpu_wrapgr_do(G1, nf, t); call cpu_wrapgr_do(Gp, nf, t)
             call gpu_wrapgr_do(Dt, Ndim, nf, t)
             !$omp target update from(Dt)
             U2 = Dt
             drt(2*nsl-t+1) = frob_rel(U2, G1); dsens(2*nsl-t+1) = frob_rel(Gp, G1)
          enddo
          dret = frob_rel(U2, G0)
          !$omp target exit data map(delete: Dt)
          deallocate(G0, G1, U2, Dt, Gp)
          write(line,'(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU wrap self-test up/down round trip flavour ', nf, &
               & ' combination ', ic, ' (', nsl, ' up wraps, then ', drt(nsl), ' <- device vs CPU; CPU vs CPU with 1e-16 perturbed start ', dsens(nsl), &
               & ' | after the down wraps: ', drt(2*nsl), ' vs ', dsens(2*nsl)
          call say(line)
          do t = 1, 2*nsl
             call acc_max(dmax(5), drt(t))
          enddo
          call acc_max(dround, dret)
          ! report: after Nwrap slices and growth
          call growth_line('wrapur ', nf, ic, dur, nsl, line); call say(line)
          call growth_line('wrapul ', nf, ic, dul, nsl, line); call say(line)
          call growth_line('wrapgr_up ', nf, ic, dgu, nsl, line); call say(line)
          call growth_line('wrapgr_do ', nf, ic, dgd, nsl, line); call say(line)
          do t = 1, nsl
             call acc_max(dmax(1), dur(t)); call acc_max(dmax(2), dul(t)); call acc_max(dmax(3), dgu(t))
             call acc_max(dmax(4), dgd(t))
          enddo
          ! the longer chains (all Ltrot slices) are reported, not part of the criterion
          do t = nsl+1, Ltrot
             call acc_max(dlong(1), dur(t)); call acc_max(dlong(2), dul(t)); call acc_max(dlong(3), dgu(t))
             call acc_max(dlong(4), dgd(t))
          enddo
       enddo
    enddo
    gpu_cb_fused = .false.; gpu_dense_left_mode = 0
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU wrap self-test sequences after Nwrap = ', nw, &
         & ' slices, max over slices, flavours, code-path combinations (device vs CPU): wrapur ', dmax(1), ', wrapul ', dmax(2), &
         & ', wrapgr_up ', dmax(3), ', wrapgr_do ', dmax(4), ' (tol ', tol
    call say(line)
    write(line,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'GPU wrap self-test sequences beyond Nwrap up to ', Ltrot, &
         & ' slices (reported, not a criterion; the round trip is as sensitive to rounding on the CPU, see the lines above): wrapur ', dlong(1), ', wrapul ', dlong(2), ', wrapgr_up ', dlong(3), &
         & ', wrapgr_do ', dlong(4), '; up/down round trip device vs CPU ', dmax(5), ', device vs start ', dround
    call say(line)
    if (.not. (dmax(1) <= tol .and. dmax(2) <= tol .and. dmax(3) <= tol .and. dmax(4) <= tol)) &
         & call selftest_fail('wrap sequences')
  end subroutine test_sequences

  Subroutine growth_line(name, nf, ic, d, nsl, line)
    Implicit none
    Character (len=*), intent(in) :: name
    Integer, intent(in) :: nf, ic, nsl
    Real (Kind=Kind(0.d0)), intent(in) :: d(0:)
    Character (len=*), intent(out) :: line
    Integer :: t
    Character (len=64) :: s1
    write(line,'(A,A,A,I0,A,I0,A,ES10.3,A,ES10.3,A)') 'GPU wrap self-test ', name, ' flavour ', nf, ' combination ', ic, &
         & ': rel diff after 1 slice ', d(1), ', after Nwrap ', d(nsl), ' (per slice, every second:'
    do t = 2, nsl, 2
       write(s1,'(ES9.2)') d(t)
       line = trim(line)//' '//trim(s1)
    enddo
    line = trim(line)//')'
  end subroutine growth_line

  Subroutine cpu_wrapur(U, nf, ntau, ntau1)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(:,:)
    Integer, intent(in) :: nf, ntau, ntau1
    Integer :: t, n
    do t = ntau + 1, ntau1
       call Hop_mod_mmthr(U, nf, t)
       do n = 1, size(Op_V,1)
          call Op_mmultR(U, Op_V(n,nf), nsigma%f(n,t), 'n', t)
       enddo
    enddo
  end subroutine cpu_wrapur

  Subroutine cpu_wrapul(U, nf, ntau, ntau1)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: U(:,:)
    Integer, intent(in) :: nf, ntau, ntau1
    Integer :: t, n
    do t = ntau1, ntau + 1, -1
       do n = size(Op_V,1), 1, -1
          call Op_mmultR(U, Op_V(n,nf), nsigma%f(n,t), 'c', t)
       enddo
       call Hop_mod_mmthlc(U, nf, t)
    enddo
  end subroutine cpu_wrapul

  Subroutine cpu_wrapgr_up(G, nf, ntau)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(:,:)
    Integer, intent(in) :: nf, ntau
    Integer :: n
    call Hop_mod_mmthr(G, nf, ntau+1)
    call Hop_mod_mmthl_m1(G, nf, ntau+1)
    do n = 1, size(Op_V,1)
       call Op_Wrapup(G, Op_V(n,nf), nsigma%f(n,ntau+1), Ndim, 1, ntau+1)
    enddo
  end subroutine cpu_wrapgr_up

  Subroutine cpu_wrapgr_do(G, nf, ntau)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(:,:)
    Integer, intent(in) :: nf, ntau
    Integer :: n
    do n = size(Op_V,1), 1, -1
       call Op_Wrapdo(G, Op_V(n,nf), nsigma%f(n,ntau), Ndim, 2, ntau)
       call Op_Wrapdo(G, Op_V(n,nf), nsigma%f(n,ntau), Ndim, 1, ntau)
    enddo
    call Hop_mod_mmthl(G, nf, ntau)
    call Hop_mod_mmthr_m1(G, nf, ntau)
  end subroutine cpu_wrapgr_do

!====================================================================
!  Timing
!====================================================================
  Subroutine gpu_wrap_timing()
    Implicit none
    Complex (Kind=Kind(0.d0)), allocatable :: H(:,:), Hs(:,:), D(:,:)
    Integer(int64) :: s
    Real (Kind=Kind(0.d0)) :: tg, tc
    Integer :: nf, ic, mode
    Character (len=512) :: line
    Character (len=64) :: lab
    Logical :: dense

    s = 777_int64
    nf = 1
    dense = hop_big(1)
    allocate(H(Ndim,Ndim), Hs(Ndim,Ndim), D(Ndim,Ndim))
    call rand_mat(Hs, s)
    !$omp target enter data map(alloc: D)
    write(line,'(A,I0,A,L1,A,I0)') 'GPU wrap timing: N = ', Ndim, ', dense hopping = ', dense, ', N_op = ', size(Op_V,1)
    call say(line)
    if (dense .and. .not. la_available) then
       call say('GPU wrap timing: dense hopping needs alf_la, skipped for this build')
       !$omp target exit data map(delete: D)
       deallocate(H, Hs, D)
       return
    endif

    call warm_up(Hs, D)
    ! B G B^-1 pair: mmthr + mmthl_m1 (cpu: same Hop_mod routines), restored by the inverse pair.
    ! ic = 1: default (cb: one launch per group; dense: left variant 2); ic = 2: cb one launch per B_T (fused), dense left variant 1
    do ic = 1, 2
       gpu_cb_fused = (ic == 2)
       if (dense) gpu_dense_left_mode = 3 - ic       ! 2 transposes + DGEMM, then 1 complex copy + ZGEMM
       call time_pair(Hs, D, tg, tc)
       if (dense) then
          write(lab,'(A,I0,A)') 'B G B^-1 (dense, left variant ', 3 - ic, ')'
       else
          if (ic == 1) lab = 'B G B^-1 (cb, one launch per group)'
          if (ic == 2) lab = 'B G B^-1 (cb, fused single launch)'
       endif
       call report_time(trim(lab), tg, tc)
    enddo
    gpu_cb_fused = .false.
    if (dense) gpu_dense_left_mode = 0
    ! single hopping applications (default modes)
    call time_hop(Hs, D, 1, 2, 'mmthr', 'mmthr_m1')
    call time_hop(Hs, D, 3, 5, 'mmthl', 'mmthl_m1')
    call time_hop(Hs, D, 4, 2, 'mmthlc', 'mmthr_m1')
    if (dense) then
       gpu_dense_left_mode = 1; call time_hop(Hs, D, 1, 2, 'mmthr[v1]', 'mmthr_m1[v1]')
       gpu_dense_left_mode = 2; call time_hop(Hs, D, 1, 2, 'mmthr[v2]', 'mmthr_m1[v2]')
       gpu_dense_left_mode = 0
    else
       gpu_cb_fused = .true.; call time_hop(Hs, D, 1, 2, 'mmthr[fused]', 'mmthr_m1[fused]')
       call time_hop(Hs, D, 3, 5, 'mmthl[fused]', 'mmthl_m1[fused]')
       gpu_cb_fused = .false.
    endif
    call time_vslice(Hs, D)
    call time_wrapgr(Hs, D)
    !$omp target exit data map(delete: D)
    deallocate(H, Hs, D)
  end subroutine gpu_wrap_timing

  ! busy device for ~0.5 s so that the clocks are up before anything is timed
  Subroutine warm_up(Hs, D)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: Hs(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(Ndim,Ndim)
    Real (Kind=Kind(0.d0)) :: t0
    D = Hs
    !$omp target update to(D)
    t0 = omp_get_wtime()
    do while (omp_get_wtime() - t0 < 0.5d0)
       call gpu_Vslice_R(D, Ndim, Ndim, 1, 1, 'n', .false.)
       call gpu_Vslice_R(D, Ndim, Ndim, 1, 1, 'n', .true.)
       call gpu_mmthr(D, Ndim, Ndim, 1); call gpu_mmthl_m1(D, Ndim, Ndim, 1)
       call gpu_mmthr_m1(D, Ndim, Ndim, 1); call gpu_mmthl(D, Ndim, Ndim, 1)
       !$omp target update from(D)       ! synchronise (and keeps the values bounded by the caller's restore)
       D = Hs
       !$omp target update to(D)
    enddo
  end subroutine warm_up

  Subroutine report_time(name, tg, tc)
    Implicit none
    Character (len=*), intent(in) :: name
    Real (Kind=Kind(0.d0)), intent(in) :: tg, tc
    Character (len=512) :: line
    write(line,'(A,A,A,I0,A,F12.2,A,F12.2,A,F9.2)') 'GPU wrap timing ', name, ' N=', Ndim, ': device us ', tg*1.d6, ', one CPU core us ', &
         & tc*1.d6, ', speedup ', tc/max(tg,1.d-12)
    call say(line)
  end subroutine report_time

  Real (Kind=Kind(0.d0)) function median5(t)
    Implicit none
    Real (Kind=Kind(0.d0)), intent(inout) :: t(5)
    Integer :: i, j
    Real (Kind=Kind(0.d0)) :: x
    do i = 2, 5
       x = t(i); j = i - 1
       do while (j >= 1)
          if (t(j) <= x) exit
          t(j+1) = t(j); j = j - 1
       enddo
       t(j+1) = x
    enddo
    median5 = t(3)
  end function median5

!> One call of op on the CPU / device, alternating op and its inverse so that the matrix stays bounded: us per call.
  Subroutine time_hop(Hs, D, opa, opb, na, nb)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: Hs(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(Ndim,Ndim)
    Integer, intent(in) :: opa, opb
    Character (len=*), intent(in) :: na, nb
    Complex (Kind=Kind(0.d0)), allocatable :: H(:,:)
    Real (Kind=Kind(0.d0)) :: t0, t1, tt(5), tg, tc, tone
    Integer :: itr, irep, nrep

    allocate(H(Ndim,Ndim))
    ! CPU: calibrate repetitions from one pair
    H = Hs
    t0 = omp_get_wtime(); call cpu_hop(opa, 1, H); call cpu_hop(opb, 1, H); t1 = omp_get_wtime()
    tone = max(t1 - t0, 1.d-7)
    nrep = max(1, min(50, int(0.05d0/tone)))
    do itr = 1, 5
       H = Hs
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call cpu_hop(opa, 1, H); call cpu_hop(opb, 1, H)
       enddo
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(2*nrep, kind(0.d0))
    enddo
    tc = median5(tt)
    ! device
    D = Hs
    !$omp target update to(D)
    call dev_hop(opa, 1, D, size(D,1), size(D,2)); call dev_hop(opb, 1, D, size(D,1), size(D,2)); call la_sync_if()
    nrep = 100
    do itr = 1, 5
       D = Hs
       !$omp target update to(D)
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call dev_hop(opa, 1, D, size(D,1), size(D,2)); call dev_hop(opb, 1, D, size(D,1), size(D,2))
       enddo
       call la_sync_if()
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(2*nrep, kind(0.d0))
    enddo
    tg = median5(tt)
    call report_time(na//'/'//nb, tg, tc)
    deallocate(H)
  end subroutine time_hop

  Subroutine la_sync_if()
    Implicit none
    if (la_up) call la_sync()
  end subroutine la_sync_if

  ! B G B^-1 then B^-1 G B (cpu: Hop_mod routines; device: gpu_ routines): us per pair member
  Subroutine time_pair(Hs, D, tg, tc)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: Hs(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(Ndim,Ndim)
    Real (Kind=Kind(0.d0)), intent(out) :: tg, tc
    Complex (Kind=Kind(0.d0)), allocatable :: H(:,:)
    Real (Kind=Kind(0.d0)) :: t0, t1, tt(5), tone
    Integer :: itr, irep, nrep

    allocate(H(Ndim,Ndim))
    H = Hs
    t0 = omp_get_wtime()
    call Hop_mod_mmthr(H,1,1); call Hop_mod_mmthl_m1(H,1,1); call Hop_mod_mmthr_m1(H,1,1); call Hop_mod_mmthl(H,1,1)
    t1 = omp_get_wtime()
    tone = max(t1 - t0, 1.d-7)
    nrep = max(1, min(50, int(0.05d0/tone)))
    do itr = 1, 5
       H = Hs
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call Hop_mod_mmthr(H,1,1); call Hop_mod_mmthl_m1(H,1,1)
          call Hop_mod_mmthr_m1(H,1,1); call Hop_mod_mmthl(H,1,1)
       enddo
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(2*nrep, kind(0.d0))
    enddo
    tc = median5(tt)
    D = Hs
    !$omp target update to(D)
    call gpu_mmthr(D,Ndim,Ndim,1); call gpu_mmthl_m1(D,Ndim,Ndim,1); call gpu_mmthr_m1(D,Ndim,Ndim,1); call gpu_mmthl(D,Ndim,Ndim,1)
    call la_sync_if()
    nrep = 100
    do itr = 1, 5
       D = Hs
       !$omp target update to(D)
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call gpu_mmthr(D,Ndim,Ndim,1); call gpu_mmthl_m1(D,Ndim,Ndim,1)
          call gpu_mmthr_m1(D,Ndim,Ndim,1); call gpu_mmthl(D,Ndim,Ndim,1)
       enddo
       call la_sync_if()
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(2*nrep, kind(0.d0))
    enddo
    tg = median5(tt)
    deallocate(H)
  end subroutine time_pair

  ! fused Op_V slice (all operators, rows 'n', and columns with alternating sign so that the matrix stays bounded)
  Subroutine time_vslice(Hs, D)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: Hs(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), allocatable :: H(:,:)
    Real (Kind=Kind(0.d0)) :: t0, t1, tt(5), tg, tc
    Integer :: itr, irep, nrep, n, t

    allocate(H(Ndim,Ndim))
    t = 1
    ! rows ('n'): Op_mmultR for n = 1..N_op; reps limited so that the repeated factors cannot overflow
    nrep = 20
    do itr = 1, 5
       H = Hs
       t0 = omp_get_wtime()
       do irep = 1, nrep
          do n = 1, size(Op_V,1);  call Op_mmultR(H, Op_V(n,1), nsigma%f(n,t), 'n', t); enddo
       enddo
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(nrep, kind(0.d0))
    enddo
    tc = median5(tt)
    D = Hs
    !$omp target update to(D)
    call gpu_Vslice_R(D, Ndim, Ndim, 1, t, 'n', .false.); call la_sync_if()
    do itr = 1, 5
       D = Hs
       !$omp target update to(D)
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call gpu_Vslice_R(D, Ndim, Ndim, 1, t, 'n', .false.)
       enddo
       call la_sync_if()                 ! target regions are synchronous; this only matters for alf_la calls
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(nrep, kind(0.d0))
    enddo
    tg = median5(tt)
    call report_time('Op_V slice rows (fused)', tg, tc)
    ! per-operator launches for comparison (rows): one launch per operator
    D = Hs
    !$omp target update to(D)
    do itr = 1, 5
       D = Hs
       !$omp target update to(D)
       t0 = omp_get_wtime()
       do irep = 1, nrep
          do n = 1, size(Op_V,1);  call gpu_Vop_R(D, Ndim, Ndim, n, 1, t, 'n'); enddo
       enddo
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(nrep, kind(0.d0))
    enddo
    tg = median5(tt)
    call report_time('Op_V slice rows (one launch per operator)', tg, tc)
    deallocate(H)
  end subroutine time_vslice

  ! WRAPGRUP + WRAPGRDO without updates (up at t, down at t = restores G): us per wrap
  Subroutine time_wrapgr(Hs, D)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(in) :: Hs(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)), allocatable :: H(:,:)
    Real (Kind=Kind(0.d0)) :: t0, t1, tt(5), tg, tc, tone
    Integer :: itr, irep, nrep
    Integer(int64) :: l0, l1
    Character (len=256) :: line

    allocate(H(Ndim,Ndim))
    H = Hs
    t0 = omp_get_wtime()
    call cpu_wrapgr_up(H, 1, 0); call cpu_wrapgr_do(H, 1, 1)
    t1 = omp_get_wtime()
    tone = max(t1 - t0, 1.d-7)
    nrep = max(1, min(50, int(0.05d0/tone)))
    do itr = 1, 5
       H = Hs
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call cpu_wrapgr_up(H, 1, 0); call cpu_wrapgr_do(H, 1, 1)
       enddo
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(2*nrep, kind(0.d0))
    enddo
    tc = median5(tt)
    D = Hs
    !$omp target update to(D)
    call gpu_wrapgr_up(D, Ndim, 1, 0); call gpu_wrapgr_do(D, Ndim, 1, 1); call la_sync_if()
    nrep = 100
    l0 = gpu_n_launch
    do itr = 1, 5
       D = Hs
       !$omp target update to(D)
       t0 = omp_get_wtime()
       do irep = 1, nrep
          call gpu_wrapgr_up(D, Ndim, 1, 0); call gpu_wrapgr_do(D, Ndim, 1, 1)
       enddo
       call la_sync_if()
       t1 = omp_get_wtime()
       tt(itr) = (t1 - t0)/real(2*nrep, kind(0.d0))
    enddo
    l1 = gpu_n_launch
    tg = median5(tt)
    call report_time('WRAPGR wrap (B G B^-1 + all Op_V wraps)', tg, tc)
    write(line,'(A,F8.1,A)') 'GPU wrap timing WRAPGR launches per wrap (B G B^-1 + fused Op_V wraps): ', &
         & real(l1-l0,kind(0.d0))/real(10*nrep,kind(0.d0)), ' (alf_la calls not included)'
    call say(line)
    deallocate(H)
  end subroutine time_wrapgr

End Module gpu_wrap_mod
#endif

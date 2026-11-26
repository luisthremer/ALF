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
!> complex matrices, built on the flat arrays of gpu_state_mod and on linalg_mod.  Compiled only with -DALF_GPU.
!> Called by the sweep: gpu_update_mod (WRAPGRUP/DO, tickets 205/208), gpu_stab_mod (WRAPUR/WRAPUL slice products, 206) and
!> gpu_tau_mod (tau pass, 207); gpu_obs_mod uses gpu_symm.
!>
!> Conventions
!>   * every matrix argument is an explicit-shape/assumed-size device array M(ld,*) that the caller has mapped; the
!>     routines never allocate on the device (scratch = gpu_state scr_c(Ndim,Ndim,2)) and use synchronous target
!>     regions only (ADR-002 amended); linalg_mod calls are ordered with them.
!>   * left routines (mmthr, mmthr_m1, mmthlc): M is Ndim x ncol, rows are sites.  Right routines (mmthl, mmthl_m1):
!>     M is nrow x Ndim (nrow <= Ndim for the dense/gather scratch), columns are sites.
!>   * the field of operator n at slice t is read on the device from nsigma_d(n,t) (the caller keeps it in sync).
!>
!> Hopping kernels (see gpu_state_mod header for groups)
!>   checkerboard (all elements nd <= 8): one kernel launch per group (launches per B_T = number of groups, 4 for the
!>     square-lattice checkerboard), threads over (element, tile of gpu_tl columns / rows); the element data (site list,
!>     matrix) are read once per tile.  (A single-launch variant with a team barrier between the groups was removed in
!>     ticket 205: slower on the A30 and wrong on the Arc B580.)
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
!> Self-test: ALF_GPU_SELFTEST=2 (gpu_wrap_hook, called from Setup_gpu after gpu_state_setup) or the test executable
!> 51-gpu-wrap; ALF_GPU_FAULT=hop|hopc|vexp|nan|order injects a fault into the device
!> copy (negative controls, the self-test must fail).
!--------------------------------------------------------------------
#ifdef ALF_GPU
Module gpu_wrap_mod

  Use iso_fortran_env, only: output_unit, error_unit, int64
  Use runtime_error_mod, only: error_generic, terminate_on_error
  Use Natural_Constants, only: Eps_machine
  Use Hamiltonian_main, only: ltrot, ndim
  Use linalg_mod, only: la_zgemm, la_zgemm_cr
  Use gpu_state_mod, only: g_d, hop_ac, hop_at, hop_aux, hop_big, hop_cval, hop_gptr, hop_gstart, hop_meta, &
       & hop_p, hop_rval, nsigma_d, scr_c, td_any, td_grp, td_n, v_exp, v_meta, v_p, v_type, vr_ent, vr_ptr

  Implicit none

  Private
  Public :: gpu_mmthr, gpu_mmthr_m1, gpu_mmthlc, gpu_mmthl, gpu_mmthl_m1, gpu_symm, gpu_symm_build, gpu_symm_td
  Public :: gpu_mmthr_d, gpu_mmthl_m1_d, gpu_mmthl_d, gpu_mmthr_m1_d, gpu_hopd_on
  Public :: gpu_Vop_R, gpu_Vop_L, gpu_Vop_wrap, gpu_Vslice_R, gpu_Vslice_L, gpu_Vwrap_all
  Public :: gpu_wrapur_slices, gpu_wrapul_slices, gpu_wrapgr_up, gpu_wrapgr_do
  Public :: gpu_dense_left_mode, gpu_tl, gpu_tr, gpu_vb, gpu_n_launch, gpu_n_lib

  Integer, Save :: gpu_dense_left_mode = 0     ! 0 auto (v1 for nd <= 400, else v2), 1 complex copy + ZGEMM, 2 transposes + DGEMM (measured, ticket 204)
  Integer(int64), Save :: gpu_n_launch = 0, gpu_n_lib = 0    ! counters: target-region launches / linalg_mod calls issued here
  Logical, Save :: fault_order = .false.
  Integer, Save :: gpu_tl = 1, gpu_tr = 2   ! columns (left) / rows (right) per thread of the per-group kernels (tuned on the A30)
  Integer, Save :: gpu_vb = 8                 ! columns/rows per thread of the fused Op_V slice kernels
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

  Integer, parameter :: HOPD_NMAX = 800
  Logical, Save :: hopd_built = .false.
  Integer, Save :: st_t = 1                    ! slice used by the self-test hop calls (cpu_hop / dev_hop)
  Complex (Kind=Kind(0.d0)), Allocatable, Save :: Hd(:,:,:,:)

Contains

!====================================================================
!  Hopping
!====================================================================

!--------------------------------------------------------------------
  ! 1 mmthr: left, mat, elements descending;  2 mmthr_m1: left, invmat, ascending;  3 mmthl: right, mat, ascending;
  ! 4 mmthlc: left, mat, ascending;  5 mmthl_m1: right, invmat, descending
  ! Time-dependent elements (td_grp, gpu_state_mod) are the LAST operators of a flavour (gpu_supported): the host applies them first in
  ! the descending sequences (mmthr, mmthl_m1) and last in the ascending ones.  t = time slice (the host argument t of Hop_mod_mm*).
  Subroutine gpu_mmthr(M, ld, ncol, nf, t)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, t
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (td_any) call td_apply(M, ld, ncol, nf, t, 1, 1, .true., .true.)
    if (ld == Ndim .and. ncol == Ndim .and. gpu_hopd_on()) then
       call hopd_left(M, nf, 1); return
    endif
    call hop_left(M, ld, ncol, nf, 1, .true.)
  end subroutine gpu_mmthr

  Subroutine gpu_mmthr_m1(M, ld, ncol, nf, t)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, t
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (ld == Ndim .and. ncol == Ndim .and. gpu_hopd_on()) then
       call hopd_left(M, nf, 4)
    else
       call hop_left(M, ld, ncol, nf, 2, .false.)
    endif
    if (td_any) call td_apply(M, ld, ncol, nf, t, 1, 2, .true., .false.)
  end subroutine gpu_mmthr_m1

  Subroutine gpu_mmthlc(M, ld, ncol, nf, t)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, t
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    call hop_left(M, ld, ncol, nf, 1, .false.)
    if (td_any) call td_apply(M, ld, ncol, nf, t, 1, 1, .true., .false.)
  end subroutine gpu_mmthlc

  Subroutine gpu_mmthl(M, ld, nrow, nf, t)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf, t
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (ld == Ndim .and. nrow == Ndim .and. gpu_hopd_on()) then
       call hopd_right(M, nf, 3)
    else
       call hop_right(M, ld, nrow, nf, 1, .false.)
    endif
    if (td_any) call td_apply(M, ld, nrow, nf, t, 1, 1, .false., .false.)
  end subroutine gpu_mmthl

  Subroutine gpu_mmthl_m1(M, ld, nrow, nf, t)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf, t
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (td_any) call td_apply(M, ld, nrow, nf, t, 1, 2, .false., .true.)
    if (ld == Ndim .and. nrow == Ndim .and. gpu_hopd_on()) then
       call hopd_right(M, nf, 2); return
    endif
    call hop_right(M, ld, nrow, nf, 2, .true.)
  end subroutine gpu_mmthl_m1

!--------------------------------------------------------------------
!> The time-dependent elements of flavour nf at slice t (kind kd: 1 = propagation, 2 = Symm half step; t = 0 means Ltrot as in the
!> host adjointaction), matrix k, from the left (M <- A M, M has n columns) or from the right (M <- M A, n rows), in ascending
!> or descending operator order.  Inactive slices apply nothing (the host skips |g|^2 <= Zero).  Host-side dispatch only.
  Subroutine td_apply(M, ld, n, nf, t, kd, k, left, desc)
    Implicit none
    Integer, intent(in) :: ld, n, nf, t, kd, k
    Logical, intent(in) :: left, desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: j, j1, j2, js, tt, g, e
    if (td_n(nf) == 0 .or. n <= 0) return
    tt = t
    if (tt == 0) tt = Ltrot
    if (tt < 1 .or. tt > Ltrot) call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
    if (desc) then
       j1 = td_n(nf); j2 = 1; js = -1
    else
       j1 = 1; j2 = td_n(nf); js = 1
    endif
    do j = j1, j2, js
       g = td_grp(tt,kd,j,nf)
       if (g == 0) cycle
       e = hop_gptr(g)
       if (left) then
          if (hop_meta(2,e) > 8) then
             call dense_left(M, ld, n, e, k)
          else
             call k_cb_left_g(M, ld, n, g, k, gpu_tl); gpu_n_launch = gpu_n_launch + 1
          endif
       else
          if (hop_meta(2,e) > 8) then
             call dense_right(M, ld, n, e, k)
          else
             call k_cb_right_g(M, ld, n, g, k, gpu_tr); gpu_n_launch = gpu_n_launch + 1
          endif
       endif
    enddo
  end subroutine td_apply

!--------------------------------------------------------------------
!> Hop_mod_Symm on the device: M <- B^{1/2}-similarity, per element (descending) left mat_1D2 then right invmat_1D2.  The elements of a
!> group commute, so a group is applied as all lefts then all rights.  M is an n x n matrix (ld = leading dimension).
  Subroutine gpu_symm(M, ld, n, nf, t1, t2)
    Implicit none
    Integer, intent(in) :: ld, n, nf, t1
    Integer, intent(in), optional :: t2
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: ga, gb, g, e, tr
    tr = t1
    if (present(t2)) tr = t2
    ! host order: elements descending, so the (last) time-dependent ones first: left A(g_t(t1)/2), right A(-g_t(t2)/2)
    if (td_any) then
       call td_apply(M, ld, n, nf, t1, 2, 3, .true., .true.)
       call td_apply(M, ld, n, nf, tr, 2, 4, .false., .true.)
    endif
    ga = hop_gstart(nf); gb = hop_gstart(nf+1) - 1
    if (gb < ga .or. n <= 0) return
    do g = gb, ga, -1
       if (.not. hop_big(nf)) then
          call k_cb_left_g(M, ld, n, g, 3, gpu_tl)
          call k_cb_right_g(M, ld, n, g, 4, gpu_tr)
          gpu_n_launch = gpu_n_launch + 2
       else
          e = hop_gptr(g)
          if (hop_meta(2,e) > 8) then
             call dense_left(M, ld, n, e, 3)
             call dense_right(M, ld, n, e, 4)
          else
             call k_cb_left_g(M, ld, n, g, 3, gpu_tl)
             call k_cb_right_g(M, ld, n, g, 4, gpu_tr)
             gpu_n_launch = gpu_n_launch + 2
          endif
       endif
    enddo
  end subroutine gpu_symm

!--------------------------------------------------------------------
!  Dense hopping products (small and medium N): every hop sequence is a fixed linear map, M -> Bl M or M -> M Br.  The four matrices
!  are the sequences applied to the identity (the same kernels), built once; applying one is a GEMM plus a copy instead of one launch
!  per checkerboard group.  Used for Ndim <= HOPD_NMAX; above that the group kernels are cheaper.
!    1 mmthr (left, mat, desc)   2 mmthl_m1 (right, invmat, desc)   3 mmthl (right, mat, asc)   4 mmthr_m1 (left, invmat, asc)
!--------------------------------------------------------------------
  Logical Function gpu_hopd_on()
    Implicit none
    gpu_hopd_on = (Ndim <= HOPD_NMAX)
  end function gpu_hopd_on

  Subroutine hopd_build()
    Implicit none
    Integer :: nf, w, i, nfl
    Complex (Kind=Kind(0.d0)), allocatable :: Id(:,:)
    if (hopd_built) return
    nfl = size(G_d,3)
    allocate(Hd(Ndim,Ndim,nfl,4))
    !$omp target enter data map(alloc: Hd)
    allocate(Id(Ndim,Ndim)); Id = (0.d0,0.d0)
    do i = 1, Ndim
       Id(i,i) = (1.d0,0.d0)
    enddo
    do nf = 1, nfl
       do w = 1, 4
          call k_set_c(Ndim*Ndim, Id, Hd(1,1,nf,w))
       enddo
       call hop_left (Hd(1,1,nf,1), Ndim, Ndim, nf, 1, .true.)
       call hop_right(Hd(1,1,nf,2), Ndim, Ndim, nf, 2, .true.)
       call hop_right(Hd(1,1,nf,3), Ndim, Ndim, nf, 1, .false.)
       call hop_left (Hd(1,1,nf,4), Ndim, Ndim, nf, 2, .false.)
    enddo
    deallocate(Id)
    hopd_built = .true.
  end subroutine hopd_build

  Subroutine hopd_left(M, nf, w)          ! M <- Hd(w) M
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(Ndim,*)
    Integer, intent(in) :: nf, w
    call hopd_build()
    call la_zgemm('N', 'N', Ndim, Ndim, Ndim, z1, Hd(1,1,nf,w), Ndim, M, Ndim, z0, scr_c(1,1,1), Ndim)
    call k_cpy_n(Ndim*Ndim, scr_c(1,1,1), M)
  end subroutine hopd_left

  Subroutine hopd_right(M, nf, w)         ! M <- M Hd(w)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(Ndim,*)
    Integer, intent(in) :: nf, w
    call hopd_build()
    call la_zgemm('N', 'N', Ndim, Ndim, Ndim, z1, M, Ndim, Hd(1,1,nf,w), Ndim, z0, scr_c(1,1,1), Ndim)
    call k_cpy_n(Ndim*Ndim, scr_c(1,1,1), M)
  end subroutine hopd_right

  Subroutine k_cpy_n(nel, a, b)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)), intent(in) :: a(1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: b(1)
    Integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nel
       b(i) = a(i)
    enddo
  end subroutine k_cpy_n

  Subroutine gpu_mmthr_d(M, nf)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(Ndim,*)
    Integer, intent(in) :: nf
    call hopd_left(M, nf, 1)
  end subroutine gpu_mmthr_d
  Subroutine gpu_mmthl_m1_d(M, nf)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(Ndim,*)
    Integer, intent(in) :: nf
    call hopd_right(M, nf, 2)
  end subroutine gpu_mmthl_m1_d
  Subroutine gpu_mmthl_d(M, nf)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(Ndim,*)
    Integer, intent(in) :: nf
    call hopd_right(M, nf, 3)
  end subroutine gpu_mmthl_d
  Subroutine gpu_mmthr_m1_d(M, nf)
    Implicit none
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(Ndim,*)
    Integer, intent(in) :: nf
    call hopd_left(M, nf, 4)
  end subroutine gpu_mmthr_m1_d

!> Symm(M) = A M B for every n x n matrix M: A (left factors) and B (right factors) are the sequences of gpu_symm applied to the
!> identity (the lefts and rights commute).  STATIC elements only: the time-dependent ones depend on the slice and are applied to the
!> data by gpu_symm_td (Symm(M) = A (A_td M B_td) B).  A and B must be device-resident n x n arrays with leading dimension ld; they are
!> set here (the identity is written by the host and uploaded).
  !> Time-dependent part of Hop_mod_Symm at slice t (t = 0 means Ltrot): M <- A_td(t) M B_td(t), in place; nothing if no element is active.
  Subroutine gpu_symm_td(M, ld, n, nf, t)
    Implicit none
    Integer, intent(in) :: ld, n, nf, t
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    if (.not. td_any) return
    call td_apply(M, ld, n, nf, t, 2, 3, .true., .true.)
    call td_apply(M, ld, n, nf, t, 2, 4, .false., .true.)
  end subroutine gpu_symm_td

  Subroutine gpu_symm_build(A, B, ld, n, nf)
    Implicit none
    Integer, intent(in) :: ld, n, nf
    Complex (Kind=Kind(0.d0)), intent(inout) :: A(ld,*), B(ld,*)
    Integer :: ga, gb, g, e, i
    Complex (Kind=Kind(0.d0)), allocatable :: Id(:,:)
    allocate(Id(ld,n)); Id = (0.d0,0.d0)
    do i = 1, n
       Id(i,i) = (1.d0,0.d0)
    enddo
    call k_set_c(ld*n, Id, A)
    call k_set_c(ld*n, Id, B)
    deallocate(Id)
    ga = hop_gstart(nf); gb = hop_gstart(nf+1) - 1
    if (gb < ga) return
    do g = gb, ga, -1
       if (.not. hop_big(nf)) then
          call k_cb_left_g(A, ld, n, g, 3, gpu_tl)
          call k_cb_right_g(B, ld, n, g, 4, gpu_tr)
       else
          e = hop_gptr(g)
          if (hop_meta(2,e) > 8) then
             call dense_left(A, ld, n, e, 3)
             call dense_right(B, ld, n, e, 4)
          else
             call k_cb_left_g(A, ld, n, g, 3, gpu_tl)
             call k_cb_right_g(B, ld, n, g, 4, gpu_tr)
          endif
       endif
    enddo
  end subroutine gpu_symm_build

  Subroutine k_set_c(nel, src, dst)
    Implicit none
    Integer, intent(in) :: nel
    Complex (Kind=Kind(0.d0)), intent(in) :: src(1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: dst(1)
    Integer :: i
    !$omp target teams distribute parallel do map(to: src(1:nel))
    do i = 1, nel
       dst(i) = src(i)
    enddo
  end subroutine k_set_c

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
       do g = first, last, step
          call k_cb_left_g(M, ld, ncol, g, k, gpu_tl)
          gpu_n_launch = gpu_n_launch + 1
       enddo
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
       do g = first, last, step
          call k_cb_right_g(M, ld, nrow, g, k, gpu_tr)
          gpu_n_launch = gpu_n_launch + 1
       enddo
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
!> GCC amdgcn: the `!GCC$ novector` lines are required, not tuning.  GCC maps OpenMP threads to wavefronts and vectorises
!> loops over the 64 lanes; vectorising these short loops (tile length tl = 1..2, nd <= 8) inlined into the group kernels
!> gave a 47-50 KB stack frame per wavefront (default GCN_STACK_SIZE 32 KB): stack overflow -> GPU memory fault, which
!> with the HIP runtime in the process is not reported, so the host waited forever (README of the AMD platform).
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
          !GCC$ novector
          do j = j1, j2
             M(r,j) = cmplx(a*real(M(r,j)), a*aimag(M(r,j)), kind(0.d0))
          enddo
       else if (nd == 2) then
          r1 = hP(po+1); r2 = hP(po+2)
          a11 = hR(mo+1); a21 = hR(mo+2); a12 = hR(mo+3); a22 = hR(mo+4)
          !GCC$ novector
          do j = j1, j2
             w1 = M(r1,j); w2 = M(r2,j)
             zr = a11*real(w1) + a12*real(w2); zi = a11*aimag(w1) + a12*aimag(w2)
             M(r1,j) = cmplx(zr, zi, kind(0.d0))
             zr = a21*real(w1) + a22*real(w2); zi = a21*aimag(w1) + a22*aimag(w2)
             M(r2,j) = cmplx(zr, zi, kind(0.d0))
          enddo
       else
          !GCC$ novector
          do j = j1, j2
             !GCC$ novector
             do l = 1, nd
                w(l) = M(hP(po+l), j)
             enddo
             !GCC$ novector
             do i = 1, nd
                a = hR(mo + i)
                zr = a*real(w(1)); zi = a*aimag(w(1))
                !GCC$ novector
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
          !GCC$ novector
          do j = j1, j2
             M(r,j) = c11 * M(r,j)
          enddo
       else if (nd == 2) then
          r1 = hP(po+1); r2 = hP(po+2)
          c11 = hC(mo+1); c21 = hC(mo+2); c12 = hC(mo+3); c22 = hC(mo+4)
          !GCC$ novector
          do j = j1, j2
             w1 = M(r1,j); w2 = M(r2,j)
             M(r1,j) = c11*w1 + c12*w2
             M(r2,j) = c21*w1 + c22*w2
          enddo
       else
          !GCC$ novector
          do j = j1, j2
             !GCC$ novector
             do l = 1, nd
                w(l) = M(hP(po+l), j)
             enddo
             !GCC$ novector
             do i = 1, nd
                z = hC(mo + i) * w(1)
                !GCC$ novector
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
          !GCC$ novector
          do i = i1, i2, istep
             M(i,c) = cmplx(a*real(M(i,c)), a*aimag(M(i,c)), kind(0.d0))
          enddo
       else if (nd == 2) then
          c1 = hP(po+1); c2 = hP(po+2)
          a11 = hR(mo+1); a21 = hR(mo+2); a12 = hR(mo+3); a22 = hR(mo+4)
          !GCC$ novector
          do i = i1, i2, istep
             w1 = M(i,c1); w2 = M(i,c2)
             zr = a11*real(w1) + a21*real(w2); zi = a11*aimag(w1) + a21*aimag(w2)
             M(i,c1) = cmplx(zr, zi, kind(0.d0))
             zr = a12*real(w1) + a22*real(w2); zi = a12*aimag(w1) + a22*aimag(w2)
             M(i,c2) = cmplx(zr, zi, kind(0.d0))
          enddo
       else
          !GCC$ novector
          do i = i1, i2, istep
             !GCC$ novector
             do l = 1, nd
                w(l) = M(i, hP(po+l))
             enddo
             !GCC$ novector
             do j = 1, nd
                a = hR(mo + 1 + (j-1)*nd)
                zr = a*real(w(1)); zi = a*aimag(w(1))
                !GCC$ novector
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
          !GCC$ novector
          do i = i1, i2, istep
             M(i,c) = M(i,c) * c11
          enddo
       else if (nd == 2) then
          c1 = hP(po+1); c2 = hP(po+2)
          c11 = hC(mo+1); c21 = hC(mo+2); c12 = hC(mo+3); c22 = hC(mo+4)
          !GCC$ novector
          do i = i1, i2, istep
             w1 = M(i,c1); w2 = M(i,c2)
             M(i,c1) = w1*c11 + w2*c21
             M(i,c2) = w1*c12 + w2*c22
          enddo
       else
          !GCC$ novector
          do i = i1, i2, istep
             !GCC$ novector
             do l = 1, nd
                w(l) = M(i, hP(po+l))
             enddo
             !GCC$ novector
             do j = 1, nd
                z = w(1) * hC(mo + 1 + (j-1)*nd)
                !GCC$ novector
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
!> One launch per group, flat: threads over (column tile, element) / (element, row tile); a thread handles
!> tl columns (left: consecutive) or rows (right: strided by the tile count), the element data are read once per thread.
  Subroutine k_cb_left_g(M, ld, ncol, g, k, tl)
    Implicit none
    Integer, intent(in) :: ld, ncol, g, k, tl
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
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
    mode = gpu_dense_left_mode
    if (mode == 0) mode = merge(1, 2, nd <= 400)
    if (hop_meta(1,e) == 1 .and. mode == 2) then
       ! (A M)^T = M^T A^T: transpose-gather, DGEMM on the interleaved view with A^T, transpose-scatter
       call k_tgather(nd, ncol, hop_P(po+1), M, ld, scr_c(1,1,1))               ! scr1(j,l) = M(P(l),j), ld = ncol
       call la_zgemm_cr(ncol, nd, nd, scr_c(1,1,1), ncol, hop_At(xo + (k-1)*nd*nd + 1), nd, scr_c(1,1,2), ncol)
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
       call la_zgemm('N', 'N', nd, ncol, nd, z1, hop_Cval(mo+1), nd, B, ldb, z0, C, Ndim)
    else
       call la_zgemm('N', 'N', nd, ncol, nd, z1, hop_Ac(xo + (k-1)*nd*nd + 1), nd, B, ldb, z0, C, Ndim)
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
       call la_zgemm('N', 'N', nrow, nd, nd, z1, B, ldb, hop_Cval(mo+1), nd, z0, C, nrow)
    else
       call la_zgemm_cr(nrow, nd, nd, B, ldb, hop_Rval(mo+1), nd, C, nrow)
    endif
  end subroutine gemm_right

  Subroutine k_gather_rows(nd, ncol, P, M, ld, W, ldw)
    Implicit none
    Integer, intent(in) :: nd, ncol, ld, ldw
    Integer, intent(in) :: P(1)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,1)
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
    Integer, intent(in) :: P(1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(ldw,1)
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
    Integer, intent(in) :: P(1)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ncol,1)
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
    Integer, intent(in) :: P(1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(ncol,1)
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
    Integer, intent(in) :: P(1)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(nrow,1)
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
    Integer, intent(in) :: P(1)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(nrow,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,1)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,1)
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
       call gpu_mmthr(U, ld, ncol, nf, t)
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
       call gpu_mmthlc(U, ld, ncol, nf, t)
    enddo
  end subroutine gpu_wrapul_slices

!> WRAPGRUP without updates: G = B(ntau+1) G B(ntau+1)^-1  (mmthr, mmthl_m1, then Op_Wrapup N_Type 1 for n = 1..N_op)
  Subroutine gpu_wrapgr_up(G, ld, nf, ntau)
    Implicit none
    Integer, intent(in) :: ld, nf, ntau
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    call gpu_mmthr(G, ld, Ndim, nf, ntau + 1)
    call gpu_mmthl_m1(G, ld, Ndim, nf, ntau + 1)
    call gpu_Vwrap_all(G, ld, nf, ntau + 1, 1)
  end subroutine gpu_wrapgr_up

!> WRAPGRDO without updates: G = B(ntau)^-1 G B(ntau)  (Op_Wrapdo N_Type 1 for n = N_op..1, then mmthl, mmthr_m1)
  Subroutine gpu_wrapgr_do(G, ld, nf, ntau)
    Implicit none
    Integer, intent(in) :: ld, nf, ntau
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    call gpu_Vwrap_all(G, ld, nf, ntau, -1)
    call gpu_mmthl(G, ld, Ndim, nf, ntau)
    call gpu_mmthr_m1(G, ld, Ndim, nf, ntau)
  end subroutine gpu_wrapgr_do
End Module gpu_wrap_mod
#endif

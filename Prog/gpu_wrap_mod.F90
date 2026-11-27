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
  Use Operator_mod, only: op_mmultl, op_mmultr, op_wrapdo, op_wrapup, operator
  Use Hamiltonian_main, only: ltrot, n_fl, ndim, nsigma, op_v
  Use QMC_runtime_var, only: get_nwrap
  Use Hop_mod, only: hop_mod_mmthl, hop_mod_mmthl_m1, hop_mod_mmthlc, hop_mod_mmthr, hop_mod_mmthr_m1, &
       & hop_mod_symm
  Use linalg_mod, only: la_zgemm, la_zgemm_cr
#if defined(ALF_LA_HIP)
  Use alf_hipk_mod, only: hk_enabled, hk_copyz, hk_scatter_rows, hk_vslice_r
#endif
  Use gpu_state_mod, only: finite_c, g_d, gpu_env_get, gpu_env_int, gpu_env_is, hop_ac, hop_at, hop_aux, hop_big, &
       & hop_cval, hop_gptr, hop_gstart, hop_meta, hop_p, hop_rval, hop_start, info_line, is_finite, nsigma_d, &
       & rand_mat, rnd, scr_c, selftest_fail, td_any, td_grp, td_n, v_exp, v_meta, v_p, v_type, vr_ent, vr_ptr, &
       & v_u, v_m, v_e, v_g, v_phi4, v_sstart, v_seg, v_hasnd, v_nmax, v_nd

  Implicit none

  Private
  Public :: gpu_mmthr, gpu_mmthr_m1, gpu_mmthlc, gpu_mmthl, gpu_mmthl_m1, gpu_symm, gpu_symm_build, gpu_symm_td
  Public :: gpu_mmthr_d, gpu_mmthl_m1_d, gpu_mmthl_d, gpu_mmthr_m1_d, gpu_hopd_on
  Public :: gpu_Vop_R, gpu_Vop_L, gpu_Vop_wrap, gpu_Vslice_R, gpu_Vslice_L, gpu_Vwrap_all, k_nd_rows
  Public :: gpu_wrapur_slices, gpu_wrapul_slices, gpu_wrapgr_up, gpu_wrapgr_do
  Public :: gpu_wrap_hook, gpu_wrap_selftest
  Public :: gpu_dense_left_mode, gpu_tl, gpu_tr, gpu_vb, gpu_n_launch, gpu_n_lib

  Integer, Save :: gpu_dense_left_mode = 0     ! 0 auto (v1 for nd <= 400, else v2), 1 complex copy + ZGEMM, 2 transposes + DGEMM (measured, ticket 204)
  Integer(int64), Save :: gpu_n_launch = 0, gpu_n_lib = 0    ! counters: target-region launches / linalg_mod calls issued here
  Logical, Save :: fault_order = .false.
  Integer, Save :: gpu_tl = 1, gpu_tr = 2   ! columns (left) / rows (right) per thread of the per-group kernels (tuned on the A30)
  Integer, Save :: gpu_vb = 8                 ! columns/rows per thread of the fused Op_V slice kernels
  Complex (Kind=Kind(0.d0)), Parameter :: z1 = (1.d0, 0.d0), z0 = (0.d0, 0.d0)

  Integer, parameter :: HOPD_NMAX = 800
  Logical, Save :: hopd_built = .false.
  Logical, Save :: hopd_opt = .false.          ! ALF_GPU_HOPD=1: dense hopping product Hd (opt-in; KNOWN_ISSUES G2)
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
!  per checkerboard group.  Used for Ndim <= HOPD_NMAX; above that the group kernels are cheaper.  Opt-in (ALF_GPU_HOPD=1): the
!  once-rounded product can make the wrapped G less accurate than the factor-by-factor path of Hop_mod (Precision Green 1.4-2.5x
!  the CPU's with Apple Accelerate, not with the reference BLAS; KNOWN_ISSUES G2); the default is the per-group path, which gives
!  the CPU's precision.
!    1 mmthr (left, mat, desc)   2 mmthl_m1 (right, invmat, desc)   3 mmthl (right, mat, asc)   4 mmthr_m1 (left, invmat, asc)
!--------------------------------------------------------------------
  Logical Function gpu_hopd_on()
    Implicit none
    gpu_hopd_on = hopd_opt .and. (Ndim <= HOPD_NMAX)
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
    Complex (Kind=Kind(0.d0)), intent(in) :: a(nel)
    Complex (Kind=Kind(0.d0)), intent(inout) :: b(nel)
    Integer :: i
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_copyz(nel, a, b)
       return
    endif
#endif
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
    Complex (Kind=Kind(0.d0)), intent(in) :: src(nel)
    Complex (Kind=Kind(0.d0)), intent(inout) :: dst(nel)
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
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do j = j1, j2
             M(r,j) = cmplx(a*real(M(r,j)), a*aimag(M(r,j)), kind(0.d0))
          enddo
       else if (nd == 2) then
          r1 = hP(po+1); r2 = hP(po+2)
          a11 = hR(mo+1); a21 = hR(mo+2); a12 = hR(mo+3); a22 = hR(mo+4)
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do j = j1, j2
             w1 = M(r1,j); w2 = M(r2,j)
             zr = a11*real(w1) + a12*real(w2); zi = a11*aimag(w1) + a12*aimag(w2)
             M(r1,j) = cmplx(zr, zi, kind(0.d0))
             zr = a21*real(w1) + a22*real(w2); zi = a21*aimag(w1) + a22*aimag(w2)
             M(r2,j) = cmplx(zr, zi, kind(0.d0))
          enddo
       else
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do j = j1, j2
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do l = 1, nd
                w(l) = M(hP(po+l), j)
             enddo
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do i = 1, nd
                a = hR(mo + i)
                zr = a*real(w(1)); zi = a*aimag(w(1))
                !GCC$ novector
#if defined(ALF_LA_HIP)
                !GCC$ novector
#endif
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
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do j = j1, j2
             M(r,j) = c11 * M(r,j)
          enddo
       else if (nd == 2) then
          r1 = hP(po+1); r2 = hP(po+2)
          c11 = hC(mo+1); c21 = hC(mo+2); c12 = hC(mo+3); c22 = hC(mo+4)
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do j = j1, j2
             w1 = M(r1,j); w2 = M(r2,j)
             M(r1,j) = c11*w1 + c12*w2
             M(r2,j) = c21*w1 + c22*w2
          enddo
       else
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do j = j1, j2
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do l = 1, nd
                w(l) = M(hP(po+l), j)
             enddo
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do i = 1, nd
                z = hC(mo + i) * w(1)
                !GCC$ novector
#if defined(ALF_LA_HIP)
                !GCC$ novector
#endif
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
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do i = i1, i2, istep
             M(i,c) = cmplx(a*real(M(i,c)), a*aimag(M(i,c)), kind(0.d0))
          enddo
       else if (nd == 2) then
          c1 = hP(po+1); c2 = hP(po+2)
          a11 = hR(mo+1); a21 = hR(mo+2); a12 = hR(mo+3); a22 = hR(mo+4)
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do i = i1, i2, istep
             w1 = M(i,c1); w2 = M(i,c2)
             zr = a11*real(w1) + a21*real(w2); zi = a11*aimag(w1) + a21*aimag(w2)
             M(i,c1) = cmplx(zr, zi, kind(0.d0))
             zr = a12*real(w1) + a22*real(w2); zi = a12*aimag(w1) + a22*aimag(w2)
             M(i,c2) = cmplx(zr, zi, kind(0.d0))
          enddo
       else
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do i = i1, i2, istep
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do l = 1, nd
                w(l) = M(i, hP(po+l))
             enddo
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do j = 1, nd
                a = hR(mo + 1 + (j-1)*nd)
                zr = a*real(w(1)); zi = a*aimag(w(1))
                !GCC$ novector
#if defined(ALF_LA_HIP)
                !GCC$ novector
#endif
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
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do i = i1, i2, istep
             M(i,c) = M(i,c) * c11
          enddo
       else if (nd == 2) then
          c1 = hP(po+1); c2 = hP(po+2)
          c11 = hC(mo+1); c21 = hC(mo+2); c12 = hC(mo+3); c22 = hC(mo+4)
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do i = i1, i2, istep
             w1 = M(i,c1); w2 = M(i,c2)
             M(i,c1) = w1*c11 + w2*c21
             M(i,c2) = w1*c12 + w2*c22
          enddo
       else
          !GCC$ novector
#if defined(ALF_LA_HIP)
          !GCC$ novector
#endif
          do i = i1, i2, istep
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do l = 1, nd
                w(l) = M(i, hP(po+l))
             enddo
             !GCC$ novector
#if defined(ALF_LA_HIP)
             !GCC$ novector
#endif
             do j = 1, nd
                z = w(1) * hC(mo + 1 + (j-1)*nd)
                !GCC$ novector
#if defined(ALF_LA_HIP)
                !GCC$ novector
#endif
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,ncol)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,Ndim)
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
    Integer, intent(in) :: P(nd)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,ncol)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ldw,ncol)
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
    Integer, intent(in) :: P(nd)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,ncol)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(ldw,ncol)
    Integer :: j, l
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_scatter_rows(nd, ncol, P, M, ld, W, ldw)
       return
    endif
#endif
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
    Integer, intent(in) :: P(nd)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,ncol)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(ncol,nd)
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
    Integer, intent(in) :: P(nd)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,ncol)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(ncol,nd)
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
    Integer, intent(in) :: P(nd)
    Complex (Kind=Kind(0.d0)), intent(in) :: M(ld,Ndim)
    Complex (Kind=Kind(0.d0)), intent(inout) :: W(nrow,nd)
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
    Integer, intent(in) :: P(nd)
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,Ndim)
    Complex (Kind=Kind(0.d0)), intent(in) :: W(nrow,nd)
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
    if (v_nd(n,nf) == 1) then                              ! P2.6.1: ZSLGEMM('L', cop, exp(g Phi O))
       call k_nd_rows(M, ld, ncol, nf, t, n, n, mexp_mode(cop), 1, 0)
       gpu_n_launch = gpu_n_launch + 1
       return
    endif
    cj = (cop == 'c' .or. cop == 'C')
    call k_vop_R(M, ld, ncol, n, nf, t, cj, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vop_R

  Subroutine k_vop_R(M, ld, ncol, n, nf, t, cj, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, ncol, n, nf, t, N1, po, eo, nt
    Logical, intent(in) :: cj
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,ncol)
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
    if (v_nd(n,nf) == 1) then                              ! P2.6.1: ZSLGEMM('r', cop, exp(sign g Phi O))
       call k_nd_cols(M, ld, nrow, nf, t, n, n, mexp_mode(cop), sign)
       gpu_n_launch = gpu_n_launch + 1
       return
    endif
    cj = (cop == 'c' .or. cop == 'C')
    call k_vop_L(M, ld, nrow, n, nf, t, cj, sign, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vop_L

  Subroutine k_vop_L(M, ld, nrow, n, nf, t, cj, sign, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nrow, n, nf, t, sign, N1, po, eo, nt
    Logical, intent(in) :: cj
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,Ndim)
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
    if (v_nd(n,nf) == 1) then
       call nd_wrap(G, ld, n, nf, t, dir, ntype)            ! P2.6.1: non-diagonal or type-4 operator
       return
    endif
    if (ntype /= 1) return
    call k_vop_wrap_rows(G, ld, Ndim, n, t, dir, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    call k_vop_wrap_cols(G, ld, Ndim, n, t, dir, V_meta(1,n,nf), V_meta(2,n,nf), V_meta(3,n,nf), V_type(n))
    gpu_n_launch = gpu_n_launch + 2
  end subroutine gpu_Vop_wrap

  Subroutine k_vop_wrap_rows(G, ld, nd, n, t, dir, N1, po, eo, nt)
    Implicit none
    Integer, intent(in) :: ld, nd, n, t, dir, N1, po, eo, nt
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,nd)
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
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,nd)
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
!> P2.6.1 (port of gpu-port ticket 401-A): Op_Wrapup (dir = +1) / Op_Wrapdo (dir = -1) of a non-table operator (non-diagonal, or
!> type 4) on G (ld x Ndim), N_Type 1 or 2, in Operator_mod's ZSLGEMM order.  Wrapup N1: columns <- G(:,P) U D(-) ('r','n'), then
!> rows <- D(+) U^H G(P,:) ('l','T');  N2: rows <- U G(P,:), columns <- G(:,P) U^H.  Wrapdo N1: rows <- U D(-) G(P,:), columns <-
!> G(:,P) D(+) U^H;  N2: columns <- G(:,P) U, rows <- U^H G(P,:).  N_Type 2 does nothing for a diagonal operator.  No skip check.
  Subroutine nd_wrap(G, ld, n, nf, t, dir, ntype)
    Implicit none
    Integer, intent(in) :: ld, n, nf, t, dir, ntype
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    if (ntype == 2 .and. V_meta(6,n,nf) == 1) return
    if (dir > 0 .and. ntype == 1) then
       call k_nd_cols(G, ld, Ndim, nf, t, n, n, 4, -1)
       call k_nd_rows(G, ld, Ndim, nf, t, n, n, 5, 1, 0)
    elseif (dir > 0) then
       call k_nd_rows(G, ld, Ndim, nf, t, n, n, 6, 1, 0)
       call k_nd_cols(G, ld, Ndim, nf, t, n, n, 7, 1)
    elseif (ntype == 1) then
       call k_nd_rows(G, ld, Ndim, nf, t, n, n, 4, -1, 0)
       call k_nd_cols(G, ld, Ndim, nf, t, n, n, 5, 1)
    else
       call k_nd_cols(G, ld, Ndim, nf, t, n, n, 6, 1)
       call k_nd_rows(G, ld, Ndim, nf, t, n, n, 7, 1, 0)
    endif
    gpu_n_launch = gpu_n_launch + 2
  end subroutine nd_wrap

  Integer function mexp_mode(cop)
    Implicit none
    Character (len=1), intent(in) :: cop
    mexp_mode = 1
    if (cop == 'T' .or. cop == 't') mexp_mode = 2
    if (cop == 'C' .or. cop == 'c') mexp_mode = 3
  end function mexp_mode

!> P2.6.1: element (a,b) of the N x N matrix WORK that ZSLGEMM applies for operator (n,nf) with field f, as the inline code of
!> k_nd_rows / k_nd_cols.  mode 1: exp(sf g Phi O) (= M_exp(sf*s) for types 1, 2; Op_exp on the fly for type 4), 2: its transpose,
!> 3: its conjugate transpose (Op_mmultL/R with cop n, T, c);  4: U(a,b) d(b), 5: d(a) conjg(U(b,a)) with d = E_exp(:, sf*s)
!> (types 1, 2) or exp(sf Phi g E) (type 4) (Op_Wrapup/Op_Wrapdo N_Type 1);  6: U, 7: U^H (N_Type 2).
!> zk(1:N) must hold exp(gp E) for modes 1-3 / d for modes 4, 5 of a type-4 operator (k_nd_factors).
!> Arrays are passed explicitly (device routine, as el_left_d): Vm = V_M, Vu = V_U, Vx = V_exp; mt = (N, eoff, uoff, moff, type).
  Complex (Kind=Kind(0.d0)) function nd_work(mt, mode, s, a, b, zk, Vm, Vu, Vx)
    Implicit none
    !$omp declare target
    Integer, intent(in) :: mt(5), mode, s, a, b
    Complex (Kind=Kind(0.d0)), intent(in) :: zk(*), Vm(*), Vu(*), Vx(*)
    Integer :: N1, eo, uo, mo, nt, ia, ib, k
    Complex (Kind=Kind(0.d0)) :: w
    N1 = mt(1); eo = mt(2); uo = mt(3); mo = mt(4); nt = mt(5)
    ia = a; ib = b
    if (mode == 2 .or. mode == 3 .or. mode == 5 .or. mode == 7) then
       ia = b; ib = a
    endif
    if (mode <= 3) then
       if (nt <= 2) then
          w = Vm(mo + (s+nt)*N1*N1 + ia + (ib-1)*N1)
       else
          ! Op_exp: Mat(I,J) = sum_n (exp(g E_n) conjg(U(J,n))) U(I,n)  (equal to its Kahan sum for N <= 2)
          w = (zk(1)*conjg(Vu(uo + ib))) * Vu(uo + ia)
          do k = 2, N1
             w = w + (zk(k)*conjg(Vu(uo + ib + (k-1)*N1))) * Vu(uo + ia + (k-1)*N1)
          enddo
       endif
    else
       w = Vu(uo + ia + (ib-1)*N1)
    endif
    if (mode == 3 .or. mode == 5 .or. mode == 7) w = conjg(w)
    if (mode == 4) then
       if (nt <= 2) then
          w = w * Vx(eo + b + (s+nt)*N1)
       else
          w = w * zk(b)
       endif
    elseif (mode == 5) then
       if (nt <= 2) then
          w = Vx(eo + a + (s+nt)*N1) * w
       else
          w = zk(a) * w
       endif
    endif
    nd_work = w
  end function nd_work

!> P2.6.1: per-operator factors of a type-4 operator for nd_work: s = sf * nint(Re f) (table index, types 1, 2), and for type 4
!> zk(k) = exp(((sf g) Phi) E_k) (modes 1-3, Op_exp argument of Op_mmultL/R) or exp(sf ((Phi g) E_k)) (modes 4, 5, Op_Wrapup/do)
!> with Phi = Fields_Phi of type 4 = Phi_st(nint(Re f),2) sqrt(1 + Im f).
!> Arrays passed explicitly (device routine): E = V_E(poff+1:), ph4 = V_phi4.
  Subroutine nd_factors(N1, nt, g, E, ph4, f, mode, sf, s, zk)
    Implicit none
    !$omp declare target
    Integer, intent(in) :: N1, nt, mode, sf
    Complex (Kind=Kind(0.d0)), intent(in) :: g, f
    Real (Kind=Kind(0.d0)), intent(in) :: E(*), ph4(-2:2)
    Integer, intent(out) :: s
    Complex (Kind=Kind(0.d0)), intent(out) :: zk(*)
    Integer :: k
    Complex (Kind=Kind(0.d0)) :: phi, x
    s = sf*nint(real(f))
    if (nt <= 2 .or. mode >= 6) return
    phi = cmplx(ph4(nint(real(f))), 0.d0, kind(0.d0)) * sqrt(cmplx(1.d0 + aimag(f), 0.d0, kind(0.d0)))
    do k = 1, N1
       if (mode <= 3) then
          zk(k) = exp(((real(sf,kind(0.d0))*g)*phi)*E(k))
       else
          x = (phi*g)*E(k)
          if (sf < 0) x = -x
          zk(k) = exp(x)
       endif
    enddo
  end subroutine nd_factors

!> P2.6.1: rows P of M (ld x ncol) <- WORK M(P,:) for every operator na..nb of flavour nf (pairwise disjoint P within the range,
!> so one launch gives the sequential result), WORK = nd_work(mode) with field sign sf.  The L/R skip flag is honoured for the
!> modes 1-3 (Op_mmultR returns early for tiny g), ignored for the wrap modes.  Sums in ZSLGEMM's order.  Thread = (operator, column).
!> tw = 1 applies WORK^T instead (rows of Y in G = G0 - X Y^T after a column operation G <- G WORK, gpu_lowrank_mod).
  Subroutine k_nd_rows(M, ld, ncol, nf, t, na, nb, mode, sf, tw)
    Implicit none
    Integer, intent(in) :: ld, ncol, nf, t, na, nb, mode, sf, tw
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,ncol)
    Integer :: j, n, i, l, N1, po, s, mt(5)
    Complex (Kind=Kind(0.d0)) :: z(V_NMAX), zk(V_NMAX), acc
    !$omp target teams distribute parallel do collapse(2) private(i, l, N1, po, s, mt, z, zk, acc)
    do j = 1, ncol
       do n = na, nb
          if (.not. (mode <= 3 .and. V_meta(4,n,nf) == 1)) then
             N1 = V_meta(1,n,nf); po = V_meta(2,n,nf)
             mt(1) = N1; mt(2) = V_meta(3,n,nf); mt(3) = V_meta(7,n,nf); mt(4) = V_meta(8,n,nf); mt(5) = V_type(n)
             call nd_factors(N1, mt(5), V_g(n,nf), V_E(po+1), V_phi4, nsigma_d(n,t), mode, sf, s, zk)
             do l = 1, N1
                z(l) = M(V_P(po+l), j)
             enddo
             do i = 1, N1
                if (tw == 0) then
                   acc = nd_work(mt, mode, s, i, 1, zk, V_M, V_U, V_exp) * z(1)
                   do l = 2, N1
                      acc = acc + nd_work(mt, mode, s, i, l, zk, V_M, V_U, V_exp) * z(l)
                   enddo
                else
                   acc = nd_work(mt, mode, s, 1, i, zk, V_M, V_U, V_exp) * z(1)
                   do l = 2, N1
                      acc = acc + nd_work(mt, mode, s, l, i, zk, V_M, V_U, V_exp) * z(l)
                   enddo
                endif
                M(V_P(po+i), j) = acc
             enddo
          endif
       enddo
    enddo
  end subroutine k_nd_rows

!> P2.6.1: columns P of M (rows 1..nrow) <- M(:,P) WORK, as k_nd_rows.  Thread = (operator, row).
  Subroutine k_nd_cols(M, ld, nrow, nf, t, na, nb, mode, sf)
    Implicit none
    Integer, intent(in) :: ld, nrow, nf, t, na, nb, mode, sf
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,Ndim)
    Integer :: j, n, i, l, N1, po, s, mt(5)
    Complex (Kind=Kind(0.d0)) :: z(V_NMAX), zk(V_NMAX), acc
    !$omp target teams distribute parallel do collapse(2) private(i, l, N1, po, s, mt, z, zk, acc)
    do j = 1, nrow
       do n = na, nb
          if (.not. (mode <= 3 .and. V_meta(4,n,nf) == 1)) then
             N1 = V_meta(1,n,nf); po = V_meta(2,n,nf)
             mt(1) = N1; mt(2) = V_meta(3,n,nf); mt(3) = V_meta(7,n,nf); mt(4) = V_meta(8,n,nf); mt(5) = V_type(n)
             call nd_factors(N1, mt(5), V_g(n,nf), V_E(po+1), V_phi4, nsigma_d(n,t), mode, sf, s, zk)
             do l = 1, N1
                z(l) = M(j, V_P(po+l))
             enddo
             do i = 1, N1
                acc = z(1) * nd_work(mt, mode, s, 1, i, zk, V_M, V_U, V_exp)
                do l = 2, N1
                   acc = acc + z(l) * nd_work(mt, mode, s, l, i, zk, V_M, V_U, V_exp)
                enddo
                M(j, V_P(po+i)) = acc
             enddo
          endif
       enddo
    enddo
  end subroutine k_nd_cols

!> P2.6.1: table operators n1..n2 of flavour nf on the rows (side 1, Op_mmultR) or columns (side 2, Op_mmultL, field sign sg) of M,
!> per element in operator order (asc or desc), skip flag honoured: the fused kernels restricted to a run of the operator order.
!> mc = number of columns of M (nother for side 1, nd for side 2)
  Subroutine k_vrun(M, ld, mc, nother, nd, nf, t, cj, sg, desc, n1, n2, side)
    Implicit none
    Integer, intent(in) :: ld, mc, nother, nd, nf, t, sg, n1, n2, side
    Logical, intent(in) :: cj, desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,mc)
    Integer :: j, r, pa, pb, q, p, n, i, sp
    Complex (Kind=Kind(0.d0)) :: x, z
    !$omp target teams distribute parallel do collapse(2) private(pa, pb, q, p, n, i, sp, x, z)
    do j = 1, nother
       do r = 1, nd
          pa = vr_ptr(r,nf); pb = vr_ptr(r+1,nf) - 1
          if (pb >= pa) then
             if (side == 1) then
                x = M(r,j)
             else
                x = M(j,r)
             endif
             do q = 0, pb - pa
                if (desc) then
                   p = pb - q
                else
                   p = pa + q
                endif
                n = vr_ent(1,p); i = vr_ent(2,p)
                if (n >= n1 .and. n <= n2 .and. V_meta(4,n,nf) /= 1) then
                   sp = sg*nint(real(nsigma_d(n,t)))
                   z = V_exp(V_meta(3,n,nf) + i + (sp+V_type(n))*V_meta(1,n,nf))
                   if (cj) z = conjg(z)
                   x = z * x
                endif
             enddo
             if (side == 1) then
                M(r,j) = x
             else
                M(j,r) = x
             endif
          endif
       enddo
    enddo
  end subroutine k_vrun

!> P2.6.1: a slice of a flavour with non-table operators (Op_mmultR side 1 / Op_mmultL side 2 for n ascending, or descending if
!> desc): the segments of build_vseg in operator order, table runs by k_vrun, non-table groups by one dense-block launch each.
  Subroutine nd_slice(M, ld, nother, nf, t, cop, sg, desc, side)
    Implicit none
    Integer, intent(in) :: ld, nother, nf, t, sg, side
    Character (len=1), intent(in) :: cop
    Logical, intent(in) :: desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,*)
    Integer :: q, k, k1, k2
    k1 = v_sstart(nf); k2 = v_sstart(nf+1) - 1
    do q = 0, k2 - k1
       k = k1 + q
       if (desc) k = k2 - q
       if (v_seg(1,k) == 0) then
          call k_vrun(M, ld, merge(nother, Ndim, side == 1), nother, Ndim, nf, t, (cop == 'c' .or. cop == 'C'), sg, desc, v_seg(2,k), v_seg(3,k), side)
       elseif (side == 1) then
          call k_nd_rows(M, ld, nother, nf, t, v_seg(2,k), v_seg(3,k), mexp_mode(cop), 1, 0)
       else
          call k_nd_cols(M, ld, nother, nf, t, v_seg(2,k), v_seg(3,k), mexp_mode(cop), sg)
       endif
       gpu_n_launch = gpu_n_launch + 1
    enddo
  end subroutine nd_slice

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
    if (v_hasnd(nf)) then
       call nd_slice(M, ld, ncol, nf, t, cop, 1, desc, 1)      ! P2.6.1: operator segments
       return
    endif
    call k_vslice_R(M, ld, ncol, Ndim, nf, t, (cop == 'c' .or. cop == 'C'), desc, gpu_vb)
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vslice_R

  Subroutine k_vslice_R(M, ld, ncol, nd, nf, t, cj, desc, vb)
    Implicit none
    Integer, intent(in) :: ld, ncol, nd, nf, t, vb
    Logical, intent(in) :: cj, desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,ncol)
    Integer :: jb, nblk, j, j1, j2, r, pa, pb, p, q, n, i, sp, nact
    Complex (Kind=Kind(0.d0)) :: x, z, zq(4)
#if defined(ALF_LA_HIP)
    if (hk_enabled()) then
       call hk_vslice_r(nd, ncol, ld, nf, t, cj, desc, size(vr_ptr,1), size(V_meta,2), size(nsigma_d,1), vr_ptr, vr_ent, &
            & V_meta, V_type, nsigma_d, V_exp, M)
       return
    endif
#endif

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
    if (v_hasnd(nf)) then
       call nd_slice(M, ld, nrow, nf, t, cop, sign, desc, 2)   ! P2.6.1: operator segments
       return
    endif
    call k_vslice_L(M, ld, nrow, Ndim, nf, t, (cop == 'c' .or. cop == 'C'), sign, desc, gpu_vb)
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vslice_L

  Subroutine k_vslice_L(M, ld, nrow, nd, nf, t, cj, sign, desc, vb)
    Implicit none
    Integer, intent(in) :: ld, nrow, nd, nf, t, sign, vb
    Logical, intent(in) :: cj, desc
    Complex (Kind=Kind(0.d0)), intent(inout) :: M(ld,Ndim)
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
    Integer :: q, n
    if (v_hasnd(nf)) then
       ! P2.6.1: the loop of Op_Wrapup / Op_Wrapdo (N_Type 1) operator by operator
       do q = 1, size(V_meta,2)
          n = q
          if (dir < 0) n = size(V_meta,2) + 1 - q
          call gpu_Vop_wrap(G, ld, n, nf, t, dir, 1)
       enddo
       return
    endif
    call k_vwrap_all(G, ld, Ndim, nf, t, dir)
    gpu_n_launch = gpu_n_launch + 1
  end subroutine gpu_Vwrap_all

  Subroutine k_vwrap_all(G, ld, nd, nf, t, dir)
    Implicit none
    Integer, intent(in) :: ld, nd, nf, t, dir
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,nd)
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
    Integer :: n
    call gpu_mmthr(G, ld, Ndim, nf, ntau + 1)
    call gpu_mmthl_m1(G, ld, Ndim, nf, ntau + 1)
    if (v_hasnd(nf)) then            ! P2.6.1: Op_Wrapup N_Type 1, then 2, operator by operator
       do n = 1, size(V_meta,2)
          call gpu_Vop_wrap(G, ld, n, nf, ntau + 1, 1, 1)
          call gpu_Vop_wrap(G, ld, n, nf, ntau + 1, 1, 2)
       enddo
       return
    endif
    call gpu_Vwrap_all(G, ld, nf, ntau + 1, 1)
  end subroutine gpu_wrapgr_up

!> WRAPGRDO without updates: G = B(ntau)^-1 G B(ntau)  (Op_Wrapdo N_Type 1 for n = N_op..1, then mmthl, mmthr_m1)
  Subroutine gpu_wrapgr_do(G, ld, nf, ntau)
    Implicit none
    Integer, intent(in) :: ld, nf, ntau
    Complex (Kind=Kind(0.d0)), intent(inout) :: G(ld,*)
    Integer :: n
    if (v_hasnd(nf)) then            ! P2.6.1: Op_Wrapdo N_Type 2, then 1, operator by operator
       do n = size(V_meta,2), 1, -1
          call gpu_Vop_wrap(G, ld, n, nf, ntau, -1, 2)
          call gpu_Vop_wrap(G, ld, n, nf, ntau, -1, 1)
       enddo
    else
       call gpu_Vwrap_all(G, ld, nf, ntau, -1)
    endif
    call gpu_mmthl(G, ld, Ndim, nf, ntau)
    call gpu_mmthr_m1(G, ld, Ndim, nf, ntau)
  end subroutine gpu_wrapgr_do

!====================================================================
!  Hook and self-test
!====================================================================
!> Called from Setup_gpu (qmc_setup_mod) after gpu_state_setup when the GPU path is active: reads the kernel tuning
!> variables (ALF_GPU_TL, _TR, _VB, _HOPD) and runs the self-test if ALF_GPU_SELFTEST=2 (the run then continues).
  Subroutine gpu_wrap_hook()
    Implicit none
    Character (len=64) :: f
    call gpu_env_int('ALF_GPU_TL', gpu_tl)
    call gpu_env_int('ALF_GPU_TR', gpu_tr)
    call gpu_env_int('ALF_GPU_VB', gpu_vb)
    hopd_opt = gpu_env_is('ALF_GPU_HOPD', '1')
    if (hopd_opt) call info_line('GPU wrap: dense hopping product Hd ON (ALF_GPU_HOPD=1; Precision Green worse, KNOWN_ISSUES G2)')
    if (gpu_env_get('ALF_GPU_FAULT', f) .and. gpu_env_is('ALF_GPU_SELFTEST', '2')) call inject_fault(trim(f))
    if (gpu_env_is('ALF_GPU_SELFTEST', '2')) then
       call gpu_wrap_selftest()
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
       hop_Rval(1) = hop_Rval(1)*(1.d0 + 1.d-6)
       !$omp target update to(hop_Rval)
    case ('hopc')
       hop_Cval(1) = hop_Cval(1) + cmplx(1.d-6, 0.d0, kind(0.d0))
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
    case (1); call Hop_mod_mmthr   (M, nf, st_t)
    case (2); call Hop_mod_mmthr_m1(M, nf, st_t)
    case (3); call Hop_mod_mmthl   (M, nf, st_t)
    case (4); call Hop_mod_mmthlc  (M, nf, st_t)
    case (5); call Hop_mod_mmthl_m1(M, nf, st_t)
    end select
  end subroutine cpu_hop

  Subroutine dev_hop(op, nf, D, ld, nn)
    Implicit none
    Integer, intent(in) :: op, nf, ld, nn
    Complex (Kind=Kind(0.d0)), intent(inout) :: D(ld,*)
    select case (op)
    case (1); call gpu_mmthr   (D, ld, nn, nf, st_t)
    case (2); call gpu_mmthr_m1(D, ld, nn, nf, st_t)
    case (3); call gpu_mmthl   (D, ld, nn, nf, st_t)
    case (4); call gpu_mmthlc  (D, ld, nn, nf, st_t)
    case (5); call gpu_mmthl_m1(D, ld, nn, nf, st_t)
    end select
  end subroutine dev_hop

  !> Time-dependent hopping elements against the host routines (Hop_mod_mm*, Hop_mod_Symm) at every slice where one is active,
  !> plus slice 1 and Ltrot (t-1 wrap of the Symm average), all five hop routines and Symm (t1 = t2 and t1 /= t2).
  Subroutine td_selftest(s, tol)
    Implicit none
    Integer(int64), intent(inout) :: s
    Real (Kind=Kind(0.d0)), intent(in) :: tol
    Integer :: nf, t, tt(3), it, op, j, ntest, nbad, t2
    Real (Kind=Kind(0.d0)) :: dmax(5), dsy, d, dctl
    Complex (Kind=Kind(0.d0)), allocatable :: A0(:,:), A1(:,:), Dt(:,:), S0(:,:,:), S1(:,:,:)
    Character (len=10), parameter :: opname(5) = (/ 'mmthr     ', 'mmthr_m1  ', 'mmthl     ', 'mmthlc    ', 'mmthl_m1  ' /)
    Character (len=256) :: line
    dmax = 0.d0; dsy = 0.d0; ntest = 0
    allocate(A0(Ndim,Ndim), A1(Ndim,Ndim), Dt(Ndim,Ndim), S0(Ndim,Ndim,N_FL), S1(Ndim,Ndim,N_FL))
    !$omp target enter data map(alloc: Dt)
    do t = 1, Ltrot
       ! slices where some element is active (propagation or Symm kind) plus the boundaries
       if (t /= 1 .and. t /= Ltrot) then
          if (.not. any(td_grp(t,:,:,:) > 0)) cycle
       endif
       do nf = 1, N_FL
          do op = 1, 5
             st_t = t
             call rand_mat(A0, s)
             A1 = A0
             call cpu_hop(op, nf, A1)
             Dt = A0
             !$omp target update to(Dt)
             call dev_hop(op, nf, Dt, Ndim, Ndim)
             !$omp target update from(Dt)
             d = frob_rel(Dt, A1)
             if (d > dmax(op)) dmax(op) = d
             ntest = ntest + 1
          enddo
          ! Symm: t1 = t2 = t, and (t, Ltrot) / (t, 0 -> Ltrot) combinations of the tau-displaced calls
          do it = 1, 3
             t2 = t
             if (it == 2) t2 = Ltrot
             if (it == 3) t2 = 1
             call rand_mat(A0, s)
             S0(:,:,1) = A0
             do j = 2, N_FL
                S0(:,:,j) = A0
             enddo
             call Hop_mod_Symm(S1, S0, t, t2)
             Dt = A0
             !$omp target update to(Dt)
             call gpu_symm(Dt, Ndim, Ndim, nf, t, t2)
             !$omp target update from(Dt)
             d = frob_rel(Dt, S1(:,:,nf))
             if (d > dsy) dsy = d
             ntest = ntest + 1
          enddo
       enddo
    enddo
    ! negative control: the device at a slice where nothing is active against the host at an active slice must differ clearly
    dctl = 0.d0
    do t = 2, Ltrot - 1
       if (all(td_grp(t,:,:,1) == 0)) then
          tt(1) = t; exit
       endif
    enddo
    do t = 1, Ltrot
       if (any(td_grp(t,1,:,1) > 0)) then
          tt(2) = t; exit
       endif
    enddo
    call rand_mat(A0, s)
    A1 = A0
    st_t = tt(2); call cpu_hop(1, 1, A1)
    Dt = A0
    !$omp target update to(Dt)
    st_t = tt(1); call dev_hop(1, 1, Dt, Ndim, Ndim)
    !$omp target update from(Dt)
    dctl = frob_rel(Dt, A1)
    st_t = 1
    !$omp target exit data map(delete: Dt)
    deallocate(A0, A1, Dt, S0, S1)
    nbad = 0
    write(line,'(A,ES10.3,A)') 'GPU wrap self-test time-dependent hopping negative control (device at an inactive slice vs host at an active one): diff ', dctl, ' (must exceed tol)'
    call say(line)
    if (.not. dctl > 10.d0*tol) call selftest_fail('time-dependent hopping: the self-test is blind to the slice')
    do op = 1, 5
       write(line,'(A,A,A,ES10.3,A,ES8.1)') 'GPU wrap self-test time-dependent hopping ', opname(op), ': max rel Frobenius diff ', dmax(op), ' (tol ', tol
       call say(line)
       if (.not. dmax(op) <= tol) nbad = nbad + 1
    enddo
    write(line,'(A,ES10.3,A,I0,A)') 'GPU wrap self-test time-dependent hopping Symm: max rel Frobenius diff ', dsy, ' (', ntest, ' cases)'
    call say(line)
    if (.not. dsy <= tol) nbad = nbad + 1
    if (nbad > 0) call selftest_fail('time-dependent hopping')
  end subroutine td_selftest

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
                   gpu_dense_left_mode = 1
                else
                   gpu_dense_left_mode = 0
                endif
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
    gpu_dense_left_mode = 0
    do op = 1, 5
       write(line,'(A,A,A,ES10.3,A,ES10.3,A,ES8.1)') 'GPU wrap self-test hopping ', opname(op), &
            & ': max rel Frobenius diff (dense-left mode 1) ', dh(op,1), ', (dense-left mode auto) ', dh(op,2), &
            & ' (tol ', tol
       call say(line)
       if (.not. (dh(op,1) <= tol .and. dh(op,2) <= tol)) call selftest_fail('hopping '//trim(opname(op)))
    enddo
    if (td_any) call td_selftest(s, tol)

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

    write(line,'(A,I0,A,I0)') 'GPU wrap self-test device launches issued ', gpu_n_launch, ', linalg_mod calls ', gpu_n_lib
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
       l0 = gpu_n_launch; call gpu_mmthr(Dt, Ndim, Ndim, 1, 1)
       write(line,'(A,I0,A)') 'GPU wrap structure: launches per B_T (mmthr): ', gpu_n_launch - l0, ' (= number of groups)'
       call say(line)
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
          nt = min(V_type(n), 2)               ! type 4: real part +-1, +-2 (Phi_st of type 2), imaginary part below
          old = nsigma_d(n,t0)
          do sp = -nt, nt
             if (sp == 0) cycle
             zsp = cmplx(real(sp,kind(0.d0)), 0.d0, kind(0.d0))
             if (V_type(n) == 4) zsp = cmplx(real(sp,kind(0.d0)), 0.3d0*real(sp,kind(0.d0)), kind(0.d0))
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
    if (any(v_hasnd)) return          ! P2.6.1: synthetic diagonal operators only; non-table flavours use the segment kernels
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
       do ic = 1, 2
          if (ic == 1) then
             gpu_dense_left_mode = 1
          else
             gpu_dense_left_mode = 0
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
    gpu_dense_left_mode = 0
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
       call Op_Wrapup(G, Op_V(n,nf), nsigma%f(n,ntau+1), Ndim, 2, ntau+1)   ! no-op for diagonal operators
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
End Module gpu_wrap_mod
#endif

! ALF_LA_HIP (GPU=amd): HIP implementations of small device kernels (alf_hip_kernels.hip.cpp) behind plain subroutines.
! libgomp-GCN launches an OpenMP target region with a parallel loop in ~150 us, a HIP kernel in ~15 us. The arguments are the host
! addresses of device-resident (mapped) arrays; the device addresses come from use_device_addr, as in alf_la_hip_mod.
#if defined(ALF_LA_HIP)
module alf_hipk_mod
  use iso_c_binding
  use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  implicit none
  private
  public :: hk_copyz, hk_blk_gather, hk_blk_scale_rows, hk_blk_scale_cols, hk_flush_gather, hk_scatter_rows
  public :: hk_corr, hk_scalars, hk_local
  public :: hk_vslice_r
  public :: hk_divrows, hk_rowdiv, hk_scale_cols, hk_tpup, hk_ns, hk_cgr_phase, hk_prec
  public :: hk_ident, hk_setd, hk_setl, hk_diagabs, hk_phase, hk_det_phase, hk_scale_row1, hk_scale_col1, hk_triu, hk_perm, hk_perm_vec, hk_ct, hk_colinv, hk_testscale, hk_colnorms, hk_permute_cols

  interface
    integer(c_int) function alfhk_copyz(n, a, b) bind(C, name='alfhk_copyz')
      import :: c_ptr, c_int, c_long
      integer(c_long), value :: n
      type(c_ptr), value :: a, b
    end function
    integer(c_int) function alfhk_blk_gather(nfl, nfa, nfb, b, nd, ld, g, bsd, gbd) bind(C, name='alfhk_blk_gather')
      import :: c_ptr, c_int
      integer(c_int), value :: nfl, nfa, nfb, b, nd, ld
      type(c_ptr), value :: g, bsd, gbd
    end function
    integer(c_int) function alfhk_blk_scale_rows(nf, fl, ld, b, nd, cbld, g, bsd, cbuf) bind(C, name='alfhk_blk_scale_rows')
      import :: c_ptr, c_int
      integer(c_int), value :: nf, fl, ld, b, nd, cbld
      type(c_ptr), value :: g, bsd, cbuf
    end function
    integer(c_int) function alfhk_blk_scale_cols(nf, fl, ld, b, nd, cbld, g, bsd, cbuf) bind(C, name='alfhk_blk_scale_cols')
      import :: c_ptr, c_int
      integer(c_int), value :: nf, fl, ld, b, nd, cbld
      type(c_ptr), value :: g, bsd, cbuf
    end function
    integer(c_int) function alfhk_flush_gather(ld, nrow, k, nf, nsl, g, s, c, qt, slot) bind(C, name='alfhk_flush_gather')
      import :: c_ptr, c_int
      integer(c_int), value :: ld, nrow, k, nf, nsl
      type(c_ptr), value :: g, s, c, qt, slot
    end function
    integer(c_int) function alfhk_scatter_rows(nd, ncol, p, m, ld, w, ldw) bind(C, name='alfhk_scatter_rows')
      import :: c_ptr, c_int
      integer(c_int), value :: nd, ncol, ld, ldw
      type(c_ptr), value :: p, m, w
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
    integer(c_int) function alfhk_local(n, norb, lnt, it, dec, wt, wr, wi, ld, inv, g, lacc1, lacc2) bind(C, name='alfhk_local')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: n, norb, lnt, it, dec, ld
      real(c_double), value :: wt, wr, wi
      type(c_ptr), value :: inv, g, lacc1, lacc2
    end function
    integer(c_int) function alfhk_divrows(n, ld, u, d) bind(C, name='alfhk_divrows')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: u, d
    end function
    integer(c_int) function alfhk_rowdiv(n, ld, u, d) bind(C, name='alfhk_rowdiv')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: u, d
    end function
    integer(c_int) function alfhk_scale_cols(n, ld, u, d) bind(C, name='alfhk_scale_cols')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: u, d
    end function
    integer(c_int) function alfhk_tpup(n, t, r, dr, dl) bind(C, name='alfhk_tpup')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: t, r, dr, dl
    end function
    integer(c_int) function alfhk_ns(n, a) bind(C, name='alfhk_ns')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: a
    end function
    integer(c_int) function alfhk_cgr_phase(n, ld, t, tau, cs, nvar, sgn, bad) bind(C, name='alfhk_cgr_phase')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld, nvar, sgn, bad
      type(c_ptr), value :: t, tau, cs
    end function
    integer(c_int) function alfhk_prec(n, a, b, r) bind(C, name='alfhk_prec')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: a, b, r
    end function
    integer(c_int) function alfhk_vslice_r(nd, ncol, ld, nf, t, cj, desc, ldptr, nop, vp, ve, vm, vt, ns, vx, m) &
         bind(C, name='alfhk_vslice_r')
      import :: c_ptr, c_int
      integer(c_int), value :: nd, ncol, ld, nf, t, cj, desc, ldptr, nop
      type(c_ptr), value :: vp, ve, vm, vt, ns, vx, m
    end function
    integer(c_int) function alfhk_ident(n, ld, U) bind(C, name='alfhk_ident')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U
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
    integer(c_int) function alfhk_diagabs(n, ld, U, D) bind(C, name='alfhk_diagabs')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, D
    end function
    integer(c_int) function alfhk_phase(n, ld, U, sgn, cj, cs) bind(C, name='alfhk_phase')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld, sgn, cj
      type(c_ptr), value :: U, cs
    end function
    integer(c_int) function alfhk_det_phase(n, ld, A, ipiv, cs) bind(C, name='alfhk_det_phase')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: A, ipiv, cs
    end function
    integer(c_int) function alfhk_scale_row1(n, ld, U, cs) bind(C, name='alfhk_scale_row1')
      import :: c_ptr, c_int
      integer(c_int), value :: n, ld
      type(c_ptr), value :: U, cs
    end function
    integer(c_int) function alfhk_scale_col1(n, U, cs) bind(C, name='alfhk_scale_col1')
      import :: c_ptr, c_int
      integer(c_int), value :: n
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
    integer(c_int) function alfhk_colinv(n, lda, A, D) bind(C, name='alfhk_colinv')
      import :: c_ptr, c_int
      integer(c_int), value :: n, lda
      type(c_ptr), value :: A, D
    end function
    integer(c_int) function alfhk_testscale(n, D, r) bind(C, name='alfhk_testscale')
      import :: c_ptr, c_int
      integer(c_int), value :: n
      type(c_ptr), value :: D, r
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

  subroutine hk_copyz(n, a, b)
    integer, intent(in) :: n
    complex(8), target :: a(*), b(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(a, b)
    rc = alfhk_copyz(int(n,c_long), c_loc(a(1)), c_loc(b(1)))
    !$omp end target data
    call chk(rc, 'copyz')
  end subroutine hk_copyz

  subroutine hk_blk_gather(nfl, nfa, nfb, b, nd, ld, g, bsd, gbd)
    integer, intent(in) :: nfl, nfa, nfb, b, nd, ld
    complex(8), target :: g(*), gbd(*)
    integer, target :: bsd(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(g, bsd, gbd)
    rc = alfhk_blk_gather(int(nfl,c_int), int(nfa,c_int), int(nfb,c_int), int(b,c_int), int(nd,c_int), int(ld,c_int), &
         c_loc(g(1)), c_loc(bsd(1)), c_loc(gbd(1)))
    !$omp end target data
    call chk(rc, 'blk_gather')
  end subroutine hk_blk_gather

  subroutine hk_blk_scale_rows(nf, fl, ld, b, nd, cbld, g, bsd, cbuf)
    integer, intent(in) :: nf, fl, ld, b, nd, cbld
    complex(8), target :: g(*), cbuf(*)
    integer, target :: bsd(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(g, bsd, cbuf)
    rc = alfhk_blk_scale_rows(int(nf,c_int), int(fl,c_int), int(ld,c_int), int(b,c_int), int(nd,c_int), int(cbld,c_int), &
         c_loc(g(1)), c_loc(bsd(1)), c_loc(cbuf(1)))
    !$omp end target data
    call chk(rc, 'blk_scale_rows')
  end subroutine hk_blk_scale_rows

  subroutine hk_blk_scale_cols(nf, fl, ld, b, nd, cbld, g, bsd, cbuf)
    integer, intent(in) :: nf, fl, ld, b, nd, cbld
    complex(8), target :: g(*), cbuf(*)
    integer, target :: bsd(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(g, bsd, cbuf)
    rc = alfhk_blk_scale_cols(int(nf,c_int), int(fl,c_int), int(ld,c_int), int(b,c_int), int(nd,c_int), int(cbld,c_int), &
         c_loc(g(1)), c_loc(bsd(1)), c_loc(cbuf(1)))
    !$omp end target data
    call chk(rc, 'blk_scale_cols')
  end subroutine hk_blk_scale_cols

  subroutine hk_flush_gather(ld, nrow, k, nf, nsl, g, s, c, qt, slot)
    integer, intent(in) :: ld, nrow, k, nf, nsl
    complex(8), target :: g(*), c(*), qt(*)
    integer, target :: s(*), slot(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(g, s, c, qt, slot)
    rc = alfhk_flush_gather(int(ld,c_int), int(nrow,c_int), int(k,c_int), int(nf,c_int), int(nsl,c_int), &
         c_loc(g(1)), c_loc(s(1)), c_loc(c(1)), c_loc(qt(1)), c_loc(slot(1)))
    !$omp end target data
    call chk(rc, 'flush_gather')
  end subroutine hk_flush_gather

  subroutine hk_scatter_rows(nd, ncol, p, m, ld, w, ldw)
    integer, intent(in) :: nd, ncol, ld, ldw
    integer, target :: p(*)
    complex(8), target :: m(*), w(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(p, m, w)
    rc = alfhk_scatter_rows(int(nd,c_int), int(ncol,c_int), c_loc(p(1)), c_loc(m(1)), int(ld,c_int), c_loc(w(1)), int(ldw,c_int))
    !$omp end target data
    call chk(rc, 'scatter_rows')
  end subroutine hk_scatter_rows

  subroutine hk_corr(kind_k, k, n, norb, nfl, ns, sq3, wr, wi, ld, g, jc, inv, acc, acc0)
    integer, intent(in) :: kind_k, k, n, norb, nfl, ld
    real(8), intent(in) :: ns, sq3, wr, wi
    complex(8), target :: g(*), acc(*), acc0(*)
    integer, target :: jc(*), inv(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(g, jc, inv, acc, acc0)
    rc = alfhk_corr(int(kind_k,c_int), int(k,c_int), int(n,c_int), int(norb,c_int), int(nfl,c_int), ns, sq3, wr, wi, &
         int(ld,c_int), c_loc(g(1)), c_loc(jc(1)), c_loc(inv(1)), c_loc(acc(1)), c_loc(acc0(1)))
    !$omp end target data
    call chk(rc, 'corr')
  end subroutine hk_corr

  subroutine hk_scalars(nk, hi, hj, hf, hc, n, norb, nfl, dec, ld, g, inv, uo, hres)
    integer, intent(in) :: nk, n, norb, nfl, dec, ld
    integer, target :: hi(*), hj(*), hf(*), inv(*)
    complex(8), target :: hc(*), g(*)
    real(8), target :: uo(*)
    real(8), intent(out) :: hres(6)
    integer(c_int) :: rc
    !$omp target data use_device_addr(hi, hj, hf, hc, g, inv, uo)
    rc = alfhk_scalars(int(nk,c_int), c_loc(hi(1)), c_loc(hj(1)), c_loc(hf(1)), c_loc(hc(1)), int(n,c_int), int(norb,c_int), &
         int(nfl,c_int), int(dec,c_int), int(ld,c_int), c_loc(g(1)), c_loc(inv(1)), c_loc(uo(1)), hres)
    !$omp end target data
    call chk(rc, 'scalars')
  end subroutine hk_scalars

  subroutine hk_local(n, norb, lnt, it, dec, wt, wr, wi, ld, inv, g, lacc1, lacc2)
    integer, intent(in) :: n, norb, lnt, it, dec, ld
    real(8), intent(in) :: wt, wr, wi
    integer, target :: inv(*)
    complex(8), target :: g(*), lacc1(*), lacc2(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(inv, g, lacc1, lacc2)
    rc = alfhk_local(int(n,c_int), int(norb,c_int), int(lnt,c_int), int(it,c_int), int(dec,c_int), wt, wr, wi, int(ld,c_int), &
         c_loc(inv(1)), c_loc(g(1)), c_loc(lacc1(1)), c_loc(lacc2(1)))
    !$omp end target data
    call chk(rc, 'local')
  end subroutine hk_local

  subroutine hk_divrows(n, ld, u, d)
    integer, intent(in) :: n, ld
    complex(8), target :: u(*), d(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(u, d)
    rc = alfhk_divrows(int(n,c_int), int(ld,c_int), c_loc(u(1)), c_loc(d(1)))
    !$omp end target data
    call chk(rc, 'divrows')
  end subroutine hk_divrows

  subroutine hk_rowdiv(n, ld, u, d)
    integer, intent(in) :: n, ld
    complex(8), target :: u(*), d(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(u, d)
    rc = alfhk_rowdiv(int(n,c_int), int(ld,c_int), c_loc(u(1)), c_loc(d(1)))
    !$omp end target data
    call chk(rc, 'rowdiv')
  end subroutine hk_rowdiv

  subroutine hk_scale_cols(n, ld, u, d)
    integer, intent(in) :: n, ld
    complex(8), target :: u(*), d(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(u, d)
    rc = alfhk_scale_cols(int(n,c_int), int(ld,c_int), c_loc(u(1)), c_loc(d(1)))
    !$omp end target data
    call chk(rc, 'scale_cols')
  end subroutine hk_scale_cols

  subroutine hk_tpup(n, t, r, dr, dl)
    integer, intent(in) :: n
    complex(8), target :: t(*), r(*), dr(*), dl(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(t, r, dr, dl)
    rc = alfhk_tpup(int(n,c_int), c_loc(t(1)), c_loc(r(1)), c_loc(dr(1)), c_loc(dl(1)))
    !$omp end target data
    call chk(rc, 'tpup')
  end subroutine hk_tpup

  subroutine hk_ns(n, a)
    integer, intent(in) :: n
    complex(8), target :: a(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(a)
    rc = alfhk_ns(int(n,c_int), c_loc(a(1)))
    !$omp end target data
    call chk(rc, 'ns')
  end subroutine hk_ns

  subroutine hk_cgr_phase(n, ld, t, tau, cs, nvar, sgn, bad)
    integer, intent(in) :: n, ld, nvar, sgn, bad
    complex(8), target :: t(*), tau(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(t, tau, cs)
    rc = alfhk_cgr_phase(int(n,c_int), int(ld,c_int), c_loc(t(1)), c_loc(tau(1)), c_loc(cs(1)), int(nvar,c_int), int(sgn,c_int), int(bad,c_int))
    !$omp end target data
    call chk(rc, 'cgr_phase')
  end subroutine hk_cgr_phase

  subroutine hk_prec(n, a, b, r)
    integer, intent(in) :: n
    complex(8), target :: a(*), b(*)
    real(8), target :: r(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(a, b, r)
    rc = alfhk_prec(int(n,c_int), c_loc(a(1)), c_loc(b(1)), c_loc(r(1)))
    !$omp end target data
    call chk(rc, 'prec')
  end subroutine hk_prec

  subroutine hk_vslice_r(nd, ncol, ld, nf, t, cj, desc, ldptr, nop, vp, ve, vm, vt, ns, vx, m)
    integer, intent(in) :: nd, ncol, ld, nf, t, ldptr, nop
    logical, intent(in) :: cj, desc
    integer, target :: vp(*), ve(*), vm(*), vt(*)
    complex(8), target :: ns(*), vx(*), m(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(vp, ve, vm, vt, ns, vx, m)
    rc = alfhk_vslice_r(int(nd,c_int), int(ncol,c_int), int(ld,c_int), int(nf,c_int), int(t,c_int), merge(1_c_int,0_c_int,cj), &
         merge(1_c_int,0_c_int,desc), int(ldptr,c_int), int(nop,c_int), c_loc(vp(1)), c_loc(ve(1)), c_loc(vm(1)), c_loc(vt(1)), &
         c_loc(ns(1)), c_loc(vx(1)), c_loc(m(1)))
    !$omp end target data
    call chk(rc, 'vslice_r')
  end subroutine hk_vslice_r

  subroutine hk_ident(n, ld, U)
    integer, intent(in) :: n, ld
    complex(8), target :: U(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U)
    rc = alfhk_ident(int(n,c_int), int(ld,c_int), c_loc(U(1)))
    !$omp end target data
    call chk(rc, 'ident')
  end subroutine hk_ident

  subroutine hk_setd(n, D)
    integer, intent(in) :: n
    complex(8), target :: D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(D)
    rc = alfhk_setd(int(n,c_int), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'setd')
  end subroutine hk_setd

  subroutine hk_setl(n, L)
    integer, intent(in) :: n
    real(8), target :: L(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(L)
    rc = alfhk_setl(int(n,c_int), c_loc(L(1)))
    !$omp end target data
    call chk(rc, 'setl')
  end subroutine hk_setl

  subroutine hk_diagabs(n, ld, U, D)
    integer, intent(in) :: n, ld
    complex(8), target :: U(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, D)
    rc = alfhk_diagabs(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'diagabs')
  end subroutine hk_diagabs

  subroutine hk_phase(n, ld, U, sgn, cj, cs)
    integer, intent(in) :: n, ld, sgn
    logical, intent(in) :: cj
    complex(8), target :: U(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, cs)
    rc = alfhk_phase(int(n,c_int), int(ld,c_int), c_loc(U(1)), int(sgn,c_int), merge(1_c_int,0_c_int,cj), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'phase')
  end subroutine hk_phase

  subroutine hk_det_phase(n, ld, A, ipiv, cs)
    integer, intent(in) :: n, ld
    complex(8), target :: A(*), cs(*)
    integer, target :: ipiv(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, ipiv, cs)
    rc = alfhk_det_phase(int(n,c_int), int(ld,c_int), c_loc(A(1)), c_loc(ipiv(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'det_phase')
  end subroutine hk_det_phase

  subroutine hk_scale_row1(n, ld, U, cs)
    integer, intent(in) :: n, ld
    complex(8), target :: U(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, cs)
    rc = alfhk_scale_row1(int(n,c_int), int(ld,c_int), c_loc(U(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'scale_row1')
  end subroutine hk_scale_row1

  subroutine hk_scale_col1(n, U, cs)
    integer, intent(in) :: n
    complex(8), target :: U(*), cs(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, cs)
    rc = alfhk_scale_col1(int(n,c_int), c_loc(U(1)), c_loc(cs(1)))
    !$omp end target data
    call chk(rc, 'scale_col1')
  end subroutine hk_scale_col1

  subroutine hk_triu(n, ld, U, ldr, R)
    integer, intent(in) :: n, ld, ldr
    complex(8), target :: U(*), R(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(U, R)
    rc = alfhk_triu(int(n,c_int), int(ld,c_int), c_loc(U(1)), int(ldr,c_int), c_loc(R(1)))
    !$omp end target data
    call chk(rc, 'triu')
  end subroutine hk_triu

  subroutine hk_perm(mode, n, X, ldx, perm, W, ldw)
    integer, intent(in) :: mode, n, ldx, ldw
    complex(8), target :: X(*), W(*)
    integer, target :: perm(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(X, perm, W)
    rc = alfhk_perm(int(mode,c_int), int(n,c_int), c_loc(X(1)), int(ldx,c_int), c_loc(perm(1)), c_loc(W(1)), int(ldw,c_int))
    !$omp end target data
    call chk(rc, 'perm')
  end subroutine hk_perm

  subroutine hk_perm_vec(n, L, perm, W)
    integer, intent(in) :: n
    real(8), target :: L(*), W(*)
    integer, target :: perm(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(L, perm, W)
    rc = alfhk_perm_vec(int(n,c_int), c_loc(L(1)), c_loc(perm(1)), c_loc(W(1)))
    !$omp end target data
    call chk(rc, 'perm_vec')
  end subroutine hk_perm_vec

  subroutine hk_ct(n, A, lda, B, ldb)
    integer, intent(in) :: n, lda, ldb
    complex(8), target :: A(*), B(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, B)
    rc = alfhk_ct(int(n,c_int), c_loc(A(1)), int(lda,c_int), c_loc(B(1)), int(ldb,c_int))
    !$omp end target data
    call chk(rc, 'ct')
  end subroutine hk_ct

  subroutine hk_colinv(n, lda, A, D)
    integer, intent(in) :: n, lda
    complex(8), target :: A(*), D(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, D)
    rc = alfhk_colinv(int(n,c_int), int(lda,c_int), c_loc(A(1)), c_loc(D(1)))
    !$omp end target data
    call chk(rc, 'colinv')
  end subroutine hk_colinv

  subroutine hk_testscale(n, D, r)
    integer, intent(in) :: n
    complex(8), target :: D(*)
    real(8), target :: r(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(D, r)
    rc = alfhk_testscale(int(n,c_int), c_loc(D(1)), c_loc(r(1)))
    !$omp end target data
    call chk(rc, 'testscale')
  end subroutine hk_testscale

  subroutine hk_colnorms(m, n, A, lda, nrm)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(*)
    real(8), target :: nrm(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, nrm)
    rc = alfhk_colnorms(int(m,c_int), int(n,c_int), c_loc(A(1)), int(lda,c_int), c_loc(nrm(1)))
    !$omp end target data
    call chk(rc, 'colnorms')
  end subroutine hk_colnorms

  subroutine hk_permute_cols(m, n, A, lda, ipvt)
    integer, intent(in) :: m, n, lda
    complex(8), target :: A(*)
    integer, target :: ipvt(*)
    integer(c_int) :: rc
    !$omp target data use_device_addr(A, ipvt)
    rc = alfhk_permute_cols(int(m,c_int), int(n,c_int), c_loc(A(1)), int(lda,c_int), c_loc(ipvt(1)))
    !$omp end target data
    call chk(rc, 'permute_cols')
  end subroutine hk_permute_cols

end module alf_hipk_mod
#endif

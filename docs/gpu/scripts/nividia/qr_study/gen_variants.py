#!/usr/bin/env python3
"""Ticket 108: insert QR_PREPIVOT=2 (two-pass) and =3 (blocked re-ranking) blocks into
Prog/cgr1_mod.F90 and Prog/cgr2_2_mod.F90 (idempotent: refuses if already inserted).
The existing 104 block (static sort) stays as variant 1 (#else branch)."""
import sys

V2 = r"""#if QR_PREPIVOT == 2
        ! Test-only (ticket 108), variant 2 (two-pass): pass 1 = static norm sort + unpivoted QR of a copy;
        ! the columns are then re-sorted by |R_jj| (residual norm in the pass-1 order), and the final unpivoted
        ! QR of the re-sorted original matrix is done by QDRP_decompose with all columns fixed (cost 2 x geqrf).
        ! IPVT has ZGEQP3 semantics (column j of A*P is column IPVT(j) of A).
        block
          Real(Kind=Kind(0.d0)), allocatable :: pp_nrm(:)
          Real(Kind=Kind(0.d0)) :: DZNRM2
          Complex(Kind=Kind(0.d0)), allocatable :: pp_tmp(:,:), pp_work(:), pp_tau(:)
          Complex(Kind=Kind(0.d0)) :: pp_z
          Integer, allocatable :: pp_perm(:), pp_perm2(:), pp_ord(:)
          Integer :: pp_i, pp_j, pp_k, pp_lw, pp_info
          Allocate(pp_nrm(@NC@), pp_tmp(@M@,@NC@), pp_perm(size(IPVT)), pp_perm2(size(IPVT)), pp_ord(@NC@), pp_tau(@NC@))
          do pp_j = 1, @NC@
             pp_nrm(pp_j) = DZNRM2(@M@, @A@(1,pp_j), 1)
             pp_perm(pp_j) = pp_j
          enddo
          do pp_j = @NC@+1, size(IPVT)
             pp_perm(pp_j) = pp_j
          enddo
          do pp_i = 2, @NC@
             pp_k = pp_perm(pp_i)
             pp_j = pp_i - 1
             do while (pp_j >= 1)
                if (pp_nrm(pp_perm(pp_j)) >= pp_nrm(pp_k)) exit
                pp_perm(pp_j+1) = pp_perm(pp_j)
                pp_j = pp_j - 1
             enddo
             pp_perm(pp_j+1) = pp_k
          enddo
          do pp_j = 1, @NC@
             pp_tmp(:,pp_j) = @A@(1:@M@,pp_perm(pp_j))
          enddo
          ! pass 1: unpivoted QR of the sorted copy (only |R_jj| is used)
          call ZGEQRF(@M@, @NC@, pp_tmp(1,1), @M@, pp_tau(1), pp_z, -1, pp_info)
          pp_lw = INT(DBLE(pp_z))
          Allocate(pp_work(pp_lw))
          call ZGEQRF(@M@, @NC@, pp_tmp(1,1), @M@, pp_tau(1), pp_work(1), pp_lw, pp_info)
          do pp_j = 1, @NC@
             pp_nrm(pp_j) = ABS(pp_tmp(pp_j,pp_j))
          enddo
          ! pass 2 order: stable descending sort of the pass-1 positions by |R_jj|
          do pp_j = 1, @NC@
             pp_ord(pp_j) = pp_j
          enddo
          do pp_i = 2, @NC@
             pp_k = pp_ord(pp_i)
             pp_j = pp_i - 1
             do while (pp_j >= 1)
                if (pp_nrm(pp_ord(pp_j)) >= pp_nrm(pp_k)) exit
                pp_ord(pp_j+1) = pp_ord(pp_j)
                pp_j = pp_j - 1
             enddo
             pp_ord(pp_j+1) = pp_k
          enddo
          pp_perm2 = pp_perm
          do pp_j = 1, @NC@
             pp_perm2(pp_j) = pp_perm(pp_ord(pp_j))
             pp_tmp(:,pp_j) = @A@(1:@M@,pp_perm2(pp_j))
          enddo
          @A@(1:@M@,1:@NC@) = pp_tmp
          IPVT = 1   ! all columns fixed: LAPACK pivoting off
          call QDRP_decompose(@M@, @NC@, @A@, @D@, IPVT, TAU, WORK, LWORK)
          IPVT = pp_perm2
          Deallocate(pp_nrm, pp_tmp, pp_perm, pp_perm2, pp_ord, pp_tau, pp_work)
        end block
#elif QR_PREPIVOT == 3
        ! Test-only (ticket 108), variant 3 (blocked re-ranking): per block of QR_PPNB columns the remaining
        ! columns are re-sorted by their current residual 2-norm, the block is factored by ZGEQRF, Q^H is applied
        ! to the trailing columns by ZUNMQR, and the trailing residual norms are recomputed exactly (DZNRM2,
        ! no downdating). Total flops = one geqrf plus O(M*N^2/QR_PPNB) for the norms. QRCP with block-level
        ! pivoting; QR_PPNB=1 would be ZGEQP3 with exact norms. D is separated as in QDRP_decompose.
        ! IPVT has ZGEQP3 semantics (column j of A*P is column IPVT(j) of A).
#if !defined(QR_PPNB)
#define QR_PPNB 32
#endif
        block
          Real(Kind=Kind(0.d0)), allocatable :: pp_nrm(:), pp_nrm2(:), pp_rwr(:)
          Real(Kind=Kind(0.d0)) :: DZNRM2, pp_x
          Complex(Kind=Kind(0.d0)), allocatable :: pp_tmp(:,:), pp_work(:)
          Complex(Kind=Kind(0.d0)) :: pp_z
          Integer, allocatable :: pp_perm(:), pp_perm2(:), pp_ord(:)
          Integer :: pp_i, pp_j, pp_j0, pp_k, pp_b, pp_r, pp_info, pp_lw, pp_nb
          pp_nb = QR_PPNB
          Allocate(pp_nrm(@NC@), pp_nrm2(@NC@), pp_tmp(@M@,@NC@), pp_perm(size(IPVT)), pp_perm2(size(IPVT)), pp_ord(@NC@))
          do pp_j = 1, @NC@
             pp_nrm(pp_j) = DZNRM2(@M@, @A@(1,pp_j), 1)
          enddo
          do pp_j = 1, size(IPVT)
             pp_perm(pp_j) = pp_j
          enddo
          pp_b = min(pp_nb, @NC@)
          call ZGEQRF(@M@, pp_b, @A@(1,1), @M@, TAU(1), pp_z, -1, pp_info)
          pp_lw = INT(DBLE(pp_z))
          call ZUNMQR('L', 'C', @M@, @NC@, pp_b, @A@(1,1), @M@, TAU(1), @A@(1,1), @M@, pp_z, -1, pp_info)
          pp_lw = max(pp_lw, INT(DBLE(pp_z)))
          Allocate(pp_work(pp_lw))
          pp_k = 1
          do while (pp_k <= @NC@)
             pp_r = @NC@ - pp_k + 1
             pp_b = min(pp_nb, pp_r)
             ! re-rank the remaining columns by current residual norm (stable, descending)
             do pp_i = 1, pp_r
                pp_ord(pp_i) = pp_k + pp_i - 1
             enddo
             do pp_i = 2, pp_r
                pp_j0 = pp_ord(pp_i)
                pp_j = pp_i - 1
                do while (pp_j >= 1)
                   if (pp_nrm(pp_ord(pp_j)) >= pp_nrm(pp_j0)) exit
                   pp_ord(pp_j+1) = pp_ord(pp_j)
                   pp_j = pp_j - 1
                enddo
                pp_ord(pp_j+1) = pp_j0
             enddo
             do pp_i = 1, pp_r
                pp_tmp(1:@M@,pp_i) = @A@(1:@M@,pp_ord(pp_i))
                pp_nrm2(pp_i) = pp_nrm(pp_ord(pp_i))
                pp_perm2(pp_i) = pp_perm(pp_ord(pp_i))
             enddo
             @A@(1:@M@,pp_k:@NC@) = pp_tmp(1:@M@,1:pp_r)
             pp_nrm(pp_k:@NC@) = pp_nrm2(1:pp_r)
             pp_perm(pp_k:@NC@) = pp_perm2(1:pp_r)
             ! factor the panel, apply Q^H to the trailing columns, recompute their residual norms exactly
             call ZGEQRF(@M@-pp_k+1, pp_b, @A@(pp_k,pp_k), @M@, TAU(pp_k), pp_work(1), pp_lw, pp_info)
             if (pp_k + pp_b <= @NC@) then
                call ZUNMQR('L', 'C', @M@-pp_k+1, @NC@-pp_k-pp_b+1, pp_b, @A@(pp_k,pp_k), @M@, TAU(pp_k), &
                     &      @A@(pp_k,pp_k+pp_b), @M@, pp_work(1), pp_lw, pp_info)
                do pp_j = pp_k + pp_b, @NC@
                   pp_nrm(pp_j) = DZNRM2(@M@-pp_k-pp_b+1, @A@(pp_k+pp_b,pp_j), 1)
                enddo
             endif
             pp_k = pp_k + pp_b
          enddo
          ! same LWORK/WORK contract as QDRP_decompose (WORK is used afterwards by ZUNMQR)
          Allocate(pp_rwr(2*@M@))
          call ZGEQP3(@M@, @NC@, @A@(1,1), @M@, IPVT, TAU(1), pp_z, -1, pp_rwr(1), pp_info)
          LWORK = INT(DBLE(pp_z))
          Allocate(WORK(LWORK))
          ! separate off D (as QDRP_decompose)
          do pp_i = 1, @NC@
             pp_x = ABS(@A@(pp_i, pp_i))
             @D@(pp_i) = pp_x
             do pp_j = pp_i, @NC@
                @A@(pp_i, pp_j) = @A@(pp_i, pp_j) / pp_x
             enddo
          enddo
          IPVT = pp_perm
          Deallocate(pp_nrm, pp_nrm2, pp_rwr, pp_tmp, pp_work, pp_perm, pp_perm2, pp_ord)
        end block
#else
"""

SITES = {
 "/Users/luis/Documents/_WORK/Coding/CLAUDE/ALF/Prog/cgr1_mod.F90":
    [dict(A="TPUP", M="N_size", NC="udvl%N_part", D="DUP")],
 "/Users/luis/Documents/_WORK/Coding/CLAUDE/ALF/Prog/cgr2_2_mod.F90":
    [dict(A="HLPB1", M="LQ2", NC="LQ2", D="D3")]*2,
}

for fn, sites in SITES.items():
    lines = open(fn).read().split("\n")
    if any("QR_PREPIVOT == 2" in l for l in lines):
        sys.exit("already inserted: " + fn)
    out = []; si = 0; i = 0; inside = False
    while i < len(lines):
        l = lines[i]
        if l.strip() == "#if defined(QR_PREPIVOT)":
            out.append(l)
            t = V2
            for k, v in sites[si].items():
                t = t.replace("@%s@" % k, v)
            out.extend(t.rstrip("\n").split("\n"))
            si += 1; inside = True
        elif inside and l.strip() == "#else":
            out.append("#endif"); out.append(l); inside = False
        else:
            out.append(l)
        i += 1
    assert si == len(sites)
    open(fn, "w").write("\n".join(out))

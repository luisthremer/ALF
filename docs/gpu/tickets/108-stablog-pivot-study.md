# 108 – STABLOG: a GPU-compatible pivot strategy for CGR/CGR2_2

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 104 (done)
- **Status:** done (2026-10-27; reviewed, key case reproduced; verdict: no variant passes unconditionally, STABLOG keeps host QRCP)

## Goal
Find a pivoting strategy for the QR in CGR and CGR2_2 that keeps STABLOG precision at upstream level (Green and
τ precision: max within 3×, mean within 1.5× of upstream ZGEQP3), and that maps onto GPU library calls (unpivoted
`geqrf` on large panels plus cheap device work). Otherwise document that STABLOG needs a host QRCP for this step,
and what that costs.

## Context
- `docs/gpu/reports/104-prepivot-qr.md` and its "Manager review": a static one-shot sort by initial column 2-norm
  works for the default stabilization but degrades STABLOG. The probable cause is that it has no norm downdating.
- The 104 implementation (`#if defined(QR_PREPIVOT)` in `Prog/cgr1_mod.F90`, `Prog/cgr2_2_mod.F90`) and the study
  scripts (`docs/gpu/scripts/nividia/qr_study/`) are the starting point.
- Candidates, cheapest first. Each is a new macro value, e.g. `-DQR_PREPIVOT=2`:
  1. **Two-pass:** static sort, unpivoted QR, then re-sort the columns by the row norms of R (or |R_ii|) and do a
     second QR of R·P₂ (N×N, upper triangular input). Cost ≈ 2× geqrf.
  2. **Blocked re-ranking:** factor the first block of b columns (b = 32…64) after the static sort, downdate the
     remaining column norms (as ZGEQP3 does, from the block's R), re-sort the rest, and repeat. Equivalent to QRCP
     with block-level instead of column-level pivoting; maps to geqrf on panels plus norm updates.
  3. For reference, true QRCP of the full matrix (= upstream) and its cost from ticket 101 (host ZGEQP3:
     26 ms at N=400, 373 ms at N=1024, one core).
- Regime: unlike 104, use parameters where **upstream** is stable: Nwrap chosen so that upstream Green max ≲ 1e-6
  (U=4: Nwrap=10; U=8: try Nwrap=5). Use L=6, 8 and 12, β = 4…40, **3 seeds per case**, thermalised (discard the first
  2 bins). Checkerboard, plus one dense case.

## Files allowed
- `Prog/cgr1_mod.F90`, `Prog/cgr2_2_mod.F90` (only inside `#if defined(QR_PREPIVOT)` blocks; macro-off code unchanged)
- `docs/gpu/scripts/nividia/qr_study/**`, `docs/gpu/reports/108-stablog-pivot.md`

## nividia rules
Testing only (develop on the laptop and rsync). Write only inside `~/alf-gpu/work/t108/` (TMPDIR inside it). Budget
8 CPU cores, runs ≤ 15 min each, no GPU, no leftovers, install nothing.

## Acceptance
- Unit tests 39/39 with each macro variant. Macro-off runs bitwise identical to upstream.
- The report has a table per variant: precision mean/max vs upstream, abort counts, E and S(π,π) agreement, and
  the extra flops relative to one geqrf. The verdict per variant uses the criteria in the Goal.
- numerics-reviewer approves the diff.

## Out of scope
GPU code, the default-stabilization path (already settled by 104), STAB1/2/3.

## Resolution (manager, 2026-10-27)
- numerics-reviewer approved the diff; the manager reproduced L8 b10 U4 w10 (s3, up and v2) bit for bit.
- Variant 1 and 3 fail; variant 2 (two-pass, 2× geqrf) meets the criteria in all class-A cases and 8/9 class-B cases, with
  one tau-max outlier (4.6×) and systematically higher means (up to 1.41×). Not an unconditional pass.
- Decision: STABLOG on the GPU uses host ZGEQP3 for CGR/CGR2_2 in the MVP (ADR-001 §4). Details and nits:
  `docs/gpu/reports/108-stablog-pivot.md`, "Manager review".

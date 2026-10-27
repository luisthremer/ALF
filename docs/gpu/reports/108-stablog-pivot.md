# 108 - STABLOG pivot strategy for CGR / CGR2_2 (CPU study)

Machine: nividia, GNU 13.3 + OpenBLAS (1 thread), MPI+HDF5 build, **STABLOG** (`configure.sh GNU MPI LOG HDF5`),
1 rank per run, Hubbard square lattice, Mz, Ltau=1, Symm=F, dtau=0.1, random start, 20 bins, n_skip=2 (first 2 bins
discarded for E and S; the `info` precision numbers cover the whole run). Checkerboard except one dense case.
Scripts: `docs/gpu/scripts/nividia/qr_study/` (`gen_variants.py`, `build108.sh`, `ctest108.sh`, `stage108.sh`,
`reseed108.sh`, `run108.sh`, `analyze108.py`, `post108.py`, `test_variants.F90`, `bench_variants.F90`).

## Implementation (only inside `#if defined(QR_PREPIVOT)`, in `Prog/cgr1_mod.F90` and both sites of `Prog/cgr2_2_mod.F90`)
`-DQR_PREPIVOT` or `=1` is the unchanged 104 block (static sort); the new blocks are `#if QR_PREPIVOT == 2 / #elif == 3`
in front of it. Macro-off code is untouched (the diff only adds lines: +456, no removed line). In all variants IPVT has
ZGEQP3 semantics (column j of A*P is column IPVT(j) of A), TAU/WORK/LWORK follow the QDRP_decompose contract, the D
separation is the one of QDRP_decompose, so the permutation-sign loop, the det(Q) phase and `solve_extended_System` are untouched.

- **Variant 2 (two-pass).** Pass 1: static norm sort + ZGEQRF of a *copy* (only |R_jj| is used). Pass 2: columns of the
  original matrix are re-sorted (stable, descending) by |R_jj| (the residual norm of each column in the pass-1 order) and
  QDRP_decompose factors them unpivoted (IPVT=1). The ticket's "second QR of R*P2" cannot be used literally: Q would then be
  a product Q1*Q2 of two reflector sets, but the downstream ZUNMQR/solve code needs one set; so pass 2 refactors A*P1*P2.
  Cost: 2 x geqrf (+ O(n^2) norms and sorts).
- **Variant 3 (blocked re-ranking, block size `QR_PPNB`, default 32).** Chosen implementation: per block (i) re-sort the
  remaining columns by their current residual 2-norm (stable), (ii) ZGEQRF of the panel, (iii) ZUNMQR('L','C') applies Q^H
  to the trailing columns (this is exactly what a blocked geqrf does anyway), (iv) the residual norms of the trailing
  columns are **recomputed exactly** with DZNRM2 on the updated rows. No norm downdating at all, so no cancellation/ZLAQPS
  tolerance logic is needed; the price is O(M N^2 / nb) flops. With nb=1 it is ZGEQP3 with exact norms.
  Flops: 1 x geqrf + (4 x sum of trailing block sizes) = about 1/(4 nb) of a geqrf: 0.8 % (nb=32), 1.6 % (16), 3.1 % (8).
- Standalone check (`test_variants.F90`, N=96, 4 graded/near-collinear random matrices, same code text): all variants
  reconstruct A*P = Q*R*D to 2e-15; variant 3 with nb=1 gives **exactly** ZGEQP3's permutation and |R_jj| (0 differing
  entries), nb=4 already differs from ZGEQP3 in 0..78 permutation entries. So the block re-ranking is a correct algorithm
  but not a QRCP replacement.

## Checks
- Unit tests (devel, noMPI, GNU, -fcheck=all, -Werror): macro off, =1, =2, =3 (nb=32) and =3 with QR_PPNB=4: **39/39 pass each**
  (build flags confirmed to contain the macro).
- Macro-off identical to upstream: build "up" (modified sources, macro off) vs "orig" (pre-104 cgr sources) on
  L6 b10 U4 w10 and L6 b10 U8 w5: all 38 HDF5 datasets bitwise equal, precision lines identical.
- Host timings (`bench_variants.F90`, one core, best of 3, complex N x N, graded columns), ratio to ZGEQP3:

| N | v1 static | v2 two-pass | v3 nb=16 | v3 nb=32 | v3 nb=64 | ZGEQP3 |
|---|---|---|---|---|---|---|
| 400 | 0.64 | 1.29 | 1.29 | 1.21 | 0.71 | 14 ms |
| 1024 | 0.77 | 1.55 | 1.72 | 1.53 | 0.90 | 186 ms |

  Flops relative to one geqrf: v1 1.00, v2 2.00, v3 1.0 + 1/(4 nb). (Ticket 101 quotes 26 ms / 373 ms for ZGEQP3 at
  N=400 / 1024, about 2x my numbers; machine load differed, only the ratios matter here.)

## Study design and what was trimmed
1. **Regime scan** (upstream only, 1 seed, short): U=4 Nwrap=10 has upstream Green max 1e-6 (L6) .. 5e-6 (L8) .. 3e-4 (L12 b20),
   U=4 Nwrap=5 at L12 reaches 3e-8..3e-7, U=8 Nwrap=5 is stable at L6 (2e-7..2e-6) and marginal at L8/L12 (1e-5), U=8 Nwrap=10
   is not stable (L6 b20: 0.13). So: U=4 Nwrap=10 (L<=8) and Nwrap=5 (L12), U=8 Nwrap=5. The regime class used below is
   taken from the measured upstream Green max over the final runs: **A** <= 1.5e-6, **B** <= 1.5e-5, **C** > 1.5e-5.
2. **Main study** (6 builds: upstream, v1, v2, v3 nb=32/16/8): 18 cases (L=6: beta 4,10,20,40 for U4 w10 and for U8 w5; L=8: beta 10,20,40 U4 w10,
   beta 10,20 U8 w5; L=12: beta 10 U4 w10, beta 10 and 20 U4 w5, beta 10 U8 w5; L6 b20 U8 w5 dense), each run sized to
   about 3 min, 20 bins.
3. **Seed bug found and fixed.** pyALF writes the *same* seeds file for every case, so the first "3 seeds" were three
   identical runs (bit-identical info). I rebuilt real seeds (`reseed108.sh`: s2/s3 use the seed on line 11/51 of the
   file) and re-ran. To stay near the time budget the re-run is **trimmed**: 14 of 18 cases (dropped: L6 b40 U4/U8, L6 b20 U8 w5
   checkerboard, L8 b40 U4 w10 -- their s1 results only; L6 b20 U8 w5 and the two b40 L6 rows show "seeds 1"), builds up/v1/v2/v3n32/v3n8
   (v3 nb=16 stays at 1 seed). The column "seeds" says how many independent seeds each row has; up/variant ratios are always
   taken on the same seeds. Total study wall time about 4 h on 8 cores (budget 3 h; the seed re-run cost about 1 h).
4. Precision numbers are `info` values: mean over all measurements and the maximum over the run (an extreme-value
   statistic, sensitive to single rare events; 3 seeds do not make it smooth).

Criteria (ticket Goal): Green **and** tau precision, max within 3x and mean within 1.5x of upstream. "x" columns below are
variant/upstream. Abort counts (ALF "Unstable matrix" / NaN) are 0 for every run of every build, including upstream.
E and S(pi,pi) (maximum of the S_zz(K) result) agree with upstream within 0.05 sigma in every row (identical seeds and
near-identical trajectories, so this is a weak check; the trajectories only change at rounding level unless a precision event occurs).

## Per-variant tables

### Variant 1 (static sort, 104)

| case | seeds | regime | up Green max | Green mean x | Green max x | tau mean x | tau max x | aborts | dE/sigma | dS/sigma | meets Goal |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L12_b10_U4_w10 | 3 | C | 3.66e-05 | 5.92 | 24.23 | 1.12 | 0.72 | 0/3 | +0.0 | +0.0 | NO |
| L12_b10_U4_w5 | 3 | A | 3.48e-08 | 1.03 | 1.43 | 0.91 | 1.00 | 0/3 | -0.0 | +0.0 | yes |
| L12_b10_U8_w5 | 3 | C | 3.72e-05 | 1.16 | 1.25 | 0.84 | 0.85 | 0/3 | +0.0 | +0.0 | yes |
| L12_b20_U4_w5 | 3 | A | 3.12e-07 | 0.99 | 0.80 | 0.86 | 1.57 | 0/3 | +0.0 | -0.0 | yes |
| L6_b10_U4_w10 | 3 | B | 2.79e-06 | 1.94 | 10.16 | 1.28 | 1.64 | 0/3 | +0.0 | +0.0 | NO |
| L6_b10_U8_w5 | 3 | A | 1.33e-06 | 1.00 | 0.89 | 1.14 | 1.00 | 0/3 | +0.0 | -0.0 | yes |
| L6_b20_U4_w10 | 3 | B | 5.86e-06 | 1.11 | 2.23 | 1.20 | 0.99 | 0/3 | -0.0 | +0.0 | yes |
| L6_b20_U8_w5 | 1 | B | 4.90e-06 | 1.00 | 0.99 | 1.09 | 1.00 | 0/1 | +0.0 | +0.0 | yes |
| L6_b20_U8_w5_dense | 3 | B | 4.99e-06 | 1.02 | 1.52 | 1.05 | 0.99 | 0/3 | +0.0 | -0.0 | yes |
| L6_b40_U4_w10 | 1 | A | 8.34e-07 | 1.01 | 1.00 | 1.15 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b40_U8_w5 | 1 | B | 2.49e-06 | 1.00 | 1.00 | 1.06 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b4_U4_w10 | 3 | A | 1.01e-06 | 1.37 | 2.71 | 1.32 | 1.02 | 0/3 | +0.0 | -0.0 | yes |
| L6_b4_U8_w5 | 3 | A | 7.15e-07 | 1.09 | 3.82 | 1.29 | 1.00 | 0/3 | +0.0 | -0.0 | NO |
| L8_b10_U4_w10 | 3 | B | 5.55e-06 | 3.95 | 7.52 | 1.29 | 0.73 | 0/3 | -0.0 | +0.0 | NO |
| L8_b10_U8_w5 | 3 | C | 2.01e-05 | 1.00 | 1.11 | 1.04 | 0.26 | 0/3 | -0.0 | +0.0 | yes |
| L8_b20_U4_w10 | 3 | B | 8.51e-06 | 2.03 | 3.84 | 1.17 | 1.00 | 0/3 | +0.0 | +0.0 | NO |
| L8_b20_U8_w5 | 3 | B | 1.10e-05 | 1.02 | 1.00 | 1.05 | 1.02 | 0/3 | -0.0 | -0.0 | yes |
| L8_b40_U4_w10 | 1 | B | 1.14e-05 | 1.06 | 1.06 | 1.12 | 1.37 | 0/1 | +0.0 | -0.0 | yes |


### Variant 2 (two-pass)

| case | seeds | regime | up Green max | Green mean x | Green max x | tau mean x | tau max x | aborts | dE/sigma | dS/sigma | meets Goal |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L12_b10_U4_w10 | 3 | C | 3.66e-05 | 1.88 | 7.99 | 1.00 | 0.15 | 0/3 | +0.0 | +0.0 | NO |
| L12_b10_U4_w5 | 3 | A | 3.48e-08 | 1.06 | 2.06 | 0.97 | 1.23 | 0/3 | -0.0 | -0.0 | yes |
| L12_b10_U8_w5 | 3 | C | 3.72e-05 | 1.24 | 1.01 | 1.03 | 1.32 | 0/3 | +0.0 | +0.0 | yes |
| L12_b20_U4_w5 | 3 | A | 3.12e-07 | 1.03 | 0.73 | 0.95 | 0.95 | 0/3 | -0.0 | +0.0 | yes |
| L6_b10_U4_w10 | 3 | B | 2.79e-06 | 1.05 | 1.16 | 1.07 | 0.75 | 0/3 | -0.0 | +0.0 | yes |
| L6_b10_U8_w5 | 3 | A | 1.33e-06 | 1.00 | 0.88 | 1.02 | 1.00 | 0/3 | -0.0 | -0.0 | yes |
| L6_b20_U4_w10 | 3 | B | 5.86e-06 | 1.02 | 0.99 | 1.05 | 0.99 | 0/3 | -0.0 | +0.0 | yes |
| L6_b20_U8_w5 | 1 | B | 4.90e-06 | 1.00 | 1.00 | 1.02 | 1.00 | 0/1 | -0.0 | -0.0 | yes |
| L6_b20_U8_w5_dense | 3 | B | 4.99e-06 | 1.01 | 1.58 | 1.00 | 0.99 | 0/3 | +0.0 | +0.0 | yes |
| L6_b40_U4_w10 | 1 | A | 8.34e-07 | 1.00 | 1.00 | 1.05 | 1.00 | 0/1 | +0.0 | +0.0 | yes |
| L6_b40_U8_w5 | 1 | B | 2.49e-06 | 0.99 | 1.00 | 1.02 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b4_U4_w10 | 3 | A | 1.01e-06 | 1.06 | 1.42 | 1.07 | 1.97 | 0/3 | +0.0 | +0.0 | yes |
| L6_b4_U8_w5 | 3 | A | 7.15e-07 | 1.00 | 1.00 | 1.04 | 1.00 | 0/3 | -0.0 | +0.0 | yes |
| L8_b10_U4_w10 | 3 | B | 5.55e-06 | 1.41 | 1.56 | 1.10 | 4.62 | 0/3 | +0.0 | +0.0 | NO |
| L8_b10_U8_w5 | 3 | C | 2.01e-05 | 0.97 | 1.00 | 0.96 | 0.19 | 0/3 | -0.0 | -0.0 | yes |
| L8_b20_U4_w10 | 3 | B | 8.51e-06 | 1.14 | 2.64 | 1.04 | 1.00 | 0/3 | -0.0 | +0.0 | yes |
| L8_b20_U8_w5 | 3 | B | 1.10e-05 | 1.00 | 1.00 | 1.01 | 1.00 | 0/3 | -0.0 | -0.0 | yes |
| L8_b40_U4_w10 | 1 | B | 1.14e-05 | 1.02 | 1.07 | 1.02 | 1.00 | 0/1 | +0.0 | -0.0 | yes |


### Variant 3, nb=32 (blocked re-ranking)

| case | seeds | regime | up Green max | Green mean x | Green max x | tau mean x | tau max x | aborts | dE/sigma | dS/sigma | meets Goal |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L12_b10_U4_w10 | 3 | C | 3.66e-05 | 5.03 | 30.86 | 0.85 | 0.46 | 0/3 | +0.0 | -0.0 | NO |
| L12_b10_U4_w5 | 3 | A | 3.48e-08 | 0.92 | 0.94 | 0.84 | 1.08 | 0/3 | -0.0 | +0.0 | yes |
| L12_b10_U8_w5 | 3 | C | 3.72e-05 | 0.83 | 0.47 | 0.64 | 0.88 | 0/3 | +0.0 | +0.0 | yes |
| L12_b20_U4_w5 | 3 | A | 3.12e-07 | 0.93 | 0.78 | 0.79 | 0.70 | 0/3 | +0.0 | +0.0 | yes |
| L6_b10_U4_w10 | 3 | B | 2.79e-06 | 1.93 | 8.88 | 1.08 | 0.69 | 0/3 | +0.0 | -0.0 | NO |
| L6_b10_U8_w5 | 3 | A | 1.33e-06 | 1.00 | 0.89 | 1.06 | 1.00 | 0/3 | +0.0 | -0.0 | yes |
| L6_b20_U4_w10 | 3 | B | 5.86e-06 | 1.10 | 1.63 | 1.07 | 0.99 | 0/3 | +0.0 | +0.0 | yes |
| L6_b20_U8_w5 | 1 | B | 4.90e-06 | 1.00 | 1.00 | 1.04 | 1.00 | 0/1 | +0.0 | +0.0 | yes |
| L6_b20_U8_w5_dense | 3 | B | 4.99e-06 | 1.03 | 1.74 | 1.02 | 0.99 | 0/3 | +0.0 | -0.0 | yes |
| L6_b40_U4_w10 | 1 | A | 8.34e-07 | 1.00 | 1.00 | 1.07 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b40_U8_w5 | 1 | B | 2.49e-06 | 1.00 | 1.00 | 1.04 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b4_U4_w10 | 3 | A | 1.01e-06 | 1.34 | 1.84 | 1.09 | 1.01 | 0/3 | +0.0 | -0.0 | yes |
| L6_b4_U8_w5 | 3 | A | 7.15e-07 | 1.09 | 3.82 | 1.11 | 1.00 | 0/3 | +0.0 | -0.0 | NO |
| L8_b10_U4_w10 | 3 | B | 5.55e-06 | 4.27 | 12.01 | 1.09 | 0.41 | 0/3 | -0.0 | +0.0 | NO |
| L8_b10_U8_w5 | 3 | C | 2.01e-05 | 1.00 | 1.11 | 0.97 | 0.26 | 0/3 | -0.0 | +0.0 | yes |
| L8_b20_U4_w10 | 3 | B | 8.51e-06 | 2.17 | 10.75 | 1.07 | 1.00 | 0/3 | +0.0 | +0.0 | NO |
| L8_b20_U8_w5 | 3 | B | 1.10e-05 | 1.02 | 1.00 | 1.02 | 1.02 | 0/3 | -0.0 | -0.0 | yes |
| L8_b40_U4_w10 | 1 | B | 1.14e-05 | 1.10 | 1.06 | 1.08 | 1.37 | 0/1 | +0.0 | -0.0 | yes |


### Variant 3, nb=8

| case | seeds | regime | up Green max | Green mean x | Green max x | tau mean x | tau max x | aborts | dE/sigma | dS/sigma | meets Goal |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L12_b10_U4_w10 | 3 | C | 3.66e-05 | 1.32 | 4.10 | 0.84 | 0.57 | 0/3 | -0.0 | +0.0 | NO |
| L12_b10_U4_w5 | 3 | A | 3.48e-08 | 0.89 | 1.06 | 0.83 | 0.94 | 0/3 | -0.0 | +0.0 | yes |
| L12_b10_U8_w5 | 3 | C | 3.72e-05 | 0.79 | 0.47 | 0.64 | 0.88 | 0/3 | +0.0 | +0.0 | yes |
| L12_b20_U4_w5 | 3 | A | 3.12e-07 | 0.91 | 0.78 | 0.79 | 0.79 | 0/3 | +0.0 | +0.0 | yes |
| L6_b10_U4_w10 | 3 | B | 2.79e-06 | 1.98 | 19.82 | 1.03 | 0.69 | 0/3 | +0.0 | +0.0 | NO |
| L6_b10_U8_w5 | 3 | A | 1.33e-06 | 1.00 | 0.88 | 1.01 | 1.00 | 0/3 | +0.0 | -0.0 | yes |
| L6_b20_U4_w10 | 3 | B | 5.86e-06 | 1.13 | 2.77 | 1.02 | 0.99 | 0/3 | +0.0 | -0.0 | yes |
| L6_b20_U8_w5 | 1 | B | 4.90e-06 | 1.00 | 0.99 | 1.00 | 1.00 | 0/1 | +0.0 | +0.0 | yes |
| L6_b20_U8_w5_dense | 3 | B | 4.99e-06 | 1.02 | 2.53 | 0.99 | 0.99 | 0/3 | +0.0 | +0.0 | yes |
| L6_b40_U4_w10 | 1 | A | 8.34e-07 | 1.02 | 1.00 | 1.04 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b40_U8_w5 | 1 | B | 2.49e-06 | 1.01 | 1.00 | 1.01 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b4_U4_w10 | 3 | A | 1.01e-06 | 1.36 | 14.26 | 1.03 | 1.02 | 0/3 | -0.0 | +0.0 | NO |
| L6_b4_U8_w5 | 3 | A | 7.15e-07 | 1.09 | 10.39 | 1.03 | 1.00 | 0/3 | -0.0 | -0.0 | NO |
| L8_b10_U4_w10 | 3 | B | 5.55e-06 | 3.36 | 15.04 | 1.02 | 0.43 | 0/3 | -0.0 | -0.0 | NO |
| L8_b10_U8_w5 | 3 | C | 2.01e-05 | 1.00 | 1.11 | 0.93 | 0.25 | 0/3 | -0.0 | -0.0 | yes |
| L8_b20_U4_w10 | 3 | B | 8.51e-06 | 1.98 | 8.83 | 1.02 | 1.00 | 0/3 | +0.0 | +0.0 | NO |
| L8_b20_U8_w5 | 3 | B | 1.10e-05 | 1.02 | 1.00 | 0.99 | 1.02 | 0/3 | -0.0 | -0.0 | yes |
| L8_b40_U4_w10 | 1 | B | 1.14e-05 | 1.11 | 1.10 | 1.06 | 1.37 | 0/1 | +0.0 | -0.0 | yes |


### Variant 3, nb=16

| case | seeds | regime | up Green max | Green mean x | Green max x | tau mean x | tau max x | aborts | dE/sigma | dS/sigma | meets Goal |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L12_b10_U4_w10 | 1 | C | 3.66e-05 | 3.07 | 6.91 | 0.84 | 0.45 | 0/1 | -0.0 | -0.0 | NO |
| L12_b10_U4_w5 | 1 | A | 3.48e-08 | 0.94 | 1.08 | 0.84 | 0.40 | 0/1 | -0.0 | +0.0 | yes |
| L12_b10_U8_w5 | 1 | C | 3.72e-05 | 0.84 | 0.60 | 0.61 | 0.29 | 0/1 | +0.0 | +0.0 | yes |
| L12_b20_U4_w5 | 1 | A | 3.12e-07 | 0.88 | 0.78 | 0.78 | 0.89 | 0/1 | -0.0 | +0.0 | yes |
| L6_b10_U4_w10 | 1 | B | 2.79e-06 | 1.45 | 3.52 | 1.06 | 1.09 | 0/1 | +0.0 | -0.0 | NO |
| L6_b10_U8_w5 | 1 | A | 1.33e-06 | 1.01 | 1.00 | 1.05 | 1.00 | 0/1 | +0.0 | -0.0 | yes |
| L6_b20_U4_w10 | 1 | B | 5.86e-06 | 1.05 | 1.00 | 1.05 | 1.00 | 0/1 | -0.0 | -0.0 | yes |
| L6_b20_U8_w5 | 1 | B | 4.90e-06 | 1.00 | 1.00 | 1.02 | 1.00 | 0/1 | +0.0 | +0.0 | yes |
| L6_b20_U8_w5_dense | 1 | B | 4.99e-06 | 1.00 | 1.04 | 1.01 | 1.01 | 0/1 | -0.0 | -0.0 | yes |
| L6_b40_U4_w10 | 1 | A | 8.34e-07 | 1.01 | 1.00 | 1.06 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b40_U8_w5 | 1 | B | 2.49e-06 | 1.01 | 1.00 | 1.02 | 1.00 | 0/1 | -0.0 | +0.0 | yes |
| L6_b4_U4_w10 | 1 | A | 1.01e-06 | 1.37 | 3.49 | 1.06 | 0.99 | 0/1 | -0.0 | -0.0 | NO |
| L6_b4_U8_w5 | 1 | A | 7.15e-07 | 1.08 | 2.12 | 1.07 | 2.10 | 0/1 | +0.0 | +0.0 | yes |
| L8_b10_U4_w10 | 1 | B | 5.55e-06 | 4.62 | 11.99 | 1.06 | 1.25 | 0/1 | +0.0 | +0.0 | NO |
| L8_b10_U8_w5 | 1 | C | 2.01e-05 | 1.00 | 1.00 | 1.03 | 1.41 | 0/1 | -0.0 | +0.0 | yes |
| L8_b20_U4_w10 | 1 | B | 8.51e-06 | 1.61 | 6.98 | 1.03 | 0.58 | 0/1 | +0.0 | -0.0 | NO |
| L8_b20_U8_w5 | 1 | B | 1.10e-05 | 1.04 | 1.00 | 1.06 | 1.09 | 0/1 | -0.0 | -0.0 | yes |
| L8_b40_U4_w10 | 1 | B | 1.14e-05 | 1.09 | 1.06 | 1.07 | 1.37 | 0/1 | +0.0 | -0.0 | yes |


Pass counts (full criterion / mean criterion only):

| variant | class A (6 cases) | class B (9) | class C (3) |
|---|---|---|---|
| 1 | 5/6, 6/6 | 6/9, 6/9 | 2/3, 2/3 |
| 2 | **6/6, 6/6** | 8/9, **9/9** | 2/3, 2/3 |
| 3 nb=32 | 5/6, 6/6 | 6/9, 6/9 | 2/3, 2/3 |
| 3 nb=16 (1 seed) | 5/6, 6/6 | 6/9, 7/9 | 2/3, 2/3 |
| 3 nb=8 | 4/6, 6/6 | 6/9, 6/9 | 2/3, 3/3 |

## Verdicts against the Goal
- **Variant 1 (static sort): fails** (as in 104). Green mean up to 3.95x (L8 b10), 5.9x (L12 b10 w10); max up to 10x (L6 b10 U4),
  24x (L12), 7.5x (L8 b10); fails 3 of 9 B cases and the class-A case L6 b4 U8 (Green max 3.8x).
- **Variant 2 (two-pass): closest, but not an unconditional pass.** Mean criterion met in all class A and B cases
  (worst: Green mean 1.41x at L8 b10 U4 w10, 1.14x L8 b20), full criterion met in all 6 class-A cases (worst max ratio
  2.06x, L12 b10 w5) and 8 of 9 class-B cases. Failures: L8 b10 U4 w10 **tau max 4.62x** (one of three seeds: upstream 2.3e-6,
  v2 1.05e-5, a single rare event; Green max there only 1.56x) and in class C (upstream max 3.7e-5, outside the regime) L12 b10 U4 w10
  (Green mean 1.88x, max 8x). Cost 2 x geqrf (host time 1.3-1.55x ZGEQP3), maps to two unpivoted geqrf calls plus a device sort/permute.
- **Variant 3 (blocked re-ranking): fails**, for every block size tried. Mean 1.3-5.0x in class B/C (L8 b10: 4.3x/3.4x for nb=32/8), max
  up to 12-31x (nb=32), 14-20x (nb=8, including class A: L6 b4 U4 14x, L6 b4 U8 10x). Smaller blocks do not help (nb=8 is
  worse than nb=32 on the class-A max); with nb=1 it would reproduce ZGEQP3 exactly but a panel width 1 is not a GPU-friendly geqrf.
  Reading: the within-block order is decided with stale norms, which is exactly where nearly parallel large columns hurt.
- **Fallback documented as the ticket asks:** if even variant 2's single outlier is not acceptable, STABLOG needs a host
  QRCP (ZGEQP3) for the N x N (CGR) and 2N x 2N (CGR2_2) factorization: 14 ms at N=400, 186 ms at N=1024 on one core here
  (ticket 101: 26 / 373 ms), plus transferring the matrix (N^2 complex = 2.6 MB at N=400, 16.8 MB at N=1024) both ways and
  keeping the downstream ZUNMQR/ZTRSM on the host or sending back the reflectors, per CGR call.

## Open points
- Class-B/C max statistics rest on 3 seeds and rare events; if variant 2 is to be adopted a longer run of the L8 b10 and
  L12 b10 U4 w10 cases (and more seeds) would tell whether the 4.6x tau-max is a tail event. Other ideas not tried:
  v3 followed by a v2-style correction; iterating v2 once more; QRCP of a small panel on the host.
- Discrepancy: pyALF 2.6 `Simulation(..., only_prep)` writes identical `seeds` files for all cases (for DISCREPANCIES.md).
  Ticket 104's "identical seeds for up/pp" is unaffected, but any "N seeds" study through pyALF prep needs explicit seeds.
- Seeds are paired between builds, so E / S differences are not an independent statistical test.

## Full aggregated table (all builds; "aborts" = failed/total seeds, "x up" = ratio to upstream on the same seeds)
| case | build | aborts | Green mean | (x up) | Green max | (x up) | tau mean | (x up) | tau max | (x up) | Phase max | E | dE/sig | S(pi,pi) | dS/sig | wall s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L12_b10_U4_w10 | up | 0/3 | 6.80e-10 | 1.00 | 3.66e-05 | 1.00 | 2.30e-11 | 1.00 | 6.78e-06 | 1.00 | 0.00e+00 | -127.3056(0.1058) | - | 7.735(0.850) | - | 169 |
| L12_b10_U4_w10 | v1 | 0/3 | 4.03e-09 | 5.92 | 8.88e-04 | 24.23 | 2.56e-11 | 1.12 | 4.90e-06 | 0.72 | 0.00e+00 | -127.3056(0.1058) | +0.0 | 7.735(0.850) | +0.0 | 167 |
| L12_b10_U4_w10 | v2 | 0/3 | 1.28e-09 | 1.88 | 2.93e-04 | 7.99 | 2.30e-11 | 1.00 | 9.97e-07 | 0.15 | 0.00e+00 | -127.3056(0.1058) | +0.0 | 7.735(0.850) | +0.0 | 179 |
| L12_b10_U4_w10 | v3n32 | 0/3 | 3.42e-09 | 5.03 | 1.13e-03 | 30.86 | 1.95e-11 | 0.85 | 3.13e-06 | 0.46 | 0.00e+00 | -127.3056(0.1058) | +0.0 | 7.735(0.850) | -0.0 | 176 |
| L12_b10_U4_w10 | v3n16 | 0/1 | 2.12e-09 | 3.07 | 2.18e-04 | 6.91 | 2.17e-11 | 0.84 | 3.03e-06 | 0.45 | 0.00e+00 | -127.3087(0.2006) | -0.0 | 6.273(1.061) | -0.0 | 185 |
| L12_b10_U4_w10 | v3n8 | 0/3 | 8.97e-10 | 1.32 | 1.50e-04 | 4.10 | 1.93e-11 | 0.84 | 3.86e-06 | 0.57 | 0.00e+00 | -127.3056(0.1058) | -0.0 | 7.735(0.850) | +0.0 | 181 |
| L12_b10_U4_w5 | up | 0/3 | 1.07e-12 | 1.00 | 3.48e-08 | 1.00 | 1.10e-13 | 1.00 | 3.67e-09 | 1.00 | 0.00e+00 | -127.3071(0.1293) | - | 6.418(0.713) | - | 212 |
| L12_b10_U4_w5 | v1 | 0/3 | 1.10e-12 | 1.03 | 4.97e-08 | 1.43 | 1.00e-13 | 0.91 | 3.66e-09 | 1.00 | 0.00e+00 | -127.3071(0.1293) | -0.0 | 6.418(0.713) | +0.0 | 208 |
| L12_b10_U4_w5 | v2 | 0/3 | 1.13e-12 | 1.06 | 7.19e-08 | 2.06 | 1.07e-13 | 0.97 | 4.51e-09 | 1.23 | 0.00e+00 | -127.3071(0.1293) | -0.0 | 6.418(0.713) | -0.0 | 231 |
| L12_b10_U4_w5 | v3n32 | 0/3 | 9.77e-13 | 0.92 | 3.28e-08 | 0.94 | 9.17e-14 | 0.84 | 3.95e-09 | 1.08 | 0.00e+00 | -127.3071(0.1293) | -0.0 | 6.418(0.713) | +0.0 | 223 |
| L12_b10_U4_w5 | v3n16 | 0/1 | 1.04e-12 | 0.94 | 3.06e-08 | 1.08 | 1.03e-13 | 0.84 | 1.28e-09 | 0.40 | 0.00e+00 | -127.3909(0.2675) | -0.0 | 5.034(0.856) | +0.0 | 236 |
| L12_b10_U4_w5 | v3n8 | 0/3 | 9.49e-13 | 0.89 | 3.69e-08 | 1.06 | 9.09e-14 | 0.83 | 3.44e-09 | 0.94 | 0.00e+00 | -127.3071(0.1293) | -0.0 | 6.418(0.713) | +0.0 | 233 |
| L12_b10_U8_w5 | up | 0/3 | 3.44e-11 | 1.00 | 3.72e-05 | 1.00 | 1.06e-11 | 1.00 | 3.68e-05 | 1.00 | 0.00e+00 | -85.0591(0.1931) | - | 2.966(0.143) | - | 213 |
| L12_b10_U8_w5 | v1 | 0/3 | 3.99e-11 | 1.16 | 4.66e-05 | 1.25 | 8.87e-12 | 0.84 | 3.13e-05 | 0.85 | 0.00e+00 | -85.0591(0.1931) | +0.0 | 2.966(0.143) | +0.0 | 208 |
| L12_b10_U8_w5 | v2 | 0/3 | 4.28e-11 | 1.24 | 3.75e-05 | 1.01 | 1.09e-11 | 1.03 | 4.84e-05 | 1.32 | 0.00e+00 | -85.0591(0.1931) | +0.0 | 2.966(0.143) | +0.0 | 228 |
| L12_b10_U8_w5 | v3n32 | 0/3 | 2.84e-11 | 0.83 | 1.74e-05 | 0.47 | 6.77e-12 | 0.64 | 3.23e-05 | 0.88 | 0.00e+00 | -85.0591(0.1931) | +0.0 | 2.966(0.143) | +0.0 | 222 |
| L12_b10_U8_w5 | v3n16 | 0/1 | 2.31e-11 | 0.84 | 5.47e-06 | 0.60 | 5.79e-12 | 0.61 | 1.07e-05 | 0.29 | 0.00e+00 | -84.7350(0.2492) | +0.0 | 3.065(0.285) | +0.0 | 234 |
| L12_b10_U8_w5 | v3n8 | 0/3 | 2.73e-11 | 0.79 | 1.75e-05 | 0.47 | 6.81e-12 | 0.64 | 3.23e-05 | 0.88 | 0.00e+00 | -85.0591(0.1931) | +0.0 | 2.966(0.143) | +0.0 | 231 |
| L12_b20_U4_w5 | up | 0/3 | 3.89e-12 | 1.00 | 3.12e-07 | 1.00 | 2.41e-13 | 1.00 | 1.28e-08 | 1.00 | 0.00e+00 | -127.2297(0.1030) | - | 4.964(0.393) | - | 171 |
| L12_b20_U4_w5 | v1 | 0/3 | 3.87e-12 | 0.99 | 2.49e-07 | 0.80 | 2.07e-13 | 0.86 | 2.00e-08 | 1.57 | 0.00e+00 | -127.2297(0.1030) | +0.0 | 4.964(0.393) | -0.0 | 168 |
| L12_b20_U4_w5 | v2 | 0/3 | 4.01e-12 | 1.03 | 2.27e-07 | 0.73 | 2.28e-13 | 0.95 | 1.21e-08 | 0.95 | 0.00e+00 | -127.2297(0.1030) | -0.0 | 4.964(0.393) | +0.0 | 186 |
| L12_b20_U4_w5 | v3n32 | 0/3 | 3.60e-12 | 0.93 | 2.42e-07 | 0.78 | 1.90e-13 | 0.79 | 8.95e-09 | 0.70 | 0.00e+00 | -127.2297(0.1030) | +0.0 | 4.964(0.393) | +0.0 | 179 |
| L12_b20_U4_w5 | v3n16 | 0/1 | 2.45e-12 | 0.88 | 2.43e-07 | 0.78 | 1.27e-13 | 0.78 | 3.19e-09 | 0.89 | 0.00e+00 | -127.1171(0.2111) | -0.0 | 4.878(0.491) | +0.0 | 176 |
| L12_b20_U4_w5 | v3n8 | 0/3 | 3.56e-12 | 0.91 | 2.42e-07 | 0.78 | 1.90e-13 | 0.79 | 1.01e-08 | 0.79 | 0.00e+00 | -127.2297(0.1030) | +0.0 | 4.964(0.393) | +0.0 | 187 |
| L6_b10_U4_w10 | up | 0/3 | 3.82e-11 | 1.00 | 2.79e-06 | 1.00 | 9.13e-12 | 1.00 | 1.49e-06 | 1.00 | 0.00e+00 | -31.6630(0.0301) | - | 4.930(0.321) | - | 92 |
| L6_b10_U4_w10 | v1 | 0/3 | 7.41e-11 | 1.94 | 2.84e-05 | 10.16 | 1.16e-11 | 1.28 | 2.44e-06 | 1.64 | 0.00e+00 | -31.6630(0.0301) | +0.0 | 4.930(0.321) | +0.0 | 92 |
| L6_b10_U4_w10 | v2 | 0/3 | 4.00e-11 | 1.05 | 3.24e-06 | 1.16 | 9.77e-12 | 1.07 | 1.11e-06 | 0.75 | 0.00e+00 | -31.6630(0.0301) | -0.0 | 4.930(0.321) | +0.0 | 97 |
| L6_b10_U4_w10 | v3n32 | 0/3 | 7.36e-11 | 1.93 | 2.48e-05 | 8.88 | 9.89e-12 | 1.08 | 1.03e-06 | 0.69 | 0.00e+00 | -31.6630(0.0301) | +0.0 | 4.930(0.321) | -0.0 | 91 |
| L6_b10_U4_w10 | v3n16 | 0/1 | 5.03e-11 | 1.45 | 3.60e-06 | 3.52 | 8.95e-12 | 1.06 | 8.00e-07 | 1.09 | 0.00e+00 | -31.7110(0.0729) | +0.0 | 5.637(0.586) | -0.0 | 96 |
| L6_b10_U4_w10 | v3n8 | 0/3 | 7.55e-11 | 1.98 | 5.53e-05 | 19.82 | 9.36e-12 | 1.03 | 1.03e-06 | 0.69 | 0.00e+00 | -31.6630(0.0301) | +0.0 | 4.930(0.321) | +0.0 | 93 |
| L6_b10_U8_w5 | up | 0/3 | 5.39e-12 | 1.00 | 1.33e-06 | 1.00 | 2.28e-12 | 1.00 | 4.66e-06 | 1.00 | 0.00e+00 | -21.2012(0.0415) | - | 6.685(0.489) | - | 101 |
| L6_b10_U8_w5 | v1 | 0/3 | 5.41e-12 | 1.00 | 1.19e-06 | 0.89 | 2.61e-12 | 1.14 | 4.67e-06 | 1.00 | 0.00e+00 | -21.2012(0.0415) | +0.0 | 6.685(0.489) | -0.0 | 101 |
| L6_b10_U8_w5 | v2 | 0/3 | 5.41e-12 | 1.00 | 1.16e-06 | 0.88 | 2.33e-12 | 1.02 | 4.65e-06 | 1.00 | 0.00e+00 | -21.2012(0.0415) | -0.0 | 6.685(0.489) | -0.0 | 110 |
| L6_b10_U8_w5 | v3n32 | 0/3 | 5.41e-12 | 1.00 | 1.19e-06 | 0.89 | 2.40e-12 | 1.06 | 4.67e-06 | 1.00 | 0.00e+00 | -21.2012(0.0415) | +0.0 | 6.685(0.489) | -0.0 | 101 |
| L6_b10_U8_w5 | v3n16 | 0/1 | 4.27e-12 | 1.01 | 7.68e-07 | 1.00 | 1.79e-12 | 1.05 | 4.26e-06 | 1.00 | 0.00e+00 | -21.4026(0.0774) | +0.0 | 13.814(1.350) | -0.0 | 106 |
| L6_b10_U8_w5 | v3n8 | 0/3 | 5.40e-12 | 1.00 | 1.17e-06 | 0.88 | 2.30e-12 | 1.01 | 4.67e-06 | 1.00 | 0.00e+00 | -21.2012(0.0415) | +0.0 | 6.685(0.489) | -0.0 | 104 |
| L6_b20_U4_w10 | up | 0/3 | 3.93e-11 | 1.00 | 5.86e-06 | 1.00 | 1.21e-11 | 1.00 | 1.17e-05 | 1.00 | 0.00e+00 | -31.7060(0.0181) | - | 5.740(0.349) | - | 94 |
| L6_b20_U4_w10 | v1 | 0/3 | 4.37e-11 | 1.11 | 1.31e-05 | 2.23 | 1.45e-11 | 1.20 | 1.16e-05 | 0.99 | 0.00e+00 | -31.7060(0.0181) | -0.0 | 5.740(0.349) | +0.0 | 93 |
| L6_b20_U4_w10 | v2 | 0/3 | 4.01e-11 | 1.02 | 5.80e-06 | 0.99 | 1.27e-11 | 1.05 | 1.16e-05 | 0.99 | 0.00e+00 | -31.7060(0.0181) | -0.0 | 5.740(0.349) | +0.0 | 99 |
| L6_b20_U4_w10 | v3n32 | 0/3 | 4.34e-11 | 1.10 | 9.56e-06 | 1.63 | 1.29e-11 | 1.07 | 1.16e-05 | 0.99 | 0.00e+00 | -31.7060(0.0181) | +0.0 | 5.740(0.349) | +0.0 | 94 |
| L6_b20_U4_w10 | v3n16 | 0/1 | 3.51e-11 | 1.05 | 1.71e-06 | 1.00 | 1.09e-11 | 1.05 | 8.74e-07 | 1.00 | 0.00e+00 | -31.6999(0.0316) | -0.0 | 7.297(0.713) | -0.0 | 93 |
| L6_b20_U4_w10 | v3n8 | 0/3 | 4.45e-11 | 1.13 | 1.62e-05 | 2.77 | 1.24e-11 | 1.02 | 1.16e-05 | 0.99 | 0.00e+00 | -31.7060(0.0181) | +0.0 | 5.740(0.349) | -0.0 | 96 |
| L6_b20_U8_w5 | up | 0/1 | 6.53e-12 | 1.00 | 4.90e-06 | 1.00 | 4.09e-12 | 1.00 | 5.35e-06 | 1.00 | 0.00e+00 | -21.2888(0.0527) | - | 3.129(0.315) | - | 108 |
| L6_b20_U8_w5 | v1 | 0/1 | 6.53e-12 | 1.00 | 4.86e-06 | 0.99 | 4.44e-12 | 1.09 | 5.35e-06 | 1.00 | 0.00e+00 | -21.2888(0.0527) | +0.0 | 3.129(0.315) | +0.0 | 104 |
| L6_b20_U8_w5 | v2 | 0/1 | 6.52e-12 | 1.00 | 4.90e-06 | 1.00 | 4.17e-12 | 1.02 | 5.35e-06 | 1.00 | 0.00e+00 | -21.2888(0.0527) | -0.0 | 3.129(0.315) | -0.0 | 114 |
| L6_b20_U8_w5 | v3n32 | 0/1 | 6.53e-12 | 1.00 | 4.89e-06 | 1.00 | 4.25e-12 | 1.04 | 5.35e-06 | 1.00 | 0.00e+00 | -21.2888(0.0527) | +0.0 | 3.129(0.315) | +0.0 | 107 |
| L6_b20_U8_w5 | v3n16 | 0/1 | 6.53e-12 | 1.00 | 4.90e-06 | 1.00 | 4.17e-12 | 1.02 | 5.35e-06 | 1.00 | 0.00e+00 | -21.2888(0.0527) | +0.0 | 3.129(0.315) | +0.0 | 107 |
| L6_b20_U8_w5 | v3n8 | 0/1 | 6.53e-12 | 1.00 | 4.87e-06 | 0.99 | 4.10e-12 | 1.00 | 5.35e-06 | 1.00 | 0.00e+00 | -21.2888(0.0527) | +0.0 | 3.129(0.315) | +0.0 | 108 |
| L6_b20_U8_w5_dense | up | 0/3 | 6.95e-12 | 1.00 | 4.99e-06 | 1.00 | 4.98e-12 | 1.00 | 1.53e-05 | 1.00 | 0.00e+00 | -21.2425(0.0552) | - | 4.573(0.194) | - | 77 |
| L6_b20_U8_w5_dense | v1 | 0/3 | 7.08e-12 | 1.02 | 7.57e-06 | 1.52 | 5.23e-12 | 1.05 | 1.52e-05 | 0.99 | 0.00e+00 | -21.2425(0.0552) | +0.0 | 4.573(0.194) | -0.0 | 76 |
| L6_b20_U8_w5_dense | v2 | 0/3 | 7.01e-12 | 1.01 | 7.87e-06 | 1.58 | 5.01e-12 | 1.00 | 1.53e-05 | 0.99 | 0.00e+00 | -21.2425(0.0552) | +0.0 | 4.573(0.194) | +0.0 | 81 |
| L6_b20_U8_w5_dense | v3n32 | 0/3 | 7.14e-12 | 1.03 | 8.68e-06 | 1.74 | 5.07e-12 | 1.02 | 1.52e-05 | 0.99 | 0.00e+00 | -21.2425(0.0552) | +0.0 | 4.573(0.194) | -0.0 | 78 |
| L6_b20_U8_w5_dense | v3n16 | 0/1 | 6.56e-12 | 1.00 | 2.43e-06 | 1.04 | 5.41e-12 | 1.01 | 1.33e-05 | 1.01 | 0.00e+00 | -21.2170(0.0826) | -0.0 | 5.033(0.441) | -0.0 | 74 |
| L6_b20_U8_w5_dense | v3n8 | 0/3 | 7.11e-12 | 1.02 | 1.26e-05 | 2.53 | 4.93e-12 | 0.99 | 1.52e-05 | 0.99 | 0.00e+00 | -21.2425(0.0552) | +0.0 | 4.573(0.194) | +0.0 | 78 |
| L6_b40_U4_w10 | up | 0/1 | 3.95e-11 | 1.00 | 8.34e-07 | 1.00 | 2.09e-11 | 1.00 | 2.50e-05 | 1.00 | 0.00e+00 | -31.7892(0.0355) | - | 4.708(0.248) | - | 104 |
| L6_b40_U4_w10 | v1 | 0/1 | 3.99e-11 | 1.01 | 8.35e-07 | 1.00 | 2.41e-11 | 1.15 | 2.50e-05 | 1.00 | 0.00e+00 | -31.7892(0.0355) | -0.0 | 4.708(0.248) | +0.0 | 94 |
| L6_b40_U4_w10 | v2 | 0/1 | 3.97e-11 | 1.00 | 8.34e-07 | 1.00 | 2.19e-11 | 1.05 | 2.50e-05 | 1.00 | 0.00e+00 | -31.7892(0.0355) | +0.0 | 4.708(0.248) | +0.0 | 107 |
| L6_b40_U4_w10 | v3n32 | 0/1 | 3.96e-11 | 1.00 | 8.35e-07 | 1.00 | 2.25e-11 | 1.07 | 2.50e-05 | 1.00 | 0.00e+00 | -31.7892(0.0355) | -0.0 | 4.708(0.248) | +0.0 | 100 |
| L6_b40_U4_w10 | v3n16 | 0/1 | 4.00e-11 | 1.01 | 8.35e-07 | 1.00 | 2.22e-11 | 1.06 | 2.50e-05 | 1.00 | 0.00e+00 | -31.7892(0.0355) | -0.0 | 4.708(0.248) | +0.0 | 105 |
| L6_b40_U4_w10 | v3n8 | 0/1 | 4.02e-11 | 1.02 | 8.35e-07 | 1.00 | 2.17e-11 | 1.04 | 2.50e-05 | 1.00 | 0.00e+00 | -31.7892(0.0355) | -0.0 | 4.708(0.248) | +0.0 | 101 |
| L6_b40_U8_w5 | up | 0/1 | 6.46e-12 | 1.00 | 2.49e-06 | 1.00 | 5.87e-12 | 1.00 | 1.30e-05 | 1.00 | 0.00e+00 | -21.2915(0.0492) | - | 3.424(0.165) | - | 112 |
| L6_b40_U8_w5 | v1 | 0/1 | 6.47e-12 | 1.00 | 2.49e-06 | 1.00 | 6.24e-12 | 1.06 | 1.30e-05 | 1.00 | 0.00e+00 | -21.2915(0.0492) | -0.0 | 3.424(0.165) | +0.0 | 113 |
| L6_b40_U8_w5 | v2 | 0/1 | 6.43e-12 | 0.99 | 2.49e-06 | 1.00 | 5.96e-12 | 1.02 | 1.30e-05 | 1.00 | 0.00e+00 | -21.2915(0.0492) | -0.0 | 3.424(0.165) | +0.0 | 113 |
| L6_b40_U8_w5 | v3n32 | 0/1 | 6.48e-12 | 1.00 | 2.49e-06 | 1.00 | 6.11e-12 | 1.04 | 1.30e-05 | 1.00 | 0.00e+00 | -21.2915(0.0492) | -0.0 | 3.424(0.165) | +0.0 | 113 |
| L6_b40_U8_w5 | v3n16 | 0/1 | 6.51e-12 | 1.01 | 2.49e-06 | 1.00 | 5.97e-12 | 1.02 | 1.30e-05 | 1.00 | 0.00e+00 | -21.2915(0.0492) | -0.0 | 3.424(0.165) | +0.0 | 113 |
| L6_b40_U8_w5 | v3n8 | 0/1 | 6.52e-12 | 1.01 | 2.49e-06 | 1.00 | 5.92e-12 | 1.01 | 1.30e-05 | 1.00 | 0.00e+00 | -21.2915(0.0492) | -0.0 | 3.424(0.165) | +0.0 | 115 |
| L6_b4_U4_w10 | up | 0/3 | 2.53e-11 | 1.00 | 1.01e-06 | 1.00 | 6.36e-12 | 1.00 | 7.71e-07 | 1.00 | 0.00e+00 | -29.6683(0.0176) | - | 3.216(0.089) | - | 92 |
| L6_b4_U4_w10 | v1 | 0/3 | 3.47e-11 | 1.37 | 2.75e-06 | 2.71 | 8.39e-12 | 1.32 | 7.86e-07 | 1.02 | 0.00e+00 | -29.6683(0.0176) | +0.0 | 3.216(0.089) | -0.0 | 91 |
| L6_b4_U4_w10 | v2 | 0/3 | 2.67e-11 | 1.06 | 1.44e-06 | 1.42 | 6.81e-12 | 1.07 | 1.52e-06 | 1.97 | 0.00e+00 | -29.6683(0.0176) | +0.0 | 3.216(0.089) | +0.0 | 96 |
| L6_b4_U4_w10 | v3n32 | 0/3 | 3.40e-11 | 1.34 | 1.87e-06 | 1.84 | 6.96e-12 | 1.09 | 7.79e-07 | 1.01 | 0.00e+00 | -29.6683(0.0176) | +0.0 | 3.216(0.089) | -0.0 | 93 |
| L6_b4_U4_w10 | v3n16 | 0/1 | 3.40e-11 | 1.37 | 3.08e-06 | 3.49 | 6.53e-12 | 1.06 | 4.29e-07 | 0.99 | 0.00e+00 | -29.6538(0.0199) | -0.0 | 3.131(0.140) | -0.0 | 94 |
| L6_b4_U4_w10 | v3n8 | 0/3 | 3.44e-11 | 1.36 | 1.45e-05 | 14.26 | 6.53e-12 | 1.03 | 7.85e-07 | 1.02 | 0.00e+00 | -29.6683(0.0176) | -0.0 | 3.216(0.089) | +0.0 | 91 |
| L6_b4_U8_w5 | up | 0/3 | 4.11e-12 | 1.00 | 7.15e-07 | 1.00 | 1.18e-12 | 1.00 | 4.75e-07 | 1.00 | 0.00e+00 | -19.8063(0.0963) | - | 7.685(0.699) | - | 99 |
| L6_b4_U8_w5 | v1 | 0/3 | 4.47e-12 | 1.09 | 2.73e-06 | 3.82 | 1.52e-12 | 1.29 | 4.74e-07 | 1.00 | 0.00e+00 | -19.8063(0.0963) | +0.0 | 7.685(0.699) | -0.0 | 97 |
| L6_b4_U8_w5 | v2 | 0/3 | 4.12e-12 | 1.00 | 7.13e-07 | 1.00 | 1.22e-12 | 1.04 | 4.75e-07 | 1.00 | 0.00e+00 | -19.8063(0.0963) | -0.0 | 7.685(0.699) | +0.0 | 105 |
| L6_b4_U8_w5 | v3n32 | 0/3 | 4.49e-12 | 1.09 | 2.73e-06 | 3.82 | 1.30e-12 | 1.11 | 4.74e-07 | 1.00 | 0.00e+00 | -19.8063(0.0963) | +0.0 | 7.685(0.699) | -0.0 | 97 |
| L6_b4_U8_w5 | v3n16 | 0/1 | 4.28e-12 | 1.08 | 1.51e-06 | 2.12 | 1.23e-12 | 1.07 | 3.55e-07 | 2.10 | 0.00e+00 | -19.9980(0.1624) | +0.0 | 8.853(1.363) | +0.0 | 101 |
| L6_b4_U8_w5 | v3n8 | 0/3 | 4.47e-12 | 1.09 | 7.43e-06 | 10.39 | 1.22e-12 | 1.03 | 4.74e-07 | 1.00 | 0.00e+00 | -19.8063(0.0963) | -0.0 | 7.685(0.699) | -0.0 | 99 |
| L8_b10_U4_w10 | up | 0/3 | 6.10e-11 | 1.00 | 5.55e-06 | 1.00 | 9.41e-12 | 1.00 | 2.28e-06 | 1.00 | 0.00e+00 | -56.4762(0.0395) | - | 8.429(0.773) | - | 120 |
| L8_b10_U4_w10 | v1 | 0/3 | 2.41e-10 | 3.95 | 4.17e-05 | 7.52 | 1.22e-11 | 1.29 | 1.66e-06 | 0.73 | 0.00e+00 | -56.4762(0.0395) | -0.0 | 8.429(0.773) | +0.0 | 122 |
| L8_b10_U4_w10 | v2 | 0/3 | 8.57e-11 | 1.41 | 8.67e-06 | 1.56 | 1.03e-11 | 1.10 | 1.05e-05 | 4.62 | 0.00e+00 | -56.4762(0.0395) | +0.0 | 8.429(0.773) | +0.0 | 135 |
| L8_b10_U4_w10 | v3n32 | 0/3 | 2.60e-10 | 4.27 | 6.66e-05 | 12.01 | 1.03e-11 | 1.09 | 9.39e-07 | 0.41 | 0.00e+00 | -56.4762(0.0395) | -0.0 | 8.429(0.773) | +0.0 | 129 |
| L8_b10_U4_w10 | v3n16 | 0/1 | 2.53e-10 | 4.62 | 6.65e-05 | 11.99 | 8.73e-12 | 1.06 | 7.64e-07 | 1.25 | 0.00e+00 | -56.4978(0.0802) | +0.0 | 12.045(1.751) | +0.0 | 123 |
| L8_b10_U4_w10 | v3n8 | 0/3 | 2.05e-10 | 3.36 | 8.34e-05 | 15.04 | 9.63e-12 | 1.02 | 9.74e-07 | 0.43 | 0.00e+00 | -56.4762(0.0395) | -0.0 | 8.429(0.773) | -0.0 | 132 |
| L8_b10_U8_w5 | up | 0/3 | 8.02e-12 | 1.00 | 2.01e-05 | 1.00 | 3.26e-12 | 1.00 | 5.60e-05 | 1.00 | 0.00e+00 | -37.7285(0.0894) | - | 10.069(0.904) | - | 146 |
| L8_b10_U8_w5 | v1 | 0/3 | 8.01e-12 | 1.00 | 2.23e-05 | 1.11 | 3.40e-12 | 1.04 | 1.44e-05 | 0.26 | 0.00e+00 | -37.7285(0.0894) | -0.0 | 10.069(0.904) | +0.0 | 146 |
| L8_b10_U8_w5 | v2 | 0/3 | 7.75e-12 | 0.97 | 2.01e-05 | 1.00 | 3.11e-12 | 0.96 | 1.05e-05 | 0.19 | 0.00e+00 | -37.7285(0.0894) | -0.0 | 10.069(0.904) | -0.0 | 164 |
| L8_b10_U8_w5 | v3n32 | 0/3 | 8.02e-12 | 1.00 | 2.23e-05 | 1.11 | 3.15e-12 | 0.97 | 1.44e-05 | 0.26 | 0.00e+00 | -37.7285(0.0894) | -0.0 | 10.069(0.904) | +0.0 | 155 |
| L8_b10_U8_w5 | v3n16 | 0/1 | 7.44e-12 | 1.00 | 4.27e-06 | 1.00 | 3.13e-12 | 1.03 | 7.78e-06 | 1.41 | 0.00e+00 | -37.7722(0.1329) | -0.0 | 7.300(1.245) | +0.0 | 157 |
| L8_b10_U8_w5 | v3n8 | 0/3 | 7.98e-12 | 1.00 | 2.23e-05 | 1.11 | 3.04e-12 | 0.93 | 1.43e-05 | 0.25 | 0.00e+00 | -37.7285(0.0894) | -0.0 | 10.069(0.904) | -0.0 | 154 |
| L8_b20_U4_w10 | up | 0/3 | 1.19e-10 | 1.00 | 8.51e-06 | 1.00 | 2.12e-11 | 1.00 | 8.05e-06 | 1.00 | 0.00e+00 | -56.6065(0.0652) | - | 4.720(0.285) | - | 130 |
| L8_b20_U4_w10 | v1 | 0/3 | 2.40e-10 | 2.03 | 3.27e-05 | 3.84 | 2.47e-11 | 1.17 | 8.05e-06 | 1.00 | 0.00e+00 | -56.6065(0.0652) | +0.0 | 4.720(0.285) | +0.0 | 128 |
| L8_b20_U4_w10 | v2 | 0/3 | 1.36e-10 | 1.14 | 2.24e-05 | 2.64 | 2.21e-11 | 1.04 | 8.05e-06 | 1.00 | 0.00e+00 | -56.6065(0.0652) | -0.0 | 4.720(0.285) | +0.0 | 138 |
| L8_b20_U4_w10 | v3n32 | 0/3 | 2.57e-10 | 2.17 | 9.14e-05 | 10.75 | 2.26e-11 | 1.07 | 8.05e-06 | 1.00 | 0.00e+00 | -56.6065(0.0652) | +0.0 | 4.720(0.285) | +0.0 | 135 |
| L8_b20_U4_w10 | v3n16 | 0/1 | 1.81e-10 | 1.61 | 3.08e-05 | 6.98 | 2.15e-11 | 1.03 | 1.49e-06 | 0.58 | 0.00e+00 | -56.5702(0.0849) | +0.0 | 5.014(0.554) | -0.0 | 134 |
| L8_b20_U4_w10 | v3n8 | 0/3 | 2.35e-10 | 1.98 | 7.51e-05 | 8.83 | 2.17e-11 | 1.02 | 8.05e-06 | 1.00 | 0.00e+00 | -56.6065(0.0652) | +0.0 | 4.720(0.285) | +0.0 | 137 |
| L8_b20_U8_w5 | up | 0/3 | 9.71e-12 | 1.00 | 1.10e-05 | 1.00 | 8.74e-12 | 1.00 | 4.18e-05 | 1.00 | 0.00e+00 | -37.7852(0.0895) | - | 4.063(0.155) | - | 151 |
| L8_b20_U8_w5 | v1 | 0/3 | 9.90e-12 | 1.02 | 1.10e-05 | 1.00 | 9.13e-12 | 1.05 | 4.24e-05 | 1.02 | 0.00e+00 | -37.7852(0.0895) | -0.0 | 4.063(0.155) | -0.0 | 150 |
| L8_b20_U8_w5 | v2 | 0/3 | 9.76e-12 | 1.00 | 1.10e-05 | 1.00 | 8.81e-12 | 1.01 | 4.17e-05 | 1.00 | 0.00e+00 | -37.7852(0.0895) | -0.0 | 4.063(0.155) | -0.0 | 168 |
| L8_b20_U8_w5 | v3n32 | 0/3 | 9.91e-12 | 1.02 | 1.10e-05 | 1.00 | 8.87e-12 | 1.02 | 4.24e-05 | 1.02 | 0.00e+00 | -37.7852(0.0895) | -0.0 | 4.063(0.155) | -0.0 | 158 |
| L8_b20_U8_w5 | v3n16 | 0/1 | 1.21e-11 | 1.04 | 1.10e-05 | 1.00 | 7.44e-12 | 1.06 | 3.99e-05 | 1.09 | 0.00e+00 | -37.8077(0.1692) | -0.0 | 4.502(0.346) | -0.0 | 162 |
| L8_b20_U8_w5 | v3n8 | 0/3 | 9.90e-12 | 1.02 | 1.10e-05 | 1.00 | 8.67e-12 | 0.99 | 4.25e-05 | 1.02 | 0.00e+00 | -37.7852(0.0895) | -0.0 | 4.063(0.155) | -0.0 | 156 |
| L8_b40_U4_w10 | up | 0/1 | 1.72e-10 | 1.00 | 1.14e-05 | 1.00 | 4.70e-11 | 1.00 | 2.30e-05 | 1.00 | 0.00e+00 | -56.5758(0.0955) | - | 3.415(0.234) | - | 128 |
| L8_b40_U4_w10 | v1 | 0/1 | 1.82e-10 | 1.06 | 1.21e-05 | 1.06 | 5.29e-11 | 1.12 | 3.15e-05 | 1.37 | 0.00e+00 | -56.5758(0.0955) | +0.0 | 3.415(0.234) | -0.0 | 128 |
| L8_b40_U4_w10 | v2 | 0/1 | 1.76e-10 | 1.02 | 1.22e-05 | 1.07 | 4.77e-11 | 1.02 | 2.30e-05 | 1.00 | 0.00e+00 | -56.5758(0.0955) | +0.0 | 3.415(0.234) | -0.0 | 136 |
| L8_b40_U4_w10 | v3n32 | 0/1 | 1.89e-10 | 1.10 | 1.21e-05 | 1.06 | 5.06e-11 | 1.08 | 3.15e-05 | 1.37 | 0.00e+00 | -56.5758(0.0955) | +0.0 | 3.415(0.234) | -0.0 | 135 |
| L8_b40_U4_w10 | v3n16 | 0/1 | 1.88e-10 | 1.09 | 1.21e-05 | 1.06 | 5.01e-11 | 1.07 | 3.15e-05 | 1.37 | 0.00e+00 | -56.5758(0.0955) | +0.0 | 3.415(0.234) | -0.0 | 134 |
| L8_b40_U4_w10 | v3n8 | 0/1 | 1.91e-10 | 1.11 | 1.25e-05 | 1.10 | 4.97e-11 | 1.06 | 3.15e-05 | 1.37 | 0.00e+00 | -56.5758(0.0955) | +0.0 | 3.415(0.234) | -0.0 | 135 |


## Manager review (2026-10-27)
- **Reproduced:** L8 b10 U4 w10 seed s3, builds `up` and `v2`, re-run by the manager from the stored parameters/seeds:
  all `Precision` lines bit-identical to the study's (incl. v2 tau max 1.05e-5 vs upstream 2.28e-6).
- **numerics-reviewer: approve.** Macro-off code unchanged (470 insertions, 0 deletions, all inside `#if defined(QR_PREPIVOT)`);
  permutation composition, single ZGEQRF-format reflector set, workspace and phase bookkeeping correct for v2 and v3.
  Nits, not fixed (test helpers only): `test_variants.F90` prints but does not assert its checks, M=96 is a multiple of
  every nb tested (no partial last block, no nb > N, no rectangular case), no zero-norm/tied columns, unseeded RNG; timings
  in `bench_variants.F90` include temporaries' allocation. These must be fixed before any of this code leaves the study.
- **How "max" is aggregated:** the ratios above take the maximum over seeds per build. Per seed the spread is larger:
  at L8 b10 U4 w10 seed s2, Green max is 8.7e-6 (v2) vs 1.6e-6 (up), 5.4x. The v2 means are systematically higher
  (1.3–1.4x in that case), i.e. close to the 1.5x limit, not noise around 1.
- **Not patched:** `udv_state_mod.F90` decompose (already pre-pivoted by STABLOG itself); the study covers CGR/CGR2_2 only.
- Ran ≈ 4 h instead of the ≲ 3 h budget (the seed re-run after the pyALF seed bug).
- **Verdict (manager):** no GPU-friendly variant meets the Goal unconditionally. For STABLOG, the GPU MVP keeps upstream
  ZGEQP3 on the host for the CGR/CGR2_2 factorization (ADR-001 §4). Variant 2 stays a test-only macro; adopting it would
  need the owner's decision plus a longer, more-seed study of the class-B cases.

# Cluster request (project 2, end of Phase 6) — 2026-10-31

Branch `feature/delayed-submatrix` (local only). Nothing pushed. No cluster access so far.
All numbers below come from logs written by the scripts named in brackets; details and commands in `CLAUDE.md`
(sections "Phase 0" … "Phase 6 results") and `docs/benchmarks/`.

## 1. What was done locally (Phases 0–6)

| Phase | Result | Evidence |
|---|---|---|
| 0 | Starting state re-verified: ctest 40/40, regress 9/9 IDENTICAL, trajectories IDENTICAL over 21 600 proposals | CLAUDE.md "Phase 0" |
| 1 | `Op_V_is_symmetric` fast path: setup at L = 20 24.5 s → 0.2 s; **L = 60 (N = 3600): 1.9 s**. Scheme / N_delay validated right after the parameter read | test 37 (270 random sets vs the dense original, 0 mismatches) |
| 2 | Design for rank-N vertices, M-form submatrix, cross-slice; independent re-derivation agrees (≤ 1.4e-14) | `docs/superpowers/specs/2026-10-31-rankN-design.md` (D1–D15), `docs/independent-derivation-rankN.md` |
| 3 | **Rank-N delayed**: all tiers pass for 14 sets (tV, Kondo-J, Spin-Peierls, Z2, Jz, projector tV, complex hopping, type-1 Ising, …) | `testsuite/lowrank/tiers.sh delayed 3` |
| 4 | **Rank-N submatrix (M-form)**: all tiers pass; trajectories identical to standard and to delayed (34 comparisons) | `tiers.sh submatrix 4` |
| 5 | Projector statistics pass; local MPI (1/2/4 ranks), threads, tempering, PARALLEL_PARAMS all pass; cross-slice prototype measured | `parallel_tier.sh --rebuild`, `docs/benchmarks/cross_slice.md` |
| 6 | **Automatic N_delay** (default `N_delay = 0`: deterministic static rule; `-1`: + warm-up timing; `≥ 1`: fixed), logged in `info` and (cluster) `data.h5` | ctest 43/43, `docs/benchmarks/auto_ndelay.md` |

Validation per tier (both schemes unless noted):
1. **Unit** (ctest 36): rank-N cases at N_delay {1, 2, 7, 32, 1000}: delayed ≤ 2.9e-12 (ratio) / 2.3e-13 (G) well-conditioned;
   ill-conditioned cases ≤ 3.6e-10 (TOL_ILL 1e-8). Mutation tests fail as they should (e.g. 142 rows without the M transform).
2. **Trajectory**: identical Metropolis decisions standard vs delayed vs submatrix for all 14 sets, also at N_delay 2 and 1000,
   with N_delay 0 and −1, and with the cross-slice key on.
3. **Statistical** (pre-registered rules): tv, kondo, jz, tv_proj, z2, doped, pv_adiab: PASS (submatrix Kondo after the
   pre-registered replicate; the first failure was an array-only failure that the standard chain with the same seed shows too).
4. **Regression**: default path byte-identical to baseline `f81b5cbc` for 13 sets (latest run at HEAD 2f98c1f2: 13/13).
5. **Stability** (paired, same trajectory): Precision Green max ratio scheme/standard ≤ 1.35 (delayed), ≤ 1.32 (submatrix);
   Phase max ≤ 1.61 / 1.54.
6. **Parallel** (local): MPI 18/18, threads (omp1 vs omp4 bit-identical), tempering 3/3, PARALLEL_PARAMS 3/3.

Mac benchmarks (M1 Pro, 1 thread):
- Hubbard L = 20 (N = 400), N_delay 32: speedup vs standard 2.2× (delayed) / 2.3× (submatrix) (`docs/benchmarks/size_scan.md`).
- Bond models at L = 12 (`docs/benchmarks/d13_bond_models.md`): delayed tV 1.59×, Kondo 2.08×, Spin-Peierls 1.30×;
  submatrix 0.97× / 1.23× / 1.21× → **delayed is recommended for bond models**, submatrix for Hubbard-type models.
- Automatic N_delay: all 20 points (Hubbard L = 12/16/20, tV L = 8/12, both schemes) within noise of the best fixed value;
  static rule worst +3.8 %, timing worst +2.4 % (`docs/benchmarks/auto_ndelay.md`).
- Cross-slice delay: no measurable gain (on/off 0.97–1.07, within noise) (`docs/benchmarks/cross_slice.md`).
- Large-N precision pre-check: see section 6.
- BLAS ceiling (Mac, indicative): see section 6.

Tolerances: **none changed** (unit 1e-10 / 1e-8 ill-conditioned; statistical rules as pre-registered in project 1).

## 2. Decisions that need your review
1. **Cross-slice prototype kept as opt-in** (`Lowrank_cross_slice = .F.` default), documented as "no measured gain, not recommended".
   This deviates from my own pre-registered rule (no gain → remove). Reason: removing it reworks validated code paths for no
   correctness or default-path benefit. → open question 00.
2. **Submatrix accuracy guard K_ND = 1 (D15)**: non-diagonal content is never kept in the M-form (found by the Kondo precision tier:
   65× worse Precision Green max with K_ND = 8). Consequence: submatrix is not faster than standard for bond vertices (D13). → open question 1.
3. **Submatrix-T (arXiv:2602.03656) not implemented** (D9); SaSb-zz covered at the unit tier only, Jz at all tiers via the test model
   `LR_Jz_Test` (D11).
4. **N_delay bounds the pending rank** (columns / pending sites), unchanged for rank-1 (D5).
5. **New default `N_delay = 0` = automatic (static rule, deterministic)**; timing only with `-1` (it breaks bit-reproducibility across
   machines). The previous default was 32.
6. **HDF5 provenance only for non-standard schemes** (standard data.h5 stays as in the baseline); compiled only against stubs here.
7. Fallbacks: vertices that are not eligible (support > 8, etc.) take the standard update after a flush (unchanged from project 1).

## 3. Open questions
- **00** Cross-slice prototype: keep (opt-in, as now) or remove?
- **0** Helma + fakeintel: ALF's `HELMA` preset links with `-static`, so `LD_PRELOAD` of fakeintel cannot act. Proposal: build the
  Helma benchmark binary without `-static` (dynamic MKL) and preload fakeintel; the jobs abort if fakeintel is requested for a static
  binary. Agree? How is fakeintel meant to be loaded on Helma (module, path, LD_PRELOAD)?
- **1** Fast submatrix for bond models: the option is a hybrid that routes non-diagonal blocks to the delayed form (not built).
  Wanted, or is "delayed for bond models" enough?

## 4. Permission request
I ask for explicit permission to start the cluster phase (Phase 7: build + validation on the test/short partitions, ≤ 1 node,
≤ 30 min, no further approval needed per your rules), and for approval of the **first campaign only: the pilot C0 + B-S on
SuperMUC-NG and H0 + B-H on Helma** (section 7). Later campaigns are requested one at a time with the pilot-calibrated estimate.

Phase 7 steps (test partitions): build all variants with the machine's toolchain (noMPI, MPI, TEMPERING, HDF5, debug, trace,
timers; Intel ifx/ifort + MKL have never been used on this branch), ctest + regress + a trajectory subset; HDF5 attribute check
(`h5dump -A -g / data.h5`, delayed vs standard); short runs of the untested combinations (HMC + low-rank, flavour reconstruction with
non-diagonal vertices, N_delay = −1 with Global_tau, cross-slice under MPI/tempering); timers overhead check; `blas_bench --probe`
and the provenance/pinning checks of the campaign generator.

## 5. Placeholders I need from you
SuperMUC-NG Phase 1:
- ssh alias (host name in your `~/.ssh/config`),
- workspace directory (I read and write only there),
- LRZ project id (`sbatch --account`),
- your default jobscript(s) (modules, EAR setting / option name, anything site-specific).

Helma:
- already known: ssh host `helma`, workspace `/hnvme/workspace/b133ae29-CLAUDE` (nothing outside it is read or written);
- jobscript / run instructions (account, partition, `I_MPI_PMI_LIBRARY` or mpirun usage, modules),
- fakeintel: library path and how to load it (see open question 0).

## 6. Local pre-checks for the cluster runs
**Large-N precision** (`docs/benchmarks/largeN_precision.md`, from `raw_largeN_precision.txt` by `precision_report.py`): Mz Hubbard
Beta 2, L = 32 (N = 1024; U = 4 and 8) and L = 40 (N = 1600; U = 4), same seed; automatic N_delay (64; submatrix 128 at L = 40) and
fixed 32 vs standard. Same acceptance in all 12 pairs (same trajectory); Precision Green max ratio 1.000, mean 0.99–1.02, Phase max
≤ 1.61 (rule ≤ 2). At U = 8 the standard run itself has Green max 8.8e-2 with Nwrap = 10, so C0 sets a shorter Nwrap at U = 8 from
the pilot.

**BLAS ceiling** (Task T2, `testsuite/lowrank/blas_bench.sh`; Mac logs `docs/benchmarks/raw_blas_bench_mac_1x{1,4}.txt`, SUMMARY lines;
quiet machine, median of 7 batches, ZGEMM as in the flush: alpha = −1, beta = 1): 1 thread, N = 400 / 1600 / 3600: ZGEMM N×N×N
243 / 258 / 252 GFLOP/s, ZGEMM N×N×64 177 / 145 / 129, ZGERU 7.6 / 8.7 / 8.2, ZGEMV 6.8 / 6.0 / 5.9, ZGEQP3 24.3 / 16.6 / 10.5
GFLOP/s; 4 threads (1 × 4): ZGEMM N×N×N 288 / 466 / 478, ZGERU 8.7 / 21.0 / 16.8. The BLAS-3 / BLAS-2 gap (1 thread: ~30–43× ZGEMM N×N×N vs ZGERU/ZGEMV on the Mac,
AMX units) is what delayed/submatrix exploit; on the clusters the same tool gives the per-layout ceilings (MPI build with a barrier
per point; `--probe` serves the MKL/fakeintel provenance check of every job).

## 7. Proposed benchmark matrix and cost
Generated by `python3 -I testsuite/lowrank/campaign/estimate_cost.py` (model and assumptions in its header; uncertainty about a
factor 2 either way; the pilot calibrates it). SuperMUC: 48 cores per node; Helma: 384 cores per node.

| campaign | jobs | node-h | core-h |
|---|---:|---:|---:|
| C0 SuperMUC pilot + large-N precision: Hubbard U=4 and U=8, Beta=4 (Ltrot 40), L=40, 60, standard + delayed/submatrix auto + fixed 32, 2x24, 1 rep | 4 | 13 | 606 |
| B-S SuperMUC BLAS ceiling (T2), layouts 1x48 2x24 4x12, N up to 8000 | 3 | 3 | 144 |
| C1 SuperMUC 2D Hubbard U=1 and U=4, Beta=2 (Ltrot 20), L=20..60, 3 schemes, auto N_delay, layouts 1x48 2x24 4x12 | 90 | 53 | 2550 |
| C1b SuperMUC 2D Hubbard U=4, Beta=4 (Ltrot 40), L=40, 60, layout 2x24 | 6 | 13 | 617 |
| C2 SuperMUC N_delay scan: L=40, 60, delayed+submatrix, fixed 16..256 (5) + auto + tune, 2x24 | 42 | 21 | 1023 |
| C3 SuperMUC 3D cubic Hubbard U=4, Beta=2, L=8, 12, 16, 20 (N up to 8000), 3 schemes, 1x48 (2x24 for L<=16) | 12 | 52 | 2520 |
| C4 SuperMUC tV (rank-2 bonds), L=20, 30, 40, standard vs delayed (+submatrix), 2x24 | 9 | 1 | 30 |
| H0 Helma pilot: Hubbard U=4, Beta=2, L=40, 60, 3 schemes, 8x48, 1 rep | 2 | 1 | 499 |
| B-H Helma BLAS ceiling (T2), layouts 8x48 16x24 2x192, fakeintel on/off, N up to 8000 | 6 | 6 | 2304 |
| H1 Helma 2D Hubbard U=4, Beta=2, L=20..60, 3 schemes, auto, layouts 8x48 16x24 2x192, fakeintel on | 45 | 15 | 5772 |
| H1b Helma fakeintel off control: L=40, 60, 8x48 | 6 | 4 | 1496 |
| H2 Helma 3D cubic L=12, 16, 20, 3 schemes, 8x48 | 9 | 45 | 17275 |
| H3 Helma N_delay scan L=40, 60, delayed+submatrix, 7 values, 8x48 | 42 | 13 | 4962 |
| **total supermuc** | | **156** | **7491** |
| **total helma** | | **84** | **32309** |

First approval requested: C0 + B-S (SuperMUC, 16 node-h / 750 core-h) and H0 + B-H (Helma, 7 node-h / 2803 core-h).

Each campaign: whole exclusive nodes, `MKL_NUM_THREADS = OMP_NUM_THREADS`, pinning checked per node, ≥ 3 repetitions with median and
spread (pilot: 1), full provenance in every log (machine, node, CPU, clock, compiler and flags, MKL version and threading layer,
fakeintel on/off, MPI, ranks × threads, model, L, β, Δτ, N_delay fixed/auto/used, commit), timers build for the time split.

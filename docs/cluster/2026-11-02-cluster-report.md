# Cluster report: delayed and submatrix updates on SuperMUC-NG (Phase 7/8, 2026-11-01/02)

Branch `feature/delayed-submatrix` (local only, nothing pushed). Machine: SuperMUC-NG Phase 1, test partition only (user rule),
2x Xeon Platinum 8174, 48 cores/node; ifort 2021.4, MKL 2020.4, Intel MPI 2019.12 (`stack/22.2.1`). Every number below comes from
logs written by the scripts named in brackets; raw logs and generated reports are under `docs/cluster/`.

## 1. Validation on SuperMUC-NG (`testsuite/lowrank/cluster/validate.sh`, jobs 5956005, 5956027, 5956034, 5956097, 5956146, 5956147)
| tier | result | log |
|---|---|---|
| unit (ctest) | 43/43 (unit-test link needs `-qopt-matmul` with ifort: `ALF_TEST_LDFLAGS`) | `logs/smuc_validate_*` |
| trajectory | 44 comparisons IDENTICAL: 11 sets x {standard vs delayed, vs submatrix} x {N_delay default, -1}, incl. HMC (acceptance 1) and N_delay = -1 with Global_tau | `logs/smuc_validate_5956005_summary.txt` |
| regression | production flags (-ipo, -fp-model fast=2): confout_0 identical for all 13 sets, only Ener_scal/Kin_scal differ at ~1e-15 relative; **both built with -no-ipo: 13/13 byte-IDENTICAL** | job 5956097 |
| debug shadow | 10/10 runs (tv, kondo, jz, z2, spin_peierls x delayed/submatrix), max relative G error <= 2e-13 | job 5956005 |
| parallel | MPI 1/2/4 ranks (half, tv), OpenMP 1 vs 4 threads, tempering, PARALLEL_PARAMS: all PASS (pre-registered replicates where an array rule failed once); mpi/half/submatrix np4 replicate: see §5 | `logs/ptier_*.txt` |
| other | HDF5 provenance attributes correct; timers build; MPI 4-rank trajectory identity incl. cross-slice | job 5956005 |

Fixes made for the Intel toolchain (all validated again locally: ctest 43/43, regress 13/13): `matmul` in the Op_V_is_symmetric fast
path -> ZGEMM (2a2dabb0); debug shadow `matmul` -> ZGEMM helper (9b763ba1). Cluster pitfalls handled in the scripts: compute nodes
without git/make; `KMP_AFFINITY=...compact` from the stack modules pinned all concurrent processes to core 0.

## 2. Pilot C0 (approved; `docs/cluster/pilot/`)
- Stage S (4 nodes): packed launch verified (one run per node, placement check OK, Intel MPI pin tables: 2 disjoint 24-core domains).
- Stage A (16 nodes, 2x24, same seed per point => same chain): paired Precision Green ratios scheme/standard 0.985–1.018, incl.
  automatic N_delay 64/128/256 at N up to 3136 (`pilot/smuc_c0_A_precision.md`). Speedups vs standard (`pilot/smuc_c0_A_report.md`):
  L=40 Beta=4: 1.61–1.86; L=48 Beta=2 U=8: delayed 1.80, submatrix 2.15; L=56 Beta=2 U=4: delayed 2.64, submatrix 2.87.
- Stage B2 (10 nodes): calibration of sweep times for the scan (`pilot/smuc_c0_B2_report.md`).

## 3. BLAS ceilings (`testsuite/lowrank/blas_bench.sh` via `cluster/blas_layouts.sbatch`, job 5956151; `blas/raw_smuc_*`)
1x48, N = 3600: ZGEMM N x N x N 2.66 TFLOP/s per node, ZGEMM N x N x 256 2.35, ZGERU 53 GFLOP/s, ZGEMV 28 GFLOP/s (BLAS-3/BLAS-2 ~50x).
2x24 and 4x12 in the raw logs (SUMMARY lines; per-rank and node sums).

## 4. N_delay scan (`campaign/descriptions/smuc_ndscan.py`, 714 runs, 3 repetitions; `ndscan/ndscan.md`, plots `ndscan/ndscan_vs_ndelay.png`, `ndscan/ndscan_abs.png`)
12 parameter points x {standard, delayed, submatrix} x N_delay {8 ... 256, auto (0), auto+timing (-1)}, plus the base point in
1x48 and 4x12; one run per node, same seed per point (same chain); 17 packed jobs (2026-11-01 14:00–17:45). `collect.py --strict`:
714/714 OK (allowed unavailable fields: job.boost_job_start, job.numa_nodes, run.boost_start — the nodes expose no boost/NUMA record).
Sweep time from the T1 timers (bins >= 2), median of 3 repetitions.

**Noise floor.** Repetition ranges are mostly < 1 %, but identical configurations measured as different runs (automatic N_delay vs
the fixed run with the same value, 28 pairs) differ by median 0.0 %, max 3.1 %: differences below ~3 % are not significant.

**Optimum.** The curves are flat near the optimum: for Hubbard L = 40 (all Ltrot, U and layouts) every N_delay from 32 to 128 is
within 3 % of the best. Smaller systems prefer smaller delayed values (L = 20: 8–32, 64–256 cost 7–10 %; L = 30: 32–64, 128 costs 8 %).
8 costs 4–27 % for submatrix (most at large N), 256 costs 4–40 % for submatrix. Best fixed values: delayed 8–128, submatrix 8–128,
mostly 32–64 (table in `ndscan.md`).
- System size (L = 20/30/40/50, Ndim 400–2500): submatrix optimum 64/64/32/128, delayed 16/32/64/16 (L=50 delayed: 16 is 2–5 %
  ahead of 32–128). Speedup best vs standard grows with N: 1.12–1.16, 1.17–1.26, 1.74–1.82, 2.86 (delayed) / 3.49 (submatrix).
- Time slices (Beta = 1/2/4, Ltrot 10/20/40) and U (1/4/8): **no effect** on the optimum (32–64 for both schemes; all within noise).
- Model / vertex rank: the optimum is **model dependent** and, for delayed, not a smooth function of Ndim: tV delayed best at 8
  (L=20, 64+ costs 18–30 %), 128 (L=30, monotonic: 32 costs 15 %) and 64 (L=40). Kondo L=20: flat (delayed 3.77–3.98 s, non-monotonic).
  Submatrix gives **no speedup for tV** (0.99–1.02 vs standard) and only 1.06 for Kondo L=20; delayed gives 1.17–1.56 (tV) and
  1.14 (Kondo). For Hubbard submatrix is 3–7 % faster than delayed at L <= 40 and 18 % at L = 50.
- Threads per rank (base point): optimum 32 (submatrix) and 8–64 (delayed, 1x48 flat), unchanged by the layout. Throughput per node
  with the best N_delay (ranks / s per sweep): 1x48 0.026, 2x24 0.060 (submatrix), 4x12 0.103 sweeps/s; standard 0.024 / 0.033 / 0.037.
  The speedup vs standard is 1.09–1.13 at 1x48, 1.74–1.82 at 2x24 and 2.72–2.80 at 4x12 (more ranks per node = less threading in
  the BLAS-3 part of standard, the delayed/submatrix advantage grows).

**Automatic N_delay (static rule, N_delay = 0).** Within 3 % of the best fixed value for 22 of 28 (point, scheme) curves. Outside:
tvL30 delayed 1.146 (rule 32, best 128), tvL40 submatrix 1.082 (flat, noisy curve without speedup), hubL50 delayed 1.078 (rule 128,
best 16; 1.054 against the fixed-128 run), kondoL20 delayed 1.055 (non-monotonic curve), hubL40 submatrix family 1.017–1.023 (rule 128,
best 32). **Timing (N_delay = -1)** is not better than the static rule here: 1.000–1.091, and its sweep times include the trial
sweeps when sweeps 2–3 fall into timed bins (1 sweep per bin at L >= 40: up to +6.6 % against the fixed run with the chosen value),
so its ratios are upper bounds for these short runs.

**Conclusion.** The static rule is good enough: the loss is <= 3 % in most cases and at most 15 % (one model point); no rule in
Ndim alone fits the tV delayed optimum (8 / 128 / 64 at L = 20/30/40), so the data do not justify changing the default. Possible
refinement (not done): submatrix cap 64 instead of 128 for Ndim <= 1600 would gain ~2 % at L = 40. The more useful finding is the
scheme choice: delayed for tV and Kondo, submatrix for Hubbard.

Caveats: (1) every run printed OpenMP affinity on each parallel region (`OMP_DISPLAY_AFFINITY`, ~1.2 M lines, 75 MB per run). The
line count is the same for all schemes of a point (e.g. 1 197 082 for standard and delayed at hubL40 2x24), so it is a constant cost per
point: it cannot move the optimum. Measured A/B (job 5957003, `cluster/affinity_overhead.sbatch`, hubL40 2x24, ABBA order,
`logs/affinity_ab_5957003/`): delayed N_delay 64 34.94 s/sweep with vs 34.39 without (+1.6 %), standard 60.96 vs 60.14 (+1.4 %);
speedup 1.745 vs 1.749, so the absolute times are ~1.5 % high and the ratios are unaffected. (2) The 51 1x48 runs were first marked
INVALID because Intel MPI 2019 truncates the pin table after ~30 cpus; the pinning check now uses the OMP affinity in that case
(f9ab32f4) and `rejudge_pinning.sh` re-judged them (51 OK; original verdict kept in run.prov).

## 5. Open items / decisions for the user
- Cross-slice prototype: parked by the user on 2026-11-01 — **reminder: decide keep (opt-in) or remove**.
- 3D cubic lattice: not in this ALF version; skipped (user decision).
- Helma: not touched (maintenance; off limits).
- Scheme choice per model (delayed for tV/Kondo, submatrix for Hubbard) is documentation only; the default is unchanged.
- Merge / push of the branch: user decision.

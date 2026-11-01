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

## 4. N_delay scan (`campaign/descriptions/smuc_ndscan.py`, 714 runs, 3 repetitions)
(filled in from `ndscan_analyze.py` when the campaign is complete)

## 5. Open items / decisions for the user
- Cross-slice prototype: parked by the user on 2026-11-01 — **reminder: decide keep (opt-in) or remove**.
- 3D cubic lattice: not in this ALF version; skipped (user decision).
- Helma: not touched (maintenance; off limits).
- Merge / push of the branch: user decision.

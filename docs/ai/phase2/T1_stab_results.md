# T1 — stabilisation candidates on the H100 (CG-ALF, 2026-11-26/27)

Owner's go: Q14 (2026-11-26), cap 4 GPU-h. Used: **0.40 GPU-h** (Helma `preempt`, 1 × H100 each):

| Job | What | Elapsed |
|---|---|---|
| 1463455 | cold-start scan; only the runs without `ALF_GPU_*` switches ran (production runs refuse switches without `ALF_GPU_DEBUG=1`) | 0:02:51 |
| 1463834 | cold-start scan, timing hook, all with `ALF_GPU_DEBUG=1` | 0:10:24 |
| 1463941 | thermalised scan (warm-up + Nwrap 10/20/30/40) | 0:10:28 |

Code: research branch `rd/cgalf-stab` (fork; 31bd7ad9 prototype, 7ad95cf6 raw logs in `docs/cgalf_T1/`), i.e. `gpu-port` plus
CholeskyQR2 in `dev_decompose` (`gpu_stab_mod.F90`, switch `ALF_GPU_CHOLQR=1`) and `la_potrf` (`alf_la_mod.F90`). Not merged anywhere.

## Candidates

- **Pre-pivoted QR** (current GPU default): column-norm pre-pivot of U·D + unpivoted `geqrf`, Q via `ormqr`, one Newton-Schulz step.
- **CholeskyQR2 after scale separation** (idea I11): pivot p from the norms of U·D; CholeskyQR2 on the *unscaled* block
  M = U(:,p) (W = M^H M, `potrf`, Q1 = M R1^-1, repeated once); R = R2 R1 diag(D(p)); then the same D separation, phase and
  V update as before. Q is explicit (no `ormqr`). If a Cholesky fails (info > 0), the pre-pivoted QR is used for that step.
  **CUDA and host backends only.** The prototype refuses the oneMKL and HIP builds. Its Q is not a product of Householder
  reflectors, so the tracked determinant phase (S1, oneMKL) is not valid on this path; CUDA takes the phase from the LU anyway.
- **Adaptive Nwrap**: measured here only as the *offline* choice of Nwrap from the scan below. No run-time controller exists.

## Setup

Production model (`/jakob_params` without the defect): Hubbard, N_FL = 2, Mz = .F., flux, U = 5.3, beta = L, dtau = 0.1, submatrix
N_delay = 32, ham_h0 = 0 (device observables), `ALF_GPU_DEBUG=1` on every compared run (no effect on numbers or time: checked
against the switch-free runs of job 1463455 to the digit).

## Table 1 — cost per stabilisation (job 1463834, `ALF_GPU_STABTIME=1`)

The decompose cost is the intercept a of the measured WRAPUR+decompose block time t(Nwrap) = a + b·Nwrap at Nwrap 5/10/20/40
(both flavours, median of 5). It is a fit, not a direct timer.

| L (N) | method | decompose a (ms) | wrap b (ms / slice) |
|---|---|---|---|
| 24 (576) | pre-pivoted QR | 6.9 | 0.13 |
| 24 (576) | CholeskyQR2 | 5.8 (−16 %) | 0.13 |
| 32 (1024) | pre-pivoted QR | 14.3 | 0.49 |
| 32 (1024) | CholeskyQR2 | 10.6 (−26 %) | 0.49 |

CGR (unchanged by the candidate): 11.8 ms (L = 24), 23.4 ms (L = 32) per call, device only.

## Table 2 — Precision Green max and time per sweep vs Nwrap

**Thermalised** (job 1463941): warm-up from a cold start with the current QR, Nwrap 10 (NBin 1; L = 24: 40 sweeps, L = 32: 25 sweeps;
script `~/cgalf/T1/job_th.sh` on Helma). Then every run started from that run's `confout_0` (as `confin_0`), 4 sweeps, same
seeds. Acceptance was identical for QR and CholQR in every pair.

| L | method | PG max Nwrap 10 | 20 | 30 | 40 | s / sweep Nwrap 10 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|
| 24 | QR | 1.6e-11 | 4.0e-11 | 1.2e-9 | 5.1e-9 | 2.83 | 2.18 | 1.99 | 1.90 |
| 24 | CholQR2 | 9.7e-12 | 4.4e-11 | 8.8e-10 | 6.8e-9 | 2.70 | 2.18 | 1.98 | 1.88 |
| 32 | QR | 6.0e-11 | 5.2e-10 | 4.6e-9 | **1.3e-8** | 7.82 | 6.26 | 5.72 | 5.43 |
| 32 | CholQR2 | 1.2e-10 | 2.5e-10 | 5.3e-9 | **1.1e-8** | 7.58 | 6.14 | 5.66 | 5.38 |

PG mean grows from ~2e-14 (Nwrap 10) to ~2e-11 (Nwrap 40), equal for both methods within 10 %. Cholesky failures in all sweeps: 0.

**Cold start** (job 1463834, 3 sweeps): PG max at Nwrap 5/10/20/40: L = 24 QR 1.2e-11 / 9.0e-12 / 8.3e-11 / 8.6e-9, CholQR2
8.0e-12 / 8.1e-12 / 6.8e-11 / 9.7e-9; L = 32 QR 3.2e-11 / 7.0e-11 / 9.6e-11 / 1.0e-8, CholQR2 6.7e-11 / 4.6e-11 / 9.9e-11 / 1.0e-8.
The thermalised configuration is up to 5× worse than the cold start at L = 32, Nwrap 20 (5.2e-10 vs 9.6e-11).

## Verdict

1. **Nwrap 20 is safe at L ≤ 32 on a thermalised configuration**: PG max ≤ 5.2e-10, a factor ~20 below the 1e-8 limit, for both
   methods. It saves ~20 % per sweep against Nwrap 10 (L = 24: 2.83 → 2.18 s, L = 32: 7.82 → 6.26 s).
   Nwrap 30 has a margin of only 2× at L = 32 (4.6e-9). Nwrap 40 breaks the limit at L = 32 (1.3e-8).
2. **CholeskyQR2 is as accurate as the pre-pivoted QR** at every Nwrap (thermalised and cold). It makes the decompose 16–26 % cheaper,
   which is 1–3 % of a sweep, because wraps and CGR dominate. Worth taking only together with other GPU-friendly stabilisation work,
   and only for the CUDA/host backends unless the determinant phase is re-derived for it.
3. **Adaptive Nwrap** pays more than any QR variant. The PG max grows with L (L = 24 → 32: ~10× at Nwrap 20), so a fixed
   Nwrap 20 needs a check at L ≥ 40. A run-time controller (raise Nwrap while PG max < 1e-9, lower it above 1e-8) is the next
   step, if wanted.

## Caveats

- Short runs (4 sweeps from one warm-up per L, one seed): PG max is an extreme value. Longer runs and more seeds can raise it.
- One model (production Hubbard). Other models (larger U, larger beta·t) need their own scan.
- The decompose cost is a fit (intercept), not a direct per-kernel timer.
- The CholQR2 path is CUDA/host only (see above). In the `ALF_GPU_STABTIME` hook it fell back to QR twice per run on the hook's
  synthetic start states; it never fell back in a sweep.

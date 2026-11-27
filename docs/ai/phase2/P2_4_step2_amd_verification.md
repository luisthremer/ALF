# P2.4 step 2, feature `gpu-amd`: device verification on the MI210 (2026-11-26/27)

Machine `amd` (AMD Instinct MI210, gfx90a), gfortran 14.2 with amdgcn offload, ROCm core-7.14, gpu-port's patched libgomp GCN
plugin, `GCN_NUM_THREADS=4`, `OMP_TARGET_OFFLOAD=MANDATORY`. Build: `source configure.sh AMDGPU noMPI GPU=amd NO-INTERACTIVE && make`.
Scripts and raw outputs: `amd:~/p2_amd/` (`gate2.sh`, `gate2eval.py`, `cmp_scal.py`, `bign.sh`, `gpchk.sh`, `prodchk.sh`).
GPU=amd is **experimental while KNOWN_ISSUES R2 is open**.

## Gate 2: default path, GPU=amd against `ALF_GPU_DISABLE=1` of the same binary, same seeds

Cases: Hubbard (Start input), Mz and SU(2), L = 6, 8, 12, 16 (and Mz L = 24, Ltau 1, submatrix), Ltau 0 and 1, delayed and submatrix
(N_delay 8), beta 5, Nwrap 10, NSweep 3 x 2 bins; builds -O3 (production flags) and -O2 (layout shaken), plus a tau lane
(`-O0 -g -fcheck=all`, L = 4, 6, 8, Mz and SU(2), Ltau 1). 72 pairs. Required: both runs rc 0, GPU path ACTIVE, `confout`
identical, no NaN/Inf, and the tolerance rule below.

**Tolerance rule** (manager decision (A), 2026-11-27). Every output file: |x - y| / max(1, |y|) <= tol over all numbers, with
tol = max(1e-11, 100 x Precision Green mean of the CPU run) (the rule of `gpu_e2e.py`). Exception, calibrated by a null comparison:
the imaginary parts of the scalar observables (`*_scal`, zero up to round-off) only need |Im x - Im y| <= 1e-11 x (sites / 16).
Real parts and every other file stay under the unchanged rule.

Null comparison behind the exception (SU(2) 16x16, Ltau 0, same seeds, `confout` identical everywhere):

| Pair | max abs Im diff, scalars | max rel Re diff, scalars | other files, max rel |
|---|---|---|---|
| GPU=amd vs CPU (-O3) | 1.03e-11 | 5.1e-15 | 6.7e-12 |
| null (i): CPU delayed vs CPU submatrix (same binary) | 2.4e-13 | 8.5e-16 | 1.7e-13 |
| null (ii): CPU -O3 (GPU=amd binary) vs GPU=none -O2 | 1.56e-10 | 6.8e-14 | 1.1e-10 |

The manager chose null (ii) (another compile and library, as GPU vs CPU) as the reference.

**Result: 70 of 72 PASS.** The two FAILs are the same case, Mz 24x24, Ltau 1, submatrix, at -O3 and -O2 (identical
signature): SpinXY_eq and SpinT_eq 1.8e-8 against tol 4.9e-9; everything else in these runs passes. Resolution of this case:

| Run, same input and seeds (all `confout` identical to the CPU run) | SpinXY_eq | SpinT_eq | Precision Green max |
|---|---|---|---|
| GPU=amd (device path) | 1.83e-8 | 1.81e-8 | 1.1e-6 |
| GPU=host (the same device algorithm on the CPU, melmac) | 6.6e-9 | 6.0e-9 | 4.4e-7 |
| CPU with `-DQR_PREPIVOT` (pre-pivoted QR on the CPU path) | 4.3e-9 | 3.3e-9 | 7.4e-7 |
| CPU GPU=none -O2 (compiler null) | 1.75e-10 | 1.46e-10 | 4.8e-7 |
| gpu-port 078d806b GPU=amd vs its own CPU path (Ltau 1, 10 sweeps x 2 bins) | 4.0e-8 | 2.1e-8 | — |

The device stabilisation's column-norm pre-pivoted QR (instead of ZGEQP3; P2.4 step 1d) moves these O(1) SpinXY/SpinT correlation
entries of the Mz chain at L = 24 by 4e-9 .. 7e-9 already on the CPU. GPU=amd is within 2.8x of GPU=host. Reading: a property of
the pivoting algorithm (all vendors), not an R2 symptom. The e2e tolerance is too tight for Mz at L = 24; a size-scaled rule for
Mz correlations is a step-3 decision.

Final tree (rebased on `Code-modularity-GPU` 75b9c881): 5 gate cases (Mz/SU(2), L = 8, 12, 16, Ltau 0/1, delayed/submatrix) PASS;
`55-linalg-hip` (both QR paths, mixed presence stops), `56-hipk-stab` (both QR paths), `53-gpu-obs` (sq1, sq2s, sq2h, bil2),
`52-gpu-stab` self-test in modes ALF_GPU_CGR_LU = 0, 1, 2 PASS; LU mode tau files within 6e-15, check mode 120 CGR calls, max
|tracked - LU| 1.1e-15. GPU=none 16/16 scenarios byte-identical to 75b9c881 (melmac); preprocessed `linalg_mod`, `gpu_*_mod`
identical to the base for GPU=none/host and for `-DALF_LA_CUDA` / `-DALF_LA_ONEMKL`.

## Large N (S1 tracked phase, hipSOLVER QR path, R1 guard)

Hubbard Mz, Ltrot 64 (beta 6.4), Nwrap 10, Ltau 0, Symm F, submatrix N_delay 32, 1 sweep:

| L (N) | rc | acceptance | Precision Phase max | Precision Green mean / max | s per sweep | host RSS |
|---|---|---|---|---|---|---|
| 56 (3136) | 0 | 0.411 | 0 | 1.8e-10 / 1.1e-6 | 328 | 13 GB |
| 60 (3600) | 0 | 0.413 | 0 | 3.0e-10 / 1.8e-6 | 775 | 17 GB |
| 64 (4096) | not run (stopped to keep amd under 85 % RAM; rerun alone) | | | | | |

The R1 guard fired at the start-up upload (identity U) in both runs, as expected.

## gpu-port and production (for the production stream)

gpu-port 078d806b, GPU=amd vs its own CPU path, same seeds, Ltau 1, NSweep 10 x 2 bins, `confout` identical everywhere:
Mz 4x4 / 6x6 / 12x12 max rel 6.6e-10 / 1.6e-10 / 3.6e-10; SU(2) 4x4 / 6x6 / 12x12 7e-13 / 4e-13 / 5e-12; Mz 24x24 see above.
Precision tau max is of the same order on the GPU and the CPU run (Mz 12x12: 2.0e-7 / 1.8e-7).
Production code (tag prod/2026-11-11-amd-L36 = b3923afe, built noMPI/no HDF5) on the production input L36_U5.25_h1_4w, 2 sweeps:
GPU vs CPU Den_eq 2.0e-9, scalar Im 3.3e-9; null CPU -O2 vs CPU -O3 Den_eq 2.2e-8, scalar Im 3.6e-8. Clean.

## Precision Green, HOPD (G2 question)

6x6, Checkerboard T, Nwrap 10, 4 seeds, same chain: Precision Green mean GPU/CPU 1.02 .. 1.28 (mean 1.12) with `ALF_GPU_HOPD=0`
(default), 1.18 .. 1.25 (mean 1.22) with `ALF_GPU_HOPD=1`. Nwrap 25 is unstable on this input for the CPU path as well (rc 3).

## Not done in this feature

Speed against gpu-port at L = 24/32 (scripts ready: `amd:~/p2_amd/spd.sh`, gpu-port build `amd:~/p2_amd/gp`); L = 64; the root
cause of R2 (GCC 15 offload is not available on amd).

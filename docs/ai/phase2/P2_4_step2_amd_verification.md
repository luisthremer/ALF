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
| 64 (4096) | 0 | 0.412 | 0 | 3.0e-10 / 2.3e-6 | 364 | 22 GB (run alone) |

The R1 guard fired at the start-up upload (identity U) in every run, as expected.

## gpu-port and production (for the production stream)

gpu-port 078d806b, GPU=amd vs its own CPU path, same seeds, Ltau 1, NSweep 10 x 2 bins, `confout` identical everywhere:
Mz 4x4 / 6x6 / 12x12 max rel 6.6e-10 / 1.6e-10 / 3.6e-10; SU(2) 4x4 / 6x6 / 12x12 7e-13 / 4e-13 / 5e-12; Mz 24x24 see above.
Precision tau max is of the same order on the GPU and the CPU run (Mz 12x12: 2.0e-7 / 1.8e-7).
Production code (tag prod/2026-11-11-amd-L36 = b3923afe, built noMPI/no HDF5) on the production input L36_U5.25_h1_4w, 2 sweeps:
GPU vs CPU Den_eq 2.0e-9, scalar Im 3.3e-9; null CPU -O2 vs CPU -O3 Den_eq 2.2e-8, scalar Im 3.6e-8. Clean.

## Precision Green, HOPD (G2 question)

6x6, Checkerboard T, Nwrap 10, 4 seeds, same chain: Precision Green mean GPU/CPU 1.02 .. 1.28 (mean 1.12) with `ALF_GPU_HOPD=0`
(default), 1.18 .. 1.25 (mean 1.22) with `ALF_GPU_HOPD=1`. Nwrap 25 is unstable on this input for the CPU path as well (rc 3).

## Speed against gpu-port (follow-up, final tree on 75b9c881+, measured 2026-11-27)

MI210 alone (one run at a time), Hubbard Mz, beta 5, dtau 0.1, Nwrap 10, Ltau 0, Checkerboard T, submatrix N_delay 32,
ham_h0 = 0 unless stated, NBin 1. Seconds per sweep in steady state = (wall(4 sweeps) - wall(2 sweeps)) / 2, same input for both
codes (`amd:~/p2_amd/spd2.sh`, results `amd:~/p2_amd/spd/results_symmF.txt`). `Symm = F` for both: with `Symm = T` the start-up
check `Op_V_is_symmetric` of the CM line (KNOWN_ISSUES P1; its fast version is LEVERAGES row 1, in review) dominates the wall time
(CM L = 24: ~230 s start-up; L = 32: 2377 s for 2 sweeps against 70 s on gpu-port).

| L (N) | gpu-port 078d806b (Hd on, its default) | CM, `ALF_GPU_HOPD=1` (Hd) | CM default (`HOPD=0`, per-factor, G2) | CM default, h0 = 1 |
|---|---|---|---|---|
| 24 (576) | 1.85 | 1.10 (1.7x faster) | 5.36 (2.9x slower) | 6.08 |
| 32 (1024) | 18.2 | 16.9 | 14.2 | 16.9 |

At L = 32 (Ndim 1024 > `HOPD_NMAX` = 800) the dense product is not used anyway; the two CM columns measure the same path and their
spread (14.2 .. 16.9 s) is the run-to-run noise: CM is 1.1-1.3x faster than gpu-port there. At L = 24 the per-factor hopping
product (the G2 default) costs 4.9x on amdgcn: every checkerboard group is one libgomp target launch (~150 us each). For AMD
production at Ndim <= 800, `ALF_GPU_HOPD=1` (Precision Green 1.2x the CPU value on this input, section above) or a HIP kernel for the
group product (next step) is the choice.

## HIP kernel for the checkerboard group product (feature `gpu-amd-hop`, 2026-11-27)

`k_cb_left_g` / `k_cb_right_g` of `gpu_wrap_mod` run as HIP kernels on GPU=amd (same arithmetic as `el_left_d` / `el_right_d`,
one HIP launch per group instead of one libgomp target region). Checks on the MI210:
- `51-gpu-wrap` self-test (device hopping products against the host `Hop_mod`, 1e-14): PASSED for checkerboard, dense and
  `ALF_GPU_HOPD=1`; faults `hop` and `nan` seen.
- Reduced gate, GPU vs `ALF_GPU_DISABLE=1` (Mz/SU(2), L = 8, 12, 16, 24, Ltau 0/1, delayed/submatrix, default HOPD = 0; rule of
  decision (A)): 28 of 32 PASS; the 4 others are Mz 24x24, SpinXY/SpinT_eq 1.8e-8, the accepted pre-pivoted-QR effect
  (KNOWN_ISSUES G7, class C), with the same values as before the change.
- Precision Green (the G2 gate): the GPU Precision Green mean is **bit-identical** before and after the change in all 25 cases
  run in both gates; the HIP kernel reproduces the OpenMP kernel exactly. GPU/CPU ratio of the Precision Green mean over the 32
  cases: median 1.25, max 1.83 (SU(2) 24x24, a case that had no earlier run).

Speed, same setup as the previous section (MI210 alone, `Symm = F`, s per sweep in steady state; `amd:~/p2_amd/spd/results_hop.txt`):

| L (N) | gpu-port 078d806b | CM default `HOPD=0`, before -> after | CM `HOPD=1` |
|---|---|---|---|
| 24 (576) | 1.78 | 5.36 -> **1.30** (4.1x; 1.4x faster than gpu-port) | 0.98 |
| 32 (1024) | 17.7 | 14.2 -> **5.1** (2.8x; 3.4x faster than gpu-port) | 3.1 |

At L = 32 the dense product is not used (Ndim > 800), so the last two columns run the same path: their spread (3.1 .. 5.1 s) is
run-to-run noise. The per-factor default (CPU-level precision, G2) is now faster than gpu-port at both sizes.

## P2.6 models on GPU=amd (feature `gpu-amd-p26`, 2026-11-27)

`p26-models` (non-diagonal, rank > 1 and type-4 operators: Kondo, tV, Spin_Peierls) was refused on ALF_LA_HIP builds. Found and
fixed on the way: R3 (HIP slice kernel read the old `V_meta` layout, every model; feature `gpu-amd-r3`, merged first) and R2
symptom 3 (libgomp zero-length map in `k_dly_gpp` / `k_dly_wrap`, delayed path of the P2.6 models). rocprofv3 on the MI210
(kernel time, OpenMP P2.6 kernels as merged): `k_nd_rows` + `k_nd_cols` 66-82 % (Kondo 12.0 of 18.0 s, tV 17.0 of 21.0 s,
Spin_Peierls 24.2 of 30.7 s), then `k_rk_gpp` / `k_rk_prep` / `k_rk_upd` (standard) and `k_gen_block` / `k_gen_accept`
(delayed); each OpenMP launch costs 120-230 us, a HIP launch ~4 us. HIP versions (same arithmetic) of `k_dly_gpp`,
`k_dly_wrap`, `k_nd_rows`, `k_nd_cols`, `k_rk_gpp` and the rank-r accept (`k_rk_prep` + `k_rk_upd`); then the refusal is removed.

Gate, GPU=amd vs `ALF_GPU_DISABLE=1` of the same binary, same seeds (`p25_run.py`, rule of `gpu_e2e.py`; `confout` identical):
- `cases_p26.py`, 16 cases (Kondo, tV, Spin_Peierls x standard / delayed 4, 16 / submatrix 4, tV tau, Hubbard 4x4 Mz
  standard / delayed / submatrix): 16/16 OK on the final tree (rebased on `Code-modularity-GPU` 13dcd1e3; preprocessed
  `alf_hipk_mod`, `gpu_lowrank_mod`, `gpu_wrap_mod`, `gpu_update_mod`, `gpu_state_mod` identical to the base for GPU=none, host,
  `-DALF_LA_CUDA`, `-DALF_LA_ONEMKL`). Max scaled diff: Kondo 2.5e-12, tV 1.2e-13, Spin_Peierls 2.1e-10
  (tol >= 2.6e-10), Hubbard 2.7e-10 (tol 6.1e-9). Also 16/16 on the intermediate trees (OpenMP P2.6 kernels with the two fixes;
  plus HIP nd).
- `cases_h100.py` big cases (Kondo 8x8 bilayer, tV 12x12, Spin_Peierls 12x12; standard / delayed 16, 32 / submatrix 16;
  NSweep 10, NBin 2), final tree before the rebase, `confout` identical everywhere, rc 0:
  - tV standard 1.0e-11, delayed 32 1.4e-11, submatrix 16 7.4e-13: OK. tV delayed 16: Pot_scal 1.67e-11 against tol 1.66e-11,
    while the CPU-only null (same binary, submatrix 16 vs delayed 16) gives 1.88e-11: tolerance-edge round-off (G9, class C).
  - Kondo, all four schemes: Dimer_tau 1.48e-11 against the floor 1e-11 (G9, as H100 1.26e-11).
  - Spin_Peierls, all four schemes: FAIL, observables up to 1.3e-7 (Part_scal, Dimer_eq, Pot_scal, SpinZ_eq), Precision Green
    mean 2.1e-10 vs 8.2e-12 on the CPU. Not AMD and not this feature: the build with the merged OpenMP P2.6 kernels shows the
    same, and so does GPU=host on melmac (no offload; Precision Green 2.6e-10 vs 9.0e-12, Part_scal 1.8e-7): known G10
    (KNOWN_ISSUES, class B, P2.6 follow-up).

Time of the whole run, MI210 alone, `cases_p26.py` inputs (s; columns: OpenMP P2.6 kernels with the two fixes, + HIP nd,
+ HIP rank-r = final; CPU = `ALF_GPU_DISABLE=1`):

| Case | OpenMP | + HIP nd | final | CPU |
|---|---|---|---|---|
| Kondo 4x4 bilayer, standard | 27.6 | 14.4 | 11.0 | 0.6 |
| Kondo, delayed 16 | 39.1 | 12.0 | 12.6 | 0.6 |
| Kondo, submatrix 4 | 28.8 | 15.0 | 11.7 | 0.7 |
| tV 4x4, standard | 30.9 | 11.2 | 6.7 | 0.3 |
| Spin_Peierls 6x6, standard | 43.5 | 16.3 | 11.4 | 0.8 |
| Spin_Peierls, delayed 16 | 66.7 | 15.4 | 14.9 | 1.1 |

2.5-4.5x faster than the OpenMP kernels. These small cases stay launch-bound and far behind the CPU, as on every vendor (H100
4.4 s vs 1.7 s); the crossover size per model is P2.6.2.

Follow-up idea (not in this feature): the remaining OpenMP launches of ~100-200 us each on GPU=amd: `k_site_wrap`,
`k_upd_prep` / `k_upd_geru` (standard update of diagonal operators), `k_gen_block` / `k_gen_accept` / `k_dly_accept`
(delayed P2.6), `k_vrun`, `k_wrap_g0`, and the small `target update` copies around the update (`rk_p`, `rk_w`, `rk_g`).

## Not done in this feature

The root cause of R2 (GCC 15 offload is not available on amd).

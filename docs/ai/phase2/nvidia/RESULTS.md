# P2.4 step 2 `gpu-nvidia`: device results (Helma, NVIDIA H100)

Measured 2026-11-26/27 by the gpu-nvidia session. Candidate trees from `git archive` of the canvas: `0b84b0e5` (= the
feature rebased onto Code-modularity-GPU 4bd27bc5, with the G5 fix; the later commits add the e2e lane tolerance and docs).
Earlier canary/linalg numbers from `c1fbafce` (same CUDA code without the G5 fix). Toolchain: nvhpc 25.7 (nvfortran),
`-mp=gpu -gpu=cc90,cuda12.9`, driver 595.91.07, H100 94 GB, partition `preempt`. Scripts in this directory (`dev.py` drives
the runs; `hbuild.sh`, `hrefs.sh`, `hstrict.sh` build; `*.sbatch` the jobs; `melmac_check.sh`, `pp2.sh` the GPU=none/host
checks). H100 use for the whole step: about 3.1 GPU-hours (sacct).

## Build (Helma login node)

```bash
source configure.sh HELMAGPU noMPI GPU=nvidia NO-INTERACTIVE && make -j16 all     # also with LOG, and MPI
source configure.sh HELMAGPU noMPI GPU=none NO-INTERACTIVE && make -j16 all       # nvfortran CPU reference
module load cmake/3.30.5     # testsuite (after configure)
```
Run environment: `OMP_TARGET_OFFLOAD=MANDATORY`. With `-Minfo=mp` every target region of the whole build generates a GPU
kernel and no "Accelerator restriction" is reported. nvfortran pitfall found on the way: a kernel that uses `merge()` on
complex values is dropped silently (CDSELECT, run-time "Failed to find device function"): use if/else.

## Correctness

| Check | Job | Result |
|---|---|---|
| 56-cuda-order: stream-order canary (RAW/WAR/WAW of cuBLAS calls against synchronous target regions; 8292 device QRs through the info slots) | 1463534 | PASS: RAW 2.8e-14, WAR 2.5e-14, WAW 0. Negative control (`-DALF_LA_TEST_BADSTREAM`, created stream): FAIL (RAW 1.07, WAR 24, WAW 0.44) as required |
| 55-cuda-linalg: device branch of every wrapper vs host LAPACK/BLAS, N = 200 | 1463616 | max 1.5e-14 (zgemm 4e-16, ztrsm 7e-16, geqrf/unmqr 3e-15, getrf+getri 1.5e-14, qrpp 1.1e-15); qrpp permutation equal to the host one with tied norms, bitwise re-run, NaN -> -101; S1 reflector-determinant phase vs LU 9e-15; mixed operands stop the run |
| Full ctest on the H100, default and STABLOG | 1465022, 1465739 | 144/146 each in 1465022; the two failures were only the check-mode byte-for-byte lane of `gpu-state-e2e` and `cgr-phase-e2e`; with the lane of the testsuite commit both pass in both builds (1465739) |
| gpu_e2e (in ctest): GPU path vs `ALF_GPU_DISABLE=1`, 17 cases, faults | 1465022 | confout byte-identical in every case; SU(2) <= 4e-14, Mz <= 2.6e-11 (Ener/Kin/Pot); faults upd, obs, tau-prop, tau-blocks seen; stab-perm, stab-det stop |
| Phase modes (gpu_e2e, cgr_phase_e2e) | 1465022 | tracked vs LU 8.6e-15 (SU(2), after the G5 fix; 1.6e-6 before), check mode phase diff <= 1.4e-15 |
| Trajectory GPU vs CPU path, same binary and seeds, Beta 10, 10 sweeps x 4 bins (12x12), 20 x 4 (6x6) (`dev.py traj`) | 1465023 | confout byte-identical in all 7 cases. Max scaled deviation: Mz 12x12 9.7e-9 (SpinXY_eq; standard, delayed 16, submatrix 16, Ltau 0 the same), SU(2) 12x12 1.4e-11, Mz 6x6 1.3e-9, SU(2) 6x6 Nwrap 25 2.6e-9. GPU=none (nvfortran) == CPU path of the CUDA binary byte for byte in all 7 |
| Same, GPU=host build on melmac (same parameters) | melmac | Mz 12x12 7.1e-8, SU(2) 1.1e-11: the deviation comes from the device-path algorithms of the common code (pre-pivoted QR, ZUNMQR), not from CUDA; CUDA is closer to the CPU path than GPU=host |
| Same, CUDA build without `-Mfprelaxed -fast` | 1465749 | Mz 12x12 9.7e-9, SU(2) 9.2e-12: the flags are not the cause |
| Replay (same binary, same seeds, twice) | 1465023 | not bitwise on the device (Mz 12x12: 6.7e-15 scaled, confout identical); bitwise on GPU=host. Property of the device run (reductions); the e2e check-mode lane has 1e-13 on device backends |
| MPI, 2 ranks on 2 H100 (rank mod number of devices) | 1465024 | each rank on its own GPU (nvidia-smi: bus 06 and A6); confout identical to the CPU path, max 4.3e-8 (SpinXY_eq), tau 2.6e-9 |
| Large-N phase test: L 56/60/64, Beta 6.4, U 5, Ltau 0, Symm F, 2 sweeps, `ALF_GPU_CGR_LU=2` | 1465023 | acceptance 0.422 / 0.425 / 0.423, no NaN, 56 CGR calls compared, max |tracked - LU| 0; wall 400 / 573 / 1109 s |

Precision Green, mean (Beta 10, 5 sweeps x 2 bins, SU(2); the G2 question):

| L, Nwrap | CPU path | CUDA (`HOPD` off, default) | CUDA `ALF_GPU_HOPD=1` |
|---|---|---|---|
| 6, 10 | 1.1e-14 | 1.6e-14 | 2.1e-14 |
| 6, 25 | 7.2e-12 | 8.9e-12 | 1.8e-11 |
| 12, 10 | 1.7e-14 | 4.6e-14 | 5.2e-14 |
| 12, 25 | 1.5e-11 | 2.1e-11 | 3.2e-11 |
| 24, 10 | 9.1e-14 | 4.0e-13 | 4.0e-13 |
| 24, 25 | 8.9e-11 | 1.5e-10 | 2.5e-10 |

The dense hopping product with cuBLAS ZGEMM raises Precision Green further, as with the host BLAS (G2): keeping it opt-in is
right for CUDA too.

## Speed (job 1465739; CPU time of the sweeps from the info file / NSweep)

L x L, Beta 32 (Ltrot 320), Nwrap 10, submatrix N_delay 32, Ltau 0, ham_h0 0, Symm F (Symm T: the start-up check P1
dominates), 3 sweeps, two repetitions each (spread < 2 %). Same parameters on all codes.

| L | this feature (tracked det) | this feature, LU (`ALF_GPU_CGR_LU=1`) | gpu-port 078d806b (LU) | rd/cgalf-s1 (tracked) |
|---|---|---|---|---|
| 24 | 3.60 s | 3.85 s | 3.79 s | 3.57 s |
| 32 | 6.91 s | 7.38 s | 7.30 s | 6.92 s |

The CM CUDA path is 5 % faster than gpu-port and as fast as rd/cgalf-s1; S1 saves 6 % (L 24 and 32) against LU on the same
binary. `ALF_GPU_HOPD=1` at L 24: 3.57 s (gpu-port uses the dense hopping product by default for N <= 800).
Start-up is slower than on gpu-port (wall time minus sweep time: L 24 ~8 s vs ~2 s, L 32 ~74 s vs ~2 s); not investigated
(candidate: host initialisation of the UDV stack before the upload, common code). Follow-up, not a CUDA issue.

## Known issues found

- **G5** (KNOWN_ISSUES): nvfortran 25.7 drops the imaginary part of a complex host scalar read in a serial `!$omp target`
  region; fixed in this feature (`map(to: ph0)` in `k_cgr_phase`). gpu-port and the production branches are not affected.
- `gpu_e2e.py` still says the `ALF_GPU_FAULT=nsig` case (stale device copy of the field) is "to be added for the vendor
  backends"; not added here (common test, all vendors).

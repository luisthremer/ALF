# P2.4 step 2 `gpu-intel`: device results (SuperMUC-NG Phase 2, Intel Max 1550)

Measured 2026-11-26/27 by the gpu-intel session on the canvas `cm/feat/gpu-intel` (code of ae9449fb; the later commits add the
SUPERMUCNG2 machine block and docs only; job 566646 repeats ctest and the gate on the final canvas rebased onto a9f25efa). Scripts in this directory; outputs in `$SCRATCH/p2_intel/<job>` on SuperMUC (user di38zuw).

## Build (login node)

```bash
export ALF_INTEL_AOT=pvc        # ahead-of-time for the Max 1550 (optional; JIT otherwise)
source configure.sh SUPERMUCNG2 noMPI GPU=intel NO-INTERACTIVE && make -j16 all
source configure.sh SUPERMUCNG2 MPI GPU=intel HDF5 NO-INTERACTIVE && make -j16      # production style
source configure.sh SUPERMUCNG2 MPI GPU=none HDF5 NO-INTERACTIVE && make -j16       # CPU build on Phase 2
```
All four variants built (ifx 2025.3.0, Intel MPI 2021.15, MKL 2025.3.0, LRZ hdf5/1.14.5-intel25). The testsuite needs a JIT
build: CMake drops the repeated `-Xopenmp-target-backend` of the AOT flags (`job_ctest.sbatch`).

Run environment (`env_smuc2.sh`): `OMP_TARGET_OFFLOAD=MANDATORY`, `ZE_FLAT_DEVICE_HIERARCHY=FLAT` (a tile is a device),
`LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576`, one tile per run by `ZE_AFFINITY_MASK=<tile>`.

## Correctness

| Check | Job | Result |
|---|---|---|
| ctest on the device (13-linalg incl. the device cases, 50-54 GPU tests, gpu-state-e2e) | 566623 | 53/53 passed; e2e: confout byte-identical active vs disabled in every case, check mode phase diff <= 1.8e-15 |
| Trajectory GPU=intel vs GPU=none (same ifx), Start parameters, 20 sweeps x 5 bins (`traj.py`) | 566621 | L6 Mz=F 3.7e-14, L6 Mz=T 1.2e-11, L12 Mz=F 1.1e-14, L12 Mz=T 6.9e-11 (worst relative deviation of all observable files and confout); confout and acceptance identical; `ALF_GPU_CGR_LU=2` phase diff <= 5.4e-15 |
| Gate after the AMD offload finding (`gate.py`): GPU path vs `ALF_GPU_DISABLE=1` and the check-mode lane, same binary and seeds; Mz and SU(2), L 6/8/12/16, Ltau 0/1, delayed/submatrix (N_delay 16); STABLOG and an -O1 build at Ltau 1; L24 once per model | 566631 | 66/66 PASS: no NaN, all ACTIVE rc 0; worst SU(2) <= 6.5e-14, Mz <= 3.6e-11 (L24 Mz 1.2e-10) against tol = 100 x Precision Green mean of the CPU path; check-mode phase diff <= 1.4e-14 |
| Negative controls of the gate | 566636 | `ALF_GPU_FAULT=upd` (standard scheme): run stops (rc 3); `ALF_GPU_FAULT=obs` (submatrix): deviation 1e-3 -> both FAIL as required. (`upd` with delayed/submatrix has no effect: KNOWN_ISSUES G6) |
| Re-check of the final canvas after the rebase onto Code-modularity-GPU a9f25efa (`job_final.sbatch`): device ctest + the gate with the default build | 566646 | ctest 54/54; gate 34/34 PASS, both negative controls FAIL |
| Determinant phase at N >~ 2950: L 56/60/64, Beta 6.4, U 5, h0 0, submatrix N_delay 32, 1 sweep, Symm = F (Symm = T stalls in Op_V_is_symmetric, KNOWN_ISSUES P1) | 566642 | acceptance 0.3145 / 0.3139 / 0.3146, no NaN; tracked det(U) and LU path identical |

## Precision Green (reference set: N_FL 2, Beta = L, U 5.3, h0 1, Nwrap 10, submatrix N_delay 128; 15 sweeps)

| L | mean | max | gpu-port doc (mean / max) |
|---|---|---|---|
| 24 | 1.0e-14 | 1.1e-11 | 1.7e-14 / 1.2e-11 |
| 28 | 2.6e-14 | 9.4e-11 | 4.6e-14 / 1.1e-10 |
| 32 | 6.3e-14 | 4.5e-10 | 1.4e-13 / 5.3e-10 |

## Speed (job 566626)

Layout: one node, 8 runs at the same time (one per tile, `ZE_AFFINITY_MASK` = tile), no `ZEX_NUMBER_OF_CCS`, `OMP_NUM_THREADS=14`,
no CPU pinning, AOT builds, staging buffer 1 MB, noMPI binaries. s/sweep = (CPU time of NBin = 3 - NBin = 1) / 10 (NSweep 5),
wall-time differences agree. Parameters: the reference production set above with N_delay 128.

| L | h0 | this line | gpu-port tip (078d806b, same layout) |
|---|---|---|---|
| 24 | 0 | 20.6 | 20.1 |
| 24 | 1 | 20.2 | 40.7 (gpu-port measures on the host when h0 /= 0; device observables with h0 came with step 1e) |
| 24 | 1, `ALF_GPU_CGR_LU=1` (LU phase) | 21.3 (tracked det(U) saves 5 %) | |
| 28 | 1 | 33.8 | |
| 32 | 1 | 51.0 | did not finish in the job |

Both binaries are about 2.5x slower than gpu-port's documented 8.25 s/sweep at L = 24 (8 ranks x 14 threads, `I_MPI_PIN_DOMAIN=auto`,
N_delay 128, MPI binary). Not a code difference (the two binaries agree at h0 = 0); to be traced in S-c (pinning, CCS, the
Intel switches).

## Mz = T deviation (investigation, P2.4 plan section 6)

Hypothesis H0: the Mz = T deviation is the rounding sensitivity of the CPU result itself, not a GPU error. Same seeds, Start parameters
(20 sweeps x 5 bins), worst relative deviation over all observable files and confout:

| comparison | L6 Mz=T | L12 Mz=T | L6 Mz=F |
|---|---|---|---|
| CPU ifx, `MKL_CBWR=AUTO` vs `COMPATIBLE` (rounding only, same algorithm; melmac) | 1.7e-12 | 5.0e-11 | 9.3e-15 |
| CPU gfortran default vs `-DQR_PREPIVOT` (the pivot rule of the device; melmac) | 9.8e-12 | 2.4e-11 | 2.9e-13 |
| Max 1550 vs CPU ifx (job 566621) | 1.2e-11 | 6.9e-11 | 3.7e-14 |
| CPU Precision Green mean | 3.2e-11 | 5.1e-11 .. 6.4e-11 | 1.2e-14 |

The device deviation is at the size of the CPU's own rounding noise (a pure change of the MKL code path already moves L12 Mz = T by
5e-11) and of the order of the Precision Green mean. H0 holds; no code error, no bug class. The tolerance rule of the checks is 100 x
the Precision Green mean of the CPU path (as in `gpu_e2e.py`). gpu-port's larger numbers (1.36e-10 at L6, 1.21e-9 at L12) come from
other parameters and its V7 switch set; its own 20-seed run reached 1.16e-9 with oneMKL alone. claude-bughunter was not needed.
(ifx vs gfortran could not be compared: the two compilers write confout in different formats.)

## Budget

About 24 GPU-h of pn93la (estimate from the job lengths; sacct was unreachable): builds on the login node; jobs 566620-566646.

## S-c: the Intel switches (feature `gpu-intel-speed`, 2026-11-27)

Switches (all oneMKL only): `ALFSH_ASYNCBLK` (asynchronous SYCL block path of the submatrix update), `ALFSH_BLKFUSE=0..5` (fused
block steps; >= 2 needs N_delay <= 64), `ALFSH_BLKPF` (prefetch of the next block, needs BLKFUSE >= 2), `ALFSH_OBSCORR` (correlation
sums as a SYCL kernel). The info file names the active values, e.g.
`GPU Intel block path (submatrix): ALFSH_ASYNCBLK=1, ALFSH_BLKFUSE=5, ALFSH_BLKPF=1` and `GPU observables ...: ALFSH_OBSCORR=1`.

| Check | Job | Result |
|---|---|---|
| Gate with ASYNCBLK=1 BLKFUSE=5 BLKPF=1 OBSCORR=1 in the GPU lanes (default and STABLOG builds; Mz/SU(2), L 6-16 + 24, Ltau 0/1, delayed/submatrix, CGR_LU=2 lane) | 566761 | 50/50 PASS, negative controls FAIL; worst deviations as with the switches off (Mz L24 1.2e-10, STABLOG Mz 3.6e-11); prefetch 30 000 hits, 0 misses |
| Fusion levels 1-5 with and without prefetch (intel box, B580, smoke) | - | all PASS (8x8 Mz T/F, Ltau 0/1) |
| GPU=none 16 scenarios vs bb4b089c, GPU=host ctest (melmac) | - | byte-identical; 143/143 |

Speed (job 566762): reference set (Beta = L, N_FL 2, U 5.3, h0 1, Nwrap 10, submatrix N_delay 32), one run per tile, 8 runs on the
node, s/sweep = (CPU time NBin 2 - NBin 1) / 3 sweeps:

| L | this line, switches off | this line, switches on | intel-batched tip, same switches | gpu-port tip |
|---|---|---|---|---|
| 24 | 32.3 | 30.0 (-7 %) | 28.5 | 46.8 |
| 32 | 77.8 | 73.7 (-5 %) | 77.1 | 156.0 |

Layout question (job 566660, L = 24, N_delay 32): one run alone with 1 / 12 pinned / 14 threads: 34.9 / 30.6 / 30.8 s; 8 runs at
once, 12 pinned threads each: 29.6-32.6 s. Pinning, thread count and node load change nothing. In the same harness this line and
intel-batched agree (above), so the factor to intel-batched's documented 3.75 s/sweep (its internal TEAM_TIMERS, mpiexec launch)
is a difference of measurement, not of code; not pursued further (open for whoever needs absolute Intel throughput).

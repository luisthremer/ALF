# Parking list (needs more than 1 node, or more than a short test job)

Nothing here has been submitted. Node-hour estimates use my measured, tuned rates (C: 24 s/sweep/chain at L=28, 42 s at L=32, 32 chains per node; noise +-15-20 %).
1 node-h = 4 GPU-h. Script names refer to `Scripts_and_Parameters_files/` and `docs/gpu/sessions/supermuc_large_prompt/jobs/`.

## P1. Production L=28 and L=32 (the task of the user's 2nd message)
- What: 10 U values (5.0, 5.1, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5, 5.6) per size, point-defect Hubbard model, beta=L, N_FL=2, Mz=F, PARALLEL_PARAMS (one set per U),
  32 walkers per node (4 per tile with `ZEX_NUMBER_OF_CCS`), `Nwrap=20`, `N_delay=64`, 1000 warm-up sweeps per chain, `NSweep=50` per bin.
- Why: the user's target (about 80k sweeps per U after warm-up).
- Command (per size; PP binary `build_pp/Prog/ALF.out` already built and verified):
  ```
  python3 Start_PP_GPU/make_pp_dirs.py RUN L28 5.0,5.1,5.2,5.25,5.3,5.35,5.4,5.45,5.5,5.6 16 50 100000 64 46.0 <seed-offset>
  # copy JobfileSuperMUC-NG-Phase2-GPU-PP.sh into RUN, fill ##...## (partition general, 5 nodes, time 47:30:00, OMPT=4), then
  sbatch JobfileSuperMUC-NG-Phase2-GPU-PP.sh      # restart after the limit: sbatch again (out_to_in_pp.sh moves confout -> confin)
  ```
- Estimated node-hours (E):
  | layout | nodes | wall | node-h |
  |---|---|---|---|
  | L=28, 16 walkers/U, 6000 net sweeps/chain | 5 | ~40 h | ~200 |
  | L=32, 16 walkers/U, 6000 net sweeps/chain | 5 | ~67 h (2 segments) | ~335 |
  | L=32, 32 walkers/U, ~3100 net sweeps/chain (fits 46 h, ~100k/U) | 10 | 46 h | ~460 |
- Depends on: user decision on layout and GPU-h; `/hppfs/work/pn93la` or a scratch retention plan for the data (`data.h5` ~180 MB per set at L=24, more at L=32).
- Open check before launching: a clean marginal timing of the final configuration (section 9 of FINDINGS.md) and equilibration length (1000 sweeps at L=24 was the reference; L=28/32 are not verified).

## P2. Full-length L=24, U=5.3 comparison with the A30 result
- What: 4 walkers, 100 bins x 50 sweeps (first 20 discarded), compare `Ener_scal` etc. with `E = -505.42 +- 0.40` (naive error over 77 bins, so compare with a binned error).
- Command: `Scripts_and_Parameters_files/JobfileSuperMUC-NG-Phase2-GPU.sh` with `ref_L24/parameters` (N_delay 32, Nwrap 10) on 1 node, 4 ranks, `--time=06:00:00`, `ZE_AFFINITY_MASK` per card. Runs on 1 node but > 1 h, hence parked.
- Estimate: 100 bins x 50 sweeps = 5000 sweeps per chain x 14.6 s = 20 h with 4 chains and no CCS (20 node-h); with 32 chains on one node (CCS, 17 s/sweep) the same wall time gives 8x the statistics. Depends on: user OK.

## P3. Two-node MPI question (analysis, reading the code; test not run)
- Communicated: per bin, `MPI_REDUCE` of the observable arrays to the group root (a few MB) and a few scalars; at start the seeds. Frequency: once per bin (here every 50 sweeps x ~25 s = ~20 min). Walkers are independent, no exchange under PARALLEL_PARAMS.
  So multi-node runs scale ~linearly and there is no reason for fewer, larger nodes except warm-up cost; more independent walkers per set = more nodes.
- Test plan (2 nodes, ~0.2 node-h): `mpiexec -n 64` L=24 with `ZEX_NUMBER_OF_CCS`, 4 sweeps, compare s/sweep/chain with 1 node (16.9 s expected) and `data.h5` identity with the 1-node 32-rank run for the first 32 ranks. Command: `sbatch --nodes=2 --ntasks-per-node=32 jobs/task0_scan.sh "64" ccs`-style job on the test partition.
- GPU-aware MPI is irrelevant here (no GPU buffers in MPI calls); `I_MPI_OFFLOAD` need not be set.

## P4. GPU-native QR/stabilisation (code work, then 1-node tests)
- Why: oneMKL `geqrf` at N=1024 costs 31-96 ms (host-assisted) against 0.45 ms for a ZGEMM; QR dominates the sweep and the GPU is idle half the time (VTune 47 % busy).
  A device-resident QR (blocked Householder in OpenMP target kernels or a CholeskyQR2 variant with a stability fall-back) could save an estimated 30-50 %.
- Test cost: < 5 node-h of 1-node jobs (t_timing, same-seed vs CPU, L=28/32 timing).
- Depends on: development time; correctness proof against the same-compiler CPU run (existing tools in `jobs/`).

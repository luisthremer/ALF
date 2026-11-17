# One-node reproducer, LRZ incident 273051 (i20r02c06s02, card 0 / tile 1)

**What it is:** ALF, a finite-temperature quantum Monte Carlo code (Fortran, OpenMP offload to the GPU via Level Zero, oneMKL), run on one node with 8 MPI ranks, one per tile. Rank r uses device r in flat hierarchy. The input is fixed: an L=44 Hubbard model with fixed seeds and 1 sweep.

**Running it:**
- From this directory: `sbatch [-w NODE] run_repro.sh [NREP]`. NREP is the number of repetitions, default 3.
- One repetition takes about 5 minutes on a healthy node (job 562002: 311-312 s).
- Toolchain: `env_intel_pvc.sh` loads intel/2025.3.0, intel-mkl/2025.3.0, intel-mpi/2021.15.0 and hdf5/1.14.5-intel25, and sets OMP_TARGET_OFFLOAD=MANDATORY.

**Expected results:**
- **On a healthy node:** every repetition ends with rc=0 and `Acceptance = 0.318327091942149`, bit-identical every time and on every node we tried. The log line reads `verdict=OK`.
- **On the faulty tile (original node, 2026-11-12):** all 4 of 4 runs aborted after about 3 minutes, with rc=3, on the rank bound to device 1. The message was `15.68 / 10.05 ... is exceeding the threshold of 10 for G difference!`. The code's own consistency check found a matrix that is wrong by O(1-10) per entry. Normal accuracy is about 1e-12. No driver or runtime error was reported. The log line reads `verdict=FAULT`.
- **Same node with tile 1 hidden** (`ZE_AFFINITY_MASK=0,2,3,4,5,6,7`): the run was clean, with the same acceptance as on other nodes.

**Files:**
- `ALF-pp-intel-5667acda.out`: the binary that was used.
- `src_5667acda.tar`: its source (git commit 5667acda of ALF with the GPU port).
- `parameters`, `seeds`: the input.
- `env_intel_pvc.sh`: the toolchain setup.
- `run_repro.sh`: the job script.
- `reference/`:
  - `info_healthy_i20r01c01s06.txt`: the result on a healthy node, job 560622.
  - `repro.562002.log`: this script on a healthy node, 2026-11-17 (both OK).
  - `run_faulty_i20r02c06s02.log` and `node.560620.log`: the abort on the faulty tile, job 560620.

**Lighter alternative, not yet shown on the faulty tile:**
- `stress/t_stress.F90` plus its build and run scripts.
- It repeats zgemm/zgetrf/zgeqrf/zunmqr/ztrsm at N=1024 on every device and checks each result bit for bit against the first.
- It catches an injected error on all 8 tiles, and showed 0 differences in 100,711 iterations on a healthy node.
- Its runs on the faulty node were cancelled before they started, so it has never been run on the faulty tile. Use the ALF reproducer above as the primary test.

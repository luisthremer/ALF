# ALF scaling benchmarks for the SuperMUC-NG Phase 2 grant proposal (2026-11)

Status 2026-11-27: version pinned, binary built. **No benchmark job submitted** (owner: "don't submit anything yet").
The plan below is a proposal and has not been approved yet.

## 1. Exact version (verified on SuperMUC, 2026-11-27)

| Item | Value |
|---|---|
| Production binary | `/hppfs/work/pn93la/di38zuw/prod_2026-11/bin/ALF-pp-intel-410c7b04.out`, md5 `d2dac44061835a957bec7402fb7db7d5` |
| Used by | `job.sh` of RUN_L24, RUN_L28, RUN_L32, RUN_L36, RUN_L24_U2, RUN_L24_U3, RUN_L28_U2, RUN_L32_U2, RUN_L32_nw5to10, RUN_L32_fast* (all in `/hppfs/work/pn93la/di38zuw/prod_2026-11/`) |
| Source | ALF commit **410c7b04** ("Unstable-matrix guard: ALF_UNSTABLE_THRESHOLD ..."), contained in tag `stp/production-mpi-intel-unstable-guard` (tip 200f32c0) |
| Build dir | `$SCRATCH/prod_2026-11/build_pp_guard` (`GIT_SHA` = 410c7b04, same md5 as the binary above) |
| Configure | `source ./configure.sh IntelLLVM PP GPU=intel HDF5 NO-INTERACTIVE` (via `jobs/build_intel.sh pp`) |
| Flags | `-std08 -cpp -O3 -fp-model=fast=2 -no-prec-div -xHost -unroll -finline-functions -no-wrap-margin -DMPI -DTEMPERING -DPARALLEL_PARAMS -DHDF5`, OpenMP offload `spir64` (JIT) |
| Toolchain | `intel/2025.3.0 intel-mpi/2021.15.0 intel-mkl/2025.3.0 hdf5/1.14.5-intel25` (`jobs/env_intel_pvc.sh`) |

Notes.
- `info` of the runs says `GPU build: git none` (the source was copied with `git archive`), so the commit is taken from the
  binary path, the md5 match and the build directory's `GIT_SHA`.
- RUN_L24, L28, L32 and L36 started on 2026-11-11/12 with `ALF-pp-intel-5667acda.out` (`job.sh.5667acda`) and switched to
  410c7b04 on 2026-11-12. Between 5667acda and 410c7b04 only `Prog/control_mod.F90` changes (the G-precision abort
  threshold can be lowered by `ALF_UNSTABLE_THRESHOLD`; the default is unchanged), so speed and numerics are the same.
- 410c7b04 vs the tag tip 200f32c0: the source differs only in the exit call of `control_read_threshold` (`alf_exit`
  instead of `stop 2`); no effect on speed.
- 410c7b04 already contains the fast `Op_V_is_symmetric` of the GPU line (early return for all-diagonal operator sets,
  per-operator local products). For this Hubbard model the check returns at once. It runs once at start-up.

Benchmark build: `$SCRATCH/grant_bench_2026-11/build_pp/Prog/ALF.out`, from `git archive 410c7b04` with
`jobs/build_intel.sh pp` (the production build script with only the paths changed), job 566780. Result: see section 4.

## 2. Proposed benchmark plan (NOT approved)

What the GCS template (V2.1-2026JAN12) asks for: scaling plots **and** tables with speedup, using the typical
parameters of the planned project, starting at the smallest possible node count and going up to the largest node count used
in production; intra-node performance discussed. The current draft (`ANTRAG_SUPERMUC/common/performance_gpu.tex`) has
three open TODOs: 5-node values at L=24 and 32, a point at 20 to 40 nodes, and L=36 measured instead of extrapolated.

All runs use the production layout (4 chains per tile, 32 ranks per node, `ZEX_NUMBER_OF_CCS` 4 per tile, 3 OpenMP
threads, N_delay 64, Nwrap 10) and the jakob_params parameters with NSweep 2, NBin 4, CPU_MAX 0. Metric: seconds per
sweep and chain from the bin completion times (bin 1 discarded), node throughput (sweeps per node-hour), and
efficiency t(1 node)/t(n nodes).

| Part | What | Cost (estimate) |
|---|---|---|
| A | Intra-node, L=24 and 32: 1 tile with 1 and 4 chains, 1 GPU (2 tiles), 2 GPUs, 1 node; 1 chain per tile on 1 node | ~2 node-h |
| B | Weak scaling, L=24 and 32: 1, 2, 4, 8, 16, 32 nodes, 64 chains per U (32 nodes = 16 U x 64, the largest run in `params_gpu.tex`) | ~19 node-h |
| C | L=36 on 1 and 4 nodes | ~1.5 node-h |
| | Total | ~22 node-h = ~90 GPU-h (~120 GPU-h with reruns) |

Strong scaling: the chains do not communicate except for the bin reduction within a U set, and one chain runs on
one tile. The speedup table for a fixed amount of statistics is therefore derived from B, with the warm-up sweeps per
chain as the part that does not parallelise, and is labelled as derived.

Open questions for the owner: scope (full / lean / minimal), partition for 8 to 32 nodes (`test` allows at most 4 nodes
and 1 h; `general` up to 96 nodes; `preempt` up to 40 nodes and can kill a job), and whether 32 nodes is the right
maximum.

## 3. Files

- `jobs/build_intel.sh`, `jobs/env_intel_pvc.sh`: build (submitted as job 566780).
- `jobs/bench.sh`: DRAFT job script for the plan above, not submitted. Still missing: `make_params.py` and the seeds file.

## 4. Build result (job 566780, 2026-11-27)

`make` rc 0 in 1.5 min; `$SCRATCH/grant_bench_2026-11/build_pp/Prog/ALF.out`, 10,855,456 bytes, md5
`4f035e04b855f70c4a4afc5d3e18ed61`. The md5 differs from the production binary (10,855,672 bytes) because absolute build
paths are compiled in; configure line, modules and `ALF_*` flags are the same apart from these paths (diff of the two
build logs with the paths replaced). `info_git` in the source tree is a tracked file of the repository and says
`localtau`; it does not describe this build (the `GIT_SHA` file does).

## 5. Builds in ANTRAG_2026-11-27 (owner request, 2026-11-27)

Branch `antrag/2026-11-27` on `fork` = commit 410c7b04 without changes (the fast `Op_V_is_symmetric` is already in the
shared `Prog/Operator_mod.F90`, used by CPU and GPU builds). Source: `git archive 410c7b04` in
`/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/src_410c7b04` (`GIT_SHA` file). Build script: `jobs/build_intel.sh`
(paths now point to this directory; `cpu-pp` variant added as in `prod_2026-11/jobs/build_intel.sh`). Jobs charged to pn93la.

| Build | Job | Configure | Binary | md5 |
|---|---|---|---|---|
| GPU, MPI (PP), Phase 2 | 566782 | `IntelLLVM PP GPU=intel HDF5 NO-INTERACTIVE` | `build_pp/Prog/ALF.out` (10,855,376 B) | `d1469afdc28d16bffbaa8d0c3237a805` |
| CPU, MPI (PP), Phase 2 host (Sapphire Rapids, `-xHost`) | 566783 | `IntelLLVM PP HDF5 NO-INTERACTIVE` | `build_cpu-pp/Prog/ALF.out` (7,877,624 B) | `c7cfe97118dbd572907a07b227e152ab` |

Both: `make` rc 0, `-DMPI -DTEMPERING -DPARALLEL_PARAMS -DHDF5`, same toolchain as production. The CPU binary has the same
size as the production CPU binary `prod_2026-11/bin/ALF-pp-cpu-410c7b04.out`. No benchmark started; waiting for the owner.

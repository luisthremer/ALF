# testsuite/lowrank/campaign: Slurm benchmark campaigns with provenance

Generates bash-only Slurm job scripts for benchmark campaigns (SuperMUC-NG, Helma, or the local test profile), captures
provenance on the compute node, and collects/reports the logs locally. Cluster side: bash (3.2-compatible), no python.
Local side (`make_campaign.py`, `collect.py`, `report.py`): python3, stdlib only.

## Workflow

1. **Build and freeze the binary** (on the cluster, in the build shell):
   ```bash
   source configure.sh SUPERMUC-NG MPI NO-INTERACTIVE && make program -j8 && \
     testsuite/lowrank/campaign/build_info.sh $WORK/alf_bin/ALF_<tag>.out
   ```
   `build_info.sh FROZEN [SOURCE]` copies `Prog/ALF.out` (never overwrites a different existing file, makes it read-only)
   and writes `FROZEN.build_info.txt`: compiler (`$ALF_FC --version`), `ALF_FLAGS_PROG`, `ALF_LIB`, `ALF_FLAGS_EXT`,
   `module list`, git commit / dirty count of `ALF_DIR`, binary sha256, `ldd` (or `otool -L`), the names of the linked
   MKL / OpenMP runtime libraries (`libs_mkl`) and the link type (`file`: static or dynamic). Every job checks the
   binary's sha256 against it and copies its KEY=VALUE lines into the job provenance (prefix `build_`). A statically
   linked binary with fakeintel = on aborts the job (LD_PRELOAD cannot act; note that the Helma preset's
   `INTELLLVMOPTFLAGS` contains `-static`).
2. **Fill the profile** `profiles/<machine>.sh`: every `[[...]]` placeholder (account, EAR switch, modules, fakeintel
   library, `BLAS_PROBE_CMD` = Task T2's `blas_bench` with tiny arguments, I_MPI_PMI_LIBRARY on Helma). The generator
   refuses a profile that still contains `[[` anywhere (profiles must not use the double-bracket test syntax) and a
   description with a `[[...]]` placeholder (regex `\[\[[^\]]*\]\]`) in any string value.
3. **Describe the campaign** in a python file defining `CAMPAIGN` (see `tests/campaign_test.py` and the docstring of
   `make_campaign.py`): profile, partition, binary (absolute path on the target), fakeintel on/off, repetitions, seed,
   `max_jobs_in_queue`, `max_attempts`, points (model keys in `make_params.py` syntax incl. `Beta`/`Dtau`, L list,
   schemes, N_delay list incl. 0 = auto and -1 = auto+tune, NSweep/NBin, layouts `RxT`, `partial_node`, `n_fl`,
   `orbitals`, `sec_per_sweep_est`, `setup_s_est` from pilot runs, optional `mem_gb_per_rank`).
4. **Generate** locally: `python3 -I testsuite/lowrank/campaign/make_campaign.py DESC.py --root DIR` (`--dry-run` only
   prints jobs, runs and the core-hour estimate). Checks: ranks x threads = cores per node unless `partial_node`; memory
   estimate (GR + work matrices + UDV stack per flavour) x ranks <= 90 % of node memory; wall estimate (x `wall_safety`
   + 5 min) <= partition limit; at least 3 repetitions for every non-local profile. Same description again: no-op; a different description under an existing id: refused.
5. **Copy** `DIR/<id>` to the cluster (rsync) into the **workspace** file system (SuperMUC `$WORK`/`$SCRATCH`
   workspace, Helma `$WORK`/workspace dirs), never under `$HOME` (`submit_all.sh` warns), and **submit** with
   `<id>/submit_all.sh [--dry-run] [--wait SECONDS]`.
6. **Copy back** and run `collect.py <id> --strict`, then `report.py <id>/collected.json --out report.md`.

## Generated campaign directory

| path | content |
|---|---|
| `campaign.json` | resolved description, profile values, run list (requested model, L, beta, dtau, scheme, N_delay, layout, repetition, seed, NSweep, NBin), sha256 of every script and job script |
| `jobs/<job>.sbatch` | one job = one point x L x layout x repetition with all schemes (standard once, each other scheme x N_delay), order rotated by the repetition; `--exclusive --no-requeue`; submitted from the campaign dir |
| `jobs/<job>/job_<slurmid>/` | `job.prov` (job provenance), raw dumps (`lscpu.txt`, `numactl.txt`, `free.txt`, `scontrol.txt`, `modules.txt`, `env.txt`, `mpi_version.txt`), `mkl_probe/`, `job.status` |
| `runs/<run>/input/` | `parameters`, `seeds` (made locally by `make_params.py`; the cluster needs no python) |
| `runs/<run>/a<k>/` | attempt k: ALF output, `run.prov`, `rank_alf_<r>.prov`, `rss_<r>.txt`, `pin_map.txt`, `omp_affinity.txt`, `stdout.log`, `stderr.log` |
| `runs/<run>/status` | `OK` / `FAILED rc=... reason` / `INVALID pinning: ...` / `NOT_RUN ...` / `RUNNING ...` (history in `status.history`), `latest` names the attempt |
| `submit_all.sh`, `campaign.env`, `jobs.list`, `scripts/` (`profile.sh`, `provenance.sh`, `rank_wrap.sh`), `logs/` | |

Run directories are never deleted: a re-run gets a new attempt directory. `submit_all.sh` skips jobs whose runs are all
OK or that are queued (squeue by name), gives up a run after `max_attempts` attempts and a job after 2 x `max_attempts`
submissions (job-level aborts), and keeps at most `max_jobs_in_queue` jobs queued; a failing `squeue` aborts it (never
read as an empty queue). Inside a job, OK runs are skipped.
A run is FAILED unless the launcher returned 0, `info` has a CPU Time line and, for delayed/submatrix, the
`Sequential scheme: <NAME>` line; it is INVALID (finished, but excluded) unless the pinning record shows exactly RANKS
pairwise disjoint domains of THREADS cpus (THREADS x threads-per-core accepted for core binding with SMT): Intel MPI
`I_MPI_DEBUG>=4` "MPI startup(): <rank> <pid> <node> {cpus}" lines or srun `--cpu-bind=...,verbose` masks
(`prov_pin_check`; no record = INVALID). The EXIT trap marks an interrupted run `FAILED rc=<rc>` (SIGTERM at the time
limit: 143) and the job's remaining runs `NOT_RUN`, except runs that an earlier job left FAILED/INVALID (kept; a note
goes to `status.history`). Job-level aborts (before any run): binary/build-info mismatch, no MKL_VERBOSE line,
fakeintel = on without AVX-512 or with a static binary, THREADS > 1 with the `sequential` MKL threading layer.

## Provenance fields

Job level (`job.prov`, collected as `job.<key>`): Slurm (`SLURM_CLUSTER_NAME`, `SLURMD_NODENAME`, `scontrol show job`:
node list, partition, exclusive = OverSubscribe=NO), CPU (`lscpu`: model, CPUs, sockets, cores per socket, threads per
core, avx512f flag), NUMA (`numactl -H`), memory (`free -m`), clock (`/sys/.../cpufreq`: governor, driver, current
frequency min/mean/max over CPUs, max frequency, boost or intel_pstate no_turbo), EAR (profile), MPI (`$MPI_VERSION_CMD`,
`I_MPI_ROOT`), `module list`, pinning environment (`I_MPI_PIN*`, `I_MPI_DEBUG`, `KMP_*`, `OMP_PLACES`, `OMP_PROC_BIND`,
`SLURM_CPU_BIND`), `OMP_NUM_THREADS = MKL_NUM_THREADS`, `MKL_DYNAMIC=FALSE`, `OMP_DISPLAY_AFFINITY`, `MKL_VERBOSE` unset
in the runs, sha256 of job script / profile / provenance.sh / rank_wrap.sh / binary / build info, all `build_*` keys,
fakeintel on/off + library path and sha256, and the MKL probe: first `MKL_VERBOSE` version line of `$BLAS_PROBE_CMD`
run with `MKL_VERBOSE=1` through the same launcher and per-rank wrapper (version, ISA, threading layer). With
fakeintel = on the job aborts unless the ISA is AVX-512.

Run level (`run.prov`, collected as `run.<key>`): requested scheme / N_delay, clock state at start and end, wall time of
the launch, launcher rc, max RSS over ranks (`/usr/bin/time -v`, BSD `-l` on the Mac, else `sacct MaxRSS`), rank hosts,
per-rank affinity before MPI_Init (`taskset -cp` in rank_wrap.sh, `rank_affinity_preinit`: Intel MPI pins inside
MPI_Init, so this is only the launcher's mask; the verified pinning is `pin_check`), fakeintel `/proc/self/maps` check and `LD_PRELOAD` (set only in `rank_wrap.sh`), pin
map (`I_MPI_DEBUG=4` "MPI startup()" lines or srun `cpu-bind` lines) and its check (`pin_check`), OpenMP affinity lines (`OMP_DISPLAY_AFFINITY`; ALF
itself has no OpenMP region, so lines appear only if MKL threads), status.

From `info` (collect.py): commit, START/FIN TIME, CPU Time, acceptance, Precision Green mean/max and Phase max, beta,
dtau, sweeps, bins, model (if printed), scheme used, N_delay mode (fixed / auto / auto+tune / auto+tune (skipped:
reason)), N_delay used, n0, Ndim. Derived: s/sweep = CPU Time / (NSweep x NBin); setup_s = launch wall time - CPU Time
(start-up, Hamiltonian setup, initial stabilisation, final output). `timers.txt` (T1, `-DALF_TIMERS`) is parsed when
present, else `timers_present = absent`. Until T1 lands, the timers path (parser, bin-1 drop, TIMER_RATE conversion,
timing-split table) is tested only on the fixture `tests/fixtures/timers.txt`.

A source the machine cannot provide is written as `unavailable(<reason>)`; `collect.py --strict` fails on any such or
"unknown" value (only `build_git_commit` may be unavailable: a binary built in a git-archive copy; info's commit is
authoritative and required). On cluster profiles `--strict` also rejects `i_mpi_root=unset` and `omp_affinity` "none
printed" with THREADS > 1. `--allow KEY` (repeatable, e.g. `--allow job.boost_job_start`) accepts a field a node
legitimately lacks; the allowed keys are recorded in collected.json.

## collect.py and report.py

`collect.py CAMPAIGN_DIR [--strict [--allow KEY ...]] [--out FILE.json]`: one row per run; fails on a missing required key, a repeated key
in a `.prov` file, a duplicate run, or a job script whose sha256 differs from `campaign.json`. FAILED runs and runs
whose info contradicts the request (scheme, N_delay mode/value, sweeps/bins, beta/dtau, layout) and INVALID runs are
excluded and listed;
runs never started are listed as missing. Warnings: build_info commit != info commit, node not exclusive.

`report.py COLLECTED.json [...] [--out FILE.md]`: refuses to mix commits, binaries or machines; groups never pool
configurations (key includes campaign, fakeintel, MKL ISA, MKL threading layer: column "config"). Tables: median
s/sweep, range and spread over repetitions with n_valid (groups with fewer than 3 are flagged `<3!`), speedup vs
standard of the same config; max RSS (marked `<=` as an upper bound when only the job-level sacct value exists); node throughput (ranks / s per sweep) per layout;
timing split (bins >= 2, slowest rank, median over repetitions, share of T_SWEEP; integer timer values are converted
with TIMER_RATE); excluded and missing runs.

## Local test

```bash
testsuite/lowrank/campaign/tests/run_local_test.sh runs/campaign_test/bin/ALF_mpi.out runs/campaign_test/bin/ALF_nompi.out
```
Uses the `local` profile with the PATH shims in `tests/shims` (`sbatch` runs the job synchronously, `srun` launches via
MPICH `mpirun` for > 1 rank and prints emulated cpu-bind masks, `squeue`, `scontrol`, `lscpu`, `numactl`, `free`,
`taskset`, `module`, `blas_probe`; failure knobs `SHIM_SRUN_FAIL_MATCH`, `SHIM_SRUN_TERM_MATCH`,
`SHIM_SRUN_PIN_OVERLAP_MATCH`, `SHIM_SQUEUE_FAIL`) and the fixtures in `tests/fixtures` (fake `/sys` cpufreq tree, `/proc/self/maps`,
MKL_VERBOSE lines for SuperMUC AVX-512 (also with the sequential layer) and Helma with/without fakeintel, pinning
records (Intel MPI table 4x12 good/overlapping, srun masks 8x48 on 384 cpus good, overlapping, wrong size), info lines in every N_delay format, a T1
`timers.txt`). Unit checks alone: `tests/test_collect.py`, `tests/test_provenance.sh SCRATCH`. All output goes to
`runs/campaign_test/` (git-excluded).

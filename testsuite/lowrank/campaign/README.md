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
   The SuperMUC profile is filled (account `pn73xu`, `PARTITIONS="test:30"`, `PARTITION_MAX_NODES="test:16"`, `--ear=off`,
   `RUNTIME_ENV_SH` = the env file that `profile_setup` sources before `module load slurm_setup`, `BLAS_PROBE_CMD`); the controller
   verifies once that sbatch accepts `--ear=off`. New profile items: `PARTITION_MAX_NODES` (`partition:nodes`), `NODE_HOST_CHECK`,
   `profile_launch_on` (packed jobs). `job.prov` records `runtime_env_sh` and its sha256.
3. **Parameter template.** `build_info.sh` also freezes `Scripts_and_Parameters_files/Start/{parameters,seeds}` of the
   build tree as `<binary>.start/`; `make_campaign.py` uses it by default when it exists locally (local profile), else
   the description's `template`, else the current checkout's template with a warning. For a cluster campaign copy the
   `.start` directory back (or extract it with `git show <commit>:Scripts_and_Parameters_files/Start/parameters`) and set
   `template`: a newer template can carry keys an older binary cannot read (e.g. `Lowrank_cross_slice`).
4. **Describe the campaign** in a python file defining `CAMPAIGN` (see `tests/campaign_test.py` and the docstring of
   `make_campaign.py`): profile, partition, binary (absolute path on the target), fakeintel on/off, repetitions, seed,
   `max_jobs_in_queue`, `max_attempts`, points (model keys in `make_params.py` syntax incl. `Beta`/`Dtau`, L list,
   schemes, N_delay list incl. 0 = auto and -1 = auto+tune, NSweep/NBin, layouts `RxT`, `partial_node`, `n_fl`,
   `orbitals`, `sec_per_sweep_est`, `setup_s_est` from pilot runs, optional `mem_gb_per_rank`).
5. **Generate** locally: `python3 -I testsuite/lowrank/campaign/make_campaign.py DESC.py --root DIR` (`--dry-run` only
   prints jobs, runs and the core-hour estimate). Checks: ranks x threads = cores per node unless `partial_node`; memory
   estimate (GR + work matrices + UDV stack per flavour) x ranks <= 90 % of node memory; wall estimate (x `wall_safety`
   + 5 min) <= partition limit; at least 3 repetitions for every non-local profile. Same description again: no-op; a different description under an existing id: refused.
6. **Copy** `DIR/<id>` to the cluster (rsync) into the **workspace** file system (SuperMUC `$WORK`/`$SCRATCH`
   workspace, Helma `$WORK`/workspace dirs), never under `$HOME` (`submit_all.sh` warns), and **submit** with
   `<id>/submit_all.sh [--dry-run] [--wait SECONDS]`.
7. **Copy back** and run `collect.py <id> --strict`, then `report.py <id>/collected.json --out report.md`.

## Pilot campaigns and truncated runs (pilot, cpu_max_h)

- `pilot` (top-level description key, bool, default False): a cluster profile normally needs >= 3 repetitions; with `pilot = True`
  exactly 1 repetition is allowed (a warning is printed), campaign.json records `"pilot": true`, collect.py copies the flag, and
  report.py marks the heading of every table (and a banner) as "pilot (1 repetition, no spread)".
- `cpu_max_h` (optional per point, float > 0): written as `CPU_MAX` into the run's parameters. ALF then sets NBin to 10^7, runs bins
  until the budget is used and stops gracefully after the bin in which the limit is reached (the CPU_MAX clock starts after the
  Hamiltonian setup). So the run uses its whole budget regardless of `nbin`: **the run's wall estimate is `cpu_max_h*3600 +
  setup_s_est`** (in the one-node jobs and in the slot packing; `nbin` and `sec_per_sweep_est` may then be omitted). The generator's
  3 x safety formula still applies on top.
- info of such a run has "Prog will stop after hours: <CPU_MAX>" and, at the end, "Effective number of bins: <n>" instead of a "Bins"
  line. collect.py takes the bins written from there: `sec_per_sweep = CPU Time / (NSweep x bins written)` (not NSweep x NBin);
  rows carry `truncated`, `bins_done`, `sweeps_done`; a run whose info disagrees with the request (CPU_MAX missing or different,
  sweeps per bin) is excluded; with a timers.txt a differing bin count is warned. report.py flags truncated runs (banner, and
  "(k truncated, n bins)" in the n_valid column of table 1).

## Packed jobs (pack_nodes)

Why: the SuperMUC-NG allocation may use only the `test` partition: at most 16 nodes x 48 cores and 30 min per job, **only one
test job runs at a time and at most 3 may be submitted (queued)**. One exclusive 1-node job per (point x L x layout x
repetition) would need ~90 sequential half-hour jobs, and at large N the schemes of one job do not fit 30 min together. A
packed job has up to 16 nodes, every node runs its own queue of runs, all nodes concurrently.

`pack_nodes` (optional top-level key of the description, int >= 1, default 1). **`pack_nodes = 1` is the one-node job
described above, unchanged** (job scripts, run inputs, `jobs.list`, `campaign.env`, `submit_all.sh` are byte-identical to
the generator before this key existed; `campaign.json` only gains the profile's `partition_max_nodes` and the new script
hashes). With `pack_nodes > 1` the generator packs RUNS, not jobs:

- **Unit of scheduling**: one run = point x L x scheme x N_delay x layout x repetition. Estimate per run =
  `setup_s_est + nsweep x nbin x sec_per_sweep_est`; a slot's wall = `ceil((sum of its runs + 120 s) x wall_safety / 60) + 5`
  min (the one-node formula, the 120 s start-up and the 5 min margin once per slot, because the slots start concurrently).
- **Slots**: per layout (ranks x threads; jobs never mix layouts) all runs are sorted by estimate (descending, then run id)
  and put into the first slot where the slot's wall still fits the partition limit (first-fit-decreasing, deterministic);
  a new slot is opened otherwise. A single run that alone exceeds the limit is refused ("a single run ... cannot be packed").
- **Jobs**: the slots, in creation order, are cut into jobs of at most `min(pack_nodes, PARTITION_MAX_NODES[partition])` slots
  (a larger `pack_nodes` is cut with a warning). `#SBATCH --nodes=<slots> --ntasks-per-node=RANKS --cpus-per-task=THREADS
  --exclusive --no-requeue`, `--time` = the longest slot. Job names `c_<id>_pack_<layout>_<NN>`. Runs of one slot are ordered by
  point, L, repetition and then the scheme order rotated by the repetition (as in the one-node jobs). All runs of one
  (point, L, repetition) are not forced onto different nodes.
- `make_campaign.py --dry-run` prints the packing: per job its nodes and `--time`, per slot its runs and the estimated wall.
- **Inside the job** (`job_template_packed.sh`): job-level checks and provenance once (binary / build info, `profile_setup`,
  Slurm job, MPI, modules, environment; `job.prov`); host list from `scontrol show hostnames "$SLURM_JOB_NODELIST"`; slot i
  (a background subshell) runs on host i through the profile function `profile_launch_on HOST RANKS THREADS CMD ...`
  (SuperMUC: `mpiexec -hosts "$host" -n "$np" -ppn "$np" "$@"`; local profile: HOST ignored; Helma: `srun --nodelist`, untested).
  The job waits for all slots and prints one status line per slot (also in `job.status`): `slot 1 host=H rc=R runs=N ok=.. failed=..
  invalid=.. not_run=..`. A failing run does not stop its slot (as before); a slot that cannot start stops only itself
  (its runs become `NOT_RUN`, earlier FAILED/INVALID are kept); the others continue. The job's rc is 0 only if every slot's is.
  SIGTERM (time limit) is forwarded to the slots: each marks its running run `FAILED rc=143` and its remaining runs `NOT_RUN`,
  and kills its launcher's process tree.
- **Node-level provenance** (`lscpu`, `numactl`, `free`, cpufreq state, hostname) is taken ON each host by `scripts/node_prov.sh`,
  started through the same launcher (`profile_launch_on HOST 1 1 ...`; chosen over `srun -w HOST` so that one mechanism, the
  profile's launcher, is used and tested for everything, and no `srun` is needed on SuperMUC where the launcher is `mpiexec`).
  The slot compares the hostname it ran on with the requested host (`NODE_HOST_CHECK=1` in the cluster profiles) and aborts
  if the launcher ignored `-hosts`. The MKL probe runs once per host with the slot's launcher and per-rank wrapper; a
  failing probe aborts that slot. Layout: `jobs/<job>/job_<id>/job.prov` (job level, once), `slot_<i>/{node.prov, lscpu.txt, numactl.txt,
  free.txt, mkl_probe/, node_prov.out, job.prov, slot.rc, abort_reason.txt}`; the slot's `job.prov` = job part + node part, and
  each run's `run.prov` has `job_prov = jobs/<job>/job_<id>/slot_<i>/job.prov` plus `slot` and `slot_host`. The run directories,
  attempts and status files are unchanged, so `collect.py` / `report.py` read packed campaigns without change (`job.node`,
  `job.cpu_model`, ... are those of the run's node).
- **Resubmission** (`submit_all.sh`, unchanged): a job with pending runs is submitted again as a whole (same `--nodes`); OK runs are
  skipped inside the slots, so slots without pending runs finish at once, but their nodes stay allocated until the job ends.
- **SuperMUC-NG test partition** (profile `supermuc.sh`: `PARTITIONS="test:30"`, `PARTITION_MAX_NODES="test:16"`): use
  `pack_nodes = 16` and **`max_jobs_in_queue = 3`**. `submit_all.sh` then keeps at most 3 jobs queued; the queue serialises the test
  jobs (one runs at a time), so a campaign of J packed jobs takes about J x 30 min wall. Run `submit_all.sh --wait 300` (or call it
  again later) until all jobs are done. A slot's estimated wall must stay <= 30 min, so the number of runs per slot is small at
  large N; the dry-run shows the slots.

## Generated campaign directory

| path | content |
|---|---|
| `campaign.json` | resolved description, profile values, run list (requested model, L, beta, dtau, scheme, N_delay, layout, repetition, seed, NSweep, NBin), sha256 of every script and job script |
| `jobs/<job>.sbatch` | one job = one point x L x layout x repetition with all schemes (standard once, each other scheme x N_delay), order rotated by the repetition; `--exclusive --no-requeue`; submitted from the campaign dir |
| `jobs/<job>/job_<slurmid>/` | `job.prov` (job provenance), raw dumps (`lscpu.txt`, `numactl.txt`, `free.txt`, `scontrol.txt`, `modules.txt`, `env.txt`, `mpi_version.txt`), `mkl_probe/`, `job.status` (packed jobs: per-slot subdirectories, see "Packed jobs") |
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
`I_MPI_DEBUG>=4` "MPI startup(): <rank> <pid> <node> <cpus>" lines (cpus as `{0,1}`, `0,1`, `0-11` or mixed, braces
optional) or srun `--cpu-bind=...,verbose` masks (`prov_pin_check`). Domains are keyed by (node, cpu), so equal cpu
numbers on different nodes are no overlap. When both records exist (Helma: srun launch with `I_MPI_DEBUG=4`), the Intel
table is judged, because it is the final pinning after MPI_Init; the srun masks are informational. No record, or an
Intel table line that cannot be parsed (quoted in the reason), = INVALID. The EXIT trap marks an interrupted run `FAILED rc=<rc>` (SIGTERM at the time
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
(start-up, Hamiltonian setup, initial stabilisation, final output). `timers.txt` (T1, `-DALF_TIMERS`, format in `testsuite/lowrank/README.md`) is parsed when present, else
`timers_present = absent`: `#` comments, header `TIMER_RATE`, `TIMER_CLOCK`, `TIMER_RANKS`, `TIMER_TICK`, per bin
`TIMER_SWEEPS`, `TIMER <bin> <id> <rank0> <max> <mean>` (seconds) for the nine exclusive categories `update flush vwrap
hop stab meas global langevin_hmc other`, `sweep` (their sum) and `bin_overhead`, and `COUNTER` lines. A second header in
the same file (a restart in the same directory) starts a new block: the last block is kept and a warning recorded.
Warnings (listed by collect.py): sweep != sum of the nine, `flush_unattributed` != 0, sum over bins of
sweep + bin_overhead off info's CPU Time by more than 1 %. The local test injects real-format fixtures
(`tests/fixtures/timers*.txt`: a real 2-bin sample from the T1 binary plus a synthetic bin 3, a restart, a non-zero
flush_unattributed) into a copy of a finished campaign; it runs no flagged binary.

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

`report.py COLLECTED.json [...] [--out FILE.md]`: timing split = the nine exclusive categories per sweep (mean over ranks,
bins >= 2), shares of `sweep`; `sweep` of the slowest rank and `bin_overhead` separately; counters per sweep
(flush_rank_mean averaged over bins); runs with flush_unattributed != 0 are flagged. It refuses to mix commits, binaries or machines; groups never pool
configurations (key includes campaign, fakeintel, MKL ISA, MKL threading layer: column "config"). Tables: median
s/sweep, range and spread over repetitions with n_valid (groups with fewer than 3 are flagged `<3!`), speedup vs
standard of the same config; max RSS (marked `<=` as an upper bound when only the job-level sacct value exists); node throughput (ranks / s per sweep) per layout;
timing split (bins >= 2, slowest rank, median over repetitions, share of T_SWEEP; integer timer values are converted
with TIMER_RATE); excluded and missing runs.

## Local test

```bash
testsuite/lowrank/campaign/tests/run_local_test.sh runs/campaign_test/bin/ALF_mpi.out runs/campaign_test/bin/ALF_nompi.out
testsuite/lowrank/campaign/tests/run_pack_test.sh   runs/campaign_test/bin/ALF_mpi.out runs/campaign_test/bin/ALF_nompi.out
```
`run_pack_test.sh` (about 40 s, at most 3 processes; also the pilot / cpu_max_h tests with `tests/campaign_pilot_test.py`) tests packed jobs with `tests/campaign_pack_test.py`: `pack_nodes = 1`
unchanged, refusals, `--dry-run` packing, the slot assignment against an independent first-fit-decreasing, the job script,
concurrent slots on fake hosts `shimnode1..3` (sbatch shim `--nodes=N`, `scontrol show hostnames` shim, srun shim host log,
`SHIM_PROBE_FAIL_HOST`), a failing run, a failing slot, collect/report, resubmission, SIGTERM, the host check.
Uses the `local` profile with the PATH shims in `tests/shims` (`sbatch` runs the job synchronously, `srun` launches via
MPICH `mpirun` for > 1 rank and prints emulated cpu-bind masks, `squeue`, `scontrol`, `lscpu`, `numactl`, `free`,
`taskset`, `module`, `blas_probe`; failure knobs `SHIM_SRUN_FAIL_MATCH`, `SHIM_SRUN_TERM_MATCH`,
`SHIM_SRUN_PIN_OVERLAP_MATCH`, `SHIM_SQUEUE_FAIL`) and the fixtures in `tests/fixtures` (fake `/sys` cpufreq tree, `/proc/self/maps`,
MKL_VERBOSE lines for SuperMUC AVX-512 (also with the sequential layer) and Helma with/without fakeintel, pinning
records (Intel MPI table 4x12 good/overlapping, srun masks 8x48 on 384 cpus good, overlapping, wrong size), info lines in every N_delay format, a T1
`timers.txt`). Unit checks alone: `tests/test_collect.py`, `tests/test_provenance.sh SCRATCH`. All output goes to
`runs/campaign_test/` (git-excluded).

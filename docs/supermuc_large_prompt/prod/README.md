# Production scripts for ALF PP runs on SuperMUC-NG Phase 2, partition `preempt`

`preempt` requeues a preempted job (`PreemptMode=REQUEUE`, `GraceTime=0`): the job is SIGKILLed at once and its batch script is
rerun from the top with the same job id (max 7 days, 40 nodes). These files make a PP (PARALLEL_PARAMS) run survive that.

## Files

| file | purpose |
|---|---|
| `job_preempt.sh` | Slurm template (placeholders `##NAME## ##ACCOUNT## ##NODES## ##RANKS_PER_NODE## ##OMPT## ##CCS## ##TARGET_BINS## ##EXECUTABLE##`). Partition preempt, 7 days, `--requeue`, `--open-mode=append`, mail END/FAIL/REQUEUE (no mail-user line). Module loads as in `JobfileSuperMUC-NG-Phase2-GPU-PP.sh`. Runs the recovery block from `recover.sh`, starts the rolling-snapshot loop in the background, then `mpiexec` (loop stopped afterwards, also by a trap on EXIT), then logs. |
| `recover.sh` | Bash functions (`rc_*`), sourced by the job script (copy it into the run directory; override with `RECOVER_SH`) or run stand-alone: `bash recover.sh RUNDIR prepare\|guard\|bins\|snap`. Also contains the snapshot loop. |
| `check_params.py RUNDIR REFFILE` | Compares `RUNDIR/parameters` and all `Temp_*/parameters` with the reference (e.g. `../ref_L24/parameters`). Only L1, L2, Beta, N_Phi, ham_U, NSweep, N_delay, CPU_MAX, NBin may differ, `&VAR_TEMP` may be added; asserts N_FL=2, ham_h0=1.d0, Mz=.F., Beta=L.d0, N_Phi=L*L/2, NSweep=10, Nwrap=10; the sets differ only in ham_U. Markdown report, exit 1 on a violation. |
| `check_seeds.py SEEDFILE...` | All seeds in all given files positive and pairwise distinct. Exit 1 otherwise. |
| `restart_test.sh` | Test-partition job (1 node, 32 ranks, 50 min) that proves the restart logic with the real binary, see below. |
| `cmp_h5dump.py A.h5 B.h5` | Numerical comparison of two HDF5 files via `h5dump -m %.17g` (no h5py): bins counted, duplicated bins found, max relative / scaled deviation, `RESULT=PASS\|FAIL`. |

## Procedure for a campaign

1. Create the run directory (login node, python3 only), e.g. L=28, 10 U values, 16 walkers per set, 320 ranks = 10 nodes:
   `python3 Start_PP_GPU/make_pp_dirs.py RUNDIR 28 5.0,5.1,5.2,5.25,5.3,5.35,5.4,5.45,5.5,5.6 16 10 100000 64 160 SEEDOFFSET`
   (`NSweep` 10; CPU_MAX in hours, a bit below the 7 days so that ALF ends cleanly at a bin boundary; NBin is ignored when CPU_MAX > 0).
   Use a different `SEEDOFFSET` for every run directory.
2. `python3 check_params.py RUNDIR ../ref_L24/parameters` (exit 0 required) and `python3 check_seeds.py RUNDIR1/seeds RUNDIR2/seeds ...` over all run directories of the campaign (exit 0 required).
3. Fill the template and copy `recover.sh` next to it (ranks = number of U values x walkers per set; nodes = ranks / 32):
   ```
   sed -e 's|##NAME##|L28_pre|' -e 's|##ACCOUNT##|pn93la|' -e 's|##NODES##|10|' -e 's|##RANKS_PER_NODE##|32|' -e 's|##OMPT##|3|' \
       -e 's|##CCS##|0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4|' -e 's|##TARGET_BINS##|400|' -e 's|##EXECUTABLE##|/path/to/alf_gpu_pp_hdf5.exe|' \
       job_preempt.sh > RUNDIR/job.sh;  cp recover.sh RUNDIR/
   ```
   `##CCS##` empty means ZEX_NUMBER_OF_CCS is not set. `TARGET_BINS` is the number of bins per set in `Temp_i/data.h5` at which the job stops by itself.
4. `cd RUNDIR; sbatch job.sh`. After a preemption nothing has to be done (automatic requeue). If a segment ends cleanly at CPU_MAX
   without reaching TARGET_BINS (no requeue then), `sbatch job.sh` again; when the target is reached the job exits at once.
5. Before the first production run: `restart_test.sh` (below) must report `RESTART TEST: PASS`.

## What the recovery block does (every start of the script, incl. requeues)

a. **Stop guard**: if every `Temp_*/data.h5` has >= TARGET_BINS bins, log and `exit 0`.
b. **RUNNING exists** (the previous segment was killed or crashed): per `Temp_i`, `h5ls -r` must read `data.h5` and every `confout_*.h5`,
   and all bin datasets (`<obs>/sign`, `/obser`, `/back`) of `data.h5` must have the same number of bins (a kill between two appends leaves them different).
   An undamaged set is not touched. A damaged set is restored from the **newest complete rolling snapshot** (`Temp_i/.snap/A|B`: `data.h5` plus its
   `confout_*.h5` copied in as `confin_*.h5`, damaged confouts removed), else from the other slot, else from `Temp_i/.segbak/` (segment start, last resort).
   The slot is re-verified before use (bins of its data.h5 = bins in its COMPLETE marker, number of confouts = walkers, all readable). Log: bins before the damage vs restored = bins lost.
   A file `Temp_i/.restoring` marks a restore in progress, so a kill during the restore is repaired again at the next start. With no usable snapshot and no backup the script stops
   (exit 3) and leaves everything for inspection. Then `RUNNING` is removed.
c. `confout_*` -> `confin_*` per `Temp_i` (refuses if plain and `.h5` confout both exist).
d. Backup of the segment start: `Temp_i/.segbak/` = `cp -p` of `data.h5` and `confin_*` (old backup cleared first, marker file `COMPLETE` written last).
e. One line per segment start and end in `segments.log` (date, `SLURM_JOB_ID`, `SLURM_RESTART_COUNT`, bins per set, recovery action, exit code).

## Rolling snapshots (limits the loss of a corrupting kill to about SNAP_INTERVAL = 30 min of one set)

`rc_snapshot_loop` runs in the background during `mpiexec` and polls every `SNAP_POLL` (30 s). For each `Temp_i` it takes a snapshot when all hold: at least `SNAP_INTERVAL` (1800 s) since the set's newest complete
snapshot (read from the COMPLETE markers, so requeues continue correctly); `data.h5` has a consistent bin count >= 1 and larger than the snapshot's; the bin is complete (as many `confout_*.h5` as walkers
= `mpi_per_parameter_set`, all readable, each with mtime >= mtime of `data.h5` minus `SNAP_CONF_SLACK`); `data.h5` untouched for `SNAP_QUIET_DATA` (10 s), the confouts for `SNAP_QUIET_CONF` (5 s).
Procedure: record size+mtime of all sources, copy into `.snap/tmp/`, re-check that the sources did not change and that the copy has the expected bin count and is readable (else discard, logged),
write `COMPLETE` (`bins=N time=T seq=S walkers=W`), remove the **older** of the slots `.snap/A` / `.snap/B` and `mv` tmp into its place: the newer slot is valid at every moment, a kill during a copy only leaves
`.snap/tmp`. Everything is logged to `snapshots.log` (set, bins, seconds taken, discards, and the reason when a set cannot be snapshotted). The settings are environment variables of the job script
(`export SNAP_INTERVAL=... ` there; the job script defaults are those above except `SNAP_CONF_SLACK=5`).

Caveat: ALF writes the `confout` of ranks other than rank 0 possibly *before* rank 0 has finished appending to `data.h5`, so the strict condition "confout newer than data.h5" can fail for every bin.
The job script therefore allows a slack of 5 s (a production bin lasts far longer, so a confout of the previous bin is still rejected). If `snapshots.log` says "... is older than data.h5", raise `SNAP_CONF_SLACK`;
if no line `slot=` appears at all, snapshots are not being taken: look at the reason lines.

**HDF5 file locking (found on SuperMUC-NG, restart test 2026-11-11):** a plain `h5ls`/`h5dump` on `data.h5` while ALF is running takes a read lock and makes ALF's next `H5Fopen(RDWR)` fail
("unable to lock the file"), which aborts the run. All reads in `recover.sh` therefore go through `rc_h5ls` (`HDF5_USE_FILE_LOCKING=FALSE`). Do the same by hand: never run bare `h5ls`/`h5dump` on a live `data.h5`
(use `HDF5_USE_FILE_LOCKING=FALSE h5ls ...`, or copy the file first).

**Disk overhead** per set: live `data.h5` + `.segbak` (1 copy) + slot A + slot B (1 copy each) + `.snap/tmp` during a copy = 4 copies of `data.h5` in steady state, 5 briefly (the confouts are small). So the 2 snapshot slots cost 2 x `data.h5` per set;
check the quota first (`du -sh Temp_*/data.h5`).

Bin count: `h5ls -r data.h5`, first dimension of the datasets `<obs>/sign`, `/obser`, `/back` (h5ls prints C order, bins are the last Fortran index,
`init_dset`/`append_dat` in `Libraries/Modules/alf_hdf5_mod.F90`). The value used is the common count; `-1` in the log = unreadable, `-2` = datasets disagree.
`Acc_Temp` is ignored.

## Duplicated bin rule (for the analysis)

ALF writes first `data.h5`, then `confout_*`. A kill in between leaves one bin more in `data.h5` than in the saved state, so after the restart that bin is computed
again: **at most one duplicated bin per kill and set**, accepted. The duplicate is the second one of two consecutive bins with identical values in all datasets.
Find them with `python3 cmp_h5dump.py data.h5 data.h5` (prints "identical consecutive bins at [...]", bin index from 0) and delete the second copy
(or the first, they are identical) before the binning analysis; `segments.log` shows the restarts. Kills during a write are not duplicates but are repaired by step b
(restore of that set from its newest snapshot; the lost bins are simply recomputed from that state with the same RNG stream).

## restart_test.sh

`sbatch --export=GPU_PP_BIN=/path/to/hdf5_pp_gpu_exe restart_test.sh` from this directory (needs `$SCRATCH/prod_2026-11/{work_restart,src,jobs}` as described in the script header;
do not use `--export=NONE`). L=8, U=5.0 and 5.6 x 16 walkers, NSweep 2, N_delay 64, CPU_MAX 0, NBin 8:
(a) uninterrupted reference; (b) same seeds, `kill -9` of the mpiexec tree after half of the bins, `recover.sh`, rerun with per-set NBin = 8 - bins done;
(c) as (b), but `Temp_0/confout_3.h5` truncated to 100 bytes, `recover.sh` must restore Temp_0 (from the segment-start backup, i.e. from its seeds); (d) with `SNAP_INTERVAL=0` and the snapshot loop: kill after half of the bins, `Temp_0/data.h5` truncated to half its size and `Temp_1/confout_3.h5` to 100 bytes, both sets must come back from a snapshot and the rerun must reproduce the reference bit for bit; the slots are checked (bins of the slot's data.h5 = COMPLETE marker); (e) as (d) but the kill lands during a snapshot copy (`RC_SNAP_TEST_DELAY`), the previous slot must stay valid. data.h5 of both sets and all final confout files
are compared with the reference by `cmp_h5dump.py`; with a duplicated bin the first 7 unique bins are compared and the final confout is skipped (its state is one bin behind). For (d)/(e) the snapshot timing can be tuned with `SNAP_POLL`, `SNAP_QUIET_DATA`, `SNAP_QUIET_CONF`, `SNAP_CONF_SLACK` (defaults 0.2/0.3/0.2/0.3 s in the test; they must be shorter than one bin). Ends with
`RESTART TEST: PASS|FAIL` and the contents of both `segments.log`. If the run finishes before the kill point is reached the check is reported as inconclusive (FAIL): increase `NSWEEP`.

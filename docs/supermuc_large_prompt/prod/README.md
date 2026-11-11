# Production scripts for ALF PP runs on SuperMUC-NG Phase 2, partition `preempt`

`preempt` requeues a preempted job (`PreemptMode=REQUEUE`, `GraceTime=0`): the job is SIGKILLed at once and its batch script is
rerun from the top with the same job id (max 7 days, 40 nodes). These files make a PP (PARALLEL_PARAMS) run survive that.

## Files

| file | purpose |
|---|---|
| `job_preempt.sh` | Slurm template (placeholders `##NAME## ##ACCOUNT## ##NODES## ##RANKS_PER_NODE## ##OMPT## ##CCS## ##TARGET_BINS## ##EXECUTABLE##`). Partition preempt, 7 days, `--requeue`, `--open-mode=append`, mail END/FAIL/REQUEUE (no mail-user line). Module loads as in `JobfileSuperMUC-NG-Phase2-GPU-PP.sh`. Runs the recovery block from `recover.sh`, then `mpiexec`, then logs. |
| `recover.sh` | Bash functions (`rc_*`), sourced by the job script (copy it into the run directory; override with `RECOVER_SH`) or run stand-alone: `bash recover.sh RUNDIR prepare\|guard\|bins`. |
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
   If all is fine, nothing is touched. Otherwise `data.h5` and `confin_*` of that `Temp_i` are restored from `Temp_i/.segbak/` and its `confout_*` are deleted
   (this loses the bins of that set from the killed segment). With no complete backup the script stops (exit 3) and leaves everything for inspection. Then `RUNNING` is removed.
c. `confout_*` -> `confin_*` per `Temp_i` (refuses if plain and `.h5` confout both exist).
d. Backup: `Temp_i/.segbak/` = `cp -p` of `data.h5` and `confin_*` (old backup cleared first, marker file `COMPLETE` written last). Needs disk space for one extra copy of `data.h5` per set.
e. One line per segment start and end in `segments.log` (date, `SLURM_JOB_ID`, `SLURM_RESTART_COUNT`, bins per set, recovery action, exit code).

Bin count: `h5ls -r data.h5`, first dimension of the datasets `<obs>/sign`, `/obser`, `/back` (h5ls prints C order, bins are the last Fortran index,
`init_dset`/`append_dat` in `Libraries/Modules/alf_hdf5_mod.F90`). The value used is the common count; `-1` in the log = unreadable, `-2` = datasets disagree.
`Acc_Temp` is ignored.

## Duplicated bin rule (for the analysis)

ALF writes first `data.h5`, then `confout_*`. A kill in between leaves one bin more in `data.h5` than in the saved state, so after the restart that bin is computed
again: **at most one duplicated bin per kill and set**, accepted. The duplicate is the second one of two consecutive bins with identical values in all datasets.
Find them with `python3 cmp_h5dump.py data.h5 data.h5` (prints "identical consecutive bins at [...]", bin index from 0) and delete the second copy
(or the first, they are identical) before the binning analysis; `segments.log` shows the restarts. Kills during a write are not duplicates but are repaired by step b
(restore of that set, whose lost bins are simply recomputed from the backup state).

## restart_test.sh

`sbatch --export=GPU_PP_BIN=/path/to/hdf5_pp_gpu_exe restart_test.sh` from this directory (needs `$SCRATCH/prod_2026-11/{work_restart,src,jobs}` as described in the script header;
do not use `--export=NONE`). L=8, U=5.0 and 5.6 x 16 walkers, NSweep 2, N_delay 64, CPU_MAX 0, NBin 8:
(a) uninterrupted reference; (b) same seeds, `kill -9` of the mpiexec tree after half of the bins, `recover.sh`, rerun with per-set NBin = 8 - bins done;
(c) as (b), but `Temp_0/confout_3.h5` truncated to 100 bytes, `recover.sh` must restore Temp_0 (which restarts from its seeds). data.h5 of both sets and all final confout files
are compared with the reference by `cmp_h5dump.py`; with a duplicated bin the first 7 unique bins are compared and the final confout is skipped (its state is one bin behind). Ends with
`RESTART TEST: PASS|FAIL` and the contents of both `segments.log`. If the run finishes before the kill point is reached the check is reported as inconclusive (FAIL): increase `NSWEEP`.

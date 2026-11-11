#!/bin/bash
#
# ALF PP production job for the SuperMUC-NG Phase 2 partition `preempt` (requeue on preemption, SIGKILL without grace time).
# Template: fill the ##...## placeholders (sed, see README.md) and sbatch it from the run directory, which must also contain recover.sh.
#   ##NAME## ##ACCOUNT## ##NODES## ##RANKS_PER_NODE## ##OMPT## ##CCS## ##TARGET_BINS## ##EXECUTABLE##
# A preempted job is rerun from the top of this script (same job id, SLURM_RESTART_COUNT incremented): the recovery block below
# repairs a killed segment (from the rolling snapshots), moves confout->confin and ALF continues from the last written configuration, appending to Temp_*/data.h5.
# Job end by CPU_MAX (clean stop at a bin boundary) does not requeue: sbatch again until every set has ##TARGET_BINS## bins
# (the stop guard then exits immediately once the target is reached).
#
#SBATCH --job-name ##NAME##
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --partition=preempt
#SBATCH --time=7-00:00:00
#SBATCH --requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=##ACCOUNT##
#SBATCH --nodes=##NODES##
#SBATCH --ntasks-per-node=##RANKS_PER_NODE##
#SBATCH --cpus-per-task=1

module load slurm_setup
module load intel/2025.3.0 intel-mpi/2021.15.0
module unload intel-mkl
module load intel-mkl/2025.3.0 hdf5/1.14.5-intel25
for tool in h5ls mpiexec; do
  command -v $tool >/dev/null 2>&1 || { echo "FATAL: $tool not found after the module block (hdf5/1.14.5-intel25, intel-mpi); aborting before touching any file" >&2; exit 2; }
done

export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=##OMPT##
# compute engines per tile (ZEX_NUMBER_OF_CCS), e.g. 0:4,1:4,...,7:4 for 4 walkers per tile; empty = not set
CCS="##CCS##"
if [ -n "$CCS" ]; then export ZEX_NUMBER_OF_CCS="$CCS"; else unset ZEX_NUMBER_OF_CCS; fi
# no I_MPI_PIN_DOMAIN (pinning every rank to one cpu slows the GPU runtime's helper threads)

TARGET_BINS=##TARGET_BINS##
cd "${SLURM_SUBMIT_DIR:-.}" || exit 2

# ---- recovery block (recover.sh: stop guard, repair of a killed segment, confout->confin, backup, segments.log) --------------------
RECOVER_SH=${RECOVER_SH:-${SLURM_SUBMIT_DIR:-.}/recover.sh}
if [ ! -f "$RECOVER_SH" ]; then echo "FATAL: $RECOVER_SH not found (copy prod/recover.sh into the run directory)" >&2; exit 2; fi
source "$RECOVER_SH"

echo "job ${SLURM_JOB_ID:-?} restart_count ${SLURM_RESTART_COUNT:-0} on $(hostname) at $(rc_ts)"
# requeue cap: a job that keeps dying (crash loop) must not be requeued forever
if [ "${SLURM_RESTART_COUNT:-0}" -gt "${MAX_RESTARTS:-100}" ]; then
  echo "$(rc_ts) | job=${SLURM_JOB_ID:-none} | restart_count=${SLURM_RESTART_COUNT:-0} | MAX_RESTARTS=${MAX_RESTARTS:-100} exceeded, exit 4" | tee -a "$RC_LOG_FILE"
  exit 4
fi
if rc_stop_guard; then
  rc_log "stop guard: every Temp_* has >= $TARGET_BINS bins, nothing to do"
  echo "$(rc_ts) | job=${SLURM_JOB_ID:-none} | restart_count=${SLURM_RESTART_COUNT:-0} | stop guard: target $TARGET_BINS bins reached, exit 0 | bins=[$(rc_bins_all)]" >> "$RC_LOG_FILE"
  exit 0
fi
rc_prepare || { echo "FATAL: recovery failed, see messages above (run directory left as it is for inspection)" >&2; exit 3; }

# ---- rolling snapshots (see recover.sh): at most ~SNAP_INTERVAL s of one set are lost when a kill corrupts its files ----------------
# SNAP_CONF_SLACK: ranks other than rank 0 may write their confout a few seconds before rank 0 has finished data.h5, so a confout
# counts as belonging to the current bin if it is at most this much older than data.h5 (a bin lasts much longer than this).
export SNAP_INTERVAL=${SNAP_INTERVAL:-1800} SNAP_POLL=${SNAP_POLL:-30} SNAP_QUIET_DATA=${SNAP_QUIET_DATA:-10} SNAP_QUIET_CONF=${SNAP_QUIET_CONF:-5}
export SNAP_CONF_SLACK=${SNAP_CONF_SLACK:-5}
rc_snapshot_loop >> snapshot_loop.out 2>&1 &
SNAP_PID=$!
trap 'kill $SNAP_PID 2>/dev/null' EXIT
sleep 3
kill -0 $SNAP_PID 2>/dev/null || { echo "$(rc_ts) | SNAPSHOT LOOP NOT RUNNING (died at start, see snapshot_loop.out): NO SNAPSHOTS WILL BE TAKEN" | tee -a snapshots.log >&2; }

mpiexec -n $SLURM_NTASKS ##EXECUTABLE##
RC_EXIT=$?
kill $SNAP_PID 2>/dev/null; wait $SNAP_PID 2>/dev/null
rc_segment_end $RC_EXIT
exit $RC_EXIT

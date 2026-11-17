#!/bin/bash
#SBATCH --job-name=repro273051
#SBATCH --output=repro.%j.log
#SBATCH --account=pn93la
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# One-node reproducer for LRZ incident 273051 (silent wrong results on one PVC tile).
# Usage: sbatch [-w NODE] run_repro.sh [NREP]      (run from this directory; NREP repetitions, default 3)
# Each repetition: 8 MPI ranks, one per tile (flat hierarchy, rank r -> device r), ALF L=44 Hubbard, fixed seeds, 1 sweep.
# Healthy node: every repetition ends with rc=0 and "Acceptance : 0.318327091942149" (bit-identical every time).
# Faulty tile:  verdict=FAULT, rc=3, "... is exceeding the threshold of 10 for G difference!" on the rank bound to that tile.
set -u
D=$(cd "$(dirname "$0")" 2>/dev/null && pwd); [ -f "$D/parameters" ] || D=$SLURM_SUBMIT_DIR
source "$D/env_intel_pvc.sh"
export OMP_NUM_THREADS=14; unset ZEX_NUMBER_OF_CCS ZE_AFFINITY_MASK
NREP=${1:-3}
echo "host $(hostname) job ${SLURM_JOB_ID:-none} $(date +%F_%T) nrep $NREP"
xpu-smi discovery 2>/dev/null | grep -E "Device ID|Device Name|PCI BDF" || true
for i in $(seq 1 "$NREP"); do
  R=$D/run_${SLURM_JOB_ID:-x}_$i; rm -rf "$R"; mkdir -p "$R"
  cp "$D/parameters" "$D/seeds" "$R/"; mkdir -p "$R/Temp_0"; cp "$D/parameters" "$R/Temp_0/"   # this build runs in Temp_0 (parameter-set layout)
  t0=$(date +%s)
  (cd "$R" && mpiexec -n 8 "$D/ALF-pp-intel-5667acda.out" > run.log 2>&1); rc=$?
  acc=$(grep -h "^  Acceptance  " "$R/Temp_0/info" 2>/dev/null | awk '{print $3}')
  bad=$(grep -h -o "Abort([0-9]*) on node [0-9]*" "$R/run.log" | head -1)
  thr=$(grep -h "threshold" "$R/run.log" | head -1)
  if [ $rc -eq 0 ] && [ "$acc" = "0.318327091942149" ]; then v=OK; elif [ -n "$thr" ]; then v=FAULT; else v=ERROR_SEE_run.log; fi
  echo "REP $i rc=$rc $(( $(date +%s)-t0 ))s acceptance=${acc:-none} verdict=$v ${bad:+| $bad (MPI rank N = device N)} ${thr:+| $thr}"
done

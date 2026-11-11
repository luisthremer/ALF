#!/bin/bash
# usage (inside a Slurm job, login-shell bash): run_pp.sh RUNDIR NGPU RANKS_PER_GPU MPS_PCT [CPU_MAX_HOURS]
# Runs the PARALLEL_PARAMS binary in RUNDIR (created by setup_pp.py) with NGPU*RANKS_PER_GPU ranks, rank r on GPU r/RPG,
# one CUDA MPS daemon per GPU (MPS_PCT = active thread percentage per client, 0 = no MPS).
# If CPU_MAX_HOURS is given, CPU_MAX in all Temp_*/parameters is set to it first (restart links use this).
# Restart: ALF resumes from confout_* and appends bins to data.h5 when RUNDIR already holds them.
RUNDIR=$1; NGPU=$2; export RPG=$3; export MPS_PCT=$4; CPUMAX=${5:-}
HERE=$(cd "$(dirname "$0")" && pwd)
export ALF_BIN=${ALF_BIN:-/hnvme/workspace/b133ae29-ALF_GPU/ALF_n257/Prog/ALF.out}
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
export MPS_DIR=${TMPDIR:-/tmp}/mps_${SLURM_JOB_ID:-$$}
cd "$RUNDIR" || exit 1
# restart: confout_* -> confin_* (what out_to_in.sh does), per Temp_i; keep the previous confin/data.h5 as .prev/.bak
# (a kill during an HDF5 write could corrupt a file; data.h5 is validated with h5ls and restored from .bak if unreadable)
H5LS=${H5LS:-/hnvme/workspace/b133ae29-ALF_GPU/hdf5_257/nvfortran_25.7-0/bin/h5ls}
for t in Temp_*; do
  if [ -f $t/data.h5 ]; then
    if $H5LS $t/data.h5 > /dev/null 2>&1; then cp -f $t/data.h5 $t/data.h5.bak
    elif [ -f $t/data.h5.bak ]; then echo "WARNING $t/data.h5 unreadable, restoring .bak" >> links.log; cp -f $t/data.h5.bak $t/data.h5; fi
  fi
  for c in $t/confout_*.h5; do
    [ -f "$c" ] || continue
    i=${c/confout/confin}
    [ -f "$i" ] && cp -f "$i" "$i.prev"
    mv -f "$c" "$i"
  done
done
[ -n "$CPUMAX" ] && sed -i -E "s/^(CPU_MAX\s*=\s*).*/\1$CPUMAX/" Temp_*/parameters
if [ "$MPS_PCT" != 0 ]; then
  for g in $(seq 0 $((NGPU-1))); do
    mkdir -p $MPS_DIR/pipe$g $MPS_DIR/log$g
    CUDA_VISIBLE_DEVICES=$g CUDA_MPS_PIPE_DIRECTORY=$MPS_DIR/pipe$g CUDA_MPS_LOG_DIRECTORY=$MPS_DIR/log$g nvidia-cuda-mps-control -d
  done
fi
nvidia-smi dmon -s um -d 10 -o T > dmon.log 2>&1 &
DMON=$!
date +%s.%N > t_start
${RUN_TIMEOUT:+timeout $RUN_TIMEOUT} mpirun --mca smsc ^knem --oversubscribe --report-bindings -n $((NGPU*RPG)) --bind-to core --map-by core "$HERE/rank_wrap.sh" > run.log 2>&1
rc=$?
date +%s.%N > t_end
echo "rc=$rc" > DONE
kill $DMON 2>/dev/null
if [ "$MPS_PCT" != 0 ]; then
  for g in $(seq 0 $((NGPU-1))); do echo quit | CUDA_MPS_PIPE_DIRECTORY=$MPS_DIR/pipe$g nvidia-cuda-mps-control; done
fi
exit $rc

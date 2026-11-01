#!/bin/bash -l
# BLAS ceiling (testsuite/lowrank/blas_bench.sh) for several layouts at once: node i of the allocation runs layout i
# (MPI binary, mpiexec pinned to that host). Submit from the ALF directory with as many nodes as layouts, e.g. on SuperMUC-NG:
#   sbatch --account=<project> --partition=test --ear=off --nodes=3 \
#          --export=ALL,ALF_MACHINE=SuperMUC-NG,ALF_ENV_SH=<workdir>/env_alf.sh,BB_LAYOUTS="1x48 2x24 4x12" \
#          testsuite/lowrank/cluster/blas_layouts.sbatch
# Optional: BB_N (default 400,1600,3600,8000), BB_REPS (default 5). Logs: runs/blas_bench/blas_bench_<host>_<layout>_<job>.txt.
#SBATCH -J alf_blas
#SBATCH -o runs/blas_layouts_%j.out
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:30:00
#SBATCH --exclusive
#SBATCH --get-user-env
set -o pipefail
hosts=($(scontrol show hostnames "$SLURM_JOB_NODELIST"))
layouts=(${BB_LAYOUTS:-1x48 2x24 4x12})
[ ${#hosts[@]} -ge ${#layouts[@]} ] || { echo "need ${#layouts[@]} nodes, got ${#hosts[@]}" >&2; exit 2; }
echo "job $SLURM_JOB_ID $(date -u +%FT%TZ): layouts ${layouts[*]} on ${hosts[*]}"
export BLAS_BENCH_MACHINE="${ALF_MACHINE:-GNU}" BLAS_BENCH_PROFILE="${BB_PROFILE:-supermuc}" BLAS_BENCH_READY_TIMEOUT=1500
pids=()
for i in "${!layouts[@]}"; do
    h="${hosts[$i]}"; lay="${layouts[$i]}"
    ( BLAS_BENCH_LAUNCH="mpiexec -hosts $h -n @RANKS@ -ppn @RANKS@" testsuite/lowrank/blas_bench.sh LAYOUT="$lay" \
          N="${BB_N:-400,1600,3600,8000}" REPS="${BB_REPS:-5}" TAG="${h}_${lay}_${SLURM_JOB_ID}" \
          >"runs/blas_layouts_${SLURM_JOB_ID}_${lay}.stdout" 2>&1
      echo "layout $lay on $h: exit $?" ) &
    pids+=($!)
done
wait "${pids[@]}"
echo "done $(date -u +%FT%TZ)"

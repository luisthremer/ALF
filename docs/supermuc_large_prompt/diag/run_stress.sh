#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# usage: sbatch run_stress.sh <mode: ccs4 (32 ranks, 4/tile, 3 thr) | one (8 ranks, 1/tile, 14 thr)> <seconds> [extra env K=V ...]
D=$SCRATCH/large_prompt/diag; source $D/env_intel_pvc.sh; cd $D
M=${1:-ccs4}; T=${2:-1200}; shift 2
if [ $M = ccs4 ]; then NR=32; export OMP_NUM_THREADS=3 ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4
else NR=8; export OMP_NUM_THREADS=14; fi
for kv in "$@"; do export "$kv"; done
env | grep -E "^(ZE|ZEX|SYCL|LIBOMPTARGET|STRESS|OMP)_" | sort
export STRESS_SEC=$T
echo "node $(hostname)"
for d in 0 1 2 3; do timeout 60 xpu-smi health -d $d 2>&1 | grep -v "^\s*$" | head -40; done
timeout 120 xpu-smi stats -d 0 -e 2>&1 | grep -i -E "ras|error|ecc|reset" | head -20
timeout 120 xpu-smi stats -d 1 -e 2>&1 | grep -i -E "ras|error|ecc|reset" | head -20
mpiexec -n $NR ./t_stress > stress_$M.$SLURM_JOB_ID.out 2>&1
grep -c ^BAD stress_$M.$SLURM_JOB_ID.out; grep "^BAD" stress_$M.$SLURM_JOB_ID.out | head -20
grep STRESS_RESULT stress_$M.$SLURM_JOB_ID.out | sort -k2n
grep -i -E "alfsh|error|fail" stress_$M.$SLURM_JOB_ID.out | head
for d in 0 1; do timeout 120 xpu-smi stats -d $d -e 2>&1 | grep -i -E "ras|error|ecc|reset" | head -20; done
dmesg 2>/dev/null | tail -5

#!/bin/bash
#SBATCH --job-name=lp_jit
#SBATCH --output=lp_jit.%j.log
#SBATCH --error=lp_jit.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Is the ~140 s start-up overhead of the L=24 GPU runs the JIT compilation of the spir64 kernels? (compare warm vs disabled persistent cache)
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_jit SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters; mkdir -p $W; cd $J
G=$L/build_gpu-txt/Prog/ALF.out
echo "host $(hostname) $(date -Is)"; ls -d ~/.cache/neo_compiler_cache $TMPDIR/neo_compiler_cache 2>&1; du -sh ~/.cache/neo_compiler_cache 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY
echo "=== warm cache (default settings), 2 sweeps, 4 ranks"
python3 run_ref_sng.py jit_warm $G 4 2 1 $J/bind_gpu_per_card.sh
echo "=== persistent cache disabled (NEO_CACHE_PERSISTENT=0)"
NEO_CACHE_PERSISTENT=0 python3 run_ref_sng.py jit_cold $G 4 2 1 $J/bind_gpu_per_card.sh
grep -h "CPU Time" $W/runs/jit_warm/info $W/runs/jit_cold/info
echo "done $(date -Is)"

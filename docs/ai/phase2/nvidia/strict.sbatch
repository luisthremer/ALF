#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH -t 00:30:00
#SBATCH -J p2nv_strict
#SBATCH -o ./strict_%j.out
module load gpu-env/2025 nvhpc/25.7 >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY
P=~/p2_nvidia
DEV_CASES=mz12-std,su2-12-std python3 -I $P/tools/dev.py traj $P/cand_0b84b0e5_strict/Prog/ALF.out $P/runs_0b84b0e5/traj_strict

#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH -c 8
#SBATCH -t 01:30:00
#SBATCH -J camp-prof
#SBATCH -o /home/hpc/b133ae/b133ae29/campaign/prof_%j.out
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader; which nsys
cd ~/camp && python3 docs/ai/phase2/campaign/profile.py nsys $PWD/Prog/ALF.out ~/campaign/prof_work ~/campaign/out h100

#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH -c 16
#SBATCH -t 04:00:00
#SBATCH -J camp-tol
#SBATCH -o /home/hpc/b133ae/b133ae29/campaign/tol_%j.out
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
cd ~/camp && python3 docs/ai/phase2/campaign/campaign.py tol --exe $PWD/Prog/ALF.out --machine h100 --work ~/campaign/tol_work --out ~/campaign/out --jobs 12 --nsweep 2 --timeout 7200

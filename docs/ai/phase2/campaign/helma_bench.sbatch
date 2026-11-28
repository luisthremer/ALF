#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH -c 4
#SBATCH -t 03:00:00
#SBATCH -J camp-bench
#SBATCH -o /home/hpc/b133ae/b133ae29/campaign/bench_%j.out
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader; hostname
cd ~/camp
python3 docs/ai/phase2/campaign/campaign.py bench --exe $PWD/Prog/ALF.out --machine h100 --work ~/campaign/bench_work --out ~/campaign/out --timeout 1800
# CPU path of the same binary (one core of this node), same grid
python3 docs/ai/phase2/campaign/campaign.py bench --cpu --exe $PWD/Prog/ALF.out --machine helma-cpu --work ~/campaign/bench_cpu_work --out ~/campaign/out --timeout 1800

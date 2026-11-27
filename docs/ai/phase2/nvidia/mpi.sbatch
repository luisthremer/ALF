#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH -t 00:30:00
#SBATCH -J p2nv_mpi
#SBATCH -o ./mpi_%j.out
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY
nvidia-smi --query-gpu=index,name --format=csv,noheader
P=~/p2_nvidia; T=$1
python3 -I $P/tools/dev.py mpi $P/cand_${T}_mpi/Prog/ALF.out $P/runs_$T/mpi
nvidia-smi --query-compute-apps=pid,gpu_uuid --format=csv,noheader

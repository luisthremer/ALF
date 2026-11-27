#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH -t 01:30:00
#SBATCH -J p2nv_dev
#SBATCH -o ./dev_%j.out
# usage: sbatch dev.sbatch <cand> <suites...>
T=$1; shift
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
P=~/p2_nvidia; W=$P/runs_$T; D=$P/tools/dev.py
for s in "$@"; do
  echo "=== $s $(date +%T)"
  case $s in
    traj)  python3 -I $D traj $P/cand_$T/Prog/ALF.out $W/traj --ref $P/cand_${T}_none/Prog/ALF.out ;;
    trajlog) python3 -I $D traj $P/cand_${T}_log/Prog/ALF.out $W/trajlog ;;
    g2)    python3 -I $D g2 $P/cand_$T/Prog/ALF.out $W/g2 ;;
    phase) python3 -I $D phase $P/cand_$T/Prog/ALF.out $W/phase ;;
    speed) python3 -I $D speed $P/cand_$T/Prog/ALF.out $W/speed ;;
    speedgp) DEV_START=$P/start_gp python3 -I $D speed $P/ref_gpu-port/Prog/ALF.out $W/speed_gp
             DEV_START=$P/start_gp python3 -I $D speed $P/ref_rd_cgalf-s1/Prog/ALF.out $W/speed_s1 ;;
  esac
done
echo "=== end $(date +%T)"

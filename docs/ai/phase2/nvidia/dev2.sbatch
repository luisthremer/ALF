#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH -t 01:30:00
#SBATCH -J p2nv_dev2
#SBATCH -o ./dev2_%j.out
T=$1
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda cmake/3.30.5 >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
P=~/p2_nvidia; W=$P/runs_$T; D=$P/tools/dev.py
for v in "" _log; do
  echo "=== e2e cand_$T$v $(date +%T)"
  ctest --test-dir $P/cand_$T$v/testsuite/tests -R "gpu-state-e2e|cgr-phase-e2e" -j2 --timeout 1500 --output-on-failure > $P/cand_$T$v/e2e2.log 2>&1
  grep -E "tests passed|Failed|byte for byte|to 1e-13|FAIL:" $P/cand_$T$v/e2e2.log | cut -c1-200
done
echo "=== g2 $(date +%T)"; python3 -I $D g2 $P/cand_$T/Prog/ALF.out $W/g2b
echo "=== speed $(date +%T)"; python3 -I $D speed $P/cand_$T/Prog/ALF.out $W/speedb
echo "=== speed gpu-port $(date +%T)"; DEV_START=$P/start_gp python3 -I $D speed $P/ref_gpu-port/Prog/ALF.out $W/speed_gp
echo "=== speed rd/cgalf-s1 $(date +%T)"; DEV_START=$P/start_gp python3 -I $D speed $P/ref_rd_cgalf-s1/Prog/ALF.out $W/speed_s1
echo "=== end $(date +%T)"

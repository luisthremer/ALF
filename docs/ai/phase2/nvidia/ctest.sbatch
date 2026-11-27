#!/bin/bash -l
#SBATCH -p preempt
#SBATCH -A b133ae
#SBATCH --gres=gpu:h100:1
#SBATCH -t 01:00:00
#SBATCH -J p2nv_ctest
#SBATCH -o ./ctest_%j.out
T=$1
module load gpu-env/2025 nvhpc/25.7 openmpi/5.0.8-nvhpc25.7-cuda cmake/3.30.5 >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
for v in "" _log; do
  echo "=== cand_$T$v"
  cd ~/p2_nvidia/cand_$T$v
  ctest --test-dir testsuite/tests -j8 --timeout 3000 --output-on-failure > ctest.log 2>&1
  grep -E "tests passed|tests failed|\*\*\*Failed|Not Run|Timeout|Exception" ctest.log | head -40
done

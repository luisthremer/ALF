#!/bin/bash
# CPU part, parallel version (replaces the cpu lines of submit_all.sh; 2026-11-27 owner: "run all simulations in parallel").
set -e
sbatch -N 4 -J nt_cpu_L16_m -p test    -t 00:40:00 multi.sh cpu 16 1:1 2:1 4:1 6:1 8:1
sbatch -N 4 -J nt_cpu_L24_m -p test    -t 00:55:00 multi.sh cpu 24 1:1 2:1 4:1 6:1 8:1
sbatch -N 4 -J nt_cpu_L32_m -p test    -t 00:55:00 multi.sh cpu 32 2:1 4:1 6:1 8:1
sbatch -N 4 -J nt_cpu_L36_m -p test    -t 00:55:00 multi.sh cpu 36 2:1 4:1 6:1 8:1
RUN_TIMEOUT=5800 sbatch -N 1 -J nt_cpu_L32t1_m -p general -t 01:40:00 multi.sh cpu 32 1:1
RUN_TIMEOUT=5800 sbatch -N 1 -J nt_cpu_L36t1_m -p general -t 01:40:00 multi.sh cpu 36 1:1

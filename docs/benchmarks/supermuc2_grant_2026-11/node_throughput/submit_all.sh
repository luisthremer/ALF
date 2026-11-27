#!/bin/bash
# Submit the node throughput test (owner OK 2026-11-27). Run on the SuperMUC login node from ANTRAG_2026-11-27/jobs/nt.
set -e
for L in 16 24 32 36; do
  for s in gpu gpub; do REPEAT=2 sbatch -J nt_${s}_L$L -p test -t 00:50:00 run_layouts.sh $s $L; done
done
sbatch -J nt_cpu_L16 -p test    -t 00:55:00 run_layouts.sh cpu 16
sbatch -J nt_cpu_L24 -p general -t 02:00:00 run_layouts.sh cpu 24
sbatch -J nt_cpu_L32 -p general -t 03:00:00 run_layouts.sh cpu 32
RUN_TIMEOUT=7200 sbatch -J nt_cpu_L36 -p general -t 04:00:00 run_layouts.sh cpu 36

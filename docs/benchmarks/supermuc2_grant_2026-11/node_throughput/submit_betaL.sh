#!/bin/bash
# beta = L series (2026-11-27, owner: as fast as possible in wall-clock time). Results in runs_betaL/, separate reports.
# Configurations: the best ones at beta = 4. NBIN 6 (5 timed sweeps) except CPU L=32/36 (NBIN 3-4).
set -e
export BETA=L RUNDIR=runs_betaL
# GPU: two 4-node test jobs, longest runs first
sbatch --export=BETA,RUNDIR -N 4 -J nb_gpu_a -p test -t 00:59:00 multi2.sh \
  gpu:36:4_3:6 gpu:32:4_3:6 gpu:36:2_7:6 gpu:32:2_7:6 gpu:24:4_3:6 gpu:24:2_7:6 gpu:16:4_3:6 gpu:16:2_7:6
sbatch --export=BETA,RUNDIR -N 4 -J nb_gpu_b -p test -t 00:59:00 multi2.sh \
  gpub:36:4_3:6 gpub:36:4_1:6 gpub:32:4_3:6 gpub:32:4_1:6 gpub:24:4_3:6 gpub:24:4_1:6 gpub:16:4_3:6 gpub:16:4_1:6
# CPU on test: L=36 t2 (3 bins), L=36 t4 (4 bins, backup), L=32 t2 (4 bins), then L=24/16 t1, t2 (6 bins)
sbatch --export=BETA,RUNDIR -N 4 -J nb_cpu_t -p test -t 01:00:00 multi2.sh \
  cpu:36:2_1:3 cpu:36:4_1:4 cpu:32:2_1:4 cpu:24:1_1:6 cpu:24:2_1:6 cpu:16:1_1:6 cpu:16:2_1:6
# CPU 1 thread at L=32/36 (best layout at beta=4): too long for test -> general, 3 bins
sbatch --export=BETA,RUNDIR -N 1 -J nb_cpu_L32t1 -p general -t 01:30:00 multi2.sh cpu:32:1_1:3
sbatch --export=BETA,RUNDIR -N 1 -J nb_cpu_L36t1 -p general -t 02:15:00 multi2.sh cpu:36:1_1:3

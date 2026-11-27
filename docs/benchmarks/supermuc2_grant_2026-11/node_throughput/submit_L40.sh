#!/bin/bash
# beta = 4 series extended to L = 40 (2026-11-27), only the fastest configuration per code at L = 16-36 (owner):
# CPU 112 x 1 thread, ALF GPU 2 proc./tile x 7 threads, Intel-opt. 4 proc./tile x 1 thread. Results in runs/.
# GPU NBIN 11 (10 timed sweeps) as before; CPU NBIN 6 (5 timed), so that ~7 min/sweep fits the 1-h test limit.
# (first version with all 5 configurations, jobs 566872/566873, cancelled before start)
set -e
export BETA=4 RUNDIR=runs
sbatch --export=BETA,RUNDIR -N 1 -J nt_L40_cpu -p test -t 00:59:00 multi2.sh cpu:40:1_1:6
sbatch --export=BETA,RUNDIR -N 1 -J nt_L40_gpu -p test -t 00:40:00 multi2.sh gpu:40:2_7:11 gpub:40:4_1:11

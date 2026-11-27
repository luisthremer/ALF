#!/bin/bash
# beta = L series L = 44..56 with the det-phase fix (GPU only), owner OK 2026-11-27. Results in runs_betaL_large/.
# NBIN 4 (3 timed sweeps), L = 56: NBIN 3 (2 timed), so that 2 proc./tile x 7 thr. (~11 min/sweep at L = 56) fits.
set -e
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
export BETA=L RUNDIR=runs_betaL_large EXE_GPU=$A/build_pp_alfgpu_detfix/Prog/ALF.out BINB=$A/bin_batched_65800b82/pp
for L in 56 52 48 44; do
  nb=4; [ $L = 56 ] && nb=3
  part=test; t=00:55:00; [ $L -ge 52 ] && { part=general; t=02:00:00; }
  sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 4 -J nb_big$L -p $part -t $t multi2.sh gpub:$L:2_7:$nb gpu:$L:2_7:$nb gpub:$L:1_14:$nb gpu:$L:1_14:$nb
done

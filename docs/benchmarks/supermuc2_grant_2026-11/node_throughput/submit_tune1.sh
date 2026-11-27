#!/bin/bash
# Tuning stage 1 (owner 2026-11-27): layouts with disjoint core slices per rank (PINMODE=excl), both GPU codes with the det-phase fix,
# beta = 4, L = 24, 40, 56, 8 layouts each. Results in runs_tune/ (tag w<w>t<t>_x). Metric: work units per node-hour.
set -e
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
export BETA=4 RUNDIR=runs_tune EXE_GPU=$A/build_pp_alfgpu_detfix/Prog/ALF.out BINB=$A/bin_batched_65800b82/pp
for L in 56 40 24; do
  nb=6; [ $L = 56 ] && nb=4
  for side in gpub gpu; do
    items=""
    for lay in 1_14 2_7 4_3 4_2 4_1 8_1 6_2 7_2; do items="$items $side:$L:$lay:$nb:x:PINMODE=excl"; done
    t=01:30:00; [ $L != 56 ] && t=00:50:00
    sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 8 -J nt_t1_${side}_$L -p general -t $t multi2.sh $items
  done
done

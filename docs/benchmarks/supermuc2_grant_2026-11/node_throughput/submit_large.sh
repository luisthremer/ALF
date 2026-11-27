#!/bin/bash
# beta = 4 series L = 44..56 with the det-phase fix (GPU only, no CPU): ALF GPU antrag/2026-11-27-detfix 90a51fcb,
# Intel-opt. intel-batched-detfix 65800b82. Layouts 2 proc./tile x 7 thr. and 1 x 14 (best at L=40). Results in runs_large/.
set -e
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
export BETA=4 RUNDIR=runs_large EXE_GPU=$A/build_pp_alfgpu_detfix/Prog/ALF.out BINB=$A/bin_batched_65800b82/pp
for L in 56 52 48 44; do
  sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 4 -J nt_big$L -p test -t 00:50:00 multi2.sh gpub:$L:2_7:11 gpu:$L:2_7:11 gpub:$L:1_14:11 gpu:$L:1_14:11
done

#!/bin/bash -l
#SBATCH -J ib_tput
#SBATCH -A pn93la
#SBATCH -p preempt
#SBATCH --qos=preempt
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 00:25:00
#SBATCH -x i20r02c06s02
#SBATCH -o %x.%j.log
#SBATCH -e %x.%j.err
# Throughput of tonight's best binary (unified build_spd_aottt, commit 819a96ec, V7 + LU-free) per node at L=24 for
# 1, 2, 4 chains per tile (8 / 16 / 32 ranks; CCS = chains per tile). Timer: TEAM_TIMERS of rank 0, 3 sweeps.
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
export W=$SCRATCH/large_prompt/batched/tput; mkdir -p $W
BIN=$SCRATCH/large_prompt/unified/build_spd_aottt/Prog/ALF.out
echo "start $(date -Is) $(hostname) $BIN"
ccs() { local c=$1 s=""; for t in 0 1 2 3 4 5 6 7; do s="$s,$t:$c"; done; echo "ZEX_NUMBER_OF_CCS=${s#,}"; }
for rep in 1 2; do
  for cpt in 1 2 4; do
    nr=$((8*cpt)); omp=$((12/cpt)); t=L24_c${cpt}_r$rep
    env $(var V7) ALFSH_QRNB=16 $(ccs $cpt) OMPT=$omp python3 $J/run_ref_sng.py $t $BIN $nr 3 1 > /dev/null
    tim $t | sed "s/^/chains_per_tile=$cpt /"
    grep -h "Precision Green" $W/runs/$t/info | head -1
  done
done
echo "end $(date -Is)"

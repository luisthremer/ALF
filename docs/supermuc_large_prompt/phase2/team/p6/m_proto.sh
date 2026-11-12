#!/bin/bash
#SBATCH --job-name=p6_proto
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
# Phase 6 prototype: wrap time per walker, default vs batched wrap (stand-in walkers), L=24, V7 switches. B = batch size, mode 1 strided/strided, 2 wide + strided, 3 separate.
P6=$L0/p6; export W=$P6/work_proto; BIN=$P6/build_aottt/Prog/ALF.out; WRAPS=$P6/w_cfg.sh
show() { grep -h -a "GPU wb" $W/runs/$1/run.log | sort -k2,2n | sed "s/^/  [$1] /"; }
ph() { grep -a "TEAM_TIMERS" $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); for(i=0;i<7;i++) s[i]+=$(NF-6+i)/n; c++} END {printf "PHASES %-4s ranks %2d per sweep [s]: sweep %.2f wrapgrup %.2f wrapgrdo %.2f stab_wrapur %.2f stab_wrapul %.2f cgr %.2f obs %.2f\n", t, c, s[0]/c,s[1]/c,s[2]/c,s[3]/c,s[4]/c,s[5]/c,s[6]/c}'; }
cs() { local tag=$1 nr=$2 omp=$3 ns=$4 ccs=$5; shift 5
  env $(var V7) $ccs "$@" OMPT=$omp python3 $J/run_ref_sng.py $tag $BIN $nr $ns 1 $WRAPS > /dev/null; show $tag; ph $tag; pg $tag; }
echo "######## C1: 8 ranks (1 per tile), one configuration per rank (rank 0 baseline; 1,2,3,4,7: B=1,2,4,8,16 mode 1; 5: B=4 mode 2; 6: B=4 mode 3)"
cs C1 8 8 2 "" WBCFG=mixed
echo "######## C2: 8 ranks all baseline (tile-to-tile spread of the baseline)"
cs C2 8 8 2 "" WBCFG=uni
echo "######## C4: 8 ranks all B=4 mode 1"
cs C4 8 8 2 "" WBCFG=uni WB_B=4 WB_M=1
echo "######## C3: 32 ranks (4/tile, CCS 4, OMP 3) all baseline: wrap time under production-like contention"
cs C3 32 3 2 "$CCS" WBCFG=uni
echo "######## C5: 8 ranks (1/tile) CCS 4 (each process gets ONE of the 4 compute slices of its tile), baseline"
cs C5 8 8 2 "$CCS" WBCFG=uni
echo "######## C6: 8 ranks (1/tile) CCS 4, all B=4 mode 1"
cs C6 8 8 2 "$CCS" WBCFG=uni WB_B=4 WB_M=1
echo "end $(date -Is)"

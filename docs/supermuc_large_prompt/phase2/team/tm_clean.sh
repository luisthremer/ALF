#!/bin/bash
#SBATCH --job-name=tm_clean
#SBATCH --output=tm_clean.%j.log
#SBATCH --error=tm_clean.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Clean before/after: in-loop sweep timer (ALF_TEAM_TIMERS: tt(0) = time inside the NSW loop, no setup), old = build_aottt
# (pre-Phase-2 source), new = aotttnew (current source); alternating, one tile, 8 host threads.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
bash build_intel.sh aotttnew > b_aotttnew.out 2>&1; ls -l $L0/build_aotttnew/Prog/ALF.out || exit 1
export W=$L0/work_clean SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
OLD=$L0/build_aottt/Prog/ALF.out; NEW=$L0/build_aotttnew/Prog/ALF.out
one() { # tag bin nsweep
  OMPT=8 python3 run_ref_sng.py $1 $2 1 $3 1 > /dev/null
  grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-12s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'
}
for r in 1 2 3 4; do one L24o$r $OLD 4; one L24n$r $NEW 4; done
export LSIZE=32
for r in 1 2; do one L32o$r $OLD 2; one L32n$r $NEW 2; done
echo "done $(date -Is)"

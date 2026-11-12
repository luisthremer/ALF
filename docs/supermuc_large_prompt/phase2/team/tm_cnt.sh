#!/bin/bash
#SBATCH --job-name=tm_cnt
#SBATCH --output=tm_cnt.%j.log
#SBATCH --error=tm_cnt.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# (1) latency per offload primitive (latbench.F90 from the earlier session, AOT build), (2) number of kernel launches and
# host<->device transfers per sweep from the LIBOMPTARGET_DEBUG log (1 vs 2 sweeps), L=24 one tile.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
AOT="-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc"
ifx -O2 -fiopenmp $AOT -qmkl -I$MKLROOT/include latbench.F90 -o latbench > b_lat.log 2>&1 || tail b_lat.log
export OMP_NUM_THREADS=8 ZE_AFFINITY_MASK=0 OMP_TARGET_OFFLOAD=MANDATORY
echo "== latbench"; ./latbench
B=$L0/build_aottim/Prog/ALF.out; W=$L0/work_cnt
for ns in 1 2; do
  D=$W/s$ns; rm -rf $D; mkdir -p $D; cd $D; cp $L0/ref_L24/seeds .
  sed -e "s/NSweep *= *[0-9]*/NSweep = $ns/" -e 's/NBin *= *[0-9]*/NBin = 1/' -e 's/CPU_MAX *= *[0-9.]*/CPU_MAX = 0.0/' $L0/ref_L24/parameters > parameters
  LIBOMPTARGET_DEBUG=1 $B 2>&1 | grep -E "^(omptarget|Target|TARGET|Libomptarget)" | sed -E 's/0x[0-9a-f]+/X/g; s/[0-9]+/N/g' | cut -c1-90 | sort | uniq -c | sort -rn > counts.txt
  echo "== s$ns: $(wc -l < counts.txt) distinct debug lines"; head -40 counts.txt
done
echo "done $(date -Is)"

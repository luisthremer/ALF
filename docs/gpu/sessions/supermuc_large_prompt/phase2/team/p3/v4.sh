#!/bin/bash
#SBATCH --job-name=p3_v4
#SBATCH --output=p3_v4.%j.log
#SBATCH --error=p3_v4.%j.err
#SBATCH --time=00:58:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# G3: old (team/build_aotttnew, branch head) vs new (p3/build_aotttnew, ALFSH_ASYNCBLK=1), staging 1 MB for both, alternating:
# L24 one tile 8 threads x4 pairs, L32 x2 pairs, node throughput (32 ranks, 4 per tile, CCS 4, OMPT=3)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p3; cd $P
export LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 OMP_TARGET_OFFLOAD=MANDATORY
OLD=$L0/build_aotttnew/Prog/ALF.out; NEW=$P/build_aotttnew/Prog/ALF.out
export W=$P/work_v4 SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters ZE_AFFINITY_MASK=0
one() { # tag bin nsweep sw
  ALFSH_ASYNCBLK=$4 OMPT=8 python3 run_ref_sng.py $1 $2 1 $3 1 > /dev/null
  grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-10s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'
  grep -h P3_TIMERS $W/runs/$1/run.log | awk -v n=$3 '{printf "    P3 per sweep: begin %5.2f (dev %5.2f)  end %5.2f  slice %5.2f\n", $(NF-6)/n,$(NF-5)/n,$(NF-4)/n,$(NF-1)/n}'
}
for r in 1 2 3 4; do one L24o$r $OLD 4 0; one L24n$r $NEW 4 1; done
export LSIZE=32
for r in 1 2; do one L32o$r $OLD 2 0; one L32n$r $NEW 2 1; done
unset LSIZE ZE_AFFINITY_MASK
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
marg() { # tag bin nranks n1 n2 env...
  local t=$1 b=$2 nr=$3 a=$4 c=$5; shift 5
  for ns in $a $c; do env "$@" python3 run_ref_sng.py ${t}_s$ns $b $nr $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${t}_s$a),$(cpu ${t}_s$c); m=(y-x)/($c-$a); print('RESULT %-16s ranks=%-3d s/sweep/chain %.2f  sweeps/s/node %.3f  sweeps/s/card %.3f'%('$t',$nr,m,$nr/m,$nr/m/4))"
  grep -h "Precision Green" $W/runs/${t}_s$c/info | head -1
}
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
marg N32_old $OLD 32 1 3 OMPT=3 $CCS ALFSH_ASYNCBLK=0; marg N32_new $NEW 32 1 3 OMPT=3 $CCS ALFSH_ASYNCBLK=1
marg N32_old2 $OLD 32 1 3 OMPT=3 $CCS ALFSH_ASYNCBLK=0; marg N32_new2 $NEW 32 1 3 OMPT=3 $CCS ALFSH_ASYNCBLK=1
echo "done $(date -Is)"

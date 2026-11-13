#!/bin/bash
#SBATCH --job-name=tm_p4eval
#SBATCH --output=tm_p4eval.%j.log
#SBATCH --error=tm_p4eval.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Lead evaluation of Phase 4 (device QR, ALFSH_DEVQR=1) with staging 1 MB: same build (p4/build_aotttnew), switch off vs on.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
export LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 OMP_TARGET_OFFLOAD=MANDATORY ALFSH_QRNB=16
export W=$L0/work_p4eval SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
B=$L0/p4/build_aotttnew/Prog/ALF.out
one() { # tag nsweep threads devqr env...
  local t=$1 ns=$2 th=$3 dq=$4; shift 4
  env "$@" ALFSH_DEVQR=$dq OMPT=$th ZE_AFFINITY_MASK=0 python3 run_ref_sng.py $t $B 1 $ns 1 > /dev/null
  grep -h TEAM_TIMERS $W/runs/$t/run.log | awk -v t=$t '{n=$(NF-7); printf "RESULT %-10s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'
  grep -h "Precision Green" $W/runs/$t/info
}
for r in 1 2 3; do one L24off$r 4 8 0; one L24on$r 4 8 1; done
for r in 1 2; do one L24t2off$r 4 2 0; one L24t2on$r 4 2 1; done
for r in 1 2; do one L32off$r 2 8 0 LSIZE=32; one L32on$r 2 8 1 LSIZE=32; done
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
node() { local t=$1 dq=$2
  for ns in 1 3; do ALFSH_DEVQR=$dq OMPT=3 ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 python3 run_ref_sng.py ${t}_s$ns $B 32 $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${t}_s1),$(cpu ${t}_s3); m=(y-x)/2; print('RESULT %-10s ranks=32 s/sweep/chain %.2f sweeps/s/node %.3f per card %.3f'%('$t',m,32/m,8/m))"
}
node N32off 0; node N32on 1; node N32off2 0; node N32on2 1
echo "done $(date -Is)"

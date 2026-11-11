#!/bin/bash
#SBATCH --job-name=tm_stage
#SBATCH --output=tm_stage.%j.log
#SBATCH --error=tm_stage.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Lead verification of LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 (Phase 3b finding): checks + timing, branch-head builds.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
SB="LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576"
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
echo "######## CHECK 1: GPU (txtnew, staging 1 MB) vs CPU same seed"
env $SB bash tm_cmp.sh txtnew stage nodefect 2>&1 | grep -E "PASS|FAIL"
echo "######## CHECK 2: txtnew staging 1 MB vs txtnew default, bitwise (GPU vs GPU)"
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
# run_cmp_gg runs 'cpu' (= A) and 'gpu' (= B) with the same env; B gets the knob via a wrapper
printf '#!/bin/bash\nexport %s\nexec "$@"\n' "$SB" > $J/sb_wrap.sh; chmod +x $J/sb_wrap.sh
for a in "S_L6_mzF L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" \
         "S_L6_mzT L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" \
         "S_L12_mzF L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"; do
  W=$L0/work_stagegg WRAP=$J/sb_wrap.sh CPUBIN=$L0/build_txtnew/Prog/ALF.out GPUBIN=$L0/build_txtnew/Prog/ALF.out CMP_TOL=0 python3 run_cmp_gg.py $a 2>&1 | grep -E "PASS|FAIL"
done
unset ZE_AFFINITY_MASK
export W=$L0/work_stage SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
B=$L0/build_aotttnew/Prog/ALF.out
one() { # tag nsweep env...
  local t=$1 ns=$2; shift 2
  env "$@" OMPT=8 ZE_AFFINITY_MASK=0 python3 run_ref_sng.py $t $B 1 $ns 1 > /dev/null
  grep -h TEAM_TIMERS $W/runs/$t/run.log | awk -v t=$t '{n=$(NF-7); printf "RESULT %-10s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'
  grep -h "Precision Green" $W/runs/$t/info
}
for r in 1 2 3; do one L24d$r 4; one L24s$r 4 $SB; done
for r in 1 2; do one L32d$r 2 LSIZE=32; one L32s$r 2 LSIZE=32 $SB; done
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
node() { # tag nranks ompt env...
  local t=$1 nr=$2 th=$3; shift 3
  for ns in 1 3; do env "$@" OMPT=$th python3 run_ref_sng.py ${t}_s$ns $L0/build_aotnew/Prog/ALF.out $nr $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${t}_s1),$(cpu ${t}_s3); m=(y-x)/2; print('RESULT %-10s ranks=%d s/sweep/chain %.2f sweeps/s/node %.3f per card %.3f'%('$t',$nr,m,$nr/m,$nr/m/4))"
}
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
node N32s 32 3 $CCS $SB
node N8s 8 8 $SB
node N16s 16 6 ZEX_NUMBER_OF_CCS=0:2,1:2,2:2,3:2,4:2,5:2,6:2,7:2 $SB
echo "done $(date -Is)"

#!/bin/bash
#SBATCH --job-name=tm_final
#SBATCH --output=tm_final.%j.log
#SBATCH --error=tm_final.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Final before/after: builds aotnew/txtnew from the current team/src; checks (new GPU vs CPU same seed; new vs old GPU bitwise);
# timing old (build_aot = branch start) vs new, alternating: L24 one tile x3, L32 one tile x2, node 32 ranks (4/tile, CCS) x1.
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
bash build_intel.sh aotnew > b_aotnew.out 2>&1 & bash build_intel.sh txtnew > b_txtnew.out 2>&1 & wait
ls -l $L0/build_aotnew/Prog/ALF.out $L0/build_txtnew/Prog/ALF.out || exit 1
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
echo "######## CHECK 1: new GPU vs CPU same seed"
bash tm_cmp.sh txtnew final nodefect 2>&1 | grep -E "PASS|FAIL"
echo "######## CHECK 2: new GPU vs old GPU (build_gpu-txt), bitwise"
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
for a in "G_L6_mzF L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" \
         "G_L12_mzF L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"; do
  W=$L0/work_finalgg CPUBIN=$L0/build_gpu-txt/Prog/ALF.out GPUBIN=$L0/build_txtnew/Prog/ALF.out CMP_TOL=0 python3 run_cmp_gg.py $a 2>&1 | grep -E "PASS|FAIL"
done
unset ZE_AFFINITY_MASK
export W=$L0/work_final SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
cpu() { grep -h "CPU Time" $W/runs/$1/info | awk '{print $NF}'; }
marg() { # tag bin nranks n1 n2 env...
  local t=$1 b=$2 nr=$3 a=$4 c=$5; shift 5
  for ns in $a $c; do env "$@" python3 run_ref_sng.py ${t}_s$ns $b $nr $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${t}_s$a),$(cpu ${t}_s$c); m=(y-x)/($c-$a); print('RESULT %-16s ranks=%-3d s/sweep/chain %.2f  sweeps/s/node %.3f  sweeps/s/card %.3f'%('$t',$nr,m,$nr/m,$nr/m/4))"
  grep -h "Precision Green" $W/runs/${t}_s$c/info | head -1
}
OLD=$L0/build_aot/Prog/ALF.out; NEW=$L0/build_aotnew/Prog/ALF.out
for r in 1 2 3; do marg L24_old_$r $OLD 1 2 6 OMPT=8 ZE_AFFINITY_MASK=0; marg L24_new_$r $NEW 1 2 6 OMPT=8 ZE_AFFINITY_MASK=0; done
for r in 1 2; do marg L32_old_$r $OLD 1 1 3 OMPT=8 ZE_AFFINITY_MASK=0 LSIZE=32; marg L32_new_$r $NEW 1 1 3 OMPT=8 ZE_AFFINITY_MASK=0 LSIZE=32; done
CCS="ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4"
marg N32_old $OLD 32 1 3 OMPT=3 $CCS; marg N32_new $NEW 32 1 3 OMPT=3 $CCS
echo "done $(date -Is)"

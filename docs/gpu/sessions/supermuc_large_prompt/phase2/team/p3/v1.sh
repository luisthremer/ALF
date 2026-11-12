#!/bin/bash
#SBATCH --job-name=p3_v1
#SBATCH --output=p3_v1.%j.log
#SBATCH --error=p3_v1.%j.err
#SBATCH --time=00:55:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14
# Phase 3b step 3: build the async-block source (aotttnew = AOT + timers, txtnew = JIT text), smoke check vs old GPU build, G1 timing pairs
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p3; J=$L0/jobs; cd $P
cp $J/run_cmp_gg.py $J/run_cmp_sng.py $J/run_ref_sng.py $J/compare.py . 2>/dev/null
rm -rf build_aotttnew build_txtnew
( bash build_intel_p3.sh aotttnew > b_aotttnew.out 2>&1 ) & ( bash build_intel_p3.sh txtnew > b_txtnew.out 2>&1 ) & wait
ls -l $P/build_aotttnew/Prog/ALF.out $P/build_txtnew/Prog/ALF.out $P/build_aotttnew/Libraries/Modules/libalfsb.so || { tail -30 b_aotttnew.out b_txtnew.out; exit 1; }
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
echo "######## smoke: L6 Mz=F, new txt build (switch ON and OFF) vs team/build_txtnew, bitwise"
for sw in 0 1; do
  ALFSH_ASYNCBLK=$sw W=$P/work_v1_$sw CPUBIN=$L0/build_txtnew/Prog/ALF.out GPUBIN=$P/build_txtnew/Prog/ALF.out CMP_TOL=0 python3 run_cmp_gg.py V1_L6_mzF_sw$sw L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "PASS|FAIL|rror|ALFSH|async"
done
echo "######## G1 timing: L=24, one tile, 8 threads; old = team/build_aotttnew, new = p3/build_aotttnew with ALFSH_ASYNCBLK=1"
export W=$P/work_v1t SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
one() { # tag bin nsweep sw
  ALFSH_ASYNCBLK=$4 OMPT=8 python3 run_ref_sng.py $1 $2 1 $3 1 > /dev/null
  grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-10s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'
  grep -h "ALFSH\|rror" $W/runs/$1/run.log | head -3
}
OLD=$L0/build_aotttnew/Prog/ALF.out; NEW=$P/build_aotttnew/Prog/ALF.out
for r in 1 2 3; do one old$r $OLD 3 0; one new$r $NEW 3 1; done
echo "done $(date -Is)"

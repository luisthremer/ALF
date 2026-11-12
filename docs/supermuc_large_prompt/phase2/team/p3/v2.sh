#!/bin/bash
#SBATCH --job-name=p3_v2
#SBATCH --output=p3_v2.%j.log
#SBATCH --error=p3_v2.%j.err
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
# G2: (a) same-seed vs CPU, (b) bitwise vs team/build_txtnew (3x), (c) planted missing wait; all with ALFSH_ASYNCBLK=1, staging 1 MB
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p3; J=$L0/jobs; cd $P
# throwaway build: the host reads the gathered block without waiting for the D2H copy
rm -rf src_plant build_txtplant; cp -r src src_plant
sed -i 's/^       call ab_drain()      ! the one unavoidable wait of a block.*$/       continue      ! PLANTED DEFECT: wait removed/' src_plant/Prog/gpu_lowrank_mod.F90
grep -n "PLANTED" src_plant/Prog/gpu_lowrank_mod.F90 || { echo "plant sed failed"; exit 1; }
SRC=$P/src_plant BLD=$P/build_txtplant bash build_intel_p3.sh txtnew > b_txtplant.out 2>&1
ls -l build_txtplant/Prog/ALF.out build_txtnew/Prog/ALF.out || exit 1
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
echo "######## (a) same-seed vs CPU build (same compiler), new GPU txt build, switch ON"
export W=$P/work_v2a CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$P/build_txtnew/Prog/ALF.out; mkdir -p $W
cd $J
Q=run_cmp_sng.py
python3 $Q L6_mzF_sub  L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
CMP_TOL=1e-10 python3 $Q L6_mzT_sub  L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
python3 $Q L12_mzF_sub L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
cd $P
echo "######## (b) bitwise vs team/build_txtnew (old GPU), CMP_TOL=0, 3 repetitions"
OLDB=$L0/build_txtnew/Prog/ALF.out; NEWB=$P/build_txtnew/Prog/ALF.out
for r in 1 2 3; do
  W=$P/work_v2b CPUBIN=$OLDB GPUBIN=$NEWB CMP_TOL=0 python3 run_cmp_gg.py B${r}_L6_mzF L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "PASS|FAIL" | sed "s/^/L6 Mz=F #$r: /"
  W=$P/work_v2b CPUBIN=$OLDB GPUBIN=$NEWB CMP_TOL=0 python3 run_cmp_gg.py B${r}_L6_mzT L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "PASS|FAIL" | sed "s/^/L6 Mz=T #$r: /"
  W=$P/work_v2b CPUBIN=$OLDB GPUBIN=$NEWB CMP_TOL=0 python3 run_cmp_gg.py B${r}_L12_mzF L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72 2>&1 | grep -E "PASS|FAIL" | sed "s/^/L12 Mz=F #$r: /"
done
W=$P/work_v2b CPUBIN=$OLDB GPUBIN=$NEWB CMP_TOL=0 python3 run_cmp_gg.py B1_L24 L=24 mz=F scheme=submatrix ndelay=32 nsweep=2 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288 2>&1 | grep -E "PASS|FAIL" | sed "s/^/L24 2 sweeps: /"
echo "switch engaged in the new-binary runs: $(grep -rl 'ALFSH_ASYNCBLK=1, submatrix' $P/work_v2b 2>/dev/null | wc -l) run logs"
echo "######## (c) planted defect (no wait before the host reads the gathered block) must FAIL"
for a in "L6 L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" "L12 L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"; do
  W=$P/work_v2c CPUBIN=$OLDB GPUBIN=$P/build_txtplant/Prog/ALF.out CMP_TOL=0 python3 run_cmp_gg.py PL_${a%% *} ${a#* } 2>&1 | grep -E "PASS|FAIL|dev" | head -3 | sed "s/^/PLANT ${a%% *}: /"
done
echo "done $(date -Is)"

#!/bin/bash
#SBATCH --job-name=ib_alf
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --exclude=i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14
# Build the intel-batched branch (txtnew variant, recipe of unified/jobs/build_spd.sh) and show that its default behaviour is unchanged:
# gg = bitwise (CMP_TOL=0) vs build_spd_txtnew (819a96ec; the oneMKL code paths are identical up to this branch's additions),
# gate = same seed vs the same-compiler CPU build (L6F/L12F 1e-13, L6T 1e-10), V7 switches.
B=$SCRATCH/large_prompt/batched; U=$SCRATCH/large_prompt/unified; L0=$SCRATCH/large_prompt/team; J=$U/jobs
source $L0/jobs/env_intel_pvc.sh
echo "start $(date -Is) $(hostname)"
BLD=$B/build_ib_txtnew; rm -rf $BLD; mkdir -p $BLD; cp -r $B/src_ib/. $BLD/; cd $BLD
source ./configure.sh IntelLLVM MPI GPU=intel NO-INTERACTIVE > configure.log 2>&1; tail -2 configure.log
export OMP_NUM_THREADS=1
( time make -j56 all ) > make.log 2>&1; tail -4 make.log; ls -l Prog/ALF.out Libraries/Modules/lib*.so
nm -D Libraries/Modules/libalfsb.so 2>/dev/null | grep -c "alfsb_.*_w" | sed 's/^/multi-walker symbols in libalfsb: /'
export W=$B/work_alf; mkdir -p $W; export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
SB=LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
V7="$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_QRNB=16"
declare -A CASE=( [L6F]="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12F]="L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
                  [L6T]="L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" )
declare -A TOL=( [L6F]=1e-13 [L12F]=1e-13 [L6T]=1e-10 )
k=0
for c in L6F L12F L6T; do
  env $V7 CMP_TOL=0 W=$W/gg_$c CPUBIN=$U/build_spd_txtnew/Prog/ALF.out GPUBIN=$BLD/Prog/ALF.out ZE_AFFINITY_MASK=$k python3 $J/run_cmp_gg.py gg_$c ${CASE[$c]} > $W/out_gg_$c.txt 2>&1 &
  k=$((k+1))
  env $V7 CMP_TOL=${TOL[$c]} W=$W/gate_$c CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$BLD/Prog/ALF.out ZE_AFFINITY_MASK=$k python3 $J/run_cmp_sng.py gate_$c ${CASE[$c]} > $W/out_gate_$c.txt 2>&1 &
  k=$((k+1))
done
wait
for f in $W/out_*.txt; do echo "#### $(basename $f)"; grep -E "rc=|PASS|FAIL|acceptance" $f | cut -c1-160; done
echo "end $(date -Is)"

#!/bin/bash
#SBATCH --job-name=p8_g1
#SBATCH --time=00:30:00
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
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
L0=$SCRATCH/large_prompt/team; P=$L0/p8/obs; cd $P/jobs
source p8lib.sh; export W=$C/work_g1; mkdir -p $W/out
export SEEDS0=$L0/src/Scripts_and_Parameters_files/Start/seeds
cp -n $L0/jobs/compare.py $P/jobs/compare.py
B=$C/build_txt/Prog/ALF.out; OLD=$L0/p5/build_txtnew/Prog/ALF.out
# one same-binary pair, OBSPAR 0 (slot "cpu") vs 1 (slot "gpu") (or plant), bitwise; $1 = tile, $2 = tag, $3 = wrapper of slot 2, rest = run_cmp8 keys
pair() { local m=$1 t=$2 w2=$3; shift 3
  env ZE_AFFINITY_MASK=$m WBIN=$B W=$W/pair CPUBIN=$P/jobs/wrap_par0.sh GPUBIN=$P/jobs/$w2 CMP_TOL=0 SEEDS=$SEEDS0 python3 $P/jobs/run_cmp8.py $t "$@" > $W/out/$t.txt 2>&1; }
dflt() { local m=$1 t=$2; shift 2   # default build (no OBSPAR) vs the old p5 binary, bitwise
  env -u ALFSH_DEVQR -u ALFSH_DEVUNMQR -u ALF_GPU_OBSPAR ZE_AFFINITY_MASK=$m W=$W/dflt CPUBIN=$OLD GPUBIN=$B CMP_TOL=0 SEEDS=$SEEDS0 python3 $P/jobs/run_cmp8.py $t "$@" > $W/out/$t.txt 2>&1; }
F6="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
F12="L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
F24="L=24 mz=F scheme=submatrix ndelay=32 nsweep=1 nbin=2 U=5.3d0 symm=T beta=24.d0 nphi=288"
( pair 0 G1_L6F wrap_par1.sh $F6; pair 0 G1_L6T wrap_par1.sh L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
  pair 0 G1_L6_nfl1 wrap_par1.sh L=6 mz=F nfl=1 nsun=2 h0=0.d0 scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 ) &
( pair 1 G1_L12F wrap_par1.sh $F12; pair 1 G1_L6_nsun4 wrap_par1.sh L=6 mz=F nsun=4 scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
  pair 1 G1_L6_bil wrap_par1.sh L=6 lat=Bilayer_square mz=T h0=0.d0 scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 ) &
( pair 2 G1_L24F wrap_par1.sh $F24; export $(var V7); pair 2 G1_L12F_V7 wrap_par1.sh $F12 ) &
( dflt 3 G1d_L6F $F6; dflt 3 G1d_L12F $F12; pair 3 G1_L6F_long wrap_par1.sh L=6 mz=F scheme=submatrix ndelay=8 nsweep=60 nbin=3 U=5.3d0 symm=T beta=6.d0 nphi=18; pair 3 G1_L12F_long wrap_par1.sh L=12 mz=F scheme=submatrix ndelay=32 nsweep=12 nbin=2 U=5.3d0 symm=T beta=12.d0 nphi=72 ) &
( pair 4 G2_L6F wrap_plant.sh $F6; pair 4 G2_L12F wrap_plant.sh $F12; pair 4 G2_L24F wrap_plant.sh $F24 ) &
wait
for f in $W/out/*.txt; do echo "## $(basename $f .txt): $(grep -E 'rc=' $f | tr '\n' ' ')"; grep -E "^(PASS|FAIL)|acceptance" $f; done
echo "end $(date -Is)"

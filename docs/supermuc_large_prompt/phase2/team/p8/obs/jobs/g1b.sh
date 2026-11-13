#!/bin/bash
#SBATCH --job-name=p8_g1b
#SBATCH --time=00:15:00
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
source p8lib.sh; export BTXT=$C/build_txt9; export W=$C/work_g1b; mkdir -p $W/out
export SEEDS0=$L0/src/Scripts_and_Parameters_files/Start/seeds
B=${BTXT:-$C/build_txt}/Prog/ALF.out
pairv() { local m=$1 v=$2 t=$3; shift 3
  env ZE_AFFINITY_MASK=$m WVAR=$v WBIN=$B W=$W/pair CPUBIN=$P/jobs/wrap_par0.sh GPUBIN=$P/jobs/wrap_var.sh CMP_TOL=0 SEEDS=$SEEDS0 python3 $P/jobs/run_cmp8.py $t "$@" > $W/out/$t.txt 2>&1; }
F6="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
F12="L=12 mz=F scheme=submatrix ndelay=32 nsweep=1 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
for v in 0 9; do ( pairv $((v%8)) $v V${v}_L6F $F6; pairv $((v%8)) $v V${v}_L12F $F12 ) & done
wait
for f in $W/out/*.txt; do echo "## $(basename $f .txt): $(grep -E '^(PASS|FAIL)' $f)"; grep -E "^(Den|Green|Spin)" $f | grep -v " 0.000e+00   0.000e+00"; done
echo "end $(date -Is)"
grep -h DBG $W/pair/runs/V9_L6F/gpu/run.log | head -40

#!/bin/bash
#SBATCH --job-name=p6_gate0
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
# default build unchanged: p5 txtnew (before the prototype) vs p6 txtnew (prototype compiled in, no switches), bitwise, V7 switches
P6=$L0/p6; export W=$P6/work_gate0
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 CPUBIN=$C/build_txtnew/Prog/ALF.out GPUBIN=$P6/build_txtnew/Prog/ALF.out CMP_TOL=0
for a in "g0_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" \
         "g0_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" \
         "g0_L24F L=24 mz=F scheme=submatrix ndelay=32 nsweep=1 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288"; do
  echo "## ${a%% *}"; env $(var V7) python3 $J/run_cmp_gg.py $a 2>&1 | grep -E "worst rel dev|PASS|FAIL|acceptance"
done
echo "end $(date -Is)"

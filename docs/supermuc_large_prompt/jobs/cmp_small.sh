#!/bin/bash
#SBATCH --job-name=lp_cmp
#SBATCH --output=lp_cmp.%j.log
#SBATCH --error=lp_cmp.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Step 2: small system, Intel GPU (1 rank, device 0) vs same-compiler ifx CPU, fixed seeds, 1 node
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
export W=$SCRATCH/large_prompt/work_cmp CPUBIN=$SCRATCH/large_prompt/build_cpu-txt/Prog/ALF.out GPUBIN=$SCRATCH/large_prompt/build_gpu-txt/Prog/ALF.out
export SEEDS=$SCRATCH/large_prompt/src/Scripts_and_Parameters_files/Start/seeds
mkdir -p $W; cd $SCRATCH/large_prompt/jobs
echo "host $(hostname) $(date -Is)"; sycl-ls 2>&1 | grep -c "level_zero:gpu"
P=run_cmp_sng.py
python3 $P L6_mzF_sub   L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
python3 $P L6_mzT_sub   L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
python3 $P L6_mzF_std   L=6 mz=F scheme=standard nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
echo "done $(date -Is)"

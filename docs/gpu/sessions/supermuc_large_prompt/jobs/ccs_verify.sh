#!/bin/bash
#SBATCH --job-name=lp_ccsv
#SBATCH --output=lp_ccsv.%j.log
#SBATCH --error=lp_ccsv.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Numerical check of the production configuration: 32 ranks (4 per tile) with ZEX_NUMBER_OF_CCS=:4 vs 32 ifx-CPU ranks, same seeds
# (text builds, L=8, Mz=F and Mz=T). Fixed-seed trajectories are only comparable with a same-compiler CPU build.
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs; cd $J
export W=$L/work_ccsv CPUBIN=$L/build_cpu-txt/Prog/ALF.out GPUBIN=${GPUBIN:-$L/build_gpu-txt/Prog/ALF.out} SEEDS=$L/ref_L24/seeds
export ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4
mkdir -p $W; echo "host $(hostname) $(date -Is) gpu=$GPUBIN"
python3 run_cmp_sng.py ccs32_mzF L=8 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=2 U=5.3d0 symm=T beta=8.d0 nphi=32 nr=32
python3 run_cmp_sng.py ccs32_mzT L=8 mz=T scheme=submatrix ndelay=32 nsweep=4 nbin=2 U=5.3d0 symm=T beta=8.d0 nphi=32 nr=32
echo "done $(date -Is)"

#!/bin/bash
#SBATCH --job-name=tm_cmp
#SBATCH --output=tm_cmp.%j.log
#SBATCH --error=tm_cmp.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Same-seed GPU (gpu-txt) vs same-compiler CPU (cpu-txt) comparison + planted defects (ALF_GPU_FAULT_STAB eps:x, ALF_GPU_NS=0).
# usage: sbatch tm_cmp.sh [GPUBUILD default gpu-txt] [TAGSUFFIX]
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J; GB=${1:-gpu-txt}; S=${2:-base}
export W=$L0/work_cmp_$S CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$L0/build_$GB/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
mkdir -p $W; echo "host $(hostname) $(date -Is) gpu=$GPUBIN"
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
P=run_cmp_sng.py
if [ "$3" != onlydefect ]; then
python3 $P L6_mzF_sub  L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
CMP_TOL=1e-10 python3 $P L6_mzT_sub  L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
python3 $P L12_mzF_sub L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
fi
if [ "$3" != nodefect ]; then
echo "=== planted defect 1: ALF_GPU_FAULT_STAB=eps:1e-12 (graded perturbation of G at every stabilisation)"
ALF_GPU_DEBUG=1 ALF_GPU_FAULT_STAB=eps:1e-12 python3 $P D1_L6_mzF_sub L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
echo "=== planted defect 2: ALF_GPU_NS=0 (no Newton-Schulz step in the device QR)"
ALF_GPU_DEBUG=1 ALF_GPU_NS=0 python3 $P D2_L12_mzF_sub L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
fi
echo "done $(date -Is)"

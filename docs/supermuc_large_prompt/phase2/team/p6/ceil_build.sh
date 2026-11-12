#!/bin/bash
#SBATCH --job-name=p6_cbuild
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
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:20:00
# build the Gate-1 ceiling benchmark against a private copy of the combined-branch AOT build (p5/build_aottt, read-only source)
echo "start $(date -Is) node $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P6=$SCRATCH/large_prompt/team/p6; P5=$SCRATCH/large_prompt/team/p5
cd $P6; [ -d build ] || cp -a $P5/build_aottt build
cd build
export ALF_HDF5_DIR=$PWD/h5ext
source ./configure.sh IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE > configure2.log 2>&1
for v in ALF_FLAGS_MODULES ALF_FLAGS_ANA ALF_FLAGS_PROG ALF_LIB; do
  eval "export $v=\"\${$v//-fopenmp-targets=spir64/-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc}\""
done
cd $P6
$ALF_FC $ALF_FLAGS_PROG -I build/Libraries/Modules -o ceil_bench ceil_bench.F90 build/Libraries/Modules/modules_90.a $ALF_LIB && ls -l ceil_bench
echo "end $(date -Is)"

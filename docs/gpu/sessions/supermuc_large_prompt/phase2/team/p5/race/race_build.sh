#!/bin/bash
#SBATCH --job-name=rc_build
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
#SBATCH --cpus-per-task=48
#SBATCH --time=00:30:00
# copy of the combined-branch AOT build (p5/build_aottt: sources identical to worktree head 2d5cc615) + t_stress linked against its modules
echo "start $(date -Is) node $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p5; R=$P/race; cd $R
rm -rf build; cp -a $P/build_aottt build; cd build
export ALF_HDF5_DIR=$PWD/h5ext
source ./configure.sh IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE > configure2.log 2>&1
for v in ALF_FLAGS_MODULES ALF_FLAGS_ANA ALF_FLAGS_PROG ALF_LIB; do
  eval "export $v=\"\${$v//-fopenmp-targets=spir64/-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc}\""
done
echo "FC=$ALF_FC"; echo "FLAGS=$ALF_FLAGS_PROG"; echo "LIB=$ALF_LIB"
cd $R
$ALF_FC $ALF_FLAGS_PROG -I build/Libraries/Modules -o t_stress t_stress.F90 build/Libraries/Modules/modules_90.a $ALF_LIB && ls -l t_stress
ls -l build/Libraries/Modules/libalfdevqr.so; md5sum build/Libraries/Modules/libalfdevqr.so $P/build_aottt/Libraries/Modules/libalfdevqr.so
echo "end $(date -Is)"

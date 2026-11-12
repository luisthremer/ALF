#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# build the production commit 5667acda (GPU=intel, MPI, no HDF5, JIT like production) and link t_stress against its modules
D=$SCRATCH/large_prompt/diag
source $D/env_intel_pvc.sh
cd $D; rm -rf build; mkdir build; tar -xf src_5667acda.tar -C build; cd build
source ./configure.sh IntelLLVM MPI GPU=intel NO-INTERACTIVE > configure.log 2>&1
export OMP_NUM_THREADS=1
make -j48 lib > make.log 2>&1 || make -j48 all >> make.log 2>&1
tail -5 make.log; ls -l Libraries/Modules/*.a
echo "FC=$ALF_FC"; echo "FLAGS=$ALF_FLAGS_PROG"; echo "LIB=$ALF_LIB"
cd $D && $ALF_FC $ALF_FLAGS_PROG -I build/Libraries/Modules -o t_stress t_stress.F90 build/Libraries/Modules/modules_90.a $ALF_LIB && ls -l t_stress

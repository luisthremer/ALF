#!/bin/bash
#SBATCH --job-name=lp_ppv
#SBATCH --output=lp_ppv.%j.log
#SBATCH --error=lp_ppv.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# PARALLEL_PARAMS build: Intel GPU vs same-compiler ifx CPU, same seeds, 3 U values x 8 walkers (24 ranks), L=8 (text builds), the
# production settings Nwrap=20 N_delay=64 and ZEX_NUMBER_OF_CCS; every Temp_<i>/ is compared with compare.py.
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs; W=$L/work_ppv; rm -rf $W; mkdir -p $W; cd $J
export ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 OMP_NUM_THREADS=1
echo "host $(hostname) $(date -Is)"
for k in gpu cpu; do
  python3 $L/src/Scripts_and_Parameters_files/Start_PP_GPU/make_pp_dirs.py $W/$k 8 5.0,5.3,5.6 8 4 2 64 0.0 7 | tail -1
  sed -i -e 's/Nwrap *= *10/Nwrap = 20/' $W/$k/parameters $W/$k/Temp_*/parameters
  grep -q "^ *Nwrap *= *20" $W/$k/Temp_1/parameters && echo "Nwrap=20 set"
done
cd $W/gpu && OMP_TARGET_OFFLOAD=MANDATORY mpiexec -n 24 $L/build_gpu-pptxt/Prog/ALF.out > run.log 2>&1; echo "gpu rc=$?"; grep -h "GPU path" Temp_0/info | head -2
cd $W/cpu && OMP_TARGET_OFFLOAD=DISABLED mpiexec -n 24 $L/build_cpu-pptxt/Prog/ALF.out > run.log 2>&1; echo "cpu rc=$?"
cd $J; for i in 0 1 2; do echo "== Temp_$i (U = $(grep -h ham_U\  $W/gpu/Temp_$i/parameters | head -1))"; python3 compare.py $W/cpu/Temp_$i $W/gpu/Temp_$i | tail -3; done
echo "done $(date -Is)"

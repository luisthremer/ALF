#!/bin/bash
#SBATCH --job-name=lp_l24
#SBATCH --output=lp_l24.%j.log
#SBATCH --error=lp_l24.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Steps 4+5: L=24, Beta=24, U=5.3 reference parameters, short (cannot equilibrate: reference needs 1000 warm-up sweeps/walker)
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_l24 SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters   # the reference run's own parameters and seeds
mkdir -p $W; cd $J
G=$L/build_gpu-txt/Prog/ALF.out; C=$L/build_cpu-txt/Prog/ALF.out; GH=$L/build_mpi/Prog/ALF.out
echo "host $(hostname) $(date -Is)"
echo "=== T1: 1 rank, 1 GPU (device 0), 10 sweeps"
OMP_TARGET_OFFLOAD=MANDATORY python3 run_ref_sng.py l24_g1 $G 1 5 2
echo "=== T2: 4 ranks, one per card, 10 sweeps each"
OMP_TARGET_OFFLOAD=MANDATORY python3 run_ref_sng.py l24_g4 $G 4 5 2 $J/bind_gpu_per_card.sh
echo "=== T3a: 4 ranks GPU one per card, 4 sweeps (same seeds as CPU)"
OMP_TARGET_OFFLOAD=MANDATORY python3 run_ref_sng.py l24_g4s $G 4 2 2 $J/bind_gpu_per_card.sh
echo "=== T3b: 4 ranks ifx CPU (1 thread each), 4 sweeps"
OMP_TARGET_OFFLOAD=DISABLED python3 run_ref_sng.py l24_c4s $C 4 2 2
python3 compare.py $W/runs/l24_c4s $W/runs/l24_g4s | tail -20
echo "=== T4: production binary (HDF5), 4 ranks one per card, 10 sweeps"
OMP_TARGET_OFFLOAD=MANDATORY python3 run_ref_sng.py l24_g4h $GH 4 5 2 $J/bind_gpu_per_card.sh
ls -l $W/runs/l24_g4h/data.h5
echo "=== D: no offload device -> MANDATORY must fail loudly"
mkdir -p $W/runs/noDev && cp $W/runs/l24_g4s/parameters $W/runs/l24_g4s/seeds $W/runs/noDev/ && cd $W/runs/noDev
ONEAPI_DEVICE_SELECTOR=opencl:cpu OMP_TARGET_OFFLOAD=MANDATORY mpiexec -n 2 $G 2>&1 | head -12; echo "rc=${PIPESTATUS[0]}"
echo "done $(date -Is)"

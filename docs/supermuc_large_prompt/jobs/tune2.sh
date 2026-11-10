#!/bin/bash
#SBATCH --job-name=lp_tune2
#SBATCH --output=lp_tune2.%j.log
#SBATCH --error=lp_tune2.%j.err
#SBATCH --time=00:55:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Speed: host threads per rank (the stabilisation's pivoted QR, ZGEQP3, runs on the host), N_delay. L=24 reference, HDF5 binary, 4 sweeps in 1 bin.
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_tune2 SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters; mkdir -p $W; cd $J
G=${1:-$L/build_mpi/Prog/ALF.out}
export OMP_TARGET_OFFLOAD=MANDATORY ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4
export I_MPI_PIN_DOMAIN=omp:compact
echo "host $(hostname) $(date -Is) bin=$G"
one() { # tag nranks threads
  OMPT=$3 python3 run_ref_sng.py t2_$1 $G $2 4 1 | tail -1
  echo "   CPU Time $(grep -h 'CPU Time' $W/runs/t2_$1/info | awk '{print $NF}') s  [ranks=$2 OMP_NUM_THREADS=$3]"
}
one r1_t1 1 1
one r1_t2 1 2
one r1_t4 1 4
one r1_t8 1 8
one r1_t16 1 16
echo "=== 32 ranks (4 per tile, CCS), host threads per rank"
one r32_t1 32 1
one r32_t3 32 3
echo "done $(date -Is)"

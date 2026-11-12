#!/bin/bash
#SBATCH --job-name=lp_tune3
#SBATCH --output=lp_tune3.%j.log
#SBATCH --error=lp_tune3.%j.err
#SBATCH --time=00:58:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Speed at the production size L=28 (Beta=L): 32 ranks (4 per tile, CCS), host threads per rank, Nwrap, N_delay. 3 sweeps (1 bin) each.
# usage: sbatch tune3.sh [L] [binary]; marginal cost from the 3-sweep run minus a 1-sweep run is skipped (start-up ~10 s): compare CPU Time only.
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt; J=$L0/jobs; cd $J
LL=${1:-28}
G=${2:-$L0/build_mpi/Prog/ALF.out}
export SEEDS=$L0/ref_L24/seeds OMP_TARGET_OFFLOAD=MANDATORY ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 I_MPI_PIN_DOMAIN=omp:compact
export W=$L0/work_tune3_L$LL; mkdir -p $W
mk() { # nwrap ndelay -> params file
  sed -e "s/L1 *= *24/L1 = $LL/" -e "s/L2 *= *24/L2 = $LL/" -e "s/N_Phi *= *288/N_Phi = $((LL*LL/2))/" -e "s/Beta *= *24.d0/Beta = $LL.d0/" \
      -e "s/Nwrap *= *10/Nwrap = $1/" -e "s/N_delay *= *32/N_delay = $2/" $L0/ref_L24/parameters > $W/params_w$1_d$2
}
one() { # tag nranks threads nwrap ndelay
  mk $4 $5; export PARAMS=$W/params_w$4_d$5
  OMPT=$3 python3 run_ref_sng.py $1 $G $2 3 1 | tail -1
  echo "   CPU Time $(grep -h 'CPU Time' $W/runs/$1/info | awk '{print $NF}') s  Precision: $(grep -h 'Precision Green' $W/runs/$1/info | head -1 | awk '{print $(NF-1), $NF}')  [L=$LL ranks=$2 threads=$3 Nwrap=$4 N_delay=$5]"
}
echo "host $(hostname) $(date -Is) bin=$G"
one a_t1_w10_d32   32 1 10 32
one b_t3_w10_d32   32 3 10 32
one c_t6_w10_d32   32 6 10 32
one d_t3_w20_d32   32 3 20 32
one e_t3_w10_d64   32 3 10 64
one f_t3_w10_d128  32 3 10 128
echo "done $(date -Is)"

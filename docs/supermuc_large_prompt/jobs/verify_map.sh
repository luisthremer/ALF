#!/bin/bash
#SBATCH --job-name=lp_map
#SBATCH --output=lp_map.%j.log
#SBATCH --error=lp_map.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Verify the rank -> tile mapping (no ZE_AFFINITY_MASK, rank mod 8) with per-tile utilisation: start 3 ranks, expect tiles 0,1,2 busy only.
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
L=$SCRATCH/large_prompt; J=$L/jobs
export W=$L/work_map SEEDS=$L/ref_L24/seeds PARAMS=$L/ref_L24/parameters; mkdir -p $W; cd $J
export OMP_TARGET_OFFLOAD=MANDATORY
echo "host $(hostname) $(date -Is)"; xpu-smi dump -h 2>&1 | grep -i -E "tile|-t," | head -5
( sleep 40; for i in 1 2 3; do echo "-- sample $i $(date +%T)"; for d in 0 1 2 3; do for t in 0 1; do printf "card %s tile %s: " $d $t; xpu-smi dump -d $d -t $t -m 0,18 -n 1 2>&1 | grep -v lspci | tail -1; done; done; sleep 5; done ) > $W/tiles.txt 2>&1 &
python3 run_ref_sng.py map3 $L/build_mpi/Prog/ALF.out 3 6 1 | tail -1
wait; cat $W/tiles.txt | cut -c1-140
echo "done $(date -Is)"

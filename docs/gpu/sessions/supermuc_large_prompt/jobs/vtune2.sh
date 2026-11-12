#!/bin/bash
#SBATCH --job-name=lp_vt2
#SBATCH --output=lp_vt2.%j.log
#SBATCH --error=lp_vt2.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# Phase 2 Part 4: inclusive (top-down) host time per routine of ONE rank, L=24 reference parameters, 3 sweeps, baseline settings
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
module load intel-oneapi-vtune/2024.2.1
L=$SCRATCH/large_prompt; W=$L/work_vtune2; rm -rf $W; mkdir -p $W/run; cd $W/run
cp $L/ref_L24/seeds seeds
sed -e 's/NSweep *= *[0-9]*/NSweep = 3/' -e 's/NBin *= *[0-9]*/NBin = 1/' -e 's/CPU_MAX *= *[0-9.]*/CPU_MAX = 0.0/' $L/ref_L24/parameters > parameters
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 ZE_AFFINITY_MASK=0
echo "host $(hostname) $(date -Is)"
vtune -collect hotspots -knob sampling-mode=sw -knob enable-stack-collection=true -result-dir $W/r_hot -- $L/build_mpi/Prog/ALF.out > $W/hot.out 2>&1; echo "rc=$?"
echo "=== top-down (inclusive)"
vtune -report top-down -result-dir $W/r_hot -call-stack-mode all -column="Total Time:Self,Total Time:Total" -limit 70 2>&1 | cut -c1-200 | head -100
echo "=== flat self"
vtune -report hotspots -result-dir $W/r_hot -limit 25 2>&1 | cut -c1-200 | head -40
echo "done $(date -Is)"

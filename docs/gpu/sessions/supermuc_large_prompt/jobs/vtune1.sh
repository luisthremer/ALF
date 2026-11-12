#!/bin/bash
#SBATCH --job-name=lp_vtune
#SBATCH --output=lp_vtune.%j.log
#SBATCH --error=lp_vtune.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Task 1: VTune GPU-offload + hotspots of ONE rank (L=24 reference parameters, 3 sweeps, HDF5 binary)
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
module load intel-oneapi-vtune/2024.2.1
L=$SCRATCH/large_prompt; J=$L/jobs; W=$L/work_vtune; rm -rf $W; mkdir -p $W/run; cd $W/run
cp $L/ref_L24/seeds seeds
sed -e 's/NSweep *= *[0-9]*/NSweep = 3/' -e 's/NBin *= *[0-9]*/NBin = 1/' -e 's/CPU_MAX *= *[0-9.]*/CPU_MAX = 0.0/' $L/ref_L24/parameters > parameters
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
echo "host $(hostname) $(date -Is)"; which vtune
vtune -collect gpu-offload -result-dir $W/r_gpuoff -- mpiexec -n 1 $L/build_mpi/Prog/ALF.out > $W/gpuoff.out 2>&1; echo "gpu-offload rc=$?"
tail -5 $W/gpuoff.out | cut -c1-200
vtune -report summary -result-dir $W/r_gpuoff 2>&1 | head -80 | cut -c1-180
echo "=== computing tasks (GPU) top"
vtune -report hotspots -group-by=computing-task -result-dir $W/r_gpuoff -limit 25 2>&1 | head -45 | cut -c1-200
echo "=== host hotspots (user-mode sampling), separate run"
vtune -collect hotspots -result-dir $W/r_hot -- mpiexec -n 1 $L/build_mpi/Prog/ALF.out > $W/hot.out 2>&1; echo "hotspots rc=$?"
vtune -report hotspots -result-dir $W/r_hot -limit 30 2>&1 | head -50 | cut -c1-200
echo "done $(date -Is)"

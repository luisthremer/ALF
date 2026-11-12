#!/bin/bash
#SBATCH --job-name=lp_mpi4
#SBATCH --output=lp_mpi4.%j.log
#SBATCH --error=lp_mpi4.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# Step 3: 4 MPI ranks (walkers) on one node; rank->device mapping; fewer devices than ranks; no device
source $SCRATCH/large_prompt/jobs/env_intel_pvc.sh
J=$SCRATCH/large_prompt/jobs
export W=$SCRATCH/large_prompt/work_mpi4 CPUBIN=$SCRATCH/large_prompt/build_cpu-txt/Prog/ALF.out GPUBIN=$SCRATCH/large_prompt/build_gpu-txt/Prog/ALF.out
export SEEDS=$SCRATCH/large_prompt/src/Scripts_and_Parameters_files/Start/seeds
mkdir -p $W; cd $J
echo "host $(hostname) $(date -Is)"; xpu-smi discovery 2>&1 | grep -c "Device Name"
mon() { # sample GPU utilisation + memory per device while a run is active
  for i in 1 2 3 4 5 6; do sleep 8; echo "-- xpu-smi sample $i"; xpu-smi dump -d -1 -m 0,18 -n 1 2>&1 | tail -10; done > $W/xpu_$1.txt 2>&1
}
P="python3 run_cmp_sng.py"
echo "=== A: default mapping (rank mod 8 devices -> devices 0..3 = 2 cards)"
mon A & MP=$!
$P mpi4_A L=8 mz=F scheme=submatrix ndelay=8 nsweep=20 nbin=3 U=5.3d0 symm=T beta=8.d0 nphi=32 nr=4
wait $MP
echo "=== B: one rank per card (ZE_AFFINITY_MASK=2*rank)"
mon B & MP=$!
WRAP=$J/bind_gpu_per_card.sh $P mpi4_B L=8 mz=F scheme=submatrix ndelay=8 nsweep=20 nbin=3 U=5.3d0 symm=T beta=8.d0 nphi=32 nr=4 only=gpu
wait $MP
python3 compare.py $W/runs/mpi4_A/cpu $W/runs/mpi4_B/gpu | tail -3
echo "=== C: only 2 devices visible (ZE_AFFINITY_MASK=0,1), 4 ranks -> expect 2 ranks per device"
ZE_AFFINITY_MASK=0,1 $P mpi4_C L=8 mz=F scheme=submatrix ndelay=8 nsweep=20 nbin=3 U=5.3d0 symm=T beta=8.d0 nphi=32 nr=4 only=gpu
grep -h "omp devices\|OpenMP device" $W/runs/mpi4_C/gpu/run.log | sort | uniq -c
python3 compare.py $W/runs/mpi4_A/cpu $W/runs/mpi4_C/gpu | tail -2
echo "=== D: no offload device (ONEAPI_DEVICE_SELECTOR=opencl:cpu), OMP_TARGET_OFFLOAD=MANDATORY"
ONEAPI_DEVICE_SELECTOR=opencl:cpu mpiexec -n 2 $GPUBIN 2>&1 | head -8
echo "done $(date -Is)"

#!/bin/bash -l
# L=8 correctness smoke test: same-compiler CPU reference (ALF_cpu, GPU=none) vs GPU build, same seeds, 1 rank each.
#SBATCH -A b133ae
#SBATCH -p preempt
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=32
#SBATCH -t 00:30:00
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
# pass GPUBIN/CPUBIN/NVHPC/OMPI via: sbatch --export=ALL,GPUBIN=...  (needs removing --export=NONE on the command line)
WS=/hnvme/workspace/b133ae29-ALF_GPU; S=$WS/ALF/docs/gpu/scripts/helma
module load gpu-env/2025 nvhpc/${NVHPC:-24.11} openmpi/${OMPI:-5.0.5-nvhpc24.11-cuda}
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
R=$WS/smoke/runs/$SLURM_JOB_ID; mkdir -p $R
for kind in cpu gpu; do
  python3 $S/setup_pp.py $R/$kind 8 10 64 0 1 5.0 --nbin 4 --master smoke
  bin=${CPUBIN:-$WS/ALF_cpu/Prog/ALF.out}; [ $kind = gpu ] && bin=${GPUBIN:-$WS/ALF/Prog/ALF.out}
  ( cd $R/$kind && t0=$(date +%s.%N); mpirun --mca smsc ^knem -n 1 $bin > run.log 2>&1; rc=$?; echo "$kind rc=$rc wall=$(echo "$(date +%s.%N) - $t0" | bc) s" )
done
/hnvme/workspace/b133ae29-CLAUDE/venv_h5/bin/python $S/cmp_h5.py $R/cpu $R/gpu
grep -h -A1 "Precision Green" $R/gpu/Temp_0/info | head -3; grep -h "GPU path\|GPU observables" $R/gpu/Temp_0/info | cut -c1-150

#!/bin/bash
#SBATCH --job-name=p8_g4
#SBATCH --time=00:50:00
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p8/obs; cd $P/jobs
source p8lib.sh; export W=$C/work_g4; mkdir -p $W; B=$C/build_aottt/Prog/ALF.out
echo "######## A: one chain (1 rank, tile 0, OMP 8), V7, alternating OBSPAR 0/1, 3 sweeps, phase timers"
for r in 1 2 3; do for p in 0 1; do BIN=$B one A_p${p}_$r V7 3 ALF_GPU_OBSPAR=$p; done; done
echo "######## B/C: per card from the CPU-time difference of 1 vs 3 sweeps (as p6/m_tput.sh), alternating OBSPAR 0/1"
run() { local nr=$1 omp=$2 ccs=$3 v=$4 p=$5 tag=$6
  for ns in 1 3; do env $(var $v) $ccs OMPT=$omp ALF_GPU_OBSPAR=$p python3 $J/run_ref_sng.py ${tag}_s$ns $B $nr $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${tag}_s1),$(cpu ${tag}_s3); m=(y-x)/2; print(\"RESULT %-12s ranks %2d OMP %2d OBSPAR $p  s/sweep/chain %.3f  sweeps/s/node %.4f  per card %.4f\"%(\"$tag\",$nr,$omp,m,$nr/m,$nr/m/4))"
  pg ${tag}_s3; tim ${tag}_s3; }
for r in 1 2; do
  for p in 0 1; do run 8 14 "" V7 $p T1_V7_p${p}_$r; done
  for p in 0 1; do run 32 3 "$CCS" V7 $p T4_V7_p${p}_$r; done
  for p in 1 0; do run 32 3 "$CCS" V2 $p T4_V2_p${p}_$r; done
done
echo "end $(date -Is)"

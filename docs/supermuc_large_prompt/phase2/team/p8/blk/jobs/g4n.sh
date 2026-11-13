#!/bin/bash
#SBATCH --job-name=p8b_g4n
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# G4n + G3: node layout (32 ranks, 4 per tile, CCS 4, OMP 3), V7, L24, old vs new alternating in one job, N_delay 64 and 32;
# sweeps/s per card from the CPU-time difference of a 1-sweep and a 3-sweep run (as p6/m_tput.sh); Precision Green of the 3-sweep runs (G3);
# old vs new run directories compared bitwise (same binary, same seeds).
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; cd $L0/p8/blk/jobs; source p8lib.sh
export W=$P8/work_g4n; mkdir -p $W
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
runn() { local m=$1 nd=$2 tag=$3
  for ns in 1 3; do env $(var V7) $CCS OMPT=3 LSIZE=24 NDELAY=$nd ALFSH_BLKPACK=$m python3 $J/run_ref_sng.py ${tag}_s$ns $BIN 32 $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${tag}_s1),$(cpu ${tag}_s3); m=(y-x)/2; print(\"RESULT %-14s ranks 32 (4/tile) nd $nd BLKPACK $m  s/sweep/chain %.3f  sweeps/s/node %.3f  per card %.3f\"%(\"$tag\",m,32/m,32/m/4))"
  pg ${tag}_s3; tim ${tag}_s3; p3 ${tag}_s3 3; }
for m in 0 1 2 2 1 0; do runn $m 64 n${m}_nd64_a$RANDOM; done
for m in 1 2 0; do runn $m 32 n${m}_nd32_a$RANDOM; done
echo "=== bitwise old vs new (3-sweep node runs, same nd): compare.py tol 0"
for nd in 64 32; do a=$(ls -d $W/runs/n0_nd${nd}_a*_s3 | head -1); b=$(ls -d $W/runs/n1_nd${nd}_a*_s3 | head -1); echo "nd=$nd $a $b"; python3 $J/compare.py $a $b --tol 0 | tail -3; done
echo "end $(date -Is)"

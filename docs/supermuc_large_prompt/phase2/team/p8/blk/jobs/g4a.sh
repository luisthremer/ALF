#!/bin/bash
#SBATCH --job-name=p8b_g4a
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
# G4a: old (ALFSH_BLKPACK unset) vs new (=1; =2 is the non-bitwise sub-switch), alternating (ABBA) in one job, V7, L24: one chain on tile 0 (in-loop phase + block timers, 3 sweeps)
# and one rank per tile (8 ranks, OMP 14, CPU-time difference 1 vs 3 sweeps, sweeps/s per card).  N_delay 32 and 64 (production).
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; cd $L0/p8/blk/jobs; source p8lib.sh
export W=$P8/work_g4a; mkdir -p $W
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
echo "=== one chain, tile 0"
for nd in 32 64; do for m in 0 1 2 2 1 0; do
  t=c${m}_nd${nd}_$RANDOM; one $t V7 3 LSIZE=24 NDELAY=$nd ALFSH_BLKPACK=$m > $W/$t.out 2>&1; echo "m=$m nd=$nd"; grep -h -E "RESULT|PG|Precision" $W/$t.out; p3 $t 3
done; done
echo "=== one rank per tile (8 ranks, OMP 14)"
runt() { local m=$1 nd=$2 tag=$3
  for ns in 1 3; do env $(var V7) OMPT=14 LSIZE=24 NDELAY=$nd ALFSH_BLKPACK=$m python3 $J/run_ref_sng.py ${tag}_s$ns $BIN 8 $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${tag}_s1),$(cpu ${tag}_s3); m=(y-x)/2; print(\"RESULT %-14s ranks 8 (1/tile) nd $nd BLKPACK $m  s/sweep/chain %.3f  sweeps/s/node %.3f  per card %.3f\"%(\"$tag\",m,8/m,2/m))"
  pg ${tag}_s3; }
for m in 0 1 2 2 1 0; do runt $m 64 t${m}_nd64_$RANDOM; done
echo "end $(date -Is)"

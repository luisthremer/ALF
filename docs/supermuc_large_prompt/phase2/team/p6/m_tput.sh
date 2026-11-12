#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
echo "start $(date -Is) $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
# Throughput vs walkers per tile (V7 = staging+async+devqr2+devunmqr), L=24, per card from CPU-time difference 1 vs 3 sweeps.
export W=$SCRATCH/large_prompt/team/p6/work_tput; B=$C/build_aottt/Prog/ALF.out
run() { local nr=$1 omp=$2 ccs=$3 tag=$4
  for ns in 1 3; do env $(var V7) $ccs OMPT=$omp python3 $J/run_ref_sng.py ${tag}_s$ns $B $nr $ns 1 > /dev/null; done
  python3 -c "x,y=$(cpu ${tag}_s1),$(cpu ${tag}_s3); m=(y-x)/2; print(\"RESULT %-6s ranks %2d (%d/tile) OMP %2d  s/sweep/chain %.2f  sweeps/s/node %.3f  per card %.3f\"%(\"$tag\",$nr,$nr//8,$omp,m,$nr/m,$nr/m/4))"
  pg ${tag}_s3; }
run 8 14 "" T1
run 16 6 "$CCS" T2
run 32 3 "$CCS" T4
run 64 1 "$CCS" T8
echo "end $(date -Is)"

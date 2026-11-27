#!/bin/bash
# P2.6.1 check on the Max 1550: the H100 case list docs/ai/phase2/p26_runs/cases_h100.py (Kondo, tV, Spin_Peierls x standard /
# delayed / submatrix, small and big, Hubbard regression) with the P2.5 runner (GPU path vs ALF_GPU_DISABLE=1, same seeds, e2e rule),
# the cases spread round-robin over the 8 tiles (one runner per tile, OMP_NUM_THREADS=1 as set by the runner). Default GPU=intel
# build (AOT; Intel switches at their defaults).
#SBATCH -J p2i_p26
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/p26_$SLURM_JOB_ID
. $S/env_smuc2.sh
mkdir -p $R
P=$B/src_intel/docs/ai/phase2
labels=$(cd $P/p26_runs && python3 -c "from cases_h100 import CASES; print(' '.join(c[0] for c in CASES))")
i=0; declare -a grp
for l in $labels; do t=$((i % 8)); grp[$t]="${grp[$t]:+${grp[$t]},}$l"; i=$((i+1)); done
for t in 0 1 2 3 4 5 6 7; do
  ( export ZE_AFFINITY_MASK=$t
    python3 $P/p25_runs/p25_run.py --exe $B/ALF-intel.out --start $B/src_intel/Scripts_and_Parameters_files/Start \
      --work $R/w$t --only ${grp[$t]} $P/p26_runs/cases_h100.py > $R/runner_$t.log 2>&1 ) &
done
wait
cat $R/w*/results.tsv | sed 's/^[^\t]*\t/max1550\t/' > $R/results_max1550.tsv
cat $R/results_max1550.tsv

#!/bin/bash -l
# A/B test: cost of OMP_DISPLAY_AFFINITY=TRUE (set by the campaign job templates) on the sweep time. One node, the smuc_ndscan
# hubL40 2x24 inputs (same parameters and seeds, same environment as the campaign profile), runs in ABBA order:
# delayed N_delay 64 on/off/off/on, standard on/off. Each run in its own directory under $OUT; timers.txt per run.
# Usage (SuperMUC-NG, from the work directory): sbatch --account=pn73xu ALF/testsuite/lowrank/cluster/affinity_overhead.sbatch
#SBATCH --job-name=aff_ab
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:30:00
#SBATCH --ear=off
#SBATCH --output=affinity_ab_%j.out
W=/hppfs/work/pn73xu/di38zuw/CLAUDE
. "$W/env_alf.sh"
C="$W/campaigns/smuc_ndscan"
BIN="$W/alf_bin/ALF_timers_mpi_78bc38ee.out"
OUT="$W/affinity_ab_${SLURM_JOB_ID}"
mkdir -p "$OUT"
unset FORT_BLOCKSIZE KMP_AFFINITY MKL_VERBOSE LD_PRELOAD
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0 I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_DOMAIN=omp:compact I_MPI_DEBUG=4
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_NUM_THREADS=24 MKL_NUM_THREADS=24 MKL_DYNAMIC=FALSE
export OMP_AFFINITY_FORMAT='OMP_AFFINITY host=%H pid=%P thread=%n affinity=%A'
host="$(hostname -s)"
i=0
for spec in delayed_nd64:TRUE delayed_nd64:FALSE delayed_nd64:FALSE delayed_nd64:TRUE standard:TRUE standard:FALSE; do
    run="${spec%%:*}"; disp="${spec#*:}"; i=$((i + 1))
    d="$OUT/${i}_${run}_display${disp}"; mkdir -p "$d"
    cp "$C/runs/hubL40_L40_2x24_r1_$run/input/parameters" "$C/runs/hubL40_L40_2x24_r1_$run/input/seeds" "$d/"
    t0=$(date +%s.%N)
    (cd "$d" && OMP_DISPLAY_AFFINITY="$disp" mpiexec -hosts "$host" -n 2 -ppn 2 "$BIN" > stdout.log 2> stderr.log); rc=$?
    t1=$(date +%s.%N)
    echo "$i $run display=$disp rc=$rc wall_s=$(echo "$t1 - $t0" | bc) stdout_bytes=$(stat -c %s "$d/stdout.log")"
done

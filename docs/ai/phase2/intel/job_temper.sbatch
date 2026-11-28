#!/bin/bash
# P2.6 tempering on the GPU path, Max 1550 (GPU=intel, Tempering build, AOT): GPU path vs ALF_GPU_DISABLE=1, same seeds, with the
# runner ../p25_runs/p25_run.py; one tile per rank (ZE_AFFINITY_MASK = local rank, Intel MPI), 2 and 4 replicas.
#SBATCH -J p2i_tp
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH -o %x.%j.out
module load slurm_setup
S=$HOME/p2_intel/scripts; B=$HOME/p2_intel/build; R=$SCRATCH/p2_intel/tp_$SLURM_JOB_ID
. $S/env_smuc2.sh; unset ZE_AFFINITY_MASK
cd $B/src_tp && . ./configure.sh SUPERMUCNG2 Tempering GPU=intel NO-INTERACTIVE > /dev/null 2>&1; cd - > /dev/null
mkdir -p $R
printf '#!/bin/bash\nexport ZE_AFFINITY_MASK=$MPI_LOCALRANKID\nexec %s\n' "$B/src_tp/Prog/ALF.out" > $R/rank_wrap.sh; chmod +x $R/rank_wrap.sh
P=$B/src_tp/docs/ai/phase2; C=$P/p26_runs/tempering/cases_tempering.py; ST=$B/src_tp/Scripts_and_Parameters_files/Start
export P25_MPIRUN=mpiexec
python3 $P/p25_runs/p25_run.py --build max1550-np2 --exe $R/rank_wrap.sh --start $ST --work $R/np2 --mpi 2 --groups 2 --timeout 1200 $C
python3 $P/p25_runs/p25_run.py --build max1550-np4 --exe $R/rank_wrap.sh --start $ST --work $R/np4 --mpi 4 --groups 4 --timeout 1200 \
  --only tp-mz-std,tp-mz-sub4,tp-mz-nex3,tp-mz-8x8-sub16,tp-neg-temper $C
grep -h "OpenMP device\|GPU Intel" $R/np2/tp-mz-std/active/Temp_*/info
cat $R/np2/results.tsv $R/np4/results.tsv > $R/results_max1550_temper.tsv

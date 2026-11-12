#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=03:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Suspect-tile test: 8 ranks, 4 on device 1 (suspect tile) and 4 on device 0 (same card, other tile), CCS 4 so the 4 ranks per tile
# run concurrently (as production), OMP 3, STRESS_SEC (default 3 h). Every iteration checked bitwise against iteration 1.
D=$SCRATCH/large_prompt/diag; source $D/env_intel_pvc.sh; cd $D
export OMP_NUM_THREADS=3 ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 STRESS_SEC=${1:-10800}
printf '#!/bin/bash\nr=${PMI_RANK:-$SLURM_PROCID}\nexport STRESS_DEV=$(( r %% 2 == 0 ? 1 : 0 ))\nexec "$@"\n' > $D/devwrap.sh; chmod +x $D/devwrap.sh
echo "node $(hostname) $(date -Is)"
mpiexec -n 8 $D/devwrap.sh ./t_stress > stress_long.$SLURM_JOB_ID.out 2>&1
grep -c ^BAD stress_long.$SLURM_JOB_ID.out; grep ^BAD stress_long.$SLURM_JOB_ID.out | head -20
grep -E "^STRESS rank|STRESS_RESULT" stress_long.$SLURM_JOB_ID.out | sort -k3n
echo "end $(date -Is)"

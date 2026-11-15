#!/bin/bash
#SBATCH --job-name=p9_cal
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:45:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Step 1 (calibration) on the dev seed sets k = 1 (original), 3, 11, 23, 37 (the p5/g5 rotations): oneMKL, V7, P0 (unedited library copy), P12/P11/P10 (planted tau)
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
cd $SCRATCH/large_prompt/team/p9/l6t/jobs; source lib.sh; export W=$P/work_cal; mkdir -p $W
md5sum $CPUBIN $GPUBIN $P/build_txt/Libraries/Modules/libalfdevqr.so $P/plant/*/libalfdevqr.so
w() { local k=$1 t=$2; mkseed $k; echo "seed set k=$k first seed $(head -1 $W/seeds/seeds_$k) tile $t"
  runc $k L6T; runc $k L6F
  for v in oneMKL V7 P0 P12 P11 P10; do rung $k L6T $v $t; cmpr $k L6T $v; done
  for v in V7 P0 P12 P11 P10; do rung $k L6F $v $t; cmpr $k L6F $v; done; }
i=0; for k in 1 3 11 23 37; do w $k $i > $W/w_$k.out 2>&1 & i=$((i+1)); done; wait
cat $W/w_*.out
grep -h ^RES $W/w_*.out > $W/res_cal.txt
python3 $J/rule.py $W/res_cal.txt V7 P0 P12 P11 P10
echo "end $(date -Is)"

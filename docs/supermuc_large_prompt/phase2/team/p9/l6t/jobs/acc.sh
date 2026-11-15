#!/bin/bash
#SBATCH --job-name=p9_acc
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Step 2 (acceptance, judged once): 20 new seed sets k = 41, 43, ..., 79 (rotations of the original seeds file; dev sets were k = 1, 3, 11, 23, 37)
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
cd $SCRATCH/large_prompt/team/p9/l6t/jobs; source lib.sh; export W=$P/work_acc; mkdir -p $W
md5sum $CPUBIN $GPUBIN $P/build_txt/Libraries/Modules/libalfdevqr.so
KS=$(seq 41 2 79)
w() { local t=$1; shift; for k in "$@"; do mkseed $k; echo "seed set k=$k first seed $(head -1 $W/seeds/seeds_$k) tile $t"
  runc $k L6T; runc $k L6F; runc $k L12F
  for v in oneMKL V7; do rung $k L6T $v $t; cmpr $k L6T $v; done
  rung $k L6F V7 $t; cmpr $k L6F V7; rung $k L12F V7 $t; cmpr $k L12F V7; done; }
NW=${NW:-5}; for t in $(seq 0 $((NW-1))); do w $t $(echo $KS | tr ' ' '\n' | awk -v t=$t -v n=$NW '(NR-1)%n==t') > $W/w_$t.out 2>&1 & done; wait
cat $W/w_*.out
grep -h ^RES $W/w_*.out > $W/res_acc.txt
python3 $J/rule.py $W/res_acc.txt V7
echo "end $(date -Is)"

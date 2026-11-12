#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --exclude=i20r02c06s02
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
P7=$SCRATCH/large_prompt/team/p7; export W=$P7/work_probe2; mkdir -p $W; cd $P7/jobs
runw() { local tag=$1 wrap=$2; shift 2; env $(var V7) "$@" OMPT=8 ZE_AFFINITY_MASK=0 python3 $J/run_ref_sng.py $tag $C/build_aottt/Prog/ALF.out 1 1 1 $wrap; }
echo "--- trivial programs with the tracing layer"
icx -fiopenmp -fopenmp-targets=spir64 t_omp.c -o t_omp 2>&1 | tail -3
icpx -fsycl t_sycl.cpp -o t_sycl 2>&1 | tail -3
for t in t_omp t_sycl; do for tl in 0 1; do echo "$t TL=$tl"; ZE_AFFINITY_MASK=0 ZE_ENABLE_TRACING_LAYER=$tl OMP_TARGET_OFFLOAD=MANDATORY ./$t 2>&1 | tail -2; done; done
echo "--- gdb on ALF with tracing layer"
which gdb-oneapi
runw G1 $P7/jobs/gdbwrap.sh ZE_ENABLE_TRACING_LAYER=1
grep -v "^ *$" $W/runs/G1/run.log | sed -n 6,80p | cut -c1-220
echo "--- advisor roofline GPU"
which advisor
runw ADV $P7/jobs/advwrap.sh
grep -v "^ *$" $W/runs/ADV/run.log | tail -30 | cut -c1-250
echo "end $(date -Is)"

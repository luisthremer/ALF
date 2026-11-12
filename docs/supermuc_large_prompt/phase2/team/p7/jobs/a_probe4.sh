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
P7=$SCRATCH/large_prompt/team/p7; export W=$P7/work_probe4; mkdir -p $W; cd $P7/jobs
runw() { local tag=$1 wrap=$2; shift 2; env $(var V7) "$@" OMPT=8 ZE_AFFINITY_MASK=0 python3 $J/run_ref_sng.py $tag $C/build_aottt/Prog/ALF.out 1 1 1 $wrap; }
ls -d /lrz/sys/intel/oneapi_2025.*/pti/latest/lib 2>&1
echo "--- T1: SYCL_UR_USE_LEVEL_ZERO_V2=1 (no tracing)"
one T1 V7 1 SYCL_UR_USE_LEVEL_ZERO_V2=1
echo "--- T1b: baseline again"
one T1b V7 1
echo "--- T2: V2 + tracing layer"
one T2 V7 1 SYCL_UR_USE_LEVEL_ZERO_V2=1 ZE_ENABLE_TRACING_LAYER=1
tail -3 $W/runs/T2/run.log | cut -c1-200
echo "--- T3: V2 + tracing layer + PTI view (2025.1.0 pti 0.11)"
runw T3 $P7/jobs/ptiwrap.sh SYCL_UR_USE_LEVEL_ZERO_V2=1 ZE_ENABLE_TRACING_LAYER=1
tim T3; grep -h "^PTI" $W/runs/T3/run.log | head -50 | cut -c1-250
tail -4 $W/runs/T3/run.log | cut -c1-200
echo "--- T4: advisor survey without mpiexec"
runw T4 /bin/true
cd $W/runs/T4; source /lrz/sys/intel/oneapi_2025.3.1/advisor/latest/advisor-vars.sh >/dev/null 2>&1
env $(var V7) OMP_NUM_THREADS=8 ZE_AFFINITY_MASK=0 timeout 400 advisor --collect=survey --profile-gpu --project-dir=$W/adv4 -- $C/build_aottt/Prog/ALF.out 2>&1 | tail -15 | cut -c1-250
echo "end $(date -Is)"

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
#SBATCH --time=00:15:00
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
P7=$SCRATCH/large_prompt/team/p7; export W=$P7/work_probe5; mkdir -p $W; cd $P7/jobs
icpx -fsycl t_sycl.cpp -o t_sycl 2>&1 | tail -3
echo "--- t_sycl default, SYCL_UR_TRACE=2"; SYCL_UR_TRACE=2 ./t_sycl 2>&1 | head -30 | cut -c1-200
echo "--- t_sycl V2=1"; SYCL_UR_USE_LEVEL_ZERO_V2=1 ./t_sycl 2>&1 | head -5
echo "--- t_sycl V2=0"; SYCL_UR_USE_LEVEL_ZERO_V2=0 ./t_sycl 2>&1 | head -5
echo "--- t_sycl ONEAPI_DEVICE_SELECTOR=level_zero:gpu"; ONEAPI_DEVICE_SELECTOR=level_zero:gpu ./t_sycl 2>&1 | head -5
echo "--- env"; env | grep -E "^(ZE|ZES|SYCL|UR|ONEAPI|NEOReadDebug|ZET|LIBOMPTARGET|OMP_|KMP)" | sort
echo "--- ALF V7 + V2=1, stderr"; one S1 V7 1 SYCL_UR_USE_LEVEL_ZERO_V2=1; grep -v "^ *$" $W/runs/S1/run.log | tail -12 | cut -c1-250
echo "--- ALF V7, SYCL_UR_TRACE=2 head (baseline)"; one S2 V7 1 SYCL_UR_TRACE=2; grep -n -i -E "adapter|error|uninit|platform" $W/runs/S2/run.log | head -20 | cut -c1-250
echo "end $(date -Is)"

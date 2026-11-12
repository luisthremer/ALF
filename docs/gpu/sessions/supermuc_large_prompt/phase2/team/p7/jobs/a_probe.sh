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
#SBATCH --time=00:20:00
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
P7=$SCRATCH/large_prompt/team/p7; export W=$P7/work_probe; mkdir -p $W
echo "--- system"; cat /proc/sys/dev/i915/perf_stream_paranoid /proc/sys/kernel/perf_event_paranoid 2>&1; ls /sys/class/drm | head; id -Gn
which sycl-ls zeinfo 2>&1; sycl-ls 2>&1 | head -12
ldconfig -p | grep -E "libmd|igdml|metrics|libze_tracing|libpti" | head
B=$C/build_aottt/Prog/ALF.out
echo "--- P0 baseline V7 1 sweep + ALFSH_TIMING"
one P0 V7 1 ALFSH_TIMING=1
grep -h ALFSH_TIMING $W/runs/P0/run.log | head -30
echo "--- P1 LIBOMPTARGET_PLUGIN_PROFILE=T"
one P1 V7 1 LIBOMPTARGET_PLUGIN_PROFILE=T
grep -v -E "^ *[0-9]+ *$" $W/runs/P1/run.log | grep -i -E "kernel|profile|plugin|total" | head -40
echo "--- P2 LIBOMPTARGET_INFO=-1"
one P2 V7 1 LIBOMPTARGET_INFO=-1
wc -l $W/runs/P2/run.log; grep -v "^ *$" $W/runs/P2/run.log | sed -n 8,60p | cut -c1-220
echo "--- P3 ZE_ENABLE_TRACING_LAYER=1 alone"
one P3 V7 1 ZE_ENABLE_TRACING_LAYER=1
tail -5 $W/runs/P3/run.log | cut -c1-200
echo "end $(date -Is)"

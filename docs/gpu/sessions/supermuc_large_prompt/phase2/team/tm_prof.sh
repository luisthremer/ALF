#!/bin/bash
#SBATCH --job-name=tm_prof
#SBATCH --output=tm_prof.%j.log
#SBATCH --error=tm_prof.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
# Own profile of ONE chain on ONE tile, L=24 reference parameters, 2 sweeps: VTune gpu-offload (device tasks, host API time)
# and hotspots with stacks (what the host computes). usage: sbatch tm_prof.sh [BUILD] [OMPT]
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
VT=${VT:-/lrz/sys/intel/oneapi_2025.3.1/vtune/latest}; source $VT/env/vars.sh >/dev/null 2>&1 || export PATH=$VT/bin64:$PATH
B=${1:-aot}; T=${2:-8}
L0=$SCRATCH/large_prompt/team; W=$L0/work_prof_${B}_t$T${TAG}; rm -rf $W; mkdir -p $W/run; cd $W/run
cp $L0/ref_L24/seeds seeds
sed -e 's/NSweep *= *[0-9]*/NSweep = ${NSW:-2}/' -e 's/NBin *= *[0-9]*/NBin = 1/' -e 's/CPU_MAX *= *[0-9.]*/CPU_MAX = 0.0/' $L0/ref_L24/parameters > parameters
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=$T ZE_AFFINITY_MASK=0
echo "host $(hostname) $(date -Is)"; which vtune; vtune --version | head -1
/usr/bin/time -f "plain wall %e s" $L0/build_$B/Prog/ALF.out > $W/plain.out 2>&1; tail -1 $W/plain.out
vtune -collect gpu-offload -result-dir $W/r_off -- $L0/build_$B/Prog/ALF.out > $W/off.out 2>&1; echo "gpu-offload rc=$?"
vtune -report summary -result-dir $W/r_off 2>&1 | cut -c1-200 | head -90
echo "=== GPU computing tasks"
vtune -report hotspots -group-by=computing-task -sort-desc="Total Time" -result-dir $W/r_off -limit 40 2>&1 | cut -c1-220 | head -60
echo "=== host tasks / API"
vtune -report hotspots -group-by=computing-task -column="Instance Count"  -result-dir $W/r_off -limit 30 2>&1 | cut -c1-220 | head -45
[ -n "$NOHOT" ] && exit 0; rm -rf run2; mkdir run2; cp parameters seeds run2/; cd run2
vtune -collect hotspots -knob sampling-mode=sw -knob enable-stack-collection=true -result-dir $W/r_hot -- $L0/build_$B/Prog/ALF.out > $W/hot.out 2>&1; echo "hotspots rc=$?"
echo "=== hotspots by module"; vtune -report hotspots -group-by=module -result-dir $W/r_hot -limit 20 2>&1 | cut -c1-200 | head -30
echo "=== hotspots flat self"; vtune -report hotspots -result-dir $W/r_hot -limit 40 2>&1 | cut -c1-200 | head -50
vtune -report top-down -result-dir $W/r_hot -call-stack-mode all -column="CPU Time:Self,CPU Time:Total" -limit 3000 > $W/topdown.txt 2>&1; wc -l $W/topdown.txt
echo "done $(date -Is)"

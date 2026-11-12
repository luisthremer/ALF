#!/bin/bash
#SBATCH --job-name=p3_v3
#SBATCH --output=p3_v3.%j.log
#SBATCH --error=p3_v3.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14
# G2 (c), second planted defect: no drain at the end-of-slice flush_all (queued block work not waited before the OpenMP hopping wraps)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p3; J=$L0/jobs; cd $P
rm -rf src_plant2 build_txtplant2; cp -r src src_plant2
python3 - <<'PY'
import os
p=os.environ['SCRATCH']+'/large_prompt/team/p3/src_plant2/Prog/gpu_lowrank_mod.F90'
s=open(p).read()
old="""    if (.not. inited .or. mode == LR_STD) return
#if defined(ALF_LA_ONEMKL)
    call ab_drain()
#endif
    if (fast) call gpu_lr_sync()"""
assert old in s
s=s.replace(old,"""    if (.not. inited .or. mode == LR_STD) return
    if (fast) call gpu_lr_sync()""")
open(p,'w').write(s)
PY
SRC=$P/src_plant2 BLD=$P/build_txtplant2 bash build_intel_p3.sh txtnew > b_txtplant2.out 2>&1
ls -l build_txtplant2/Prog/ALF.out || exit 1
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
OLDB=$L0/build_txtnew/Prog/ALF.out
for r in 1 2 3; do for a in "L6 L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" "L12 L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"; do
  W=$P/work_v3 CPUBIN=$OLDB GPUBIN=$P/build_txtplant2/Prog/ALF.out CMP_TOL=0 python3 run_cmp_gg.py PL2_${a%% *}_$r ${a#* } 2>&1 | grep -E "PASS|FAIL|deviation|Terminating" | head -3 | sed "s/^/PLANT2 ${a%% *} #$r: /"
done; done
echo "done $(date -Is)"

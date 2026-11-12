#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; C=$L0/comb; cd $C/jobs
# Gate after the devqr race fix: same seed vs CPU with V4 and V3, single tile L24/L32 timing + Precision Green V2/V4
source c_lib.sh; export W=$C/work_gone
CMPA="W=$C/work_gcmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py"
for v in V4 V3; do
  echo "######## CHECK $v vs CPU, same seed"
  env $(var $v) CMP_TOL=1e-13 bash -c "$CMPA ${v}_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18; $CMPA ${v}_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" 2>&1 | grep -E "PASS|FAIL"
  env $(var $v) CMP_TOL=1e-10 bash -c "$CMPA ${v}_L6T L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" 2>&1 | grep -E "PASS|FAIL"
done
for r in 1 2; do for v in V2 V4; do one L24_${v}_$r $v 4; done; done
for v in V2 V4; do one L32_${v}_1 $v 2 LSIZE=32; done
echo "done $(date -Is)"

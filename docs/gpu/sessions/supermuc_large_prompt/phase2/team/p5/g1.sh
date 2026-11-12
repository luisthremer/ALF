#!/bin/bash
#SBATCH --job-name=p5_g1
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p5; cd $P/jobs
source p5lib.sh; export W=$C/work_g1
echo "host $(hostname) $(date -Is)"
echo "######## default build unchanged: p5 txtnew vs comb txtnew, no switches, bitwise (CMP_TOL=0)"
for a in "g_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" "g_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"; do
  env -u ALFSH_DEVQR -u ALFSH_DEVUNMQR W=$W/gg CPUBIN=$L0/comb/build_txtnew/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out CMP_TOL=0 SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py $a 2>&1 | grep -E "PASS|FAIL"
done
echo "######## Gate C same seed vs CPU, V7 (async + staging + devqr 2 + devunmqr)"
CMPA="W=$W/gcmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py"
for v in V7 V6; do
  env $(var $v) CMP_TOL=1e-13 bash -c "$CMPA ${v}_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18; $CMPA ${v}_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" 2>&1 | grep -E "PASS|FAIL"
  env $(var $v) CMP_TOL=1e-10 bash -c "$CMPA ${v}_L6T L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" 2>&1 | grep -E "PASS|FAIL"
done
echo "######## Gate B inside ALF (ALFSH_TIMING, no ASYNCBLK): V3 = devqr 1, V5 = devqr 2, V6 = devqr 2 + devunmqr; marginal per call from 2 vs 4 sweeps"
for L in 24 32; do for v in V3 V5 V6; do
  for ns in 2 4; do t=B_L${L}_${v}_s$ns; if [ $L = 32 ]; then one $t $v $ns LSIZE=32 ALFSH_TIMING=1 >/dev/null; else one $t $v $ns ALFSH_TIMING=1 >/dev/null; fi; done
  echo "### L=$L $v"; pc B_L${L}_${v}_s2 B_L${L}_${v}_s4; pg B_L${L}_${v}_s4; tim B_L${L}_${v}_s4
done; done
echo "done $(date -Is)"

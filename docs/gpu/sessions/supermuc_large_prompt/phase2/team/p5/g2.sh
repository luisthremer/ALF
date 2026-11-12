#!/bin/bash
#SBATCH --job-name=p5_g2
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
source p5lib.sh; export W=$C/work_g2
echo "host $(hostname) $(date -Is)"
CMPA="W=$W/gcmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py"
echo "######## planted defects, same seed vs CPU, V7 (must FAIL)"
for pl in tau conj; do
  env $(var V7) ALFSH_DEVQR_LIB=$C/plant/$pl/libalfdevqr.so CMP_TOL=1e-13 bash -c "$CMPA P${pl}_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18; $CMPA P${pl}_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" 2>&1 | grep -E "PASS|FAIL" | sed "s/^/[plant $pl] /"
done
echo "######## node config (32 ranks, 4/tile, CCS 4, OMPT=3), Precision Green of all ranks, 3 sweeps"
for r in 1 2; do node N_V7_$r V7 3; pg N_V7_$r; done
node N_V4_1 V4 3; pg N_V4_1
for pl in race_tau race_ex race_v tau; do
  node NP_$pl V7 3 ALFSH_DEVQR_LIB=$C/plant/$pl/libalfdevqr.so; echo "[plant $pl]"; pg NP_$pl
done
echo "done $(date -Is)"

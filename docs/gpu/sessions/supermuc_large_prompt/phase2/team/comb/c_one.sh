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
# single chain L=24 (one tile, 8 threads), 5 variants x 3 rounds, + checks with the switches on
source c_lib.sh; export W=$C/work_one
echo "######## CHECK V4 (both on) vs CPU, same seed"
env $(var V4) CMP_TOL=1e-13 bash -c "W=$C/work_cmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py c_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18; W=$C/work_cmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py c_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" 2>&1 | grep -E "PASS|FAIL"
env $(var V4) CMP_TOL=1e-10 bash -c "W=$C/work_cmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py c_L6T L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" 2>&1 | grep -E "PASS|FAIL"
echo "######## CHECK V2 (async only, combined build) vs team/build_txtnew (branch head), bitwise"
printf "#!/bin/bash\nexport $SB ALFSH_ASYNCBLK=1\nexec \"\$@\"\n" > $C/jobs/v2wrap.sh; chmod +x $C/jobs/v2wrap.sh
for a in "g_L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" "g_L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"; do
  W=$C/work_gg WRAP=$C/jobs/v2wrap.sh CPUBIN=$L0/build_txtnew/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out CMP_TOL=0 SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py $a 2>&1 | grep -E "PASS|FAIL"
done
for r in 1 2 3; do for v in V0 V1 V2 V3 V4; do one L24_${v}_$r $v 4; done; done
echo "done $(date -Is)"

# L6T same seed vs CPU for the components: V1 oneMKL, V3 devqr 1, devqr 2 only, devunmqr only, both (what is the Mz=T noise level and which kernel moves it?)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p5; cd $P/jobs
source p5lib.sh; export W=$C/work_g4
CMPA="W=$W/gcmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py"
for cfg in "ALFSH_DEVQR=0 ALFSH_DEVUNMQR=0" "ALFSH_DEVQR=1 ALFSH_DEVUNMQR=0" "ALFSH_DEVQR=2 ALFSH_DEVUNMQR=0" "ALFSH_DEVQR=0 ALFSH_DEVUNMQR=1" "ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1" "ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_QRNB=32"; do
  tg=$(echo $cfg | tr -d ' =A-Z_'); echo "### $cfg"
  env $SB $cfg CMP_TOL=1e-10 bash -c "$CMPA L6T_$tg L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" 2>&1 | grep -E "PASS|FAIL"
done
echo "done"

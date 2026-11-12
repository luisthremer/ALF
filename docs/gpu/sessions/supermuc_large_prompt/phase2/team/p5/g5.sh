# L6T noise level vs the seed: 4 other seed sets (rotated seed file), oneMKL vs devqr2+devunmqr vs CPU, same seed each
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p5; cd $P/jobs
source p5lib.sh; export W=$C/work_g5; mkdir -p $W
S0=$L0/src/Scripts_and_Parameters_files/Start/seeds
for k in 3 11 23 37; do
  (tail -n +$k $S0; head -n $((k-1)) $S0) > $W/seeds_$k
  for cfg in "ALFSH_DEVQR=0 ALFSH_DEVUNMQR=0" "ALFSH_DEVQR=1 ALFSH_DEVUNMQR=0" "ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1"; do
    tg=$(echo $cfg | tr -d ' =A-Z_'); echo -n "seedrot $k  $cfg : "
    env $SB $cfg CMP_TOL=1e-10 bash -c "W=$W/gcmp CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$C/build_txtnew/Prog/ALF.out SEEDS=$W/seeds_$k ZE_AFFINITY_MASK=0 python3 $J/run_cmp_sng.py S${k}_$tg L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" 2>&1 | grep -E "PASS|FAIL"
  done
done
echo done

# shared helpers (sourced) for the Phase 5 jobs
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; C=$L0/p5; PJ=$C/jobs
SB=LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
# variants (all with the 1 MB staging buffer): V1 oneMKL; V3 phase-4 device QR; V5 phase-5 geqrf; V6 = V5 + device unmqr; V2 async only; V4 async + V3; V7 async + V6
var() { case $1 in V1) echo "$SB ALFSH_ASYNCBLK=0 ALFSH_DEVQR=0";; V2) echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=0";;
  V3) echo "$SB ALFSH_ASYNCBLK=0 ALFSH_DEVQR=1";; V4) echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=1";;
  V5) echo "$SB ALFSH_ASYNCBLK=0 ALFSH_DEVQR=2";; V6) echo "$SB ALFSH_ASYNCBLK=0 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1";;
  V7) echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1";; esac; }
export ALFSH_QRNB=16 OMP_TARGET_OFFLOAD=MANDATORY SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
tim() { grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-14s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'; }
one() { local t=$1 v=$2 ns=$3; shift 3
  env $(var $v) "$@" OMPT=8 ZE_AFFINITY_MASK=0 python3 $J/run_ref_sng.py $t ${BIN:-$C/build_aottt/Prog/ALF.out} 1 $ns 1 > /dev/null; tim $t; grep -h "Precision Green" $W/runs/$t/info; }

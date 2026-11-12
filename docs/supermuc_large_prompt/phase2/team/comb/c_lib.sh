# shared helpers (sourced): variants V0 nostage/both off, V1 both off, V2 async, V3 devqr, V4 both
var() { case $1 in V0) echo "ALFSH_ASYNCBLK=0 ALFSH_DEVQR=0";; V1) echo "$SB ALFSH_ASYNCBLK=0 ALFSH_DEVQR=0";; V2) echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=0";;
  V3) echo "$SB ALFSH_ASYNCBLK=0 ALFSH_DEVQR=1";; V4) echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=1";; esac; }
SB=LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
export ALFSH_QRNB=16 OMP_TARGET_OFFLOAD=MANDATORY SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
tim() { grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-14s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'; }
one() { local t=$1 v=$2 ns=$3; shift 3
  env $(var $v) "$@" OMPT=8 ZE_AFFINITY_MASK=0 python3 $J/run_ref_sng.py $t $C/build_aottt/Prog/ALF.out 1 $ns 1 > /dev/null; tim $t; grep -h "Precision Green" $W/runs/$t/info; }

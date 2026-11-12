# shared helpers (sourced) of the Phase-7 (kernel-level GPU measurement) jobs
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh        # var(), CCS, C (= p5), J (= team/jobs), SEEDS, PARAMS
P7=$SCRATCH/large_prompt/team/p7; PJ7=$P7/jobs; BIN7=$C/build_aottt/Prog/ALF.out
# one chain on tile 0:  r1 TAG VARIANT NSWEEP [env ...]      (wrapped: stdout of the rank in runs/TAG/rank0.log)
r1() { local t=$1 v=$2 ns=$3; shift 3
  env $(var $v) "$@" OMPT=8 ZE_AFFINITY_MASK=0 python3 $PJ7/run_p7.py $t $BIN7 1 $ns 1 $PJ7/capw.sh; }
# node: 32 ranks, 4 per tile, CCS 4, 3 OpenMP threads:  r32 TAG VARIANT NSWEEP [env ...]
r32() { local t=$1 v=$2 ns=$3; shift 3
  env $(var $v) "$@" $CCS OMPT=3 python3 $PJ7/run_p7.py $t $BIN7 32 $ns 1 $PJ7/capw.sh; }
# phase times per sweep from the TEAM_TIMERS line of rank file(s)
tm() { grep -h TEAM_TIMERS $W/runs/$1/rank0.log | awk -v t=$1 '{n=$(NF-7); printf "RESULT %-14s sweep %6.2f  upd %5.2f  stab %5.2f  cgr %5.2f  obs %5.2f s/sweep\n", t, $(NF-6)/n, ($(NF-5)+$(NF-4))/n, ($(NF-3)+$(NF-2))/n, $(NF-1)/n, $NF/n}'; }

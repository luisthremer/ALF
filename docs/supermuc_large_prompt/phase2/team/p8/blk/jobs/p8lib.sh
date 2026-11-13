# shared helpers (sourced) for the p8/blk jobs: p5lib (variants V1..V7, one, node, tim, pg, cpu) + paths of the p8/blk builds
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh       # sets L0, J, C, var(), one(), node(), tim(), pg(), cpu(), CCS; exports OMP_TARGET_OFFLOAD, ALFSH_QRNB=16, SEEDS, PARAMS
P8=$L0/p8/blk; PB=$P8/jobs
B_AOT=$P8/build_aottxt/Prog/ALF.out        # AOT for PVC + phase timers + text observables (all gates)
B_TXT=$P8/build_txtnew/Prog/ALF.out        # JIT text build, same flags as the p5 txtnew (default-vs-previous gate)
B_PREV=$L0/p5/build_txtnew/Prog/ALF.out    # previous binary
SEEDS_CMP=$L0/src/Scripts_and_Parameters_files/Start/seeds
export BIN=$B_AOT
# block timers of the run: blk_begin / blk_end seconds per sweep (P3_TIMERS line of main.F90: blocks, blk_begin, begin_dev, blk_end, end_la_sync, end_dev, slice_level)
p3() { grep -h P3_TIMERS $W/runs/$1/run.log | awk -v t=$1 -v ns=$2 '{printf "RESULT %-14s blocks/sweep %d  blk_begin %.3f  begin_dev %.3f  blk_end %.3f  end_dev %.3f  slice %.3f s/sweep\n", $(NF-7)/ns, $(NF-6)/ns, $(NF-5)/ns, $(NF-4)/ns, $(NF-2)/ns, $(NF-1)/ns}' | sed "s/RESULT /RESULT $1 /"; }

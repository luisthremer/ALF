# ALF builds from the worktree copy p5/src_alf: aottt (AOT for PVC + phase timers; timing/PG gates) and txtnew (text observables: same-seed gates)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p5; cd $P/jobs; export SRC=$P/src_alf
BLD=$P/build_aottt bash build_comb.sh aottt > b_aottt.out 2>&1 &
BLD=$P/build_txtnew bash build_comb.sh txtnew > b_txtnew.out 2>&1 &
wait
ls -l $P/build_aottt/Prog/ALF.out $P/build_txtnew/Prog/ALF.out $P/build_*/Libraries/Modules/lib*.so
tail -3 b_aottt.out b_txtnew.out

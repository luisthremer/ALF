bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "error" -A3 | head
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
./bin/qr2_test correct 16 2>&1 | grep -E "FAIL|failed"
for nbs in 16 24 32; do for sm in 256 384; do echo "NBS=$nbs NBSM=$sm"; ALFDEVQR_NBS=$nbs ALFDEVQR_NBSM=$sm ./bin/qr2_test time 16 2>&1 | grep -E "TIME|rror|xception"; done; done

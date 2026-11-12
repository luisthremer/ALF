bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E " error"
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
echo "#### final kernels, no plant"; ./bin/qr2_test correct 16 2>&1 | grep -E "FAIL|failed"
for p in 1 2 3; do echo "#### DEVQR_PLANT=$p (must FAIL)"; DEVQR_PLANT=$p ./bin/qr2_test correct 16 2>&1 | grep -E "FAIL|failed" | head -6; done
./bin/qr2_test time 16 2>&1 | grep TIME

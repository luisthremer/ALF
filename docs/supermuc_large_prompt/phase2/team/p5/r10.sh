bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "error"
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
./bin/qr2_test apply 16 2>&1 | grep -E "APPLY|rror|xception"

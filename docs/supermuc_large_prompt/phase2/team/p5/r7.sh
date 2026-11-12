bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "error|real|spilled" | sed -e 's/_ZTSZZN8alfdevqr//' -e 's/EvRN4sycl.*compiled/ compiled/' | cut -c1-140
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
./bin/qr2_test correct 16 2>&1 | grep -E "FAIL|failed"
./bin/qr2_test panel 16 2>&1 | grep -E "jb=16|rror|xception"
./bin/qr2_test time 16 2>&1 | grep -E "TIME|rror|xception"

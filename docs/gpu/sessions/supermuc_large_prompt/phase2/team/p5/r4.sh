bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "error|spilled|real" | sed -e 's/_ZTSZZN8alfdevqr//' | cut -c1-120
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
./bin/qr2_test correct 16 2>&1 | grep -E "FAIL|failed"
./bin/qr2_test time 16 2>&1 | grep -E "TIME|rror|xception"
for g in 2 3; do echo "#### G=$g"; ALFDEVQR_G=$g ./bin/qr2_test panel 16 2>&1 | grep -E "jb=16|rror|xception"; done
bash $SCRATCH/large_prompt/team/p5/jobs/t_abl.sh

bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "error|spilled|real" | sed -e 's/_ZTSZZN8alfdevqr//' | cut -c1-110
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
for g in 1 2 3; do echo "#### G=$g"; ALFDEVQR_G=$g ./bin/qr2_test panel 16 2>&1 | grep "jb=16"; done
./bin/qr2_test time 16 | grep TIME

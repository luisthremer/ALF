bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "error|spilled" -A3 | sed -e 's/_ZTSZZN8alfdevqr//' -e 's/EvRN4sycl.*compiled/ compiled/' | cut -c1-140| head
source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
./bin/qr2_test correct 16 2>&1 | grep -E "FAIL|failed"
for anv in 2 4; do echo "ANV=$anv"; ALFDEVQR_ANV=$anv ALFDEVQR_NBS=32 ./bin/qr2_test time 16 2>&1 | grep -E "TIME|rror|xception"; done
ALFDEVQR_ANV=4 ./bin/qr2_test apply 16 2>&1 | grep -E "M= 576|M=1024"

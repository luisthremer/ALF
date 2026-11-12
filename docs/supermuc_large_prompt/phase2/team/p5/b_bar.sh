source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
eval icpx -fsycl -O2 -std=c++17 $AOT $P/jobs/bar_test.cpp -o $P/bin/bar_test 2>&1 | grep -E "error" -A4
./bin/bar_test

# build the standalone test (JIT and AOT for pvc)
source $SCRATCH/large_prompt/team/p5/jobs/env.sh
cd $P; mkdir -p bin
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
( time eval icpx -fsycl -O2 -std=c++17 $AOT -I$P/src -qmkl=parallel $P/jobs/qr2_test.cpp -o $P/bin/qr2_test ) 2>&1 | tail -30
ls -l $P/bin

source $SCRATCH/large_prompt/team/p5/jobs/env.sh; cd $P
ls -l bin
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
eval icpx -fsycl -O2 -std=c++17 $AOT -DABL_A -I$P/src -qmkl=parallel $P/jobs/qr2_test.cpp -o $P/bin/qr2_abl_A 2>&1 | grep -E "error" -A5 | head -30
ls -l bin

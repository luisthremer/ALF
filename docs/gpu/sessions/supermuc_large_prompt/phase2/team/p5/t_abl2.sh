source $SCRATCH/large_prompt/team/p5/jobs/env.sh
cd $P; mkdir -p bin
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
for v in Q S QS; do
  D=""; for ch in $(echo $v | sed 's/./& /g'); do D="$D -DABL_$ch"; done
  eval icpx -fsycl -O2 -std=c++17 $AOT $D -I$P/src -qmkl=parallel $P/jobs/qr2_test.cpp -o $P/bin/qr2_abl2_$v > /dev/null 2>&1 &
done; wait
for v in Q S QS; do echo "#### apply ablation $v"; ./bin/qr2_abl2_$v apply 16 2>&1 | grep -E "M= 576|M= 128 nc=512"; done

# ablation of the panel step (timing only; results wrong): A no scalar chain, B no lane reductions, C no work-group barrier, ABC all
source $SCRATCH/large_prompt/team/p5/jobs/env.sh
cd $P; mkdir -p bin
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
for v in A B C ABC; do
  D=""; for ch in $(echo $v | sed 's/./& /g'); do D="$D -DABL_$ch"; done
  eval icpx -fsycl -O2 -std=c++17 $AOT $D -I$P/src -qmkl=parallel $P/jobs/qr2_test.cpp -o $P/bin/qr2_abl_$v > /dev/null 2>&1 &
done; wait
for v in A B C ABC; do echo "#### ablation $v"; for g in 1 2; do ALFDEVQR_G=$g ./bin/qr2_abl_$v panel 16 2>&1 | grep -E "M=  64 jb=16|M= 256 jb=16|M= 576 jb=16" | sed "s/^/G=$g /"; done; done

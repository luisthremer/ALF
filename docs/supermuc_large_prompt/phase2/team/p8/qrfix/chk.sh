# syntax-only compile of the test with every header copy (login node, light)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
Q=$SCRATCH/large_prompt/team/p8/qrfix; cd $Q
for v in old new; do for d in "" "-DQR2_TESTDELAY"; do echo "== $v $d"; icpx -fsycl -O0 -std=c++17 -fsyntax-only $d -I$Q/hdr_$v -I$Q/src_head/Libraries/Modules -qmkl=parallel qr2_test.cpp 2>&1 | grep -E "error" -A3 | head -30; done; done
echo chk done

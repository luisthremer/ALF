# planted-defect libraries (copies of the final sources with one edit each), AOT for PVC
source $SCRATCH/large_prompt/team/p5/jobs/env.sh
cd $P; mkdir -p plant; cd plant
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
mk() { # name sed-expr file
  d=$1; rm -rf $d; mkdir $d; cp $P/src/alf_la_sycl_devqr.cpp $P/src/alf_la_sycl_devqr.hpp $P/src/alf_la_sycl_qr2.hpp $d/
  sed -i "$2" $d/$3; echo "$d: $(diff $P/src/$3 $d/$3 | grep -c '^>') lines changed"
  ( cd $d && eval icpx -fsycl -O2 -std=c++17 -fPIC -shared $AOT -I$MKLROOT/include alf_la_sycl_devqr.cpp -o libalfdevqr.so -qmkl -lmkl_sycl_blas > build.log 2>&1 ) &
}
# 1 tau(0) of every panel * (1+1e-12)
mk tau 's|tau\[2 \* i\] = t.r; tau\[2 \* i + 1\] = t.i; tb|tau[2 * i] = t.r * (i == 0 ? 1.0 + 1e-12 : 1.0); tau[2 * i + 1] = t.i * (i == 0 ? 1.0 + 1e-12 : 1.0); tb|' alf_la_sycl_qr2.hpp
# 2 wrong conj in the device unmqr (Left: conj flag flipped)
mk conj 's|(int)n, c, c, (const double \*)A|(int)n, c, !c, (const double *)A|' alf_la_sycl_qr2.hpp
# 3 race: ONE tau slot in the panel instead of 2 (owner of step i+1 overwrites tau_i while others still read it)
mk race_tau 's|tb\[2 \* buf|tb[2 * 0|g' alf_la_sycl_qr2.hpp
# 4 race: ONE partial-dot exchange slot in the apply kernel (G > 1) instead of 2 parities
mk race_ex 's|(s \& 1) \* CPW \* G|0 * CPW * G|g' alf_la_sycl_qr2.hpp
# 5 race: ONE v buffer in the panel instead of 2 (owner of step i+1 overwrites v_i while others still read it)
mk race_v 's|const int vo = buf \* MP;|const int vo = 0;|' alf_la_sycl_qr2.hpp
wait
ls -l */libalfdevqr.so; grep -c error */build.log

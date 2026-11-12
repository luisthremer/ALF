source $SCRATCH/large_prompt/team/p5/jobs/env.sh
cd $P/plant; AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
mk() { d=$1; rm -rf $d; mkdir $d; cp $P/src/alf_la_sycl_devqr.cpp $P/src/alf_la_sycl_devqr.hpp $P/src/alf_la_sycl_qr2.hpp $d/
  python3 - $d "$2" "$3" <<'PY'
import sys
d,a,b=sys.argv[1:4]; f=d+"/alf_la_sycl_qr2.hpp"; s=open(f).read(); assert a in s; open(f,"w").write(s.replace(a,b,1))
PY
  ( cd $d && eval icpx -fsycl -O2 -std=c++17 -fPIC -shared $AOT -I$MKLROOT/include alf_la_sycl_devqr.cpp -o libalfdevqr.so -qmkl -lmkl_sycl_blas > build.log 2>&1 ) & }
# the apply kernel's chunk-end barrier: the next buffer is written while slower sub-groups still read the other one
mk nobar "        if (cc + 1 < nch) store(buf ^ 1);
        it.barrier(sycl::access::fence_space::local_space);" "        if (cc + 1 < nch) store(buf ^ 1);"
# early reuse of the staging buffer: the chunk after next is stored into the buffer still being read (single buffer)
mk onebuf "        if (cc + 1 < nch) store(buf ^ 1);" "        if (cc + 1 < nch) store(buf);"
wait; ls -l nobar/libalfdevqr.so onebuf/libalfdevqr.so

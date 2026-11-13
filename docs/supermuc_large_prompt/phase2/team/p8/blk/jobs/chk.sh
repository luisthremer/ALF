source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh >/dev/null 2>&1
cd $SCRATCH/large_prompt/team/p8/blk; mkdir -p chk; cd chk
time icpx -fsycl -fPIC -shared -O2 -std=c++17 -I$MKLROOT/include ../src_alf/Libraries/Modules/alf_la_sycl_blk.cpp -o libchk.so -qmkl -lmkl_sycl_blas 2>&1 | tail -30
ls -l libchk.so

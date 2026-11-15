#!/bin/bash
#SBATCH --job-name=p9_build
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# GPU text build (txtnew recipe = configure IntelLLVM MPI GPU=intel, JIT, text observables) of git archive HEAD, plus planted device-QR libraries
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID HEAD $(cat $SCRATCH/large_prompt/team/p9/l6t/src_head/GIT_HEAD)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p9/l6t; cd $P/jobs; export SRC=$P/src_head
BLD=$P/build_txt bash $SCRATCH/large_prompt/team/p5/jobs/build_comb.sh txtnew > b_txt.out 2>&1 &
# planted libraries: copies of src_head's device-QR sources, one edit each, compiled exactly like the Makefile's libalfdevqr.so (ALF_SYCL_BLK_FLAGS empty = JIT)
# P0 = no edit (control: must reproduce V7 bitwise); Pnn = tau(0) of every panel times (1 + 1e-nn) (the p5 'tau' plant, other sizes)
mkdir -p $P/plant
mk() { d=$P/plant/$1; rm -rf $d; mkdir -p $d; cp $SRC/Libraries/Modules/alf_la_sycl_devqr.cpp $SRC/Libraries/Modules/alf_la_sycl_devqr.hpp $SRC/Libraries/Modules/alf_la_sycl_qr2.hpp $d/
  if [ -n "$2" ]; then python3 - $d/alf_la_sycl_qr2.hpp "$2" <<'PY'
import sys
f,e=sys.argv[1:3]; a="tau[2 * i] = t.r; tau[2 * i + 1] = t.i; tb"
b="tau[2 * i] = t.r * (i == 0 ? 1.0 + %s : 1.0); tau[2 * i + 1] = t.i * (i == 0 ? 1.0 + %s : 1.0); tb"%(e,e)
s=open(f).read(); assert s.count(a)==1; open(f,"w").write(s.replace(a,b,1))
PY
  fi
  echo "plant $1: $(diff $SRC/Libraries/Modules/alf_la_sycl_qr2.hpp $d/alf_la_sycl_qr2.hpp | grep -c '^>') line(s) changed"
  ( cd $d && icpx -fsycl -O2 -std=c++17 -fPIC -shared -I$MKLROOT/include alf_la_sycl_devqr.cpp -o libalfdevqr.so -qmkl -lmkl_sycl_blas > build.log 2>&1; echo "plant $1 rc=$?" ) &
}
mk P0 ""; mk P12 1e-12; mk P11 1e-11; mk P10 1e-10
wait
diff $P/plant/P10/alf_la_sycl_qr2.hpp $SRC/Libraries/Modules/alf_la_sycl_qr2.hpp
ls -l $P/build_txt/Prog/ALF.out $P/build_txt/Libraries/Modules/lib*.so $P/plant/*/libalfdevqr.so
tail -3 b_txt.out
md5sum $P/build_txt/Libraries/Modules/libalfdevqr.so $P/plant/*/libalfdevqr.so
echo "end $(date -Is)"

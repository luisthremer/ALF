#!/bin/bash
# Phase 0 step 4: CPU hotspot profiles with gprof (perf and VTune are blocked on nividia:
# perf_event_paranoid=4, yama ptrace_scope=1).
# usage: profile.sh build [cg]                 -> bin/ALF-gnu-PROF.out  (or bin/ALF-gnu-PROFCG.out with "cg")
#        profile.sh run BIN CONFIG NSWEEP NBIN [pmp] -> runs/prof/<BIN>_<CONFIG>: 1 MPI rank, 1 thread,
#                                                   gprof.txt, or stacks.txt with "pmp" (pmp.py)
# PROF:   reference flags (-O3 ...) + -g, objects NOT instrumented, linked with -pg and *static* OpenBLAS
#         (BLAS+LAPACK) so the 100 Hz sampling histogram also covers BLAS/LAPACK code. Unbiased self times.
# PROFCG: objects compiled with -pg as well: exact call counts and call graph, but mcount overhead
#         inflates small, frequently called routines. Use it for counts, PROF for percentages.
# CONFIG names from testsuite/gpu_reference/configs.py.
set -euo pipefail
. ~/alf-gpu/env/gnu.sh
src=~/alf-gpu/src/ALF-gnu
case $1 in
build)
  tag=PROF; ext="-g"; [ "${2:-}" = cg ] && { tag=PROFCG; ext="-g -pg"; }
  cd $src; make clean >/dev/null 2>&1
  export ALF_FLAGS_EXT="$ext"
  set +eu; . ./configure.sh GNU MPI HDF5 NO-INTERACTIVE; set -eu
  export ALF_LIB=${ALF_LIB/-llapack -lblas/\/usr\/lib\/x86_64-linux-gnu\/openblas-pthread\/libopenblas.a -lpthread}
  echo "ALF_FLAGS_PROG=$ALF_FLAGS_PROG"; echo "ALF_LIB=$ALF_LIB"
  make -j6 all > ~/alf-gpu/logs/build_gnu_$tag.log 2>&1
  rm -f Prog/ALF.out; make -C Prog ALF.out ALF_FC="$ALF_FC -pg" >> ~/alf-gpu/logs/build_gnu_$tag.log 2>&1  # link-only -pg
  cp Prog/ALF.out ~/alf-gpu/bin/ALF-gnu-$tag.out && echo BUILD_OK $tag ;;
run)
  bin=$2; cfg=$3; out=~/alf-gpu/runs/prof; d=$out/${bin}_$cfg
  tmp=$out/.prep_$bin.$$; rm -rf $d $tmp; cd ~/alf-gpu/gpu_reference
  python prep.py $tmp $cfg --alf-dir $src --nsweep $4 --nbin $5 >/dev/null
  mv $tmp/$cfg $d; rm -rf $tmp
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  cd $d; start=$(date +%s)
  if [ "${6:-}" = pmp ]; then  # stack sampling under gdb (see pmp.py) instead of gprof
    PMP_OUT=$d/stacks.txt mpirun -n 1 /opt/intel/oneapi/debugger/latest/opt/debugger/bin/gdb-oneapi \
      -batch -x ~/alf-gpu/pmp.py --args ~/alf-gpu/bin/ALF-gnu-$bin.out > run.log 2>&1
  else
    mpirun -n 1 ~/alf-gpu/bin/ALF-gnu-$bin.out > run.log 2>&1
    gprof -b ~/alf-gpu/bin/ALF-gnu-$bin.out gmon.out > gprof.txt
  fi
  echo "wall_s=$(( $(date +%s) - start ))" > timing.txt
  echo "${bin}_$cfg $(cat timing.txt)" ;;
esac

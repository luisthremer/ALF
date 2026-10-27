#!/bin/bash
# Ticket 202, run on nividia: t_order negative control. Builds module + test helpers + t_order standalone twice from
# ~/alf-gpu/work/t202/src (see t202_sync.sh): normal (must PASS) and with -DALF_LA_TEST_BADSTREAM (must FAIL on RAW, WAR, WAW).
# Run alone, not under ctest -j. usage: bash t202_negctl.sh
W=$HOME/alf-gpu/work/t202; export TMPDIR=$W/tmp; mkdir -p $TMPDIR; . ~/alf-gpu/env/nvhpc.sh
export OMP_TARGET_OFFLOAD=MANDATORY
nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader
for V in normal badstream; do
  D=$W/exp_$V; rm -rf $D; mkdir -p $D; cd $D || exit 1
  S=$W/src; cp $S/Libraries/Modules/{alf_la_mod.F90,runtime_error_mod.F90} $S/testsuite/alf_la.tests/{la_test_util.F90,t_order.F90} .
  X=""; [ $V = badstream ] && X="-DALF_LA_TEST_BADSTREAM"
  F="-Mpreprocess -O3 -Mfprelaxed -fast -mp=gpu -gpu=cc80 -DALF_GPU -DALF_LA_CUDA $X"
  pgfortran $F -c runtime_error_mod.F90 alf_la_mod.F90 2>&1 | grep -v W-0155
  pgfortran $F -Kieee -c la_test_util.F90 2>&1
  pgfortran $F -c t_order.F90 2>&1
  pgfortran $F -cudalib=cublas,cusolver *.o -o t_order -llapack -lblas 2>&1
  echo "=== t_order $V (alone)"; ./t_order 2>&1 | grep -v ieee; echo "exit code ${PIPESTATUS[0]}"
  cd $W; rm -rf $D
done

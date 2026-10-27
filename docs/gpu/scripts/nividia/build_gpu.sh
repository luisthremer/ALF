#!/bin/bash
# Ticket 201: build a GPU target from a tree.  usage: build_gpu.sh TREE nvidia|intel [HDF5] [ctest]
# nvidia -> PGI noMPI GPU=nvidia ; intel -> IntelLLVM noMPI GPU=intel.  Logs in TREE/../logs. Cores: make -j4.
TREE=$1; T=$2; shift 2
case $T in nvidia) . ~/alf-gpu/env/nvhpc.sh; M=PGI;; intel) . ~/alf-gpu/env/intel.sh; M=IntelLLVM;; *) echo bad target; exit 1;; esac
H=""; CT=0; for a in "$@"; do [ $a = HDF5 ] && H=HDF5; [ $a = ctest ] && CT=1; done
L=$TREE/../logs; mkdir -p $L; cd $TREE || exit 1
make clean >/dev/null 2>&1
. ./configure.sh $M noMPI $H GPU=$T NO-INTERACTIVE > $L/conf_$T$H.log 2>&1 || { echo "configure failed"; exit 2; }
env | grep '^ALF_' | sort > $L/env_$T$H.txt
/usr/bin/time -f "make wall %e s" make -j4 all > $L/make_$T$H.log 2>&1 || { echo "make failed"; tail -20 $L/make_$T$H.log; exit 3; }
tail -1 $L/make_$T$H.log
ldd Prog/ALF.out | grep -Ei 'cublas|cusolver|omp|mkl|sycl' 
if [ $CT = 1 ]; then
  ( cd testsuite && rm -rf tests && mkdir tests && cd tests && cmake .. >/dev/null && cmake --build . -j4 > $L/cmake_$T$H.log 2>&1 \
    && ctest -j4 --output-on-failure > $L/ctest_$T$H.log 2>&1; tail -4 $L/ctest_$T$H.log )
fi

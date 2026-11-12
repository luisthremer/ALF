#!/bin/bash
#SBATCH --job-name=tm_build
#SBATCH --output=lp_build.%j.log
#SBATCH --error=lp_build.%j.err
#SBATCH --time=00:55:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
# usage: sbatch build_intel.sh <variant: mpi|nompi|cpu-mpi|cpu-nompi|gpu-txt|cpu-txt|gpu-ser|cpu-ser> [ctest]
V=${1:-mpi}
SRC=$SCRATCH/large_prompt/team/src
BLD=$SCRATCH/large_prompt/team/build_$V
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
module -t list 2>&1; echo MKLROOT=$MKLROOT
rm -rf $BLD; mkdir -p $BLD; cp -r $SRC/. $BLD/; cd $BLD
case $V in
  mpi)       ARGS="IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE" ;;
  nompi)     ARGS="IntelLLVM noMPI GPU=intel HDF5 NO-INTERACTIVE" ;;
  cpu-mpi)   ARGS="IntelLLVM MPI HDF5 NO-INTERACTIVE" ;;
  cpu-nompi) ARGS="IntelLLVM noMPI HDF5 NO-INTERACTIVE" ;;
  txtnew|gpu-txt)   ARGS="IntelLLVM MPI GPU=intel NO-INTERACTIVE" ;;   # text observables (no HDF5): same-seed comparisons
  cpu-txt)   ARGS="IntelLLVM MPI NO-INTERACTIVE" ;;
  aotttnew|aotnew|aottt|aottim|aot)       ARGS="IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE" ;;   # like mpi, but ahead-of-time compiled for PVC (spir64_gen)
  pp)        ARGS="IntelLLVM PP GPU=intel HDF5 NO-INTERACTIVE" ;;      # PARALLEL_PARAMS production binary (merged -mpi-pp)
  gpu-pptxt) ARGS="IntelLLVM PP GPU=intel NO-INTERACTIVE" ;;
  cpu-pptxt) ARGS="IntelLLVM PP NO-INTERACTIVE" ;;
  pp-aot)    ARGS="IntelLLVM PP GPU=intel HDF5 NO-INTERACTIVE" ;;
  gpu-ser)   ARGS="IntelLLVM noMPI GPU=intel NO-INTERACTIVE" ;;  # serial text build: repo ctest suite
  cpu-ser)   ARGS="IntelLLVM noMPI NO-INTERACTIVE" ;;
esac
# no internet: point configure.sh at the LRZ HDF5 (directory name = compiler version string, as configure computes it)
ifx -o $BLD/gcv $BLD/get_compiler_version.F90; CV=$($BLD/gcv | sed 's/[ ,()]/_/g'); rm -f $BLD/gcv
H5=$(dirname $(dirname $(which h5fc)))
# LRZ HDF5 uses libhdf5_hl_fortran.a (configure.sh expects libhdf5hl_fortran.a) and prints CMake-style 'm;dl': build a thin shim dir
HD=$BLD/h5ext/$CV; mkdir -p $HD/lib $HD/bin
ln -s $H5/include $HD/include
for f in $H5/lib/*.a; do ln -s $f $HD/lib/$(basename $f); done
ln -s $H5/lib/libhdf5_hl_fortran.a $HD/lib/libhdf5hl_fortran.a
printf '#!/bin/bash\n%s/bin/h5fc "$@" | sed -e "s/Extra libraries:.*/Extra libraries: -lm -ldl -lz/"\n' $H5 > $HD/bin/h5fc; chmod +x $HD/bin/h5fc
export ALF_HDF5_DIR=$BLD/h5ext
echo "configure $ARGS"; date -Is
source ./configure.sh $ARGS > configure.log 2>&1; tail -12 configure.log   # never pipe a sourced script
# LRZ HDF5 has separate Fortran C-stub libs (libhdf5_f90cstub.a, libhdf5_hl_f90cstub.a); they must follow libhdf5_fortran.a
case "$ARGS" in *HDF5*) export ALF_LIB="$ALF_LIB $H5/lib/libhdf5_hl_f90cstub.a $H5/lib/libhdf5_hl.a $H5/lib/libhdf5_f90cstub.a $H5/lib/libhdf5.a -lz -ldl -lm";; esac
case "$V" in aot|aotttnew|aotnew|aottim|aottt|pp-aot) AOTB=1;; *) AOTB=0;; esac
if [ "$AOTB" = 1 ]; then   # AOT for PVC instead of JIT: spir64 -> spir64_gen + ocloc device pvc (ocloc: /usr/bin/ocloc)
  for v in ALF_FLAGS_MODULES ALF_FLAGS_ANA ALF_FLAGS_PROG ALF_LIB; do
    eval "export $v=\"\${$v//-fopenmp-targets=spir64/-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc}\""   # no quotes: the Makefiles re-split the flags
  done
  echo "AOT flags: $ALF_FLAGS_PROG" | cut -c1-600
fi
case "$V" in aottt|aotttnew) export ALF_FLAGS_PROG="$ALF_FLAGS_PROG -DALF_TEAM_TIMERS";; esac   # opt-in phase timers (main.F90)
env | grep ^ALF_ | cut -c1-400
export OMP_NUM_THREADS=1
( time make -j48 all ) 2>&1 | tail -60
date -Is; ls -l Prog/ALF*.out 2>&1
if [ "$2" = ctest ]; then   # repo unit/regression suite (testsuite/), serial builds only
  module load cmake/3.30.0
  export OMP_NUM_THREADS=1
  [ "${V#gpu}" != "$V" ] && export OMP_TARGET_OFFLOAD=MANDATORY
  cd testsuite && rm -rf tests && mkdir tests && cd tests && cmake .. > cmake_conf.log 2>&1 && cmake --build . -j16 > cmake_build.log 2>&1
  ( time ctest -j2 --output-on-failure ) > ctest.log 2>&1; tail -25 ctest.log
fi

#!/bin/bash -l
# builds txt (JIT, gates) and aottt (AOT PVC + HDF5 + TEAM_TIMERS) of $M/src, plus planted libs. Run on the login node (compile only).
M=$SCRATCH/mgr_corr5; L0=$SCRATCH/large_prompt/team
source $L0/jobs/env_intel_pvc.sh
echo "start $(date -Is) $(hostname)"
bld() { V=$1; ARGS=$2; BLD=$M/build_$V
  rm -rf $BLD; mkdir -p $BLD; cp -r $M/src/. $BLD/; cd $BLD
  ifx -o $BLD/gcv $BLD/get_compiler_version.F90; CV=$($BLD/gcv | sed 's/[ ,()]/_/g'); rm -f $BLD/gcv
  H5=$(dirname $(dirname $(which h5fc)))
  HD=$BLD/h5ext/$CV; mkdir -p $HD/lib $HD/bin
  ln -s $H5/include $HD/include
  for f in $H5/lib/*.a; do ln -s $f $HD/lib/$(basename $f); done
  ln -s $H5/lib/libhdf5_hl_fortran.a $HD/lib/libhdf5hl_fortran.a
  printf '#!/bin/bash\n%s/bin/h5fc "$@" | sed -e "s/Extra libraries:.*/Extra libraries: -lm -ldl -lz/"\n' $H5 > $HD/bin/h5fc; chmod +x $HD/bin/h5fc
  export ALF_HDF5_DIR=$BLD/h5ext
  echo "[$V] configure $ARGS"
  source ./configure.sh $ARGS > configure.log 2>&1; tail -3 configure.log
  case "$ARGS" in *HDF5*) export ALF_LIB="$ALF_LIB $H5/lib/libhdf5_hl_f90cstub.a $H5/lib/libhdf5_hl.a $H5/lib/libhdf5_f90cstub.a $H5/lib/libhdf5.a -lz -ldl -lm";; esac
  if [ $V = aottt ]; then
    export ALF_SYCL_BLK_FLAGS='-fsycl-targets=spir64_gen -Xs "-device pvc"'
    for v in ALF_FLAGS_MODULES ALF_FLAGS_ANA ALF_FLAGS_PROG ALF_LIB; do
      eval "export $v=\"\${$v//-fopenmp-targets=spir64/-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc}\""
    done
    export ALF_FLAGS_PROG="$ALF_FLAGS_PROG -DALF_TEAM_TIMERS"
  fi
  env | grep ^ALF_ | cut -c1-300 > env.log
  export OMP_NUM_THREADS=1
  ( time make -j16 all ) > make.log 2>&1
  tail -5 make.log; echo "[$V] done $(date -Is)"; ls -l Prog/ALF*.out Libraries/Modules/lib*.so 2>&1
}
( bld txt "IntelLLVM MPI GPU=intel NO-INTERACTIVE" ) > $M/jobs/bld_txt.log 2>&1 &
( bld aottt "IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE" ) > $M/jobs/bld_aottt.log 2>&1 &
wait
NEW=$M/build_txt/Prog/ALF.out
[ -x $NEW ] && [ -x $M/build_aottt/Prog/ALF.out ] || { echo "BUILD FAILED"; grep -i -B2 -A5 "error" $M/build_txt/make.log $M/build_aottt/make.log | head -60; exit 1; }
# planted-defect libraries (same flags as the Makefile, txt: no AOT flags)
mkdir -p $M/wrap
for n in 1 2 3 4 5 6 7 8 9; do mkdir -p $M/plant$n; done; cd $M/build_txt/Libraries/Modules
for n in 1 2 3 4 5 6 7 8 9; do
  icpx -fsycl -fPIC -shared -O2 -std=c++17 -DALFSB_PLANT_FUSE=$n -I$MKLROOT/include alf_la_sycl_blk.cpp -o $M/plant$n/libalfsb.so -qmkl -lmkl_sycl_blas 2>&1 | tail -3
done
ls -l $M/plant*/
mk() { printf '#!/bin/bash\n%s\nexec %s "$@"\n' "$2" "$3" > $M/wrap/$1; chmod +x $M/wrap/$1; }
mk off0.sh "unset ALFSH_BLKFUSE ALFSH_BLKPF ALFSH_OBSCORR" $NEW
mk off.sh "unset ALFSH_BLKFUSE ALFSH_BLKPF; export ALFSH_OBSCORR=1" $NEW
for n in 1 2 3 4; do mk on$n.sh "unset ALFSH_BLKPF; export ALFSH_BLKFUSE=$n ALFSH_OBSCORR=1" $NEW; done
mk on5.sh "export ALFSH_BLKFUSE=2 ALFSH_BLKPF=1 ALFSH_OBSCORR=1" $NEW
mk on6.sh "export ALFSH_BLKFUSE=4 ALFSH_BLKPF=1 ALFSH_OBSCORR=1" $NEW
mk on7.sh "export ALFSH_BLKFUSE=5 ALFSH_BLKPF=0 ALFSH_OBSCORR=1" $NEW
mk on8.sh "export ALFSH_BLKFUSE=5 ALFSH_BLKPF=1 ALFSH_OBSCORR=1" $NEW
for n in 1 2 3 4 5 6 7 8 9; do
  case $n in 1) lvl=2;pf=0;; 2) lvl=4;pf=0;; 3) lvl=1;pf=0;; 4) lvl=2;pf=0;; 5) lvl=4;pf=0;; 6) lvl=2;pf=1;; 7) lvl=2;pf=1;; 8) lvl=5;pf=1;; 9) lvl=5;pf=1;; esac
  mk plant$n.sh "export ALFSH_BLKFUSE=$lvl ALFSH_BLKPF=$pf ALFSH_OBSCORR=1 LD_LIBRARY_PATH=$M/plant$n:\$LD_LIBRARY_PATH" $NEW
done
echo "end $(date -Is)"

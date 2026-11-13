#!/bin/bash
#SBATCH --job-name=mc_gate
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --exclude=i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14
# intel-batched: observables corr as SYCL kernel (ALFSH_OBSCORR=1). Builds txt (JIT, gates) and aottt (timing) of $M/src concurrently,
# a planted-defect libalfsb, then the gates: off vs unmodified (bitwise), on vs off, on/off vs same-compiler CPU, planted.
M=$SCRATCH/mgr_corr; U=$SCRATCH/large_prompt/unified; L0=$SCRATCH/large_prompt/team; J=$U/jobs
REF=$SCRATCH/large_prompt/batched/build_ib_txtnew/Prog/ALF.out   # unmodified source (Prog identical to HEAD), same recipe
source $L0/jobs/env_intel_pvc.sh
echo "start $(date -Is) $(hostname)"; module -t list 2>&1 | tr '\n' ' '; echo
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
  ( time make -j32 all ) > make.log 2>&1
  tail -5 make.log; echo "[$V] done $(date -Is)"; ls -l Prog/ALF*.out Libraries/Modules/lib*.so 2>&1
}
( bld txt "IntelLLVM MPI GPU=intel NO-INTERACTIVE" ) > $M/jobs/bld_txt.log 2>&1 &
( bld aottt "IntelLLVM MPI GPU=intel HDF5 NO-INTERACTIVE" ) > $M/jobs/bld_aottt.log 2>&1 &
wait
cat $M/jobs/bld_txt.log $M/jobs/bld_aottt.log
NEW=$M/build_txt/Prog/ALF.out
[ -x $NEW ] && [ -x $M/build_aottt/Prog/ALF.out ] || { echo "BUILD FAILED"; grep -i "error" $M/build_txt/make.log $M/build_aottt/make.log | head -20; exit 1; }
nm -D $M/build_txt/Libraries/Modules/libalfsb.so | grep -c alfsb_corr | sed 's/^/alfsb_corr symbols in libalfsb: /'
# planted-defect library (two wrong indices, -DALFSB_PLANT_CORR), same flags as the Makefile (txt: no AOT flags)
mkdir -p $M/plant; cd $M/build_txt/Libraries/Modules
icpx -fsycl -fPIC -shared -O2 -std=c++17 -DALFSB_PLANT_CORR -I$MKLROOT/include alf_la_sycl_blk.cpp -o $M/plant/libalfsb.so -qmkl -lmkl_sycl_blas 2>&1 | tail -3
ls -l $M/plant
# wrappers (the switch differs per run, the binary is the same)
mkdir -p $M/wrap
mk() { printf '#!/bin/bash\n%s\nexec %s "$@"\n' "$2" "$3" > $M/wrap/$1; chmod +x $M/wrap/$1; }
mk off.sh "unset ALFSH_OBSCORR" $NEW
mk on.sh "export ALFSH_OBSCORR=1" $NEW
mk plant.sh "export ALFSH_OBSCORR=1 LD_LIBRARY_PATH=$M/plant:\$LD_LIBRARY_PATH" $NEW
export W0=$M/work; rm -rf $W0; mkdir -p $W0; export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
SB=LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
V7="$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_QRNB=16"
declare -A CASE=( [L6F]="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12F]="L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
                  [L6T]="L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12T]="L=12 mz=T scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" )
declare -A TOL=( [L6F]=1e-13 [L12F]=1e-13 [L6T]=1e-10 [L12T]=1e-10 )
CPU=$L0/build_cpu-txt/Prog/ALF.out
TASKS=()
add() { TASKS+=("$*"); }   # name driver cpubin gpubin tol case
for c in L6F L12F L6T L12T; do
  add a_$c gg $REF $NEW 0 $c          # switch off vs unmodified build, bitwise
  add b_$c gg $M/wrap/off.sh $M/wrap/on.sh 1e-12 $c   # same binary, switch on vs off
  add con_$c sng $CPU $M/wrap/on.sh ${TOL[$c]} $c     # on vs same-compiler CPU
  add coff_$c sng $CPU $M/wrap/off.sh ${TOL[$c]} $c   # off vs same-compiler CPU
done
for c in L6F L6T; do
  add p_$c gg $M/wrap/off.sh $M/wrap/plant.sh 1e-12 $c
  add pc_$c sng $CPU $M/wrap/plant.sh ${TOL[$c]} $c
done
stream() { local s=$1 i=0 t
  for t in "${TASKS[@]}"; do
    if [ $((i % 8)) -eq $s ]; then set -- $t; local name=$1 drv=$2 cb=$3 gb=$4 tol=$5 cs=$6
      env $V7 CMP_TOL=$tol W=$W0/$name CPUBIN=$cb GPUBIN=$gb ZE_AFFINITY_MASK=$s python3 $J/run_cmp_$drv.py $name ${CASE[$cs]} > $W0/out_$name.txt 2>&1
    fi
    i=$((i+1)); set -- $s
  done; }
for s in 0 1 2 3 4 5 6 7; do stream $s & done
wait
for t in "${TASKS[@]}"; do set -- $t; f=$W0/out_$1.txt; echo "#### $1 ($2: $3 vs $4, tol $5)"; grep -E "rc=|PASS|FAIL" $f | cut -c1-160; done
echo "---- on-runs that announced the switch (gpu run.log):"; grep -l "ALFSH_OBSCORR=1" $W0/*/runs/*/gpu/run.log | wc -l
echo "---- per-file deviation, b_L12T and p_L6T:"; for n in b_L12T p_L6T p_L6F; do echo "## $n"; grep -vE "rc=" $W0/out_$n.txt | head -30; done
echo "end $(date -Is)"

#!/bin/bash
# GPU=none / GPU=host identity checks for cm/dev/gpu-nvidia (P2.4 step 2), melmac

R=~/work/ALF_P2/wt-gpu-nvidia; C=~/work/ALF_P2/nv_check; cd $C
BASE=${BASE:-889f4633}
COMMITS=$(cd $R && git rev-list --reverse $BASE..HEAD)
FILES="Libraries/Modules/linalg_mod.F90 Prog/gpu_state_mod.F90 Prog/gpu_wrap_mod.F90 Prog/gpu_update_mod.F90 Prog/gpu_lowrank_mod.F90 Prog/gpu_stab_mod.F90 Prog/gpu_obs_mod.F90 Prog/gpu_tau_mod.F90"
pp() { # $1 tree $2 defines -> concatenated preprocessed output
  for f in $FILES; do gfortran -E -cpp -P $2 $1/$f 2>/dev/null | grep -v "^[[:space:]]*$"; done | md5sum | cut -d' ' -f1; }
rm -rf trees; mkdir trees
for c in $BASE $COMMITS; do mkdir trees/$c; (cd $R && git archive $c) | tar x -C trees/$c; done
echo "== preprocessed identity (GPU=none, GPU=host; default/STABLOG; noMPI/MPI)"
for d in "" "-DSTABLOG" "-DMPI" "-DALF_GPU" "-DALF_GPU -DSTABLOG" "-DALF_GPU -DMPI"; do
  b=$(pp trees/$BASE "$d")
  for c in $COMMITS; do x=$(pp trees/$c "$d"); [ "$x" = "$b" ] && r=same || r=DIFF; echo "  [$d] ${c:0:8} $r"; done
done
echo "== configure output GPU=none/host unchanged (flags)"
for c in $BASE $COMMITS; do
  (cd trees/$c && for g in none host; do (. ./configure.sh GNU noMPI GPU=$g NO-INTERACTIVE >/dev/null 2>&1; echo "$g|$ALF_FLAGS_PROG|$ALF_FLAGS_MODULES|$ALF_LIB|${ALF_GPU:-}") ; done | sed "s#$C/trees/$c#T#g" | md5sum | cut -c1-8 | sed "s/^/  ${c:0:8} /")
done
echo "== every commit builds GPU=none and GPU=host"
for c in $COMMITS; do for g in none host; do
  rm -rf bld; cp -r trees/$c bld
  (cd bld && . ./configure.sh GNU noMPI GPU=$g NO-INTERACTIVE >/dev/null 2>&1 && make -j16 all > build.log 2>&1) && r=ok || r=FAIL
  echo "  ${c:0:8} GPU=$g $r"
done; done
TIP=$(echo $COMMITS | awk '{print $NF}')
echo "== tip ${TIP:0:8}: 16 scenarios byte-identical vs base (GPU=none)"
for t in $BASE $TIP; do rm -rf ref_$t; cp -r trees/$t ref_$t; (cd ref_$t && . ./configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && make -j16 all > build.log 2>&1); done
python3 -I bitwise_vs_ref.py $C/trees/$TIP $C/ref_$BASE/Prog/ALF.out $C/ref_$TIP/Prog/ALF.out $C/work16 2>&1 | tail -20
echo "== tip ctest GPU=none and GPU=host"
for g in none host; do
  rm -rf ct_$g; cp -r trees/$TIP ct_$g
  (cd ct_$g && . ./configure.sh GNU noMPI GPU=$g NO-INTERACTIVE >/dev/null 2>&1 && make -j16 all > build.log 2>&1 && cmake -S testsuite -B testsuite/tests >/dev/null && cmake --build testsuite/tests -j16 >/dev/null 2>&1 && ctest --test-dir testsuite/tests -j16 --timeout 6000 2>&1 | tail -4) | sed "s/^/  GPU=$g /"
done
echo DONE

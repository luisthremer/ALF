#!/bin/bash
#SBATCH --job-name=qf_q1
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:45:00
# Tasks 1 + 2 (standalone): build the test with the header copies (AOT for PVC), run on the 8 tiles in parallel.
#  race_apply : unmqr2 with tau=0 patterns, old/new header, delay no/yes;  Gate A vs the NEW header (clean + plants 1/2/3)
#  race_panel : geqrf2, plants race_tau / race_v / race_ex (apply) with the delay, real header with the same delay
echo "start $(date -Is) node $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
export SYCL_UR_USE_LEVEL_ZERO_V2=1 MKL_NUM_THREADS=6 OMP_NUM_THREADS=6
unset OMP_TARGET_OFFLOAD
Q=$SCRATCH/large_prompt/team/p8/qrfix; cd $Q; mkdir -p bin out
echo "GIT_HEAD of the source tree: $(cat src_head/GIT_HEAD)"
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
bld() { # name hdr [extra flags]
  ( eval icpx -fsycl -O2 -std=c++17 $AOT $3 -I$Q/hdr_$2 -I$Q/src_head/Libraries/Modules -qmkl=parallel qr2_test.cpp -o bin/qr2_$1 > out/build_$1.log 2>&1; echo "build $1 rc=$?" ) &
}
t0=$(date +%s)
bld old_plain old; bld new_plain new
bld old_d old -DQR2_TESTDELAY; bld new_d new -DQR2_TESTDELAY; bld tau_d tau -DQR2_TESTDELAY; bld v_d v -DQR2_TESTDELAY; bld ex_d ex -DQR2_TESTDELAY
wait; echo "builds done in $(( $(date +%s) - t0 )) s"; ls -l bin | tail -8; grep -l -E "error" out/build_*.log
md5sum hdr_*/alf_la_sycl_qr2.hpp src_head/Libraries/Modules/alf_la_sycl_qr2.hpp

slot=0
run() { # name mode [env...]   (runs on tile slot%8, 8 in parallel)
  local name=$1 mode=$2 bin=$3; shift 3
  local tile=$((slot % 8)); slot=$((slot + 1))
  ( s=$(date +%s); env ZE_AFFINITY_MASK=$tile "$@" ./bin/qr2_$bin $mode > out/$name.log 2>&1; echo "rc=$? wall=$(( $(date +%s) - s )) s" >> out/$name.log ) &
  [ $((slot % 8)) -eq 0 ] && wait
}
# ---- race_apply (Left/Right x N/C, G=2 and G=4): header old/new, delay 0 (reps 20) and delays (reps 3)
run A_old_plain   race_apply old_plain RACE_REPS=20
run A_new_plain   race_apply new_plain RACE_REPS=20
run A_old_d0      race_apply old_d     RACE_REPS=20 QR2_DELAY=0
run A_new_d0      race_apply new_d     RACE_REPS=20 QR2_DELAY=0
for d in 300 3000 30000; do
  run A_old_d$d   race_apply old_d     RACE_REPS=3 QR2_DELAY=$d
  run A_new_d$d   race_apply new_d     RACE_REPS=3 QR2_DELAY=$d
done
run A_ex_d3000    race_apply ex_d      RACE_REPS=3 QR2_DELAY=3000   # plant race_ex (single exchange slot), sensitivity control
run A_ex_d0       race_apply ex_d      RACE_REPS=20 QR2_DELAY=0
# ---- Gate A vs the NEW header: clean + plants 1/2/3
run GA_clean      correct new_plain
for p in 1 2 3; do run GA_plant$p correct new_plain DEVQR_PLANT=$p; done
wait
# ---- race_panel: G=1 (default for M<=256, forced for the larger M)
for b in new_d tau_d v_d; do
  run P_${b}_d0     race_panel $b ALFDEVQR_G=1 RACE_REPS=10 QR2_PDELAY=0
  for d in 500 5000; do run P_${b}_d$d race_panel $b ALFDEVQR_G=1 RACE_REPS=3 QR2_PDELAY=$d; done
done
for b in new_d tau_d v_d; do run P_${b}_Gdef_d5000 race_panel $b RACE_REPS=3 QR2_PDELAY=5000; done   # default G (G=2 for M>256): barrier in the step, plants not racy there
run P_old_d5000   race_panel old_d ALFDEVQR_G=1 RACE_REPS=3 QR2_PDELAY=5000
wait
echo "end $(date -Is)"
for f in out/A_*.log out/P_*.log out/GA_*.log; do echo "--- $(basename $f .log): $(grep -E '^race_apply:|^race_panel:|^correctness:' $f | tr '\n' ' ') $(tail -1 $f)"; done

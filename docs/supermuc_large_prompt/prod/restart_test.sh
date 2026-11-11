#!/bin/bash
#
# Restart / preemption test for the PP production workflow (partition test, 1 node, 32 ranks = 2 U values x 16 walkers, L=8).
#   (a) reference: one uninterrupted run of NBIN bins
#   (b) same seeds: run, kill -9 the whole mpiexec process tree after about half of the bins, leave RUNNING, recover.sh, rerun the rest
#   (c) like (b), but one confout_*.h5 of Temp_0 is truncated to 100 bytes before recovery -> recover.sh must restore Temp_0 from .segbak
# Compared with cmp_h5dump.py (data.h5 of both sets, and every final confout_*.h5 where no duplicated bin is involved).
#
# Submit from the directory that holds recover.sh and cmp_h5dump.py (this prod/ directory):
#   sbatch --export=GPU_PP_BIN=/path/to/alf_gpu_pp_hdf5.exe restart_test.sh
# Optional env (add to --export=): NBIN (8), NSWEEP (2), WORKDIR ($SCRATCH/prod_2026-11/work_restart),
#   SRCDIR ($SCRATCH/prod_2026-11/src/Scripts_and_Parameters_files/Start_PP_GPU), PROD_DIR (dir of recover.sh / cmp_h5dump.py),
#   ENV_FILE ($SCRATCH/prod_2026-11/jobs/env_intel_pvc.sh), FAIL_TIMEOUT_S (900, wait for the kill point).
# (No --export=NONE line on purpose: it would drop GPU_PP_BIN; the environment file loads the modules itself.)
#
#SBATCH --job-name restart_test
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --get-user-env
#SBATCH --output=restart_test.%j.log
#SBATCH --error=restart_test.%j.log

set -u
SCRATCH=${SCRATCH:?SCRATCH not set}
ENV_FILE=${ENV_FILE:-$SCRATCH/prod_2026-11/jobs/env_intel_pvc.sh}
set +u; source "$ENV_FILE"; set -u
export ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 OMP_NUM_THREADS=3

GPU_PP_BIN=${GPU_PP_BIN:?set GPU_PP_BIN to the HDF5 PP GPU executable (sbatch --export=GPU_PP_BIN=... restart_test.sh)}
[ -x "$GPU_PP_BIN" ] || { echo "FATAL: $GPU_PP_BIN not executable"; exit 2; }
W=${WORKDIR:-$SCRATCH/prod_2026-11/work_restart}
SRC=${SRCDIR:-$SCRATCH/prod_2026-11/src/Scripts_and_Parameters_files/Start_PP_GPU}
PROD=${PROD_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}}
NBIN=${NBIN:-8}; NSWEEP=${NSWEEP:-2}; L=8; ULIST=5.0,5.6; NPER=16; NTASKS=${SLURM_NTASKS:-32}
TIMEOUT=${FAIL_TIMEOUT_S:-900}
for f in recover.sh cmp_h5dump.py; do [ -f "$PROD/$f" ] || { echo "FATAL: $PROD/$f missing"; exit 2; }; done
[ -f "$SRC/make_pp_dirs.py" ] || { echo "FATAL: $SRC/make_pp_dirs.py missing"; exit 2; }
if ! command -v h5ls >/dev/null || ! command -v h5dump >/dev/null; then echo "FATAL: h5ls/h5dump not found (hdf5 module)"; exit 2; fi
source "$PROD/recover.sh"
echo "restart_test: bin=$GPU_PP_BIN work=$W NBIN=$NBIN NSWEEP=$NSWEEP ranks=$NTASKS host=$(hostname) $(date)"

SUMMARY=""; NFAIL=0
res() {   # res NAME STATUS(PASS|FAIL|INFO) DETAIL
  SUMMARY="$SUMMARY
$(printf '%-5s %-50s %s' "$2" "$1" "$3")"
  [ "$2" = FAIL ] && NFAIL=$((NFAIL + 1))
  echo "[check] $2 $1 $3"
  return 0
}
pf() { [ "$1" -eq 0 ] && echo PASS || echo FAIL; }

mkdirs() {   # mkdirs NAME : fresh run dir; same seed offset (0) in every variant -> identical seeds
  local d="${W:?}/${1:?}"
  rm -rf "$d"
  python3 "$SRC/make_pp_dirs.py" "$d" $L $ULIST $NPER "$NSWEEP" "$NBIN" 64 0.0 0 || { echo "FATAL: make_pp_dirs failed"; exit 2; }
}
set_nbin() {  # set_nbin DIR SET N : NBin in Temp_<SET>/parameters
  local f="$1/Temp_$2/parameters"
  sed "s/^\([[:space:]]*NBin[[:space:]]*=\).*/\1 $3/" "$f" > "$f.new" && mv "$f.new" "$f"
}
run_alf() {   # run_alf DIR TAG : mpiexec in DIR (foreground), logs the segment end to segments.log
  ( cd "$1" && mpiexec -n "$NTASKS" "$GPU_PP_BIN" > "alf_$2.log" 2>&1; rc=$?; rc_segment_end $rc >/dev/null; exit $rc )
}
kill_tree() { local c; for c in $(pgrep -P "$1" 2>/dev/null); do kill_tree "$c"; done; echo "$1"; }
min_bins() {   # smallest bin count over the sets of DIR ($1); -1 if some file cannot be read right now (being written)
  local d m=1000000 n
  for d in $(cd "$1" && rc_sets); do
    n=$(cd "$1" && rc_bins "$d/data.h5"); if [ "$n" -lt 0 ]; then echo -1; return; fi
    if [ "$n" -lt "$m" ]; then m=$n; fi
  done
  echo $m
}

# start ALF in the background, kill -9 everything as soon as every set has >= NBIN/2 bins; leaves RUNNING behind
killed_run() {   # killed_run DIR
  local D=$1 half=$((NBIN / 2)) pid t0=$SECONDS n pids
  ( cd "$D" && rc_prepare ) || { echo "FATAL: rc_prepare failed"; return 1; }
  ( cd "$D" && exec mpiexec -n "$NTASKS" "$GPU_PP_BIN" > alf_seg1_killed.log 2>&1 ) &
  pid=$!
  while :; do
    n=$(min_bins "$D")
    if [ "$n" -ge "$half" ]; then break; fi
    kill -0 $pid 2>/dev/null || { echo "run ended before the kill point (bins=$n)"; return 2; }
    if [ $((SECONDS - t0)) -gt "$TIMEOUT" ]; then echo "timeout waiting for $half bins"; kill -9 $pid; return 3; fi
    sleep 0.1
  done
  pids="$(kill_tree $pid) $(pgrep -f "$GPU_PP_BIN" 2>/dev/null | tr '\n' ' ')"
  kill -9 $pids 2>/dev/null
  wait $pid 2>/dev/null
  sleep 3
  if pgrep -f "$GPU_PP_BIN" >/dev/null 2>&1; then pkill -9 -f "$GPU_PP_BIN"; sleep 2; fi
  echo "killed with kill -9, bins per set: $(cd "$D" && rc_bins_all)"
  [ -e "$D/RUNNING" ] || { echo "RUNNING is gone after the kill (run had finished?)"; return 2; }
  return 0
}

# recover.sh (stand-alone), then rerun the remaining bins (per set NBIN - bins in data.h5)
recover_and_finish() {   # recover_and_finish DIR TAG
  local D=$1 tag=$2 i n rem
  bash "$PROD/recover.sh" "$D" prepare || { res "$tag recover.sh prepare" FAIL "returned non-zero"; return 1; }
  if [ -e "$D/RUNNING" ]; then res "$tag RUNNING removed by recover.sh" FAIL "still present"; else res "$tag RUNNING removed by recover.sh" PASS ""; fi
  res "$tag recovery action" INFO "$(tail -n 1 "$D/segments.log" | sed 's/.*recovery: //')"
  for i in $(cd "$D" && rc_sets | sed 's/Temp_//'); do
    n=$(cd "$D" && rc_bins "Temp_$i/data.h5"); rem=$((NBIN - n))
    if [ $rem -lt 1 ]; then res "$tag Temp_$i remaining bins" INFO "bins=$n >= NBIN, running 1 bin"; rem=1; fi
    set_nbin "$D" "$i" "$rem"
    echo "$tag Temp_$i has $n bins -> NBin=$rem"
  done
  run_alf "$D" seg2; res "$tag second segment exit code" "$(pf $?)" ""
  if [ -e "$D/RUNNING" ]; then res "$tag RUNNING after clean end" FAIL "still present"; else res "$tag RUNNING after clean end" PASS ""; fi
}

compare_runs() {   # compare_runs TAG DIR
  local tag=$1 D=$2 i k out dups r bad
  for i in 0 1; do
    out=$(python3 "$PROD/cmp_h5dump.py" "$W/ref/Temp_$i/data.h5" "$D/Temp_$i/data.h5"); r=$?
    echo "$out"
    dups=$(echo "$out" | sed -n 's/^DUPS=//p')
    res "$tag Temp_$i data.h5 vs reference" "$(pf $r)" "$(echo "$out" | grep -E '^(NBINS|MAXREL)' | tr '\n' ' ') duplicated bins: ${dups:-?}"
    if [ "${dups:-1}" = 0 ]; then
      bad=0
      for k in $(seq 0 $((NPER - 1))); do
        python3 "$PROD/cmp_h5dump.py" "$W/ref/Temp_$i/confout_$k.h5" "$D/Temp_$i/confout_$k.h5" >/dev/null 2>&1 || { bad=$((bad + 1)); echo "confout_$k differs"; }
      done
      res "$tag Temp_$i final confout_*.h5 vs reference" "$(pf $bad)" "$NPER files, $bad differ"
    else
      res "$tag Temp_$i final confout_*.h5 vs reference" INFO "skipped: one duplicated bin, so the final state is one bin behind the reference (expected)"
    fi
  done
}

mkdir -p "$W"
# ---------------------------------------------------------------- (a) reference
echo "=== (a) reference run: $NBIN bins uninterrupted"
mkdirs ref
( cd "$W/ref" && rc_prepare ) && run_alf "$W/ref" seg1; res "(a) reference run exit code" "$(pf $?)" ""
for i in 0 1; do
  n=$(cd "$W/ref" && rc_bins "Temp_$i/data.h5")
  if [ "$n" = "$NBIN" ]; then res "(a) Temp_$i has $NBIN bins" PASS ""; else res "(a) Temp_$i has $NBIN bins" FAIL "bins=$n"; fi
done
if [ -e "$W/ref/RUNNING" ]; then res "(a) RUNNING removed at clean end" FAIL ""; else res "(a) RUNNING removed at clean end" PASS ""; fi
if [ "$NFAIL" -gt 0 ]; then echo "reference run failed, aborting"; echo "$SUMMARY"; echo "RESTART TEST: FAIL"; exit 1; fi

# ---------------------------------------------------------------- (b) kill, recover, rerun
echo "=== (b) kill -9 after ~half of the bins, recover, finish"
mkdirs kill
cmp -s "$W/ref/seeds" "$W/kill/seeds"; res "(b) same seeds as reference" "$(pf $?)" ""
killed_run "$W/kill"; krc=$?
if [ $krc -ne 0 ]; then res "(b) kill during run" FAIL "not killed cleanly (code $krc): inconclusive, try a larger NSWEEP"; else
  res "(b) kill during run" PASS "bins per set at kill: $(cd "$W/kill" && rc_bins_all)"
  recover_and_finish "$W/kill" "(b)"
  compare_runs "(b)" "$W/kill"
fi

# ---------------------------------------------------------------- (c) kill + corrupt confout
echo "=== (c) kill -9, truncate Temp_0/confout_3.h5 to 100 bytes, recover (must restore Temp_0 from .segbak), finish"
mkdirs corrupt
cmp -s "$W/ref/seeds" "$W/corrupt/seeds"; res "(c) same seeds as reference" "$(pf $?)" ""
killed_run "$W/corrupt"; krc=$?
if [ $krc -ne 0 ]; then res "(c) kill during run" FAIL "not killed cleanly (code $krc): inconclusive"; else
  res "(c) kill during run" PASS "bins per set at kill: $(cd "$W/corrupt" && rc_bins_all)"
  cf="$W/corrupt/Temp_0/confout_3.h5"
  if [ -f "$cf" ]; then
    head -c 100 "$cf" > "$cf.tmp" && mv "$cf.tmp" "$cf"; res "(c) truncated Temp_0/confout_3.h5" INFO "now $(wc -c < "$cf") bytes"
  else res "(c) truncate Temp_0/confout_3.h5" FAIL "file missing"; fi
  recover_and_finish "$W/corrupt" "(c)"
  if grep -q "restored Temp_0" "$W/corrupt/segments.log"; then res "(c) recover.sh restored Temp_0 from .segbak" PASS ""; else res "(c) recover.sh restored Temp_0 from .segbak" FAIL "no 'restored Temp_0' in segments.log"; fi
  compare_runs "(c)" "$W/corrupt"
fi

echo; echo "================= segments.log (b: kill) ================="; cat "$W/kill/segments.log" 2>/dev/null
echo "================= segments.log (c: corrupt) =============="; cat "$W/corrupt/segments.log" 2>/dev/null
echo; echo "================= SUMMARY ================="; echo "$SUMMARY"; echo
if [ "$NFAIL" -eq 0 ]; then echo "RESTART TEST: PASS"; exit 0; else echo "RESTART TEST: FAIL ($NFAIL check(s) failed)"; exit 1; fi

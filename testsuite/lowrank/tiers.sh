#!/usr/bin/env bash
# Validation tiers 2-5 for one low-rank scheme (project 2, Tasks 4 and 6). Run from anywhere; long (~30-60 min):
# start it in the background.
#
# Usage:  testsuite/lowrank/tiers.sh SCHEME SEED [PHASE ...]
#   SCHEME  delayed | submatrix
#   SEED    seed of the candidate statistical chains (3 for delayed, 4 for submatrix; standard references: 1, 2)
#   PHASE   any of: traj debug stat regress (default: all, in this order)
# Output: one summary line per check on stdout; logs in runs/tiers/<SCHEME>/.
#   traj    compare_trajectory.sh <set> standard SCHEME for every set of SETS and TRAJ_ONLY_SETS, plus tv/kondo/jz at
#           N_delay=2 and 1000, and for submatrix also delayed vs submatrix for every set.
#   debug   builds an -DALF_LOWRANK_DEBUG binary into builds/ALF_debug.out (in place, then restores the normal build) and
#           runs tv kondo z2 spin_peierls jz tv_proj tv_flux with it (NSweep=3 NBin=2); prints the shadow summary.
#   stat    stat_runs.sh SET SCHEME SEED for tv kondo jz tv_proj z2 doped (4 in parallel), then stat_compare.py vs
#           standard_s1, and the paired stability ratios (Precision Green/Phase max and mean, SCHEME_sSEED / standard_sSEED;
#           the standard chain with seed SEED is run as well if it does not exist: same seed => same trajectory).
#   regress regress_vs_baseline.sh and ctest.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
SCHEME="${1:?scheme}"; SEED="${2:?seed}"; shift 2
PHASES="${*:-traj debug stat regress}"
LOG="$ALFDIR/runs/tiers/$SCHEME"; mkdir -p "$LOG"
cd "$ALFDIR"
source "$HERE/sets.sh"
setnames() { for e in "${SETS[@]}" "${TRAJ_ONLY_SETS[@]}"; do printf '%s\n' "${e%%|*}"; done; }
build_env() { source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1; }

if [[ " $PHASES " == *" traj "* ]]; then
  first="--rebuild"
  for s in $(setnames); do
    extra=""; [ "$s" = z2gt ] && extra="N_delay=5"
    r=$(testsuite/lowrank/compare_trajectory.sh $first "$s" standard "$SCHEME" $extra >"$LOG/traj_$s.log" 2>&1; echo $?)
    first=""
    echo "traj $s standard-$SCHEME exit=$r $(grep -E 'RESULT|DIVERGED' "$LOG/traj_$s.log" | head -1)"
  done
  for s in tv kondo jz; do for nd in 2 1000; do
    r=$(testsuite/lowrank/compare_trajectory.sh "$s" standard "$SCHEME" N_delay=$nd >"$LOG/traj_${s}_nd$nd.log" 2>&1; echo $?)
    echo "traj $s standard-$SCHEME N_delay=$nd exit=$r $(grep -E 'RESULT|DIVERGED' "$LOG/traj_${s}_nd$nd.log" | head -1)"
  done; done
  if [ "$SCHEME" = submatrix ]; then
    for s in $(setnames); do
      r=$(testsuite/lowrank/compare_trajectory.sh "$s" delayed submatrix >"$LOG/traj_${s}_dl.log" 2>&1; echo $?)
      echo "traj $s delayed-submatrix exit=$r $(grep -E 'RESULT|DIVERGED' "$LOG/traj_${s}_dl.log" | head -1)"
    done
  fi
fi

if [[ " $PHASES " == *" debug "* ]]; then
  ( build_env && make -C Prog clean >/dev/null 2>&1 && ALF_FLAGS_EXT="-DALF_LOWRANK_DEBUG" make program -j8 >"$LOG/debug_build.log" 2>&1 \
      && mkdir -p builds && cp Prog/ALF.out builds/ALF_debug.out ) || echo "debug build FAILED (see $LOG/debug_build.log)"
  ( build_env && make -C Prog clean >/dev/null 2>&1 && make program -j8 >"$LOG/normal_build.log" 2>&1 ) || echo "normal rebuild FAILED"
  for s in tv kondo z2 spin_peierls jz tv_proj tv_flux; do
    keys=""; for e in "${SETS[@]}" "${TRAJ_ONLY_SETS[@]}"; do IFS='|' read -r n common _ <<<"$e"; [ "$n" = "$s" ] && keys="$common"; done
    d="$ALFDIR/runs/tiers/$SCHEME/debug_$s"; /bin/rm -rf "$d"
    read -r -a keyarr <<<"NSweep=3 NBin=2 $keys Sequential_scheme=\"$SCHEME\""
    python3 -I "$HERE/make_params.py" "$d" "${keyarr[@]}" >/dev/null
    r=$(cd "$d" && OMP_NUM_THREADS=1 "$ALFDIR/builds/ALF_debug.out" >stdout.log 2>&1; echo $?)
    echo "debug $s exit=$r $(tr '\n' ' ' < "$d/lowrank_debug_summary.txt" 2>/dev/null)"
  done
fi

if [[ " $PHASES " == *" stat "* ]]; then
  STATS="tv kondo jz tv_proj z2 doped"
  # candidate chains, and the standard chain with the same seed (paired stability) where it does not exist yet
  jobs=""; for s in $STATS; do jobs="$jobs $s:$SCHEME"; [ -e "runs/stat/$s/standard_s$SEED/info" ] || jobs="$jobs $s:standard"; done
  printf '%s\n' $jobs | xargs -P 4 -I{} sh -c 'j={}; testsuite/lowrank/stat_runs.sh ${j%%:*} ${j##*:} '"$SEED"' >"'"$LOG"'/stat_${j%%:*}_${j##*:}.log" 2>&1 || echo "stat run $j FAILED"'
  for s in $STATS; do
    v=$(python3 -I "$HERE/stat_compare.py" "runs/stat/$s/standard_s1" "runs/stat/$s/${SCHEME}_s$SEED" 2>&1 | tee "$LOG/statcmp_$s.log" | grep -E "VERDICT" | head -1)
    echo "stat $s ${SCHEME}_s$SEED vs standard_s1: $v"
    if [ -e "runs/stat/$s/standard_s$SEED/info" ]; then
      python3 -I - "runs/stat/$s/${SCHEME}_s$SEED/info" "runs/stat/$s/standard_s$SEED/info" <<'PY'
import re, sys
def prec(f):
    t = open(f).read()
    g = re.search(r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)", t)
    p = re.search(r"Precision Phase,\s*Max\s*:\s*(\S+)", t)
    return [float(g.group(1)), float(g.group(2)), float(p.group(1)) if p else float("nan")]
a, b = prec(sys.argv[1]), prec(sys.argv[2])
print("  paired stability cand/std: Green mean %.2f max %.2f, Phase max %.2f" % tuple(x / y if y else float("nan") for x, y in zip(a, b)))
PY
    else
      echo "  (no standard_s$SEED: paired stability skipped)"
    fi
  done
fi

if [[ " $PHASES " == *" regress "* ]]; then
  testsuite/lowrank/regress_vs_baseline.sh >"$LOG/regress.log" 2>&1
  echo "regress: $(grep -c ' IDENTICAL' "$LOG/regress.log") IDENTICAL, $(grep -c ' DIFFERENT' "$LOG/regress.log") DIFFERENT"
  ( build_env && cd testsuite/tests && cmake .. >/dev/null && make -j8 >/dev/null 2>&1 && ctest -j8 2>&1 | grep -E "tests passed|tests failed" )
fi

#!/bin/bash
# recover.sh -- restart / preemption recovery for ALF PARALLEL_PARAMS (PP) runs. Sourced by job_preempt.sh, or run stand-alone for tests.
#
# Run-directory layout (cwd = run directory): Temp_<i>/data.h5, Temp_<i>/confin_<rank>.h5, Temp_<i>/confout_<rank>.h5, RUNNING
# ALF per bin: (1) append to Temp_<i>/data.h5, (2) write Temp_<i>/confout_<rank>.h5. At start it reads confin_<rank>.h5 if present.
# A SIGKILL (Slurm preemption, GraceTime=0) can therefore leave
#   - RUNNING behind (ALF would abort on the next start),
#   - a half-written data.h5 or confout (corrupt, or data.h5 whose datasets have different numbers of bins),
#   - one bin more in data.h5 than in the confout state (-> that bin is computed twice after the restart: ONE DUPLICATED BIN, accepted).
#
# Functions (all prefixed rc_, safe under `set -u`, no `set -e` needed; they never exit the calling shell, they return non-zero):
#   rc_require_h5     abort-guard: h5ls must be available (otherwise every file would look corrupt)
#   rc_bins FILE      print number of bins of a data.h5: N (all datasets agree), -1 unreadable, -2 datasets disagree (kill mid-write), 0 missing
#   rc_stop_guard     return 0 if every Temp_*/data.h5 has >= $TARGET_BINS bins (nothing left to do), else 1
#   rc_recover        step b: if RUNNING exists: check files, restore Temp_i from Temp_i/.segbak where damaged, remove RUNNING
#   rc_out_to_in      step c: confout_* -> confin_* in every Temp_i (refuses if plain and .h5 confout coexist)
#   rc_backup         step d: Temp_i/.segbak/ := data.h5 + confin_*  (cp -p; old backup cleared first; COMPLETE marker last)
#   rc_log_segment    step e: append a line to segments.log (date, job id, restart count, bins per set, recovery action)
#   rc_prepare        b + c + d + e in this order; return non-zero = do NOT start ALF
#   rc_segment_end RC log exit code and bins after the run to segments.log
# Stand-alone:  bash recover.sh RUNDIR prepare|guard|bins   (env TARGET_BINS for guard; exit status of guard: 10 = target reached)

RC_ACTION=${RC_ACTION:-none}
RC_LOG_FILE=${RC_LOG_FILE:-segments.log}

rc_ts()  { date '+%Y-%m-%d %H:%M:%S'; }
rc_log() { echo "[recover $(rc_ts)] $*"; }

# Temp_<number> directories of the current directory, numerically sorted
rc_sets() {
  local d
  for d in $(ls -d Temp_[0-9]* 2>/dev/null | sort -t_ -k2 -n); do
    case "$d" in Temp_*[!0-9]*) continue ;; esac
    [ -d "$d" ] && echo "$d"
  done
  return 0
}

rc_require_h5() {
  if ! command -v h5ls >/dev/null 2>&1; then
    rc_log "FATAL: h5ls not found (module load hdf5/1.14.5-intel25); refusing to judge any file"
    return 1
  fi
  return 0
}

# Fills RC_MIN / RC_MAX with the smallest / largest number of bins over all bin datasets (<group>/sign, /obser, /back; first dimension in
# h5ls = C order = last Fortran dimension = bin index, see init_dset/append_dat in Libraries/Modules/alf_hdf5_mod.F90).
# Return: 0 ok, 1 file unreadable by h5ls, 2 file missing.
rc_info() {
  RC_MIN=0; RC_MAX=0
  local f=$1 out res
  [ -f "$f" ] || return 2
  out=$(h5ls -r "$f" 2>/dev/null) || return 1
  res=$(printf '%s\n' "$out" | awk '
    /Dataset/ && match($0, /\{[0-9]+/) {
      name = $1
      if (name ~ /\/(sign|obser|back)$/ && name !~ /Acc_Temp/) {
        n = substr($0, RSTART + 1, RLENGTH - 1) + 0
        if (!c || n < mn) mn = n
        if (!c || n > mx) mx = n
        c++
      }
    }
    END { if (c) print mn, mx; else print 0, 0 }')
  set -- $res
  RC_MIN=${1:-0}; RC_MAX=${2:-0}
  return 0
}

rc_bins() {
  local rc
  rc_info "$1"; rc=$?
  if   [ $rc -eq 2 ]; then echo 0
  elif [ $rc -ne 0 ]; then echo -1
  elif [ "$RC_MIN" != "$RC_MAX" ]; then echo -2
  else echo "$RC_MIN"; fi
}

# "0:12 1:12 2:-1 ..." for the log
rc_bins_all() {
  local d s=""
  for d in $(rc_sets); do s="$s ${d#Temp_}:$(rc_bins "$d/data.h5")"; done
  echo "${s# }"
}

# a) stop guard: 0 = every set has >= TARGET_BINS valid bins (also needs at least one Temp_ dir)
rc_stop_guard() {
  local target=${TARGET_BINS:-0} d n any=0
  [ "$target" -gt 0 ] 2>/dev/null || return 1
  rc_require_h5 || return 1
  for d in $(rc_sets); do
    any=1
    n=$(rc_bins "$d/data.h5")
    [ "$n" -ge "$target" ] 2>/dev/null || return 1
  done
  [ $any -eq 1 ]
}

rc_restore_set() {   # $1 = Temp_i ; restore data.h5 + confin_* from .segbak, drop confout_*
  local d=$1 b=$1/.segbak f
  if [ ! -f "$b/COMPLETE" ]; then
    rc_log "FATAL: $d is damaged but $b/COMPLETE is missing (no usable backup)"
    return 1
  fi
  rm -f "$d"/confout_* "$d"/confin_* "$d"/data.h5 || return 1
  if [ -f "$b/data.h5" ]; then cp -p "$b/data.h5" "$d/data.h5" || return 1; fi
  for f in "$b"/confin_*; do
    [ -e "$f" ] || continue
    cp -p "$f" "$d/" || return 1
  done
  return 0
}

# b) recovery after a killed segment
rc_recover() {
  RC_ACTION=none
  [ -e RUNNING ] || return 0
  rc_require_h5 || return 1
  RC_ACTION="RUNNING found (previous segment killed);"
  local d f bad rc nkept=0 nrest=0
  for d in $(rc_sets); do
    bad=""
    if [ -f "$d/data.h5" ]; then
      rc_info "$d/data.h5"; rc=$?
      if   [ $rc -ne 0 ]; then bad="data.h5 unreadable"
      elif [ "$RC_MIN" != "$RC_MAX" ]; then bad="data.h5 datasets disagree (bins $RC_MIN..$RC_MAX)"
      fi
    fi
    for f in "$d"/confout_*.h5; do
      [ -e "$f" ] || continue
      h5ls -r "$f" >/dev/null 2>&1 || bad="${bad:+$bad; }${f##*/} unreadable"
    done
    if [ -n "$bad" ]; then
      rc_log "$d: DAMAGED ($bad) -> restoring from $d/.segbak"
      rc_restore_set "$d" || { RC_ACTION="$RC_ACTION restore of $d FAILED;"; return 1; }
      RC_ACTION="$RC_ACTION restored $d ($bad);"
      nrest=$((nrest + 1))
    else
      nkept=$((nkept + 1))
    fi
  done
  rc_log "previous segment was killed: $nkept set(s) kept as they are, $nrest restored from backup"
  RC_ACTION="$RC_ACTION kept $nkept, restored $nrest, RUNNING removed"
  rm -f RUNNING
  return 0
}

# c) confout -> confin per Temp_i (same as Start_PP_GPU/out_to_in_pp.sh, plus a check of every file)
rc_out_to_in() {
  local d f nh np
  for d in $(rc_sets); do
    nh=0; np=0
    for f in "$d"/confout_*; do
      [ -e "$f" ] || continue
      case "$f" in *.h5) nh=1 ;; *) np=1 ;; esac
    done
    if [ $nh -eq 1 ] && [ $np -eq 1 ]; then
      rc_log "FATAL: $d: both plain text and HDF5 confout present"
      return 1
    fi
    for f in "$d"/confout_*; do
      [ -e "$f" ] || continue
      mv -f "$f" "$d/confin_${f##*/confout_}" || return 1
    done
  done
  return 0
}

# d) backup of the state a new segment starts from
rc_backup() {
  local d b f
  for d in $(rc_sets); do
    b=$d/.segbak
    rm -rf "$b" || return 1
    mkdir "$b" || return 1
    if [ -f "$d/data.h5" ]; then cp -p "$d/data.h5" "$b/" || return 1; fi
    for f in "$d"/confin_*; do
      [ -e "$f" ] || continue
      cp -p "$f" "$b/" || return 1
    done
    : > "$b/COMPLETE"
  done
  return 0
}

# e) one line per segment start
rc_log_segment() {
  local line
  line="$(rc_ts) | job=${SLURM_JOB_ID:-none} | restart_count=${SLURM_RESTART_COUNT:-0} | start | bins_before=[$(rc_bins_all)] | recovery: $RC_ACTION"
  echo "$line" >> "$RC_LOG_FILE"
  rc_log "$line"
}

rc_segment_end() {   # $1 = exit code of mpiexec
  local line
  line="$(rc_ts) | job=${SLURM_JOB_ID:-none} | restart_count=${SLURM_RESTART_COUNT:-0} | end   | exit_code=${1:-?} | bins_after=[$(rc_bins_all)]"
  echo "$line" >> "$RC_LOG_FILE"
  rc_log "$line"
}

rc_prepare() {
  rc_recover   || { rc_log "recovery failed, not starting ALF"; return 1; }
  rc_out_to_in || { rc_log "out_to_in failed, not starting ALF"; return 1; }
  rc_backup    || { rc_log "backup failed, not starting ALF"; return 1; }
  rc_log_segment
  return 0
}

# ---- stand-alone use ------------------------------------------------------------------------------------------------------------
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  set -u
  [ $# -ge 1 ] || { echo "usage: bash recover.sh RUNDIR [prepare|guard|bins]   (env TARGET_BINS for guard)" >&2; exit 2; }
  cd "$1" || exit 2
  case "${2:-prepare}" in
    prepare) rc_prepare; exit $? ;;
    guard)   if rc_stop_guard; then echo "DONE: all sets have >= ${TARGET_BINS:-0} bins"; exit 10; else echo "CONTINUE"; exit 0; fi ;;
    bins)    rc_bins_all; exit 0 ;;
    *)       echo "unknown action $2" >&2; exit 2 ;;
  esac
fi

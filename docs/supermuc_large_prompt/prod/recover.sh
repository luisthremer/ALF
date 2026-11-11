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
#   rc_recover        step b: if RUNNING exists: check files; a damaged Temp_i is restored from its newest complete rolling snapshot
#                     (Temp_i/.snap/A|B), else the other slot, else Temp_i/.segbak (segment start); undamaged sets are not touched; remove RUNNING
#   rc_out_to_in      step c: confout_* -> confin_* in every Temp_i (refuses if plain and .h5 confout coexist)
#   rc_backup         step d: Temp_i/.segbak/ := data.h5 + confin_*  (cp -p; old backup cleared first; COMPLETE marker last)
#   rc_log_segment    step e: append a line to segments.log (date, job id, restart count, bins per set, recovery action)
#   rc_prepare        b + c + d + e in this order; return non-zero = do NOT start ALF
#   rc_segment_end RC log exit code and bins after the run to segments.log
#   rc_snapshot_loop  background loop (started by job_preempt.sh right before mpiexec): rolling snapshots of every set, see below
#   rc_snap_try SET   one snapshot attempt for one set (used by the loop; stand-alone: bash recover.sh RUNDIR snap)
# Rolling snapshots (Temp_i/.snap/A and .snap/B, each: data.h5 + confout_*.h5 + COMPLETE "bins=N time=T seq=S walkers=W"):
#   a set gets a snapshot when, at the same time, (env in brackets)
#     - >= SNAP_INTERVAL s (1800) since its newest complete snapshot (state is read from the COMPLETE markers, so a requeue continues correctly),
#     - data.h5 has a consistent bin count >= 1 larger than that of the newest snapshot,
#     - the bin is complete: as many confout_*.h5 as walkers (mpi_per_parameter_set), all readable, mtime >= data.h5 mtime - SNAP_CONF_SLACK (0 s),
#     - data.h5 untouched for >= SNAP_QUIET_DATA s (10), the confouts for >= SNAP_QUIET_CONF s (5).
#   Procedure: sizes+mtimes recorded, copy into .snap/tmp, sources re-checked (unchanged) and the copy verified (bin count, readability),
#   COMPLETE written, then the OLDER slot is replaced by mv (the newer slot stays valid at every moment). Log: snapshots.log.
#   The loop polls every SNAP_POLL s (30). RC_SNAP_TEST_DELAY (s) sleeps after the first copied file (tests: kill during a copy).
# Stand-alone:  bash recover.sh RUNDIR prepare|guard|bins|snap   (env TARGET_BINS for guard; exit status of guard: 10 = target reached)

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

# IMPORTANT: every read of an HDF5 file that ALF may have open goes through this wrapper. With HDF5's default file locking a concurrent
# h5ls holds a read lock and makes ALF's next H5Fopen(RDWR) of data.h5 fail ("unable to lock the file") -> ALF aborts (seen on SuperMUC-NG).
# Without locking, a reader cannot disturb ALF; a read of a file in the middle of a write simply fails or looks inconsistent and is retried / discarded.
rc_h5ls() { HDF5_USE_FILE_LOCKING=FALSE h5ls "$@"; }

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
  out=$(rc_h5ls -r "$f" 2>/dev/null) || return 1
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

# ---- rolling snapshots -----------------------------------------------------------------------------------------------------------
rc_slot_meta() {   # $1 = slot dir; sets M_BINS M_TIME M_SEQ M_W from its COMPLETE marker; return 1 if there is no valid marker
  M_BINS=""; M_TIME=""; M_SEQ=""; M_W=""
  local line kv
  [ -f "$1/COMPLETE" ] || return 1
  read -r line < "$1/COMPLETE" || true
  for kv in $line; do
    case "$kv" in
      bins=*) M_BINS=${kv#bins=} ;; time=*) M_TIME=${kv#time=} ;; seq=*) M_SEQ=${kv#seq=} ;; walkers=*) M_W=${kv#walkers=} ;;
    esac
  done
  [ -n "$M_BINS" ] && [ -n "$M_TIME" ] && [ -n "$M_SEQ" ] && [ -n "$M_W" ]
}

rc_snap_order() {   # $1 = Temp_i ; valid slots (marker present), newest first, one per line
  local sa=-1 sb=-1
  if rc_slot_meta "$1/.snap/A"; then sa=$M_SEQ; fi
  if rc_slot_meta "$1/.snap/B"; then sb=$M_SEQ; fi
  if [ "$sa" -ge "$sb" ]; then
    [ "$sa" -ge 0 ] && echo A
    [ "$sb" -ge 0 ] && echo B
  else
    echo B
    [ "$sa" -ge 0 ] && echo A
  fi
  return 0
}

rc_snap_target() {   # $1 = Temp_i ; slot to overwrite: an invalid one, else the older one
  local sa sb
  rc_slot_meta "$1/.snap/A" || { echo A; return 0; }; sa=$M_SEQ
  rc_slot_meta "$1/.snap/B" || { echo B; return 0; }; sb=$M_SEQ
  if [ "$sa" -le "$sb" ]; then echo A; else echo B; fi
}

rc_mtime() { stat -c %.9Y "$1" 2>/dev/null || stat -f %Fm "$1" 2>/dev/null; }   # GNU || BSD, seconds with fraction
rc_size()  { stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null; }
rc_now_f() { local t=".rc_now.$$"; : > "$t"; rc_mtime "$t"; rm -f "$t"; }          # "now" on the file system clock
rc_ge()    { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 >= b + 0) }'; }          # a >= b (floats)

rc_walkers() {   # mpi_per_parameter_set from ./parameters
  sed -n 's/^[[:space:]]*[Mm][Pp][Ii]_[Pp][Ee][Rr]_[Pp][Aa][Rr][Aa][Mm][Ee][Tt][Ee][Rr]_[Ss][Ee][Tt][[:space:]]*=[[:space:]]*\([0-9][0-9]*\).*/\1/p' parameters 2>/dev/null | head -n 1
}

rc_sigs() {   # $1 = Temp_i ; "path size:mtime" of data.h5 and all confout_*.h5 (change detector)
  local f
  for f in "$1/data.h5" "$1"/confout_*.h5; do
    [ -e "$f" ] && echo "$f $(rc_size "$f"):$(rc_mtime "$f")"
  done
  return 0
}

rc_snap_log() { echo "$(rc_ts) | $*" >> "${SNAP_LOG_FILE:-snapshots.log}"; }

rc_snap_why() {   # log a blocking reason only when it changes (per set, state in Temp_i/.snap/why)
  mkdir -p "$1/.snap"
  if [ "$(cat "$1/.snap/why" 2>/dev/null)" != "$2" ]; then echo "$2" > "$1/.snap/why"; rc_snap_log "$1 | not snapshotting: $2"; fi
}

# one snapshot attempt for one set; always returns 0
rc_snap_try() {
  local d=$1 w first last_t="" last_b=0 now nowf dm cm n nconf f confs="" sig1 sig2 tmp slot seq t0 qd qc slack iv off omin="" omax=""
  iv=${SNAP_INTERVAL:-1800}; qd=${SNAP_QUIET_DATA:-10}; qc=${SNAP_QUIET_CONF:-5}; slack=${SNAP_CONF_SLACK:-0}
  [ -f "$d/data.h5" ] || return 0
  w=$(rc_walkers)
  if [ -z "$w" ]; then rc_snap_why "$d" "cannot read mpi_per_parameter_set from ./parameters"; return 0; fi
  first=$(rc_snap_order "$d" | head -n 1)
  if [ -n "$first" ] && rc_slot_meta "$d/.snap/$first"; then last_t=$M_TIME; last_b=$M_BINS; fi
  now=$(date +%s)
  if [ -n "$last_t" ] && [ $((now - last_t)) -lt "$iv" ]; then return 0; fi
  # cheap file system checks first
  nconf=0
  for f in "$d"/confout_*.h5; do [ -e "$f" ] && { nconf=$((nconf + 1)); confs="$confs $f"; }; done
  if [ "$nconf" -ne "$w" ]; then
    [ "$nconf" -gt 0 ] && rc_snap_why "$d" "number of confout files differs from the number of walkers ($w): bin not complete"
    return 0
  fi
  nowf=$(rc_now_f); dm=$(rc_mtime "$d/data.h5")
  rc_ge "$(awk -v a="$nowf" -v b="$dm" 'BEGIN { print a - b }')" "$qd" || return 0
  for f in $confs; do
    cm=$(rc_mtime "$f")
    rc_ge "$(awk -v a="$nowf" -v b="$cm" 'BEGIN { print a - b }')" "$qc" || return 0
    off=$(awk -v a="$cm" -v b="$dm" 'BEGIN { printf "%.4f", a - b }')     # confout mtime minus data.h5 mtime (diagnostic)
    omin=$(awk -v x="$omin" -v y="$off" 'BEGIN { print (x == "" || y + 0 < x + 0) ? y : x }')
    omax=$(awk -v x="$omax" -v y="$off" 'BEGIN { print (x == "" || y + 0 > x + 0) ? y : x }')
    if ! rc_ge "$cm" "$(awk -v a="$dm" -v b="$slack" 'BEGIN { print a - b }')"; then
      rc_snap_why "$d" "${f##*/} is older than data.h5 (bin not complete, or ranks write confout before rank 0 finishes data.h5: raise SNAP_CONF_SLACK)"
      return 0
    fi
  done
  n=$(rc_bins "$d/data.h5")
  if [ "$n" -lt 1 ]; then rc_snap_why "$d" "data.h5 bin count is $n (unreadable or datasets disagree)"; return 0; fi
  [ "$n" -gt "$last_b" ] || return 0
  for f in $confs; do rc_h5ls -r "$f" >/dev/null 2>&1 || { rc_snap_why "$d" "${f##*/} unreadable"; return 0; }; done
  # copy
  mkdir -p "$d/.snap" || return 0
  tmp=$d/.snap/tmp; rm -rf "$tmp"; mkdir "$tmp" || return 0
  t0=$SECONDS
  sig1=$(rc_sigs "$d")
  cp -p "$d/data.h5" "$tmp/" || { rm -rf "$tmp"; return 0; }
  if [ -n "${RC_SNAP_TEST_DELAY:-}" ]; then sleep "$RC_SNAP_TEST_DELAY"; fi
  for f in $confs; do cp -p "$f" "$tmp/" || { rm -rf "$tmp"; return 0; }; done
  sig2=$(rc_sigs "$d")
  if [ "$sig1" != "$sig2" ]; then rc_snap_log "$d | bins=$n | DISCARDED: sources changed during the copy"; rm -rf "$tmp"; return 0; fi
  if [ "$(rc_bins "$tmp/data.h5")" != "$n" ]; then rc_snap_log "$d | bins=$n | DISCARDED: copied data.h5 does not have $n bins"; rm -rf "$tmp"; return 0; fi
  for f in "$tmp"/confout_*.h5; do
    rc_h5ls -r "$f" >/dev/null 2>&1 || { rc_snap_log "$d | bins=$n | DISCARDED: copy of ${f##*/} unreadable"; rm -rf "$tmp"; return 0; }
  done
  seq=0
  for slot in A B; do if rc_slot_meta "$d/.snap/$slot" && [ "$M_SEQ" -gt "$seq" ]; then seq=$M_SEQ; fi; done
  seq=$((seq + 1))
  echo "bins=$n time=$(date +%s) seq=$seq walkers=$w" > "$tmp/COMPLETE.new" && mv "$tmp/COMPLETE.new" "$tmp/COMPLETE" || { rm -rf "$tmp"; return 0; }
  slot=$(rc_snap_target "$d")
  rm -rf "$d/.snap/$slot"
  mv "$tmp" "$d/.snap/$slot" || return 0
  rm -f "$d/.snap/why"
  rc_snap_log "$d | bins=$n | slot=$slot seq=$seq | took $((SECONDS - t0)) s | confout-minus-data.h5 mtime: min $omin max $omax s"
  return 0
}

rc_snapshot_loop() {
  local d sp
  rc_require_h5 || return 1
  trap '[ -n "${sp:-}" ] && kill "$sp" 2>/dev/null; exit 0' TERM INT
  rc_snap_log "snapshot loop started: interval=${SNAP_INTERVAL:-1800}s poll=${SNAP_POLL:-30}s quiet=${SNAP_QUIET_DATA:-10}/${SNAP_QUIET_CONF:-5}s slack=${SNAP_CONF_SLACK:-0}s pid=$$"
  while :; do
    for d in $(rc_sets); do rc_snap_try "$d"; done
    sleep "${SNAP_POLL:-30}" & sp=$!
    wait "$sp"
  done
}

# ---- restore ---------------------------------------------------------------------------------------------------------------------
rc_restore_slot() {   # $1 = Temp_i, $2 = slot dir ; return 0 if restored from it
  local d=$1 s=$2 f b n=0
  rc_slot_meta "$s" || return 1
  [ -f "$s/data.h5" ] || return 1
  b=$(rc_bins "$s/data.h5")
  [ "$b" = "$M_BINS" ] || { rc_log "$s: data.h5 has $b bins, marker says $M_BINS"; return 1; }
  for f in "$s"/confout_*.h5; do
    [ -e "$f" ] || continue
    rc_h5ls -r "$f" >/dev/null 2>&1 || { rc_log "$s: ${f##*/} unreadable"; return 1; }
    n=$((n + 1))
  done
  [ "$n" = "$M_W" ] || { rc_log "$s: $n confout files, marker says $M_W"; return 1; }
  rm -f "$d"/confout_* "$d"/confin_* "$d/data.h5" || return 1
  for f in "$s"/confout_*.h5; do cp -p "$f" "$d/confin_${f##*/confout_}" || return 1; done
  cp -p "$s/data.h5" "$d/data.h5" || return 1      # data.h5 last
  RC_RBINS=$b
  return 0
}

# $1 = Temp_i, $2 = bins before the damage ("?" if unreadable). Sets RC_RWHAT (what was used) and RC_RBINS (bins restored).
rc_restore_set() {
  local d=$1 before=${2:-?} b=$1/.segbak f s
  RC_RWHAT=""; RC_RBINS=0
  : > "$d/.restoring"           # present until the restore is complete: a kill in between makes the next start restore again
  for s in $(rc_snap_order "$d"); do
    if rc_restore_slot "$d" "$d/.snap/$s"; then RC_RWHAT="snapshot $s"; break; fi
    rc_log "$d: snapshot slot $s not usable"
  done
  if [ -z "$RC_RWHAT" ]; then
    if [ ! -f "$b/COMPLETE" ]; then rc_log "FATAL: $d is damaged but there is neither a usable snapshot nor $b/COMPLETE"; return 1; fi
    rm -f "$d"/confout_* "$d"/confin_* "$d/data.h5" || return 1
    for f in "$b"/confin_*; do [ -e "$f" ] || continue; cp -p "$f" "$d/" || return 1; done
    if [ -f "$b/data.h5" ]; then cp -p "$b/data.h5" "$d/data.h5" || return 1; RC_RBINS=$(rc_bins "$d/data.h5"); fi
    RC_RWHAT="segment-start backup"
  fi
  rm -f "$d/.restoring"
  rc_log "$d: restored from $RC_RWHAT: $RC_RBINS bins (before the damage: $before) -> bins lost: $( [ "$before" = "?" ] && echo "unknown (data.h5 unreadable)" || echo $((before - RC_RBINS)) )"
  return 0
}

# b) recovery after a killed segment
rc_recover() {
  RC_ACTION=none
  [ -e RUNNING ] || return 0
  rc_require_h5 || return 1
  RC_ACTION="RUNNING found (previous segment killed);"
  local d f bad rc before nkept=0 nrest=0
  for d in $(rc_sets); do
    bad=""; before="?"
    [ -e "$d/.restoring" ] && bad="interrupted restore"
    if [ -f "$d/data.h5" ]; then
      rc_info "$d/data.h5"; rc=$?
      if   [ $rc -ne 0 ]; then bad="${bad:+$bad; }data.h5 unreadable"
      elif [ "$RC_MIN" != "$RC_MAX" ]; then bad="${bad:+$bad; }data.h5 datasets disagree (bins $RC_MIN..$RC_MAX)"; before=$RC_MIN
      else before=$RC_MIN
      fi
    fi
    for f in "$d"/confout_*.h5; do
      [ -e "$f" ] || continue
      rc_h5ls -r "$f" >/dev/null 2>&1 || bad="${bad:+$bad; }${f##*/} unreadable"
    done
    if [ -n "$bad" ]; then
      rc_log "$d: DAMAGED ($bad)"
      rc_restore_set "$d" "$before" || { RC_ACTION="$RC_ACTION restore of $d FAILED;"; return 1; }
      RC_ACTION="$RC_ACTION restored $d from $RC_RWHAT, bins $before -> $RC_RBINS ($bad);"
      nrest=$((nrest + 1))
    else
      nkept=$((nkept + 1))
    fi
  done
  rc_log "previous segment was killed: $nkept set(s) kept as they are, $nrest restored"
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
  [ $# -ge 1 ] || { echo "usage: bash recover.sh RUNDIR [prepare|guard|bins|snap]   (env TARGET_BINS for guard)" >&2; exit 2; }
  cd "$1" || exit 2
  case "${2:-prepare}" in
    prepare) rc_prepare; exit $? ;;
    guard)   if rc_stop_guard; then echo "DONE: all sets have >= ${TARGET_BINS:-0} bins"; exit 10; else echo "CONTINUE"; exit 0; fi ;;
    bins)    rc_bins_all; exit 0 ;;
    snap)    for d in $(rc_sets); do rc_snap_try "$d"; done; exit 0 ;;
    *)       echo "unknown action $2" >&2; exit 2 ;;
  esac
fi

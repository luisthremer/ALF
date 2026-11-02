#!/bin/bash
# Ticket 207-3 (development only, dev seeds), nividia: run one binary on one prepared config (single process, 1 thread).
#   usage: [SET=devN] [TAG=dir] [SETROOT=dir] t207_3_run.sh NAME CONFIG [env assignments...]
#   SETROOT: where the prepared seed sets are read (default $W/runs; e.g. the previous agent's t207-3/runs, read only)
#   -> t207-3/runs/TAG/CONFIG (refuses to overwrite): run.log, stamp.txt (binary sha256, source tree hash, env), provenance.json.
#   GPU binaries (nv_* except nv_none*) run inside flock ~/alf-gpu/gpu.lock; inside the lock the script waits (up to 60 min, polling)
#   until no compute process is on the GPU, so only one GPU process runs at a time.  No fixed run-time limit (decision 16); a 2 h guard
#   timeout only stops a hung run.  Same records as t207_d5_run.sh (the QREF pipeline reads them).
#   MINFREE_MB (optional, fix round t207-3c, for the dump-writing beta = 4 runs): inside the lock, once the GPU is free, the run starts
#   only if / has at least MINFREE_MB + 1600 * (other processes waiting on gpu.lock) MB free (a waiting run of another agent has
#   done its own disk check before the lock and will write ~1.6 GB of dumps after this run); else the lock is released, the script
#   sleeps 300 s and tries again, for at most DISKWAIT_S (default 28800) s, then gives up (rc=8, logged, nothing run).  Fix round
#   t207-3e: "disk too low" is signalled by the flag file $D/.diskwait (no in-band exit code); the waiter count in the check is the larger
#   of (a) the flock-command waiters seen by pgrep and (b) the processes named flock that have gpu.lock open (lsof -t, own user; this
#   sees `flock 9` waiters of a chain that holds the lock on a file descriptor), without our own ancestors.  Limit: lsof sees only
#   processes of the same user (other users' waiters are invisible; gpu.lock is used by the agents' common account only), and a waiter
#   that is not a process named flock (e.g. a python fcntl.flock caller) is not counted.
#   T207_3_LOCKED=1: the caller (t207_3_chain.sh with HOLDLOCK=1) already holds gpu.lock on its own file descriptor through the run and
#   its processing; no second flock (it would deadlock), waiting runs are not counted (they start only after the thinning), and a
#   failed disk check does NOT wait inside the held lock (t207-3e): the script leaves the flag $D/.diskwait and exits 0; the caller
#   (t207_3_chain.sh rp) releases the lock, sleeps, takes it again and calls the script again (it then reuses the run directory $D,
#   which holds only the copied inputs and stamp.txt, and removes the flag).  DISKWAIT_S is the give-up bound (rc=8) as without LOCKED;
#   a caller passes DISKWAIT_S=0 for its last attempt.
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
export TMPDIR=$W/tmp
N=$1; C=$2; shift 2; ENVS="$*"
case $N in gnu*) . ~/alf-gpu/env/gnu.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-dev1}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
# provenance.json of the run in the current directory: prov RC START END
prov() {
~/alf-gpu/venv/bin/python3 -I - "$W" "$N" "$S" "$C" "$T" "$1" "$2" "$3" "$ENVS" <<'PY'
import hashlib, json, os, socket, subprocess, sys
W, N, S, C, T, rc, s, e, env = sys.argv[1:10]
def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else None
def rd(p):
    return open(p).read().strip() if os.path.exists(p) else None
gpu = subprocess.run(['nvidia-smi', '--query-gpu=name,driver_version', '--format=csv,noheader'], capture_output=True, text=True).stdout.strip()
json.dump(dict(ticket='207-3 (development evidence, dev seeds; not acceptance)', tag=T, build=N, seed_set=S, config=C, rc=int(rc),
               wall_s=round(float(e) - float(s), 1), binary_sha256=sha(os.path.join(W, 'bin', 'ALF_%s.out' % N)),
               srctree=rd(os.path.join(W, 'bin', 'ALF_%s.srctree' % N)), srcstate=rd(os.path.join(W, 'bin', 'ALF_%s.srcstate' % N)),
               env=env.split(), seeds_sha256=sha('seeds'), parameters_sha256=sha('parameters'), run_log_sha256=sha('run.log'),
               info_sha256=sha('info'), data_h5_sha256=sha('data.h5'), host=socket.gethostname(), gpu=gpu),
          open('provenance.json', 'w'), indent=1)
PY
}
D="$W/runs/$T/$C"; RESUME=0; [ -e "$D/.diskwait" ] && RESUME=1      # resume after a disk wait (flag left by an earlier call)
# refusal: while resuming (flag $D/.diskwait set by an earlier call) the flag is dropped and, unless a provenance.json exists (never
# overwritten), a provenance.json with rc 1 is written, so the caller (t207_3_chain.sh rp) sees the failure instead of retrying forever (t207-3f)
bad() { echo "refusing: $*"
        if [ $RESUME = 1 ]; then rm -f "$D/.diskwait"
          [ -e "$D/provenance.json" ] || { local t; t=$(date +%s.%N); (cd "$D" && prov 1 $t $t); }
        fi; exit 1; }
[ -e "$D" ] && [ $RESUME = 0 ] && bad "$D exists"
[ -x "$W/bin/ALF_$N.out" ] || bad "binary $N missing"
SR=${SETROOT:-$W/runs}; [ -d "$SR/$S/$C" ] || bad "$SR/$S/$C missing"
if [ $RESUME = 1 ]; then   # a resumed directory must hold only the copied inputs and stamp.txt, with the source set's seeds and parameters
  for f in run.log data.h5 info provenance.json; do [ -e "$D/$f" ] && bad "cannot resume, $D/$f exists"; done
  for f in seeds parameters; do
    [ "$(sha256sum < "$D/$f" 2>/dev/null)" = "$(sha256sum < "$SR/$S/$C/$f" 2>/dev/null)" ] || bad "cannot resume, $f differs from $SR/$S/$C/$f"; done
fi
[ $RESUME = 0 ] && { mkdir -p "$W/runs/$T"; cp -r "$SR/$S/$C" "$D"; }
cd "$D" || exit 1; mkdir -p qref; rm -f .diskwait   # dumps of ALF_GPU_TAUDIAG_QREF
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
{ [ $RESUME = 1 ] && echo "resumed after a disk wait $(date -Is)"; echo "binary $(cat $W/bin/ALF_$N.out.sha256)"; echo "srctree $(cat $W/bin/ALF_$N.srctree 2>/dev/null)"; head -2 $W/bin/ALF_$N.srcstate 2>/dev/null; echo "set $S config $C setroot $SR";
  echo "env $*"; echo "seeds $(sha256sum seeds | cut -c1-16) parameters $(sha256sum parameters | cut -c1-16)"; date -Is; } >> stamp.txt
LOCK=(); GPU=0
case $N in nv_none*|gnu*) ;; *) export OMP_TARGET_OFFLOAD=MANDATORY; LOCK=(flock "$HOME/alf-gpu/gpu.lock"); GPU=1 ;; esac
OWN=1; DSLEEP=300; [ "$T207_3_LOCKED" = 1 ] && { LOCK=(); OWN=0; }
# wall time measured after the lock is taken and the GPU is free (GPU minutes = wall of the run)
dw0=$(date +%s)
while :; do
"${LOCK[@]}" bash -c '
  if [ "$1" = 1 ]; then
    for i in $(seq 1 360); do [ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" -eq 0 ] && break; sleep 10; done
    if [ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" -ne 0 ]; then
      echo 9 > .rc; date +%s.%N > .t0; cp .t0 .t1; echo "GPU busy for 60 min, not started" > run.log; exit 0
    fi
  fi
  if [ -n "$2" ]; then
    nw=$(( $(pgrep -u "$USER" -x flock -a | grep -c "alf-gpu/gpu[.]lock") - 1 )); [ $nw -lt 0 ] && nw=0
    if command -v lsof > /dev/null; then   # fd-style waiters (flock 9) of a chain; our own ancestors (the flock holding the lock) excluded
      anc=" "; q=$$; while [ "$q" -gt 1 ] 2>/dev/null; do anc="$anc$q "; q=$(ps -o ppid= -p $q | tr -d " "); done
      n2=0; for q in $(lsof -t "$HOME/alf-gpu/gpu.lock" 2>/dev/null); do
        case "$anc" in *" $q "*) continue;; esac
        [ "$(ps -o comm= -p $q 2>/dev/null)" = flock ] && n2=$((n2 + 1)); done
      [ $n2 -gt $nw ] && nw=$n2
    fi
    [ "$3" = 0 ] && nw=0     # the caller holds the lock through this run and its processing: a waiting run starts only after the thinning
    fm=$(df -B1M --output=avail / | tail -1 | tr -d " "); need=$(( $2 + 1600 * nw ))
    echo "disk check in the lock $(date -Is): / $fm MB free, need $need MB ($2 + 1600 x $nw waiting on gpu.lock)" >> stamp.txt
    if [ "$fm" -lt "$need" ]; then : > .diskwait; exit 0; fi
  fi
  date +%s.%N > .t0; timeout -k 30 7200 "$0" > run.log 2>&1; echo $? > .rc; date +%s.%N > .t1' $W/bin/ALF_$N.out $GPU "$MINFREE_MB" $OWN
  [ -f .diskwait ] || break
  if [ $(( $(date +%s) - dw0 )) -ge "${DISKWAIT_S:-28800}" ]; then
    rm -f .diskwait; echo 8 > .rc; date +%s.%N > .t0; cp .t0 .t1; echo "not started: / below the disk floor for ${DISKWAIT_S:-28800} s" > run.log; break
  fi
  # LOCKED: the caller holds the lock; leave the flag, nothing is run, the caller gives the lock up while waiting (see the header)
  [ "$T207_3_LOCKED" = 1 ] && { echo "disk below the floor: flag left in $D/.diskwait, caller to release the lock and retry"; exit 0; }
  rm -f .diskwait; sleep $DSLEEP
done
s=$(cat .t0); e=$(cat .t1); rc=$(cat .rc); rm -f .t0 .t1 .rc
printf '%s %s %s %s rc=%s wall=%.1fs bin %s env [%s]\n' "$T" $N $S $C $rc $(echo "$e - $s" | bc) "$(cut -c1-16 $W/bin/ALF_$N.out.sha256)" "$*" | tee -a stamp.txt
prov $rc $s $e

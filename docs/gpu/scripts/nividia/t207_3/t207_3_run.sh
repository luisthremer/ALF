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
#   sleeps 300 s and tries again, for at most DISKWAIT_S (default 28800) s, then gives up (rc=8, logged, nothing run).
#   T207_3_LOCKED=1: the caller (t207_3_chain.sh with HOLDLOCK=1) already holds gpu.lock on its own file descriptor through the run and
#   its processing; no second flock (it would deadlock), waiting runs are not counted (they start only after the thinning), and a
#   failed disk check is retried every 30 s inside the held lock (other agents' processing goes on and thins their dumps).
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
case $N in gnu*) . ~/alf-gpu/env/gnu.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-dev1}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
[ -x "$W/bin/ALF_$N.out" ] || { echo "refusing: binary $N missing"; exit 1; }
SR=${SETROOT:-$W/runs}; [ -d "$SR/$S/$C" ] || { echo "refusing: $SR/$S/$C missing"; exit 1; }
mkdir -p "$W/runs/$T"; cp -r "$SR/$S/$C" "$D"; cd "$D" || exit 1; mkdir -p qref   # dumps of ALF_GPU_TAUDIAG_QREF
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
{ echo "binary $(cat $W/bin/ALF_$N.out.sha256)"; echo "srctree $(cat $W/bin/ALF_$N.srctree 2>/dev/null)"; head -2 $W/bin/ALF_$N.srcstate 2>/dev/null; echo "set $S config $C setroot $SR";
  echo "env $*"; echo "seeds $(sha256sum seeds | cut -c1-16) parameters $(sha256sum parameters | cut -c1-16)"; date -Is; } > stamp.txt
LOCK=(); GPU=0
case $N in nv_none*|gnu*) ;; *) export OMP_TARGET_OFFLOAD=MANDATORY; LOCK=(flock "$HOME/alf-gpu/gpu.lock"); GPU=1 ;; esac
OWN=1; DSLEEP=300; [ "$T207_3_LOCKED" = 1 ] && { LOCK=(); OWN=0; DSLEEP=30; }
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
    [ "$3" = 0 ] && nw=0     # the caller holds the lock through this run and its processing: a waiting run starts only after the thinning
    fm=$(df -B1M --output=avail / | tail -1 | tr -d " "); need=$(( $2 + 1600 * nw ))
    echo "disk check in the lock $(date -Is): / $fm MB free, need $need MB ($2 + 1600 x $nw waiting on gpu.lock)" >> stamp.txt
    if [ "$fm" -lt "$need" ]; then echo 75 > .rc; exit 0; fi
  fi
  date +%s.%N > .t0; timeout -k 30 7200 "$0" > run.log 2>&1; echo $? > .rc; date +%s.%N > .t1' $W/bin/ALF_$N.out $GPU "$MINFREE_MB" $OWN
  [ "$(cat .rc)" = 75 ] || break
  rm -f .rc
  if [ $(( $(date +%s) - dw0 )) -ge "${DISKWAIT_S:-28800}" ]; then
    echo 8 > .rc; date +%s.%N > .t0; cp .t0 .t1; echo "not started: / below the disk floor for ${DISKWAIT_S:-28800} s" > run.log; break
  fi
  sleep $DSLEEP
done
s=$(cat .t0); e=$(cat .t1); rc=$(cat .rc); rm -f .t0 .t1 .rc
printf '%s %s %s %s rc=%s wall=%.1fs bin %s env [%s]\n' "$T" $N $S $C $rc $(echo "$e - $s" | bc) "$(cut -c1-16 $W/bin/ALF_$N.out.sha256)" "$*" | tee -a stamp.txt
~/alf-gpu/venv/bin/python3 -I - "$W" "$N" "$S" "$C" "$T" "$rc" "$s" "$e" "$*" <<'PY'
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

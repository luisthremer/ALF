#!/bin/bash
# Ticket 207-4 (development only, dev seeds; copy of t207_3_run.sh with W = t207-4 and a disk guard in the lock), nividia: run one binary on one prepared config (single process, 1 thread).
#   usage: [SET=devN] [TAG=dir] t207_4_run.sh NAME CONFIG [env assignments...]
#   -> t207-4/runs/TAG/CONFIG (refuses to overwrite): run.log, stamp.txt (binary sha256, source tree hash, env), provenance.json.
#   GPU binaries (nv_* except nv_none*) run inside flock ~/alf-gpu/gpu.lock; inside the lock the script waits (up to 60 min, polling)
#   until no compute process is on the GPU, so only one GPU process runs at a time.  No fixed run-time limit (decision 16); a 2 h guard
#   timeout only stops a hung run.  Same records as t207_d5_run.sh (the QREF pipeline reads them).
#   Disk guard (pattern of 207-3c b2c71200): inside the GPU lock, once the GPU is free, the run starts only if / has at least DISK_NEED_KB
#   (default 6500000 = 5 GB floor + one unthinned beta = 10 dump set, measured 1.40 GB, rounded to 1.5 GB) free, measured fresh; else it waits (lock released and
#   retaken every 300 s; with T207_4_LOCKED=1 the lock stays held and the check is retried every 30 s) for at most DISKWAIT_S (28800) s,
#   then gives up (rc=8, nothing run).  T207_4_LOCKED=1: the caller (t207_4_chain.sh) holds gpu.lock on its own fd through this run and
#   its processing; no second flock (it would deadlock).
W=$HOME/alf-gpu/work/t207-4; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
case $N in gnu*) . ~/alf-gpu/env/gnu.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-dev1}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
[ -x "$W/bin/ALF_$N.out" ] || { echo "refusing: binary $N missing"; exit 1; }
NEED=${DISK_NEED_KB:-6500000}
mkdir -p "$W/runs/$T"; cp -r "$W/runs/$S/$C" "$D"; cd "$D" || exit 1; mkdir -p qref   # dumps of ALF_GPU_TAUDIAG_QREF
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
{ echo "binary $(cat $W/bin/ALF_$N.out.sha256)"; echo "srctree $(cat $W/bin/ALF_$N.srctree 2>/dev/null)"; head -2 $W/bin/ALF_$N.srcstate 2>/dev/null; echo "set $S config $C";
  echo "env $*"; echo "seeds $(sha256sum seeds | cut -c1-16) parameters $(sha256sum parameters | cut -c1-16)"; date -Is; } > stamp.txt
LOCK=(); GPU=0
case $N in nv_none*|gnu*) ;; *) export OMP_TARGET_OFFLOAD=MANDATORY; LOCK=(flock "$HOME/alf-gpu/gpu.lock"); GPU=1 ;; esac
DSLEEP=300; [ "$T207_4_LOCKED" = 1 ] && { LOCK=(); DSLEEP=30; }
# wall time measured after the lock is taken, the GPU is free and the disk check passed (GPU minutes = wall of the run)
dw0=$(date +%s)
while :; do
"${LOCK[@]}" bash -c '
  if [ "$1" = 1 ]; then
    for i in $(seq 1 360); do [ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" -eq 0 ] && break; sleep 10; done
    if [ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" -ne 0 ]; then
      echo 9 > .rc; date +%s.%N > .t0; cp .t0 .t1; echo "GPU busy for 60 min, not started" > run.log; exit 0
    fi
  fi
  fk=$(df -k --output=avail / | tail -1 | tr -d " ")
  echo "disk check in the lock $(date -Is): / $fk kB free, need $2 kB" >> stamp.txt
  if [ "$fk" -lt "$2" ]; then echo 75 > .rc; exit 0; fi
  date +%s.%N > .t0; timeout -k 30 7200 "$0" > run.log 2>&1; echo $? > .rc; date +%s.%N > .t1' $W/bin/ALF_$N.out $GPU $NEED
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
json.dump(dict(ticket='207-4 (development evidence, dev seeds; not acceptance)', tag=T, build=N, seed_set=S, config=C, rc=int(rc),
               wall_s=round(float(e) - float(s), 1), binary_sha256=sha(os.path.join(W, 'bin', 'ALF_%s.out' % N)),
               srctree=rd(os.path.join(W, 'bin', 'ALF_%s.srctree' % N)), srcstate=rd(os.path.join(W, 'bin', 'ALF_%s.srcstate' % N)),
               env=env.split(), seeds_sha256=sha('seeds'), parameters_sha256=sha('parameters'), run_log_sha256=sha('run.log'),
               info_sha256=sha('info'), data_h5_sha256=sha('data.h5'), host=socket.gethostname(), gpu=gpu),
          open('provenance.json', 'w'), indent=1)
PY

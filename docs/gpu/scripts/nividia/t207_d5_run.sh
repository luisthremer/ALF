#!/bin/bash
# Ticket 207 (dev5: accuracy test against a quad reference, development only), nividia: run one binary on one prepared config
# (single process, 1 thread).   usage: [SET=devN] [TAG=dir] t207_d5_run.sh NAME CONFIG [env assignments...]
#   -> dev5/runs/TAG/CONFIG (refuses to overwrite), run.log, stamp.txt (binary sha256, source tree hash, env); prints one summary line.
#   GPU binaries (nv_* except nv_none*) run under flock ~/alf-gpu/gpu.lock after a check that no compute process is on the GPU.
#   No fixed run-time limit (decision 16); a 2 h guard timeout only stops a hung run.
W=$HOME/alf-gpu/work/t207/dev5; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
case $N in gnu*) . ~/alf-gpu/env/gnu.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-dev1}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
mkdir -p "$W/runs/$T"; cp -r "$W/runs/$S/$C" "$D"; cd "$D" || exit 1; mkdir -p qref   # dumps of ALF_GPU_TAUDIAG_QREF
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
{ echo "binary $(cat $W/bin/ALF_$N.out.sha256)"; echo "srctree $(cat $W/bin/ALF_$N.srctree 2>/dev/null)"; head -2 $W/bin/ALF_$N.srcstate 2>/dev/null; echo "set $S config $C";
  echo "env $*"; echo "seeds $(sha256sum seeds | cut -c1-16) parameters $(sha256sum parameters | cut -c1-16)"; date -Is; } > stamp.txt
LOCK=()
case $N in nv_none*|gnu*) ;; *)
  export OMP_TARGET_OFFLOAD=MANDATORY; LOCK=(flock "$HOME/alf-gpu/gpu.lock")
  busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
  [ "$busy" -eq 0 ] || { sleep 120; busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l); }
  [ "$busy" -eq 0 ] || { echo "GPU busy, not starting $N $C"; exit 9; } ;;
esac
# wall time measured after the lock is taken (GPU minutes = wall of the locked section)
"${LOCK[@]}" bash -c 'date +%s.%N > .t0; timeout -k 30 7200 "$0" > run.log 2>&1; echo $? > .rc; date +%s.%N > .t1' $W/bin/ALF_$N.out
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
json.dump(dict(ticket='207 dev5 (development evidence, dev seeds; not acceptance)', tag=T, build=N, seed_set=S, config=C, rc=int(rc),
               wall_s=round(float(e) - float(s), 1), binary_sha256=sha(os.path.join(W, 'bin', 'ALF_%s.out' % N)),
               srctree=rd(os.path.join(W, 'bin', 'ALF_%s.srctree' % N)), srcstate=rd(os.path.join(W, 'bin', 'ALF_%s.srcstate' % N)),
               env=env.split(), seeds_sha256=sha('seeds'), parameters_sha256=sha('parameters'), run_log_sha256=sha('run.log'),
               info_sha256=sha('info'), data_h5_sha256=sha('data.h5'), host=socket.gethostname(), gpu=gpu),
          open('provenance.json', 'w'), indent=1)
PY

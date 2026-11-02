#!/bin/bash
# Ticket 302c-2 (development only, dev seeds), on intel: run one binary on one prepared config (single process, 1 thread), with the
# records the QREF pipeline reads (as t207_3_run.sh / t207_d5_run.sh on nividia): run.log, stamp.txt (binary sha256, source tree hash,
# env), provenance.json.   usage: [SET=devN] [TAG=dir] t302c2_run.sh NAME CONFIG [ALF_GPU_* assignments...]
#   -> t302c-2/runs/TAG/CONFIG (refuses to overwrite).  Only ALF_GPU_* assignments are accepted as arguments (they are the run's env
#   in stamp.txt and provenance.json, which t207_d5_qparse.py compares with the program's own 'GPU debug env' line).
#   GPU builds (in_* except in_none*) get the mandatory B580 environment (ONEAPI_DEVICE_SELECTOR=level_zero:0,
#   OMP_TARGET_OFFLOAD=MANDATORY, FORT_BUFFERED=false) and run inside flock ~/alf-gpu/gpu.lock.  Load check before every run
#   (CLAUDE.md, intel is shared): uptime and the top processes are logged; if a process of another user uses > 50 % CPU the script
#   waits (polling every 60 s, up to 60 min) and then refuses.  No fixed run-time limit (decision 16); a 2 h guard timeout only stops a
#   hung run.
W=${T302C2_W:-$HOME/alf-gpu/work/t302c-2}   # round 2: T302C2_W=$HOME/alf-gpu/work/t302c-2b (only these two)
case "$W" in "$HOME/alf-gpu/work/t302c-2"|"$HOME/alf-gpu/work/t302c-2b") ;; *) echo "refusing work dir $W"; exit 1;; esac
:; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
. ~/alf-gpu/env/intel.sh 2>/dev/null
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-dev1}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
for kv in "$@"; do case "$kv" in ALF_GPU_*=*) ;; *) echo "refusing: argument $kv is not an ALF_GPU_* assignment"; exit 1;; esac; done
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
[ -x "$W/bin/ALF_$N.out" ] || { echo "refusing: binary $N missing"; exit 1; }
[ "$(sha256sum "$W/bin/ALF_$N.out" | cut -d' ' -f1)" = "$(cut -d' ' -f1 "$W/bin/ALF_$N.out.sha256")" ] || { echo "refusing: binary $N differs from its .sha256"; exit 1; }
busy() { ps -eo user=,pcpu=,comm= | awk -v me="$(id -un)" '$1 != me && $1 != "root" && $2 > 50 {print}'; }
for i in $(seq 1 60); do [ -z "$(busy)" ] && break; echo "load check: other users busy, waiting: $(busy | head -3 | tr '\n' ';')"; sleep 60; done
[ -z "$(busy)" ] || { echo "refusing: other users busy for 60 min: $(busy | head -3 | tr '\n' ';')"; exit 9; }
mkdir -p "$W/runs/$T"; cp -r "$W/runs/$S/$C" "$D"; cd "$D" || exit 1; mkdir -p qref   # dumps of ALF_GPU_TAUDIAG_QREF
{ uptime; ps -eo user,pcpu,comm --sort=-pcpu | head -6; } > loadcheck.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
{ echo "binary $(cat $W/bin/ALF_$N.out.sha256)"; echo "srctree $(cat $W/bin/ALF_$N.srctree 2>/dev/null)"; head -2 $W/bin/ALF_$N.srcstate 2>/dev/null; echo "set $S config $C";
  echo "env $*"; echo "seeds $(sha256sum seeds | cut -c1-16) parameters $(sha256sum parameters | cut -c1-16)"; date -Is; } > stamp.txt
LOCK=()
case $N in in_none*|gnu*) ;; *) export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY FORT_BUFFERED=false; LOCK=(flock "$HOME/alf-gpu/gpu.lock") ;; esac
tw=$(date +%s.%N)   # round 2: time before the lock, so a wait for ~/alf-gpu/gpu.lock (another agent's GPU run) is recorded (lockwait)
"${LOCK[@]}" bash -c 'date +%s.%N > .t0; timeout -k 30 7200 "$0" > run.log 2>&1; echo $? > .rc; date +%s.%N > .t1' $W/bin/ALF_$N.out
s=$(cat .t0); e=$(cat .t1); rc=$(cat .rc); rm -f .t0 .t1 .rc
printf '%s %s %s %s rc=%s wall=%.1fs lockwait=%.1fs bin %s env [%s]\n' "$T" $N $S $C $rc $(echo "$e - $s" | bc) $(echo "$s - $tw" | bc) "$(cut -c1-16 $W/bin/ALF_$N.out.sha256)" "$*" | tee -a stamp.txt
~/alf-gpu/venv/bin/python3 -I - "$W" "$N" "$S" "$C" "$T" "$rc" "$s" "$e" "$*" "$tw" <<'PY'
import hashlib, json, os, socket, subprocess, sys
W, N, S, C, T, rc, s, e, env, tw = sys.argv[1:11]
def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else None
def rd(p):
    return open(p).read().strip() if os.path.exists(p) else None
pk = subprocess.run(['dpkg-query', '-W', 'intel-opencl-icd', 'libze-intel-gpu1', 'libze1', 'libigc2'], capture_output=True, text=True).stdout.split('\n')
gpu = 'Intel Arc B580 (ONEAPI_DEVICE_SELECTOR=%s); %s' % (os.environ.get('ONEAPI_DEVICE_SELECTOR', 'unset: CPU build'), '; '.join(x for x in pk if x))
json.dump(dict(ticket='302c-2 (development evidence, dev seeds; not acceptance)', tag=T, build=N, seed_set=S, config=C, rc=int(rc),
               wall_s=round(float(e) - float(s), 1), lock_wait_s=round(float(s) - float(tw), 1), binary_sha256=sha(os.path.join(W, 'bin', 'ALF_%s.out' % N)),
               srctree=rd(os.path.join(W, 'bin', 'ALF_%s.srctree' % N)), srcstate=rd(os.path.join(W, 'bin', 'ALF_%s.srcstate' % N)),
               env=env.split(), seeds_sha256=sha('seeds'), parameters_sha256=sha('parameters'), run_log_sha256=sha('run.log'),
               info_sha256=sha('info'), data_h5_sha256=sha('data.h5'), host=socket.gethostname(), gpu=gpu),
          open('provenance.json', 'w'), indent=1)
PY

#!/bin/bash
# Ticket 401 (development): run one binary on one prepared config (single process, 1 thread), on nividia or intel.
#   usage: [TAG=dir] run.sh BIN CONFIG [env assignments...]   -> runs/TAG(default BIN)/CONFIG with run.log; prints rc and wall time
# GPU runs on nividia take ~/alf-gpu/gpu.lock (flock) when GPULOCK=1.  On intel: B580 only (level_zero:0).
W=$HOME/alf-gpu/work/${T401:-t401b}; export TMPDIR=$W/tmp
B=$1; C=$2; shift 2
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$B}; safe "$B"; safe "$C"; safe "$T"; : "${HOME:?}"
if [ -d /opt/nvidia ]; then . ~/alf-gpu/env/nvhpc.sh; else . ~/alf-gpu/env/intel.sh; export ONEAPI_DEVICE_SELECTOR=level_zero:0; fi
D="$W/runs/$T/$C"; rm -rf -- "$D"; mkdir -p "$W/runs/$T"; cp -r "$W/runs/base/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_TARGET_OFFLOAD=MANDATORY
case $B in cpu*) unset OMP_TARGET_OFFLOAD;; esac
for kv in "$@"; do export "$kv"; done
s=$(date +%s.%N)
if [ "$GPULOCK" = 1 ]; then flock ~/alf-gpu/gpu.lock $W/bin/ALF_$B.out > run.log 2>&1; rc=$?; else $W/bin/ALF_$B.out > run.log 2>&1; rc=$?; fi
e=$(date +%s.%N)
# 401-A: provenance stamp of this run (binary and source-tree hashes, host, env assignments, rc)
python3 - "$W" "$B" "$C" "$rc" "$@" > provenance.json <<'PYEOF'
import hashlib, json, os, socket, subprocess, sys, time
W, B, C, rc = sys.argv[1:5]
def h(p):
    d = hashlib.sha256()
    with open(p, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            d.update(blk)
    return d.hexdigest()
src = os.path.join(W, 'src'); files = []
for root, _, fs in os.walk(src):
    for fn in fs:
        if fn.endswith(('.F90', '.f90', '.sh', '.py', '.json', '.mk')) or fn in ('Makefile', 'configure.sh'):
            files.append(os.path.join(root, fn))
d = hashlib.sha256()
for p in sorted(files):
    d.update(os.path.relpath(p, src).encode()); d.update(h(p).encode())
print(json.dumps({"binary": B, "binary_sha256": h(os.path.join(W, 'bin', 'ALF_%s.out' % B)), "config": C, "rc": int(rc),
                  "src_tree_sha256": d.hexdigest(), "src_files": len(files), "host": socket.gethostname(),
                  "time": time.strftime('%Y-%m-%dT%H:%M:%S'), "env": sys.argv[5:], "cwd": os.getcwd()}, indent=1))
PYEOF
printf '%s %s rc=%s wall=%.1fs\n' $B $C $rc $(echo "$e - $s" | bc)

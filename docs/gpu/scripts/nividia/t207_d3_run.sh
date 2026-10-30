#!/bin/bash
# Ticket 207 (dev3: diagnosis of the default-path paired spread, development only), nividia: run one binary on one prepared config
# (single process, 1 thread).   usage: [SET=devN] [TAG=dir] t207_d3_run.sh NAME CONFIG [env assignments...]
#   -> dev3/runs/TAG/CONFIG (refuses to overwrite), run.log, stamp.txt (binary sha256, source tree hash, env); prints one summary line.
#   GPU binaries (nv_* except nv_none*) run under flock ~/alf-gpu/gpu.lock after a check that no compute process is on the GPU.
#   No fixed run-time limit (decision 16); a 2 h guard timeout only stops a hung run.
W=$HOME/alf-gpu/work/t207/dev3; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
case $N in gnu*) . ~/alf-gpu/env/gnu.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-dev1}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
mkdir -p "$W/runs/$T"; cp -r "$W/runs/$S/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
{ echo "binary $(cat $W/bin/ALF_$N.out.sha256)"; echo "srctree $(cat $W/bin/ALF_$N.srctree 2>/dev/null)"; echo "set $S config $C";
  echo "env $*"; echo "seeds $(sha256sum seeds | cut -c1-16) parameters $(sha256sum parameters | cut -c1-16)"; date -Is; } > stamp.txt
LOCK=()
case $N in nv_none*|gnu*) ;; *)
  export OMP_TARGET_OFFLOAD=MANDATORY; LOCK=(flock "$HOME/alf-gpu/gpu.lock")
  busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
  [ "$busy" -eq 0 ] || { sleep 120; busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l); }
  [ "$busy" -eq 0 ] || { echo "GPU busy, not starting $N $C"; exit 9; } ;;
esac
s=$(date +%s.%N); "${LOCK[@]}" timeout -k 30 7200 $W/bin/ALF_$N.out > run.log 2>&1; rc=$?
e=$(date +%s.%N)
printf '%s %s %s %s rc=%s wall=%.1fs bin %s env [%s]\n' "$T" $N $S $C $rc $(echo "$e - $s" | bc) "$(cut -c1-16 $W/bin/ALF_$N.out.sha256)" "$*" | tee -a stamp.txt

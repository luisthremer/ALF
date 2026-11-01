#!/bin/bash
# Ticket 207 (rebase on 212, development only), nividia: run one binary on one prepared config (single process, 1 thread).
#   usage: [SET=base|devN] [TAG=dir] t207_d2_run.sh NAME CONFIG [env assignments...]
#   -> dev2/runs/TAG(default NAME)/CONFIG (refuses to overwrite), run.log; prints rc, wall time, binary sha256 prefix, env.
#   GPU binaries (nv_* names) run under flock ~/alf-gpu/gpu.lock after a check that no compute process is on the GPU.
W=$HOME/alf-gpu/work/t207/dev2; export TMPDIR=$W/tmp
N=$1; C=$2; shift 2
case $N in gnu*|base_gnu*) . ~/alf-gpu/env/gnu.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
safe() { case "$1" in ""|*..*|*/*) echo "refusing unsafe name: [$1]" >&2; exit 1;; esac; }
T=${TAG:-$N}; S=${SET:-base}; safe "$N"; safe "$C"; safe "$T"; safe "$S"; : "${HOME:?}"
D="$W/runs/$T/$C"; [ -e "$D" ] && { echo "refusing: $D exists"; exit 1; }
mkdir -p "$W/runs/$T"; cp -r "$W/runs/$S/$C" "$D"; cd "$D" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
LOCK=()
case $N in nv_none*|gnu*|base_gnu*) ;; *)
  export OMP_TARGET_OFFLOAD=MANDATORY; LOCK=(flock "$HOME/alf-gpu/gpu.lock")
  busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
  [ "$busy" -eq 0 ] || { sleep 60; busy=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l); }
  [ "$busy" -eq 0 ] || { echo "GPU busy, not starting $N $C"; exit 9; } ;;
esac
s=$(date +%s.%N); "${LOCK[@]}" timeout -k 30 900 $W/bin/ALF_$N.out > run.log 2>&1; rc=$?
e=$(date +%s.%N)
printf '%s %s %s rc=%s wall=%.1fs bin %s env [%s]\n' "$T" $N $C $rc $(echo "$e - $s" | bc) "$(cut -c1-16 $W/bin/ALF_$N.out.sha256)" "$*"

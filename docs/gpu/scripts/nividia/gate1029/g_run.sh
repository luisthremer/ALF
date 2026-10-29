#!/bin/bash
# Gate 2026-10-29, nividia: one stamped run under `env -i` with an explicit whitelist.
#   usage: g_run.sh RUNDIR BIN CONFIG [VAR=value ...]
#     RUNDIR must not exist (never overwrite evidence); CONFIG is copied from ~/alf-gpu/work/t207/runs/base/CONFIG.
#     Only the listed VAR=value are added to the environment; ALF_GPU_* variables must be listed explicitly.
#   Writes RUNDIR/provenance.txt (source sha, binary sha256, toolchain, GPU, full env, input hashes, times, exit code).
#   Exit code: the program's exit code (124 on timeout); 2 on a usage or setup error.
set -o pipefail
D=$1; BIN=$2; C=$3; shift 3
[ -e "$D" ] && { echo "refusing: $D exists"; exit 2; }
[ -x "$BIN" ] && [ -f "$BIN.build.txt" ] || { echo "no stamped binary $BIN"; exit 2; }
SRCCFG=$HOME/alf-gpu/work/t207/runs/base/$C; [ -d "$SRCCFG" ] || { echo "no config $C"; exit 2; }
mkdir -p "$(dirname "$D")" && cp -r "$SRCCFG" "$D" || exit 2
ENVX=(); for kv in "$@"; do case "$kv" in *=*) ENVX+=("$kv");; *) echo "bad env item $kv"; exit 2;; esac; done
P=$D/provenance.txt
{ echo "run_dir $D"; echo "config $C"; cat "$BIN.build.txt"; echo "binary $BIN"
  echo "binary_sha256_now $(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "gpu $(nvidia-smi --query-gpu=name,uuid,driver_version,compute_mode --format=csv,noheader)"
  echo "gpu_apps_before $(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)"
  echo "parameters_sha256 $(sha256sum "$D/parameters" | cut -d' ' -f1)"; echo "seeds_sha256 $(sha256sum "$D/seeds" | cut -d' ' -f1)"
  echo "extra_env ${ENVX[*]:-none}"; echo "start $(date -Is)"; } > "$P"
cd "$D" || exit 2
env -i HOME="$HOME" USER="$USER" LANG=C PATH=/usr/bin:/bin "${ENVX[@]}" bash --noprofile --norc -c '
  . "$HOME/alf-gpu/env/nvhpc.sh" > /dev/null 2>&1
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  env | sort > env.txt
  exec timeout 900 "$0" > run.log 2>&1' "$BIN"
rc=$?
{ echo "end $(date -Is)"; echo "exit_code $rc"
  echo "gpu_apps_after $(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)"
  echo "env_sha256 $(sha256sum env.txt | cut -d' ' -f1)"
  echo "alf_gpu_env $(grep '^ALF_GPU_' env.txt | tr '\n' ' ')"
  for f in info data.h5 run.log; do [ -f $f ] && echo "out_sha256 $f $(sha256sum $f | cut -d' ' -f1)"; done; } >> "$P"
exit $rc

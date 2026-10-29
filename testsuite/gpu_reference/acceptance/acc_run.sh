#!/bin/bash
# Acceptance harness, MACHINE side: one run of a stamped binary under `env -i` with an explicit whitelist.
#   usage: acc_run.sh RUNDIR BINARY CONFIGDIR [--allow-debug] [VAR=value ...]
#     RUNDIR    must not exist (created atomically; evidence is never overwritten);
#     CONFIGDIR a prepared ALF run directory (parameters, seeds), copied into RUNDIR;
#     VAR=value extra variables. Every ALF_GPU_* / ALFSH_* variable needs --allow-debug (given by the driver only for
#               runs the acceptance file declares allow_debug); other variables only from the whitelist below.
#   Environment: env -i + the build tool's env file + OMP/BLAS threads = 1 + TMPDIR (inside the run's tmp/) + for
#   GPU builds OMP_TARGET_OFFLOAD=MANDATORY (no silent host fallback). Timeout ACC_TIMEOUT (default 900 s), then
#   SIGKILL after 30 s; rc 124 or 137 = timeout. With ACC_GPU_LOCK set, GPU builds run under flock on that file.
#   Writes RUNDIR/provenance.json.
# Exit code: the program's exit code; 2 on a usage or setup error (nothing run).
set -o pipefail
export LC_ALL=C
D=$1; BIN=$2; CFG=$3; shift 3
ALLOW=0; [ "${1:-}" = --allow-debug ] && { ALLOW=1; shift; }
BJ=${BIN%.out}.build.json; [ -x "$BIN" ] && [ -f "$BJ" ] || { echo "no stamped binary $BIN"; exit 2; }
[ -f "$CFG/parameters" ] && [ -f "$CFG/seeds" ] || { echo "no prepared config in $CFG"; exit 2; }
ENVX=()
for kv in "$@"; do
  case "$kv" in
    ALF_GPU_*=*|ALFSH_*=*) [ $ALLOW = 1 ] || { echo "refusing debug/fault variable without --allow-debug: $kv"; exit 2; }; ENVX+=("$kv");;
    CUDA_MPS_PIPE_DIRECTORY=*|CUDA_MPS_LOG_DIRECTORY=*|CUBLAS_WORKSPACE_CONFIG=*|ONEAPI_DEVICE_SELECTOR=*|ZE_AFFINITY_MASK=*) ENVX+=("$kv");;
    *=*) echo "refusing variable outside the whitelist: ${kv%%=*}"; exit 2;;
    *) echo "bad argument $kv"; exit 2;;
  esac
done
read -r TOOL GPUB < <(python3 -c "import json,sys; b=json.load(open(sys.argv[1])); print(b['tool'], 1 if ('GPU=nvidia' in b['configure'] or 'GPU=intel' in b['configure']) else 0)" "$BJ") || { echo "bad build record"; exit 2; }
ENVF=$HOME/alf-gpu/env/$TOOL.sh; [ -f "$ENVF" ] || { echo "no env file $ENVF"; exit 2; }
mkdir -p "$(dirname "$D")" && mkdir "$D" || { echo "refusing: $D exists or cannot be created"; exit 2; }
cp -r "$CFG/." "$D/" || exit 2
mkdir -p "$D/tmp"
cd "$D" || exit 2
GPUCK=-1
if command -v nvidia-smi > /dev/null; then
  nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader > gpu_before.txt 2>&1; GPUCK=$?
else : > gpu_before.txt; fi
OFF=""; [ "$GPUB" = 1 ] && OFF="OMP_TARGET_OFFLOAD=MANDATORY"
T=${ACC_TIMEOUT:-900}
LOCK=(); [ "$GPUB" = 1 ] && [ -n "${ACC_GPU_LOCK:-}" ] && LOCK=(flock "$ACC_GPU_LOCK")
date -Is > start.txt
"${LOCK[@]}" env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$D/tmp" $OFF "${ENVX[@]}" \
  bash --noprofile --norc -c '
  . "$1" > env_file.log 2>&1 || { echo "env file failed: $1" >&2; exit 2; }
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  env | sort > env.txt
  exec timeout -k 30 "$2" "$0" > run.log 2>&1' "$BIN" "$ENVF" "$T"
rc=$?
date -Is > end.txt
GPUCK2=-1
if command -v nvidia-smi > /dev/null; then
  nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader > gpu_after.txt 2>&1; GPUCK2=$?
else : > gpu_after.txt; fi
python3 - "$D" "$BIN" "$BJ" "$rc" "$ALLOW" "$GPUCK" "$GPUCK2" "$ENVF" "${ENVX[@]}" <<'PY'
import json, sys, os, hashlib
d, b, bj, rc, allow, ck, ck2, envf = sys.argv[1:9]; extra = sys.argv[9:]
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.isfile(p) else None
env = dict(l.rstrip("\n").split("=", 1) for l in open(os.path.join(d, "env.txt")) if "=" in l) if os.path.isfile(os.path.join(d, "env.txt")) else {}
ds = {}
try:
    import h5py
    if os.path.isfile(os.path.join(d, "data.h5")):
        with h5py.File(os.path.join(d, "data.h5"), "r") as f:
            f.visititems(lambda n, o: ds.__setitem__(n, hashlib.sha256(o[()].tobytes()).hexdigest())
                         if isinstance(o, h5py.Dataset) else None)
except Exception as e:
    ds = {"_error": f"{type(e).__name__}: {e}"}
rc = int(rc)
p = {"schema": "alf-acc-prov/1", "run_dir": d, "binary": b, "binary_sha256": h(b), "build": json.load(open(bj)),
     "inputs": {"parameters_sha256": h(os.path.join(d, "parameters")), "seeds_sha256": h(os.path.join(d, "seeds"))},
     "extra_env": extra, "allow_debug": allow == "1", "env_file": envf, "env_file_sha256": h(envf),
     "env": env, "debug_env_present": sorted(k for k in env if k.startswith(("ALF_GPU_", "ALFSH_"))),
     "gpu_check_rc": int(ck), "gpu_check_after_rc": int(ck2),
     "gpu_before": open(os.path.join(d, "gpu_before.txt")).read().splitlines(),
     "gpu_after": open(os.path.join(d, "gpu_after.txt")).read().splitlines(),
     "start": open(os.path.join(d, "start.txt")).read().strip(), "end": open(os.path.join(d, "end.txt")).read().strip(),
     "exit_code": rc, "timed_out": rc in (124, 137),
     "outputs": {f: h(os.path.join(d, f)) for f in ("info", "data.h5", "run.log")}, "datasets_sha256": ds}
json.dump(p, open(os.path.join(d, "provenance.json"), "w"), indent=1)
PY
exit $rc

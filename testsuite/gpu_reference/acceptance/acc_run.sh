#!/bin/bash
# Acceptance harness, MACHINE side: one run of a stamped binary under `env -i` with an explicit whitelist.
#   usage: acc_run.sh RUNDIR BINARY CONFIGDIR [--allow-debug] [VAR=value ...]
#     RUNDIR    must not exist (evidence is never overwritten);
#     CONFIGDIR a prepared ALF run directory (parameters, seeds), copied into RUNDIR;
#     VAR=value the only extra variables in the environment (plus the env file of the build's tool).
#   Debug and fault variables (ALF_GPU_FAULT*, _TAUSWAP, _TAU_NS, _NS, _STABCHECK, _STABMIX, _TAUDIAG*, _NOUPDATE,
#   _SELFTEST, _UPDTEST, _EXIT, _TAUTIME, _STABTIME, _DISABLE) are refused unless --allow-debug is given; the acceptance
#   driver passes it only for runs that the acceptance file declares as self-test or control runs.
#   Writes RUNDIR/provenance.json; timeout 900 s (exit 124).
# Exit code: the program's exit code; 2 on a usage or setup error (nothing run).
set -o pipefail
D=$1; BIN=$2; CFG=$3; shift 3
ALLOW=0; [ "${1:-}" = --allow-debug ] && { ALLOW=1; shift; }
[ -e "$D" ] && { echo "refusing: $D exists"; exit 2; }
BJ=${BIN%.out}.build.json; [ -x "$BIN" ] && [ -f "$BJ" ] || { echo "no stamped binary $BIN"; exit 2; }
[ -f "$CFG/parameters" ] && [ -f "$CFG/seeds" ] || { echo "no prepared config in $CFG"; exit 2; }
ENVX=()
for kv in "$@"; do
  case "$kv" in
    ALF_GPU_FAULT*=*|ALF_GPU_TAUSWAP=*|ALF_GPU_TAU_NS=*|ALF_GPU_NS=*|ALF_GPU_STABCHECK=*|ALF_GPU_STABMIX=*|ALF_GPU_TAUDIAG*=*|\
    ALF_GPU_NOUPDATE=*|ALF_GPU_SELFTEST=*|ALF_GPU_UPDTEST=*|ALF_GPU_EXIT=*|ALF_GPU_TAUTIME=*|ALF_GPU_STABTIME=*|ALF_GPU_DISABLE=*)
      [ $ALLOW = 1 ] || { echo "refusing debug/fault variable without --allow-debug: $kv"; exit 2; }; ENVX+=("$kv");;
    *=*) ENVX+=("$kv");;
    *) echo "bad argument $kv"; exit 2;;
  esac
done
TOOL=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['tool'])" "$BJ") || exit 2
mkdir -p "$(dirname "$D")" && cp -r "$CFG" "$D" || exit 2
cd "$D" || exit 2
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader > gpu_before.txt 2>/dev/null
date -Is > start.txt
env -i HOME="$HOME" USER="$USER" LANG=C PATH=/usr/bin:/bin "${ENVX[@]}" bash --noprofile --norc -c '
  . "$HOME/alf-gpu/env/$1.sh" > /dev/null 2>&1
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  env | sort > env.txt
  exec timeout 900 "$0" > run.log 2>&1' "$BIN" "$TOOL"
rc=$?
date -Is > end.txt
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader > gpu_after.txt 2>/dev/null
python3 - "$D" "$BIN" "$BJ" "$rc" "$ALLOW" "${ENVX[@]}" <<'PY'
import json, sys, os, hashlib, glob
d, b, bj, rc, allow = sys.argv[1:6]; extra = sys.argv[6:]
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.isfile(p) else None
env = dict(l.rstrip("\n").split("=", 1) for l in open(os.path.join(d, "env.txt")) if "=" in l)
out = {}
try:
    import h5py
    if os.path.isfile(os.path.join(d, "data.h5")):
        with h5py.File(os.path.join(d, "data.h5"), "r") as f:
            f.visititems(lambda n, o: out.__setitem__(n, hashlib.sha256(o[()].tobytes()).hexdigest())
                         if isinstance(o, h5py.Dataset) else None)
except ImportError:
    out = {"_error": "h5py not available"}
p = {"schema": "alf-acc-prov/1", "run_dir": d, "binary": b, "binary_sha256": h(b), "build": json.load(open(bj)),
     "inputs": {"parameters_sha256": h(os.path.join(d, "parameters")), "seeds_sha256": h(os.path.join(d, "seeds"))},
     "extra_env": extra, "allow_debug": allow == "1",
     "env": {k: v for k, v in env.items() if not k.startswith(("LS_COLORS",))},
     "debug_env_present": sorted(k for k in env if k.startswith("ALF_GPU_")),
     "gpu_before": open(os.path.join(d, "gpu_before.txt")).read().splitlines(),
     "gpu_after": open(os.path.join(d, "gpu_after.txt")).read().splitlines(),
     "start": open(os.path.join(d, "start.txt")).read().strip(), "end": open(os.path.join(d, "end.txt")).read().strip(),
     "exit_code": int(rc), "timed_out": int(rc) == 124,
     "outputs": {f: h(os.path.join(d, f)) for f in ("info", "data.h5", "run.log")}, "datasets_sha256": out}
json.dump(p, open(os.path.join(d, "provenance.json"), "w"), indent=1)
PY
exit $rc

#!/bin/bash
# Acceptance harness, MACHINE side: build one variant of a shipped candidate.
#   usage: acc_build.sh TICKET SHA12 NAME TOOL [configure args...]     (ALF_FLAGS_EXT from the caller is recorded)
#   SRC  = ~/alf-gpu/work/t<TICKET>/cand/<SHA12>/src (or base_src if NAME starts with "base_")
#   TOOL = gnu | nvhpc | intel (the env file in ~/alf-gpu/env/), sourced inside env -i (no caller environment)
#   ALF_HDF5_DIR (default ~/alf-gpu/hdf5) must hold a prebuilt HDF5 for the compiler (never downloaded here).
#   The tarball of the source is re-hashed against .alf_candidate before building.
#   Result: cand/<SHA12>/bin/<NAME>.out (read-only) and <NAME>.build.json. Refuses to overwrite an existing binary.
# Exit code: 0 only if every check, configure and make succeed.
set -o pipefail
export LC_ALL=C
TICKET=$1; S12=$2; NAME=$3; TOOL=$4; shift 4
case "$TICKET$S12$NAME" in *..*|*/*|"") echo "unsafe name"; exit 2;; esac
C=$HOME/alf-gpu/work/t$TICKET/cand/$S12; [ -f $C/.alf_candidate ] || { echo "no shipped candidate $C"; exit 2; }
SRC=src; case "$NAME" in base_*) SRC=base_src;; esac
[ -d "$C/$SRC" ] || { echo "no source $C/$SRC"; exit 2; }
key=tar_sha256; [ $SRC = base_src ] && key=base_tar_sha256
want=$(awk -v k=$key '$1==k {print $2}' $C/.alf_candidate); have=$(sha256sum $C/$SRC.tar | cut -d' ' -f1)
[ -n "$want" ] && [ "$want" = "$have" ] || { echo "TAR HASH MISMATCH for $SRC ($have vs $want)"; exit 2; }
ENVF=$HOME/alf-gpu/env/$TOOL.sh; case $TOOL in gnu|nvhpc|intel) ;; *) echo "bad tool"; exit 2;; esac
[ -f $ENVF ] || { echo "no env file $ENVF"; exit 2; }
H5=${ALF_HDF5_DIR:-$HOME/alf-gpu/hdf5}; [ -n "$(ls -A "$H5" 2>/dev/null)" ] || { echo "no prebuilt HDF5 in $H5 (see README)"; exit 2; }
export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR $C/bin $C/logs
B=$C/bin/$NAME.out; [ -e "$B" ] && { echo "refusing: $B exists (never overwrite a stamped binary)"; exit 2; }
T=$C/tree_$NAME; mkdir "$T" || { echo "refusing: $T exists"; exit 2; }
cp -a "$C/$SRC/." "$T/" && chmod -R u+w "$T" || exit 2
env -i HOME="$HOME" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$TMPDIR" ALF_HDF5_DIR="$H5" \
    ALF_FLAGS_EXT="${ALF_FLAGS_EXT:-}" LOGC=$C/logs/conf_$NAME.log LOGE=$C/logs/env_$NAME.txt LOGM=$C/logs/make_$NAME.log \
    bash --noprofile --norc -c '
  . "$1" > /dev/null 2>&1 || { echo "env file failed"; exit 2; }; shift
  cd "$0" || exit 2
  . ./configure.sh "$@" NO-INTERACTIVE > "$LOGC" 2>&1 || { echo "configure failed: $LOGC"; exit 3; }
  env | grep "^ALF_" | sort > "$LOGE"
  make -j4 program > "$LOGM" 2>&1 || { echo "make failed: $LOGM"; exit 4; }' "$T" "$ENVF" "$@"
rc=$?; [ $rc = 0 ] || exit $rc
cp "$T/Prog/ALF.out" "$B" && chmod a-w "$B" || exit 5
python3 - "$C" "$NAME" "$B" "$TOOL" "$ENVF" "$SRC" "$*" <<'PY'
import json, sys, os, hashlib, subprocess, datetime
c, name, b, tool, envf, src, conf = sys.argv[1:8]
cand = dict(l.split(None, 1) for l in open(os.path.join(c, ".alf_candidate")).read().splitlines() if " " in l)
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
alf = dict(l.rstrip("\n").split("=", 1) for l in open(os.path.join(c, "logs", f"env_{name}.txt")) if "=" in l)
sh = lambda cmd: subprocess.run(["bash", "-c", f". {envf} >/dev/null 2>&1; {cmd}"], capture_output=True, text=True).stdout.strip()
d = {"schema": "alf-acc-build/1", "name": name, "binary": b, "binary_sha256": h(b), "source": src,
     "candidate": cand, "tool": tool, "configure": conf, "ALF_FLAGS_EXT": os.environ.get("ALF_FLAGS_EXT", ""),
     "ALF_FC": alf.get("ALF_FC", ""), "ALF_FLAGS_PROG": alf.get("ALF_FLAGS_PROG", ""), "ALF_LIB": alf.get("ALF_LIB", ""),
     "env_file_sha256": h(envf), "make_log_sha256": h(os.path.join(c, "logs", f"make_{name}.log")),
     "toolchain": {"nvfortran": sh("nvfortran --version 2>/dev/null | sed -n 2p"), "gfortran": sh("gfortran --version 2>/dev/null | head -1"),
                   "ifx": sh("ifx --version 2>/dev/null | head -1"), "nvcc": sh("nvcc --version 2>/dev/null | tail -1"),
                   "driver": sh("nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null")},
     "built": datetime.datetime.now().astimezone().isoformat()}
json.dump(d, open(b[:-4] + ".build.json", "w"), indent=1)
print("built", name, d["binary_sha256"][:16])
PY

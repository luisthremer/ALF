#!/bin/bash
# Acceptance harness, MACHINE side: build one variant of a shipped candidate.
#   usage: acc_build.sh TICKET SHA12 NAME TOOL [configure args...]     (ALF_FLAGS_EXT from the caller is recorded)
#   SRC  = ~/alf-gpu/work/t<TICKET>/cand/<SHA12>/src (or base_src if NAME starts with "base_")
#   TOOL = gnu | nvhpc | intel (the env file in ~/alf-gpu/env/)
#   Result: cand/<SHA12>/bin/<NAME>.out and <NAME>.build.json (provenance). Refuses to overwrite an existing binary.
# Exit code: 0 only if configure and make succeed.
set -o pipefail
TICKET=$1; S12=$2; NAME=$3; TOOL=$4; shift 4
case "$TICKET$S12$NAME" in *..*|*/*|"") echo "unsafe name"; exit 2;; esac
C=$HOME/alf-gpu/work/t$TICKET/cand/$S12; [ -f $C/.alf_candidate ] || { echo "no shipped candidate $C"; exit 2; }
SRC=$C/src; case "$NAME" in base_*) SRC=$C/base_src;; esac
[ -d "$SRC" ] || { echo "no source $SRC"; exit 2; }
export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR $C/bin $C/logs
case $TOOL in gnu|nvhpc|intel) . ~/alf-gpu/env/$TOOL.sh > /dev/null 2>&1;; *) echo "bad tool"; exit 2;; esac
[ -n "${ALF_HDF5_DIR:-}" ] || export ALF_HDF5_DIR=$HOME/alf-gpu/hdf5   # prebuilt HDF5 per compiler (read only)
[ -d "$ALF_HDF5_DIR" ] || { echo "no prebuilt HDF5 in $ALF_HDF5_DIR (see README)"; exit 2; }
B=$C/bin/$NAME.out; [ -e "$B" ] && { echo "refusing: $B exists (never overwrite a stamped binary)"; exit 2; }
T=$C/tree_$NAME; [ -e "$T" ] && { echo "refusing: $T exists"; exit 2; }
cp -a "$SRC/." "$T/" 2>/dev/null || { mkdir -p "$T" && cp -a "$SRC/." "$T/"; }; chmod -R u+w "$T"; cd "$T" || exit 2
. ./configure.sh "$@" NO-INTERACTIVE > $C/logs/conf_$NAME.log 2>&1 || { echo "configure failed: $C/logs/conf_$NAME.log"; exit 3; }
make -j4 program > $C/logs/make_$NAME.log 2>&1 || { echo "make failed: $C/logs/make_$NAME.log"; exit 4; }
cp Prog/ALF.out "$B" && chmod a-w "$B" || exit 5
python3 - "$C" "$NAME" "$B" "$TOOL" "$*" <<'PY'
import json, sys, os, hashlib, subprocess, datetime
c, name, b, tool, conf = sys.argv[1:6]
cand = dict(l.split(None, 1) for l in open(os.path.join(c, ".alf_candidate")).read().splitlines() if " " in l)
sh = lambda cmd: subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
d = {"schema": "alf-acc-build/1", "name": name, "binary": b, "binary_sha256": h(b),
     "source": "base_src" if name.startswith("base_") else "src",
     "candidate": cand, "tool": tool, "configure": conf, "ALF_FLAGS_EXT": os.environ.get("ALF_FLAGS_EXT", ""),
     "ALF_FC": os.environ.get("ALF_FC", ""), "ALF_FLAGS_PROG": os.environ.get("ALF_FLAGS_PROG", ""),
     "ALF_LIB": os.environ.get("ALF_LIB", ""),
     "make_log_sha256": h(os.path.join(c, "logs", f"make_{name}.log")),
     "toolchain": {"nvfortran": sh("nvfortran --version 2>/dev/null | sed -n 2p"), "gfortran": sh("gfortran --version 2>/dev/null | head -1"),
                   "ifx": sh("ifx --version 2>/dev/null | head -1"), "nvcc": sh("nvcc --version 2>/dev/null | tail -1"),
                   "driver": sh("nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null")},
     "built": datetime.datetime.now().astimezone().isoformat()}
json.dump(d, open(b.replace(".out", ".build.json"), "w"), indent=1)
print("built", name, d["binary_sha256"][:16])
PY

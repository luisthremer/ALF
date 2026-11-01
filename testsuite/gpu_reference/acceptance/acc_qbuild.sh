#!/bin/bash
# Acceptance harness, MACHINE side (ticket 207-5): build the offline quad-reference program (docs/gpu/scripts/nividia/t207_d5_qref.F90) and
# the GNU CPU-oracle objects it links, from the hash-checked tarball of a shipped candidate, with provenance.
#   usage: [ALF_FLAGS_EXT=-DQR_PREPIVOT] acc_qbuild.sh TICKET WHO SHA12 def|log "CONFIGURE" QFLAGS
#     e.g.  ALF_FLAGS_EXT=-DQR_PREPIVOT acc_qbuild.sh 207 mgr 0123456789ab def "GNU noMPI" -DQR_PREPIVOT
#           acc_qbuild.sh 207 mgr 0123456789ab log "GNU noMPI LOG" -DSTABLOG
#   Work tree: cand/<WHO>/<SHA12>/qw (= D5W of the 207-2 scripts): src/ is extracted from src.tar after re-hashing it against
#   .alf_candidate; the 207-2 scripts t207_d5_build.sh (GNU oracle build, tree qw/tree_gnu_<pp|log>) and t207_d5_qbuild.sh run from that
#   extracted tree under env -i.  Result: qw/bin/qref_<pp|log>.out (+ .sha256, .src.sha256) and cand/.../bin/qref_<path>.qbuild.json
#   (tar hash, candidate, sha256 of every 207-2 script used, of the program and its source, compiler version, log hashes).
#   Refuses to overwrite.
#   ACC_QBUILD_HOME (default $HOME): the HOME the 207-2 scripts see (they read $HOME/alf-gpu/env/gnu.sh and $HOME/alf-gpu/{hdf5,venv});
#   only for a machine without ~/alf-gpu/env/gnu.sh (intel: a directory under the ticket dir with that layout); recorded in the json.
set -o pipefail
export LC_ALL=C
TICKET=$1; WHO=$2; S12=$3; P=$4; CONF=$5; QF=$6
case "$TICKET$WHO$S12" in *..*|*/*|"") echo "unsafe name"; exit 2;; esac
case $P in def) TR=gnu_pp; OUT=qref_pp;; log) TR=gnu_log; OUT=qref_log;; *) echo "path must be def or log"; exit 2;; esac
C=$HOME/alf-gpu/work/t$TICKET/cand/$WHO/$S12; [ -f $C/.alf_candidate ] || { echo "no shipped candidate $C"; exit 2; }
want=$(awk '$1=="tar_sha256" {print $2}' $C/.alf_candidate); have=$(sha256sum $C/src.tar | cut -d' ' -f1)
[ -n "$want" ] && [ "$want" = "$have" ] || { echo "TAR HASH MISMATCH for src ($have vs $want)"; exit 2; }
case "$(printf %s "$CONF" | tr a-z A-Z)" in *GPU=*) echo "the oracle is a CPU build: no GPU= token"; exit 2;; esac
QH=${ACC_QBUILD_HOME:-$HOME}
W=$C/qw; export TMPDIR=$HOME/alf-gpu/work/t$TICKET/tmp; mkdir -p $TMPDIR $C/bin $C/logs $W/logs $W/bin $W/tmp
[ -e "$C/bin/qref_$P.qbuild.json" ] && { echo "refusing: qref_$P.qbuild.json exists"; exit 2; }
if [ ! -d $W/src ]; then mkdir $W/src && tar -xf $C/src.tar -C $W/src || exit 2; fi   # from the tarball just re-hashed
S=$W/src/docs/gpu/scripts/nividia
cp $S/t207_d5_qref.F90 $W/ || exit 2
run() { env -i HOME="$QH" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin TMPDIR="$W/tmp" D5W="$W" J=4 ALF_FLAGS_EXT="${ALF_FLAGS_EXT:-}" "$@"; }
# 1. GNU CPU oracle (the 207-2 script sources the machine's gnu.sh itself); 2. the offline program against its objects
run bash $S/t207_d5_build.sh src $TR gnu $CONF > $C/logs/qbuild_oracle_$P.txt 2>&1 || { echo "oracle build failed"; tail -20 $C/logs/qbuild_oracle_$P.txt; exit 3; }
run bash $S/t207_d5_qbuild.sh $TR $OUT $QF > $C/logs/qbuild_prog_$P.txt 2>&1 || { echo "offline program build failed"; tail -20 $C/logs/qbuild_prog_$P.txt; exit 4; }
python3 - "$C" "$W" "$P" "$OUT" "$TR" "$CONF" "$QF" "$have" "$QH" "${ALF_FLAGS_EXT:-}" <<'PY'
import hashlib, json, os, subprocess, sys, datetime
c, w, p, out, tr, conf, qf, tarh, qh, ext = sys.argv[1:11]
h = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
s = os.path.join(w, "src", "docs", "gpu", "scripts", "nividia")
tools = {n: h(os.path.join(s, n)) for n in sorted(f for f in os.listdir(s) if f.startswith("t207_d5_"))}
b = os.path.join(w, "bin", out + ".out")
cand = dict(l.split(None, 1) for l in open(os.path.join(c, ".alf_candidate")).read().splitlines() if " " in l)
ev = {l.split("=", 1)[0]: l.rstrip("\n").split("=", 1)[1] for l in open(os.path.join(w, "logs", "env_" + tr + ".txt")) if "=" in l}
sidecar = open(b + ".sha256").read().split()[0]
src_side = open(os.path.join(w, "bin", out + ".src.sha256")).read().split()[0]
assert sidecar == h(b) and src_side == h(os.path.join(w, "t207_d5_qref.F90")) == tools["t207_d5_qref.F90"]
gf = subprocess.run(["bash", "-c", ". %s/alf-gpu/env/gnu.sh >/dev/null 2>&1; gfortran --version | head -1" % qh], capture_output=True, text=True).stdout.strip()
d = {"schema": "alf-acc-qbuild/1", "path": p, "program": b, "prog_sha256": h(b), "tar_sha256_built": tarh, "candidate": cand,
     "configure": conf, "ALF_FLAGS_EXT": ext, "qflags": qf, "oracle_tree": tr, "ALF_FLAGS_PROG": ev.get("ALF_FLAGS_PROG", ""),
     "ALF_LIB": ev.get("ALF_LIB", ""), "qref_src_sha256": h(os.path.join(w, "t207_d5_qref.F90")), "tools_sha256": tools,
     "gfortran": gf, "home": qh, "env_file_sha256": h(os.path.join(qh, "alf-gpu", "env", "gnu.sh")),
     "logs_sha256": {f: h(os.path.join(c, "logs", f)) for f in ("qbuild_oracle_%s.txt" % p, "qbuild_prog_%s.txt" % p)},
     "built": datetime.datetime.now().astimezone().isoformat()}
json.dump(d, open(os.path.join(c, "bin", "qref_%s.qbuild.json" % p), "w"), indent=1)
print("built", out, d["prog_sha256"][:16], gf)
PY

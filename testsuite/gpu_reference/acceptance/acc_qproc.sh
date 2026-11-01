#!/bin/bash
# Acceptance harness, MACHINE side (ticket 207-5): process the QREF dumps of one finished run before the next run starts.
#   usage: acc_qproc.sh TICKET WHO SHA12 RUNID def|log QS NULCODES HOSTIN KEEP
#   RUNID = cand/<WHO>/<SHA12>/runs/<RUNID>, written by acc_run.sh --qref (provenance.json of schema alf-acc-prov/1, qref/ with the dumps).
#   Bridge to the 207-2 pipeline, which reads records in the layout of its own runner: runs/<RUNID>/qproc/ gets copies of parameters,
#   info, run.log, a stamp.txt and a provenance.json in that layout, both generated from the run's harness provenance (verdict.py checks
#   every field of them against the harness provenance again), and the dumps (qref/ is moved, same file system).  Then
#   t207_d5_proc.sh (from the shipped candidate tree, under env -i, D5W = cand/.../qw, D5RUNS = cand/.../runs) does, in its own order:
#   spec -> manifest -> program -> strict check -> qref_proc.json -> thinning (every KEEP-th dump by the gcd keep rule).
#   qproc/acc_qproc.json records the arguments, the status, the qbuild record, the sha256 of the scripts used, the interpreter of the parser
#   steps (D5_PYTHON = ACC_PYTHON, default ~/alf-gpu/venv/bin/python3; the driver passes its own interpreter) and the machine (qref_gate.py
#   machine: host name, sha256 of /etc/machine-id, the libopenblas the offline program loads, with its sha256), which verdict.py requires to be
#   the machine it runs the gate on (the gate's exact re-processing depends on that library).
#   Record format (207-5 decision): the 207-2 parser (t207_d5_qparse.derive_run, reviewed) reads its own runner's records (rc, env token list,
#   run_log_sha256, info_sha256, seeds_sha256, parameters_sha256, binary_sha256, stamp.txt); this bridge writes exactly those, generated from the
#   harness provenance (alf-acc-prov/1), instead of teaching the reviewed parser a second format; verdict.py re-checks every bridged field
#   against the harness provenance and the attempt-log stamps before it calls the gate.
#   Exit status: 0 only if proc.sh exited 0.  Nothing is deleted before the strict check passed (proc.sh); on a failure every dump stays.
set -o pipefail
export LC_ALL=C
TICKET=$1; WHO=$2; S12=$3; RID=$4; P=$5; QS=$6; NUL=$7; HI=$8; KEEP=$9
case "$TICKET$WHO$S12$RID" in *..*|*/*|"") echo "unsafe name"; exit 2;; esac
[ $# = 9 ] || { echo "usage: acc_qproc.sh TICKET WHO SHA12 RUNID def|log QS NULCODES HOSTIN KEEP"; exit 2; }
C=$HOME/alf-gpu/work/t$TICKET/cand/$WHO/$S12; R=$C/runs/$RID; W=$C/qw; S=$W/src/docs/gpu/scripts/nividia
QB=$C/bin/qref_$P.qbuild.json
[ -f $R/provenance.json ] && [ -d $R/qref ] && [ -f $QB ] && [ -d $S ] || { echo "missing provenance/qref/qbuild record/scripts for $RID"; exit 2; }
[ -e $R/qproc ] && { echo "refusing: $R/qproc exists"; exit 2; }
QH=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['home'])" $QB) || exit 2
Q=$R/qproc; mkdir $Q || exit 2
python3 - "$R" "$Q" <<'PY' || { echo "bridge failed"; exit 2; }
import hashlib, json, os, shutil, sys
r, q = sys.argv[1:3]
h = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
p = json.load(open(os.path.join(r, "provenance.json")))
assert p.get("schema") == "alf-acc-prov/1" and p["exit_code"] == 0 and not p.get("timed_out"), "run did not end normally"
for f in ("parameters", "info", "run.log"):
    shutil.copy2(os.path.join(r, f), os.path.join(q, f))
tok = [t for t in p["extra_env"] if t.startswith("ALF_GPU_")]
assert all(t.startswith("ALF_GPU_") for t in p["extra_env"]), "extra variables other than ALF_GPU_* in a QREF run"
envtok = sorted("%s=%s" % (k, v) for k, v in p["env"].items() if k.startswith("ALF_GPU_"))
assert envtok == sorted(tok), "the recorded env has ALF_GPU_* %s, the declared extra env %s" % (envtok, sorted(tok))
seeds, par = p["inputs"]["seeds_sha256"], p["inputs"]["parameters_sha256"]
assert par == h(os.path.join(q, "parameters")) and p["outputs"]["info"] == h(os.path.join(q, "info")) and p["outputs"]["run.log"] == h(os.path.join(q, "run.log"))
with open(os.path.join(q, "stamp.txt"), "w") as f:       # the lines t207_d5_qparse.py reads
    f.write("binary %s\nenv %s\nseeds %s parameters %s\n" % (p["binary_sha256"], " ".join(tok), seeds[:16], par[:16]))
json.dump({"rc": p["exit_code"], "env": tok, "seeds_sha256": seeds, "parameters_sha256": par, "run_log_sha256": p["outputs"]["run.log"],
           "info_sha256": p["outputs"]["info"], "binary_sha256": p["binary_sha256"], "generated_from": "provenance.json (alf-acc-prov/1) by acc_qproc.sh"},
          open(os.path.join(q, "provenance.json"), "w"), indent=1)
PY
mv $R/qref $Q/qref || exit 2
BASEF=$(df -BG --output=avail $C | tail -1 | tr -dc 0-9)
PYQ=${ACC_PYTHON:-$HOME/alf-gpu/venv/bin/python3}
PROG=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['program'])" $QB) || exit 2
MACH=$(python3 "$(dirname "$0")/qref_gate.py" machine "$PROG") || { echo "machine record failed"; exit 2; }
s=$(date +%s)
env -i HOME="$QH" USER="$USER" LANG=C LC_ALL=C PATH=/usr/bin:/bin D5W="$W" D5RUNS="$C/runs" D5_PYTHON="$PYQ" \
  bash $S/t207_d5_proc.sh $Q $P $QS $NUL $HI $KEEP > $Q/proc.log 2>&1
rc=$?
python3 - "$Q" "$QB" "$S" "$rc" "$((`date +%s`-s))" "$P $QS $NUL $HI $KEEP" "$BASEF" "$PYQ" "$MACH" <<'PY'
import hashlib, json, os, sys
q, qb, s, rc, wall, args, free, pyq, mach = sys.argv[1:10]
h = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
P, QS, NUL, HI, KEEP = args.split()
json.dump({"schema": "alf-acc-qproc/1", "args": {"path": P, "qs": int(QS), "nul": NUL, "hostin": int(HI), "keep": int(KEEP)}, "proc_rc": int(rc),
           "wall_s": int(wall), "free_GB_before": int(free or -1), "qbuild_sha256": h(qb), "python": pyq, "machine": json.loads(mach),
           "tools_sha256": {n: h(os.path.join(s, n)) for n in sorted(f for f in os.listdir(s) if f.startswith("t207_d5_"))},
           "proc_log_sha256": h(os.path.join(q, "proc.log"))}, open(os.path.join(q, "acc_qproc.json"), "w"), indent=1)
PY
tail -3 $Q/proc.log
exit $rc

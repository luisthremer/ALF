#!/bin/bash
# Ticket 207 dev5 / 207-2 (accuracy test against a quad reference, development only), nividia, CPU only: run the offline quad-reference
# program on the dumps of one run directory, check its output strictly, and only then thin the dumps.
#   usage: [D5W=workdir] [D5RUNS=runsdir] [DUMP_COVER=K] [ALLOW_V1=1] t207_d5_proc.sh RUNDIR def|log QS NULCODES HOSTIN KEEP
#   RUNDIR must be under $D5RUNS (default $D5W/runs; default $D5W = ~/alf-gpu/work/t207/dev5) and holds the run's records (parameters, info, run.log,
#   stamp.txt, provenance.json) and qref/q*.bin; DUMP_COVER=K (default 1): the directory holds only every K-th dump of the run
#   (a directory thinned earlier; dumps may be symbolic links to them); ALLOW_V1=1 accepts version 1 dumps (smoke runs, no B, Bi).  Steps, in this order, any failure stops the script:
#     1. t207_d5_qparse.py spec -> qref_spec.json: what the run's own records say must be in the directory and in qref.txt (the
#        dump files present must be exactly the expected ones);
#     2. qref_manifest.txt: sha256 of every dump (written before the program runs; an existing manifest must equal the fresh one);
#     3. the program -> qref.txt (stdout), qref_err.txt (stderr), exit status;
#     4. t207_d5_qparse.py check -> qref_check.txt (every dump processed, every member present, every number finite and in the
#        format, ...); 5. qref_proc.json (program hashes, arguments, status, wall time, check result);
#     6. only if the program exited 0 and the check passed and KEEP > 0: t207_d5_qparse.py thin deletes the dumps that are not every
#        KEEP-th of the covered ones, only files listed in the manifest whose sha256 still matches and that are not symbolic links
#        (list in qref_thinned.txt).  Kept: records, manifest, full qref.txt, every KEEP-th dump.  Exit status: 0 only if all passed.
#   Program hash: the sha256 of the binary is computed and must equal the one in its .sha256 file, and the .src.sha256 file must name the
#   source t207_d5_qref.F90 next to this script (else refusal); the binary hash is recorded in qref_proc.json ('prog'; 'prog_src' = the
#   verified source hash); t207_d5_agg.py needs the record for every run and t207_d5_gate.py requires the value of its spec (prog_sha, src_sha).  KEEP > 0 must satisfy the keep rule of
#   qparse (gcd(cover*KEEP, P) = 1); the kept dumps allow pointwise re-checking only, the recomputable result is the full qref.txt plus the manifest.
#   qref_proc.json also records qref_sha256 (the sha256 of qref.txt after the check): agg and the gate refuse a qref.txt that differs from it;
#   and manifest_sha256 (fix round 4: the sha256 of qref_manifest.txt), which the gate requires to equal the hash of the manifest it parses.
#   These records sit in the run directory and can be rewritten together with what they describe: they detect accidents, not forgery.  The
#   gate re-runs the program on the kept dumps (authenticity of the numbers of those dumps); binding qref.txt, qref_proc.json and the manifest
#   to the run beyond that is the harness's job (207-5: stamps of all of them in the attempt log outside the run directory).
#   The gate (t207_d5_gate.py) refuses a run processed with DUMP_COVER != 1 (unless its spec is marked test_only) or ALLOW_V1 (always),
#   and (fix round 5) a run whose QS, NULCODES, HOSTIN or KEEP differ from its spec or that keeps fewer dumps than the spec's min_reprocessed;
#   it re-computes the kept dumps and demands the qref.txt lines of them exactly (timing tq aside).
#   4 OpenMP threads, OpenBLAS 1 thread.
#   Environment: ~/alf-gpu/env/gnu.sh is sourced only if it exists (it does on nividia; the intel machine has no such file, there the
#   caller's environment is used as it is; the program binary itself needs no LD_LIBRARY_PATH on nividia, see program_env in the gate).
W=${D5W:-$HOME/alf-gpu/work/t207/dev5}; export TMPDIR=$W/tmp
if [ -f ~/alf-gpu/env/gnu.sh ]; then . ~/alf-gpu/env/gnu.sh > /dev/null; fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R=$1; P=$2; QS=$3; NUL=$4; HI=$5; KEEP=$6; COVER=${DUMP_COVER:-1}; V1=""; [ -n "$ALLOW_V1" ] && V1=--allow-v1
PY=~/alf-gpu/venv/bin/python3
[ $# -eq 6 ] || { echo "usage: t207_d5_proc.sh RUNDIR def|log QS NULCODES HOSTIN KEEP"; exit 1; }
RUNS=${D5RUNS:-$W/runs}                                 # the only tree this script works in (the self-test uses its own)
case "$R" in $RUNS/*) ;; *) echo "refusing: $R not under $RUNS"; exit 1;; esac
case "$R" in *..*) echo "refusing: $R"; exit 1;; esac
case $P in def) B=$W/bin/qref_pp.out;; log) B=$W/bin/qref_log.out;; *) echo "bad path"; exit 1;; esac
case "$KEEP" in ''|*[!0-9]*) echo "KEEP must be a non-negative integer"; exit 1;; esac
case "$COVER" in ''|*[!0-9]*|0) echo "DUMP_COVER must be a positive integer"; exit 1;; esac
[ -x "$B" ] || { echo "refusing: $B missing"; exit 1; }
# the program that runs is the one its .sha256 file names (review fix 7): compute, compare, record
PROGSHA=$(sha256sum "$B" | cut -d' ' -f1); WANTSHA=$(cut -d' ' -f1 "$B.sha256" 2>/dev/null)
[ -n "$PROGSHA" ] && [ "$PROGSHA" = "$WANTSHA" ] || { echo "refusing: sha256 of $B is '$PROGSHA', $B.sha256 says '$WANTSHA'"; exit 1; }
# prog_src in qref_proc.json is only recorded if the .src.sha256 file names the source that sits next to this script (fix round 2)
SRCSHA=$(cut -d' ' -f1 "${B%.out}.src.sha256" 2>/dev/null); HERESRC=$(sha256sum "$HERE/t207_d5_qref.F90" 2>/dev/null | cut -d' ' -f1)
[ -n "$SRCSHA" ] && [ "$SRCSHA" = "$HERESRC" ] || { echo "refusing: ${B%.out}.src.sha256 says '$SRCSHA', $HERE/t207_d5_qref.F90 has '$HERESRC' (the binary was not built from this source)"; exit 1; }
cd "$R" || exit 1
for f in qref.txt qref_spec.json qref_proc.json qref_check.txt qref_thinned.txt; do [ -e $f ] && { echo "refusing: $R/$f exists"; exit 1; }; done
$PY -I $HERE/t207_d5_qparse.py spec "$R" --path $P --qs "$QS" --nul "$NUL" --hostin "$HI" --cover "$COVER" --keep "$KEEP" || { echo "proc FAILED: spec"; exit 2; }
( cd qref && ls q*.bin ) | sed 's|^|qref/|' | LC_ALL=C sort > qref_list.txt
if [ -e qref_manifest.txt ]; then
  sha256sum $(cat qref_list.txt) > qref_manifest.fresh.txt || exit 2
  cmp -s qref_manifest.txt qref_manifest.fresh.txt || { echo "proc FAILED: existing qref_manifest.txt differs from the fresh one"; exit 2; }
  rm -f qref_manifest.fresh.txt
else
  sha256sum $(cat qref_list.txt) > qref_manifest.txt || exit 2
fi
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
s=$(date +%s)
$B "$QS" "$NUL" "$HI" < qref_list.txt > qref.txt 2> qref_err.txt; rc=$?
e=$(date +%s)
$PY -I $HERE/t207_d5_qparse.py check "$R" $V1 > qref_check.txt 2>&1; crc=$?
$PY -I - "$R" "$B" "$P" "$QS" "$NUL" "$HI" "$COVER" "$KEEP" $rc $crc $((e-s)) "$PROGSHA" <<'PY'
import hashlib, json, sys, os
R, B, P, QS, NUL, HI, COVER, KEEP, rc, crc, wall, progsha = sys.argv[1:13]
def first(p):
    return open(p).read().split()[0]
json.dump(dict(ticket='207-2 pipeline (development evidence, dev seeds)', run=R, prog=progsha,
               prog_src=first(B.replace('.out', '.src.sha256')), args=dict(path=P, qs=int(QS), nul=NUL, hostin=int(HI), cover=int(COVER), keep=int(KEEP)),
               rc=int(rc), check_rc=int(crc), qref_sha256=hashlib.sha256(open('qref.txt', 'rb').read()).hexdigest(),
               manifest_sha256=hashlib.sha256(open('qref_manifest.txt', 'rb').read()).hexdigest(), check=open('qref_check.txt').read().strip(), wall_s=int(wall),
               lines=sum(1 for l in open('qref.txt') if l.startswith('QREF file')),
               dumps=sum(1 for l in open('qref_list.txt'))), open('qref_proc.json', 'w'), indent=1)
PY
[ $? -eq 0 ] && [ -s qref_proc.json ] || { echo "proc FAILED: could not write qref_proc.json"; exit 2; }
echo "proc $R rc=$rc check_rc=$crc wall=$((e-s))s: $(tail -1 qref_check.txt)"
if [ $rc -ne 0 ] || [ $crc -ne 0 ]; then
  echo "proc FAILED: program exit status $rc, check status $crc: every dump is kept"; exit 3
fi
if [ "$KEEP" -gt 0 ]; then
  $PY -I $HERE/t207_d5_qparse.py thin "$R" $V1 || { echo "proc FAILED: thin (nothing deleted)"; exit 4; }
  $PY -I $HERE/t207_d5_qparse.py check "$R" --dumps $V1 >> qref_check.txt 2>&1 || { echo "proc FAILED: dump directory after thinning"; exit 5; }
fi
exit 0

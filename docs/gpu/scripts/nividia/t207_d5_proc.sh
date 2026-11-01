#!/bin/bash
# Ticket 207 dev5 (accuracy test against a quad reference, development only), nividia, CPU only: run the offline quad-reference
# program on the dumps of one run directory.   usage: t207_d5_proc.sh RUNDIR def|log QS NULCODES HOSTIN KEEP
#   -> RUNDIR/qref_manifest.txt (sha256 of every dump, written before anything else), RUNDIR/qref.txt (one QREF line per dump),
#      RUNDIR/qref_proc.json (program sha256, arguments, wall time).  Afterwards the dumps are removed except every KEEP-th file
#      (intermediate debug data: 2-5 MB each; the manifest keeps their hashes).  4 OpenMP threads, OpenBLAS 1 thread.
W=$HOME/alf-gpu/work/t207/dev5; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh > /dev/null
R=$1; P=$2; QS=$3; NUL=$4; HI=$5; KEEP=$6
case "$R" in $W/runs/*) ;; *) echo "refusing: $R not under dev5/runs"; exit 1;; esac
case "$R" in *..*) echo "refusing: $R"; exit 1;; esac
case $P in def) B=$W/bin/qref_pp.out;; log) B=$W/bin/qref_log.out;; *) echo "bad path"; exit 1;; esac
cd "$R" || exit 1
[ -e qref.txt ] && { echo "refusing: $R/qref.txt exists"; exit 1; }
ls qref/q*.bin > qref_list.txt 2>/dev/null || { echo "no dumps in $R"; exit 1; }
[ -e qref_manifest.txt ] || sha256sum qref/q*.bin > qref_manifest.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
s=$(date +%s)
$B "$QS" "$NUL" "$HI" < qref_list.txt > qref.txt 2> qref_err.txt; rc=$?
e=$(date +%s)
~/alf-gpu/venv/bin/python3 -I -c "import json,sys; json.dump(dict(ticket='207 dev5 (development evidence, dev seeds)', run=sys.argv[1], prog=open(sys.argv[2]+'.sha256').read().split()[0], prog_src=open(sys.argv[2].replace('.out','.src.sha256')).read().split()[0], args=sys.argv[3:6], rc=int(sys.argv[6]), wall_s=int(sys.argv[7]), lines=sum(1 for l in open('qref.txt') if l.startswith('QREF')), dumps=sum(1 for l in open('qref_list.txt'))), open('qref_proc.json','w'), indent=1)" "$R" "$B" "$QS" "$NUL" "$HI" $rc $((e-s))
echo "proc $R rc=$rc wall=$((e-s))s lines=$(grep -c '^QREF' qref.txt) dumps=$(wc -l < qref_list.txt)"
if [ $rc -eq 0 ] && [ "$(grep -c '^QREF' qref.txt)" -eq "$(wc -l < qref_list.txt)" ] && [ "${KEEP:-0}" -gt 0 ]; then
  i=0
  while read -r f; do
    i=$((i+1))
    case "$f" in qref/q[0-9]*.bin) ;; *) continue;; esac
    [ $(( (i-1) % KEEP )) -eq 0 ] || rm -f -- "$f"
  done < qref_list.txt
fi

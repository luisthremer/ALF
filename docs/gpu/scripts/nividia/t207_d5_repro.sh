#!/bin/bash
# Ticket 207-2 (nividia, CPU only): re-run the new pipeline (t207_d5_proc.sh: spec, manifest, program, strict check) on the dumps that dev5
# kept, without touching dev5: for every dev5 run directory a new directory $D5W/runs/<tag>/hub12x12_cb_p with copies of the run's
# records (parameters info run.log stamp.txt provenance.json) and qref/ holding symbolic links to dev5's kept dumps; DUMP_COVER = the
# spacing of the kept dumps; KEEP = 0 (nothing is deleted).  The recorded results are compared by t207_d5_repro.py.
#   usage: D5W=~/alf-gpu/work/t207-2 t207_d5_repro.sh TAG [TAG...]       TAG in m_def_dev1..4 m_log_dev1..4 c_def_dev1 c_log_dev1
#                                                                          d_def_dev1 d_log_dev1 s_def s_log c_def_dev1_q20 c_log_dev1_q20
SRC=$HOME/alf-gpu/work/t207/dev5                       # read only
W=${D5W:?set D5W (the work dir of the ticket, never dev5 itself)}
[ "$W" != "$SRC" ] || { echo "refusing: D5W is dev5"; exit 1; }
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NUL=0,4,100004,1,100001,2,1608,4808,5,100005; NULL=0,4,100004,1,100001,2,5,100005
for T in "$@"; do
  case $T in
    m_def_dev1) P=def; QS=20; NU=$NUL; HI=0; CV=20;;
    m_log_dev1) P=log; QS=20; NU=$NULL; HI=0; CV=20;;
    m_def_dev[234]) P=def; QS=20; NU=$NUL; HI=0; CV=40;;
    m_log_dev[234]) P=log; QS=20; NU=$NULL; HI=0; CV=40;;
    c_def_dev1) P=def; QS=0; NU=$NUL; HI=0; CV=40;;
    c_log_dev1) P=log; QS=0; NU=$NULL; HI=0; CV=40;;
    # the same control dumps processed with qs 20 like the main runs (the gate takes one qs per gate run; not compared with dev5, whose
    # record of these runs has qs 0)
    c_def_dev1_q20) P=def; QS=20; NU=$NUL; HI=0; CV=40;;
    c_log_dev1_q20) P=log; QS=20; NU=$NULL; HI=0; CV=40;;
    d_def_dev1) P=def; QS=0; NU=$NUL; HI=1; CV=10;;
    d_log_dev1) P=log; QS=0; NU=$NULL; HI=1; CV=10;;
    s_def) P=def; QS=1; NU=$NUL; HI=0; CV=1;;
    s_log) P=log; QS=1; NU=$NULL; HI=0; CV=1;;
    *) echo "unknown tag $T"; exit 1;;
  esac
  S=$SRC/runs/${T%_q20}/hub12x12_cb_p; D=$W/runs/$T/hub12x12_cb_p
  [ -d $S/qref ] || { echo "no dumps in $S"; exit 1; }
  [ -e $D ] && { echo "refusing: $D exists"; exit 1; }
  mkdir -p $D/qref || exit 1
  for f in parameters info run.log stamp.txt provenance.json; do cp -p $S/$f $D/$f || exit 1; done
  for f in $S/qref/q*.bin; do ln -s $f $D/qref/$(basename $f); done
  echo "== $T path $P QS $QS cover $CV: $(ls $D/qref | wc -l) dumps"
  V1=""; case $T in s_def|s_log) V1=1;; esac          # smoke runs: dump version 1 (older binaries, no B, Bi)
  ALLOW_V1=$V1 DUMP_COVER=$CV D5W=$W bash $HERE/t207_d5_proc.sh $D $P $QS $NU $HI 0
done

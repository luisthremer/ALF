#!/bin/bash
# Ticket 302c-2 (development only), on intel: all calibration tables, generated (t302c2_tables.py -> t207_d5_agg.py of 5d876291) from the
# processed runs; output in tab/ (agg_*.json/.txt per measure and reference list, tables_<tag>.{json,txt}).
#   Reference lists (agg --ref LIST syntax; 'nvh' is the pipeline's name of the in-run host member, here the ifx/MKL host CGR2_2 'ifh'):
#     ens  = the frozen list of 207-2 with ifh in place of nvh: default cpu,nvh,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005; LOG
#            without n1608, n4808
#     noh  = the same without the in-run host member (sensitivity of S_ens to that member)
#   Tags: b<beta>_<path>; runs: main m<beta>_<path>_devN, ctl c<beta>_<path>_dev1 where present; --prog-sha = sha256 of the binary128 program of the path (bin/qref_{pp,log}.out.sha256).
cd ~/alf-gpu/work/t302c-2 || exit 1; export TMPDIR=$PWD/tmp
P=~/alf-gpu/venv/bin/python3; mkdir -p tab
RD=cpu,nvh,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005; RL=cpu,nvh,n4,n100004,n1,n100001,n2,n5,n100005
for b in 4 10; do
  c=hub12x12_cb_p; [ $b = 10 ] && c=hub12x12_cb_b10_p
  for p in def log; do
    args=()
    for s in dev1 dev2 dev3 dev4; do f=runs/m${b}_${p}_$s/$c/qref.txt; [ -s $f ] && args+=("main:$s=$f"); done
    f=runs/c${b}_${p}_dev1/$c/qref.txt; [ -s $f ] && args+=("ctl:ctl1=$f")
    [ ${#args[@]} -gt 0 ] || continue
    r=$RD; q=bin/qref_pp.out; [ $p = log ] && { r=$RL; q=bin/qref_log.out; }
    $P -I t302c2_tables.py --prog-sha $(cut -d' ' -f1 $q.sha256) tab b${b}_$p ens=$r noh=${r/nvh,/} -- "${args[@]}" > tab/tables_b${b}_$p.log 2>&1; echo "tables b${b}_$p rc=$? (${#args[@]} runs)"
  done
done

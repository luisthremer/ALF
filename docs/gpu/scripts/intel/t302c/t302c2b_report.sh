#!/bin/bash
# Ticket 302c-2 round 2 (development only), on intel: all tables, generated from the processed runs of t302c-2b/runs, per beta and path:
#   tab/gates_b<beta>_<path>.{json,txt}   t302c2b_gates.py (approved 4118a1de reader + agg.per_run; gate-style pooled / each_run / by_stratum
#                                         of S_ens, S_ens_q95, S_ens_max; separation ranges; NVIDIA-proposed gates at their bounds)
#   tab/cal_b<beta>_<path>.{json,txt}     t207_3_cal.py unchanged (207-3's generator of the NVIDIA numbers; flavour subsets, strata, pairs),
#                                         with the strict reader of 4118a1de (T207_3_QPDIR = the work dir)
#   Runs: main m<beta>_<path>_dev*, faults f<beta>_<path>_* (only processed ones, qref_proc.json present).
cd ~/alf-gpu/work/t302c-2b || exit 1; export TMPDIR=$PWD/tmp
P=~/alf-gpu/venv/bin/python3; mkdir -p tab
for b in 4 10; do
  c=hub12x12_cb_p; [ $b = 10 ] && c=hub12x12_cb_b10_p
  for p in def log; do
    g=(); k=()
    for d in runs/m${b}_${p}_dev*/$c; do [ -f $d/qref_proc.json ] && { l=$(basename $(dirname $d)); g+=("main:$l=$d"); k+=("main:$l=$d"); }; done
    for d in runs/f${b}_${p}_*/$c; do [ -f $d/qref_proc.json ] && { l=$(basename $(dirname $d)); g+=("flt:$l=$d"); k+=("flt:$l=$d"); }; done
    [ ${#g[@]} -gt 0 ] || continue
    q=bin/qref_pp.out; [ $p = log ] && q=bin/qref_log.out
    $P -I t302c2b_gates.py tab/gates_b${b}_$p.json b${b}_$p $p $(cut -d' ' -f1 $q.sha256) "${g[@]}" > tab/gates_b${b}_$p.txt 2>&1
    echo "gates b${b}_$p rc=$? (${#g[@]} runs)"
    rm -f tab/cal_b${b}_$p.json
    T207_3_STRATA=1 T207_3_QPDIR=$PWD $P -I t207_3_cal.py tab/cal_b${b}_$p.json $p "${k[@]}" > tab/cal_b${b}_$p.txt 2>&1
    echo "cal b${b}_$p rc=$?"
  done
done

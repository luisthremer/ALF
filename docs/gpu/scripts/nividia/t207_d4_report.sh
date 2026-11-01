#!/bin/bash
# Ticket 207 dev4 (development only), nividia: statistics of the main runs (t207_d4_stat.py per run, t207_d4_agg.py per path over the
# seed sets), the control stand-in check on dev1 (t207_d4_ctlcmp.py against the dev2 paired / control runs), the data.h5 identity of
# the dev1 runs with the dev2 paired runs, and the GPU minutes of all dev4 runs.  Output: dev4/logs/report_*.txt / *.json.
W=$HOME/alf-gpu/work/t207/dev4; cd $W || exit 1
P=~/alf-gpu/venv/bin/python3; D2=$HOME/alf-gpu/work/t207/dev2/runs
for p in def log; do
  args=()
  for s in dev1 dev2 dev3 dev4; do
    R=runs/m_${p}_$s/hub12x12_cb_p
    $P -I t207_d4_stat.py $R/run.log 331 $R/stat.json > logs/report_stat_${p}_$s.txt || echo "stat failed $p $s"
    $P -I t207_diagstat.py $R/run.log 331 --expect-points 3200 --expect-used 560 --check 0.05 1.25 2 > logs/report_diagstat_${p}_$s.txt
    echo "diagstat rc=$?" >> logs/report_diagstat_${p}_$s.txt
    args+=("$s=$R/stat.json")
  done
  $P -I t207_d4_agg.py logs/report_agg_$p.json "${args[@]}" > logs/report_agg_$p.txt
done
{
  echo "== default, dev1: control stand-in V[8] vs production ALF_GPU_TAU_NS=0 run (dev2 t0c1, same seeds)"
  $P -I t207_d4_ctlcmp.py runs/m_def_dev1/hub12x12_cb_p/run.log $D2/t0c1_pd_def_b4/hub12x12_cb_p/run.log $D2/t0c1_ctl_pd_tauns0_def/hub12x12_cb_p/run.log 331
  echo "== LOG, dev1: control stand-in V[8] vs production ALF_GPU_TAUSWAP=8 run (dev2 t0c1, same seeds)"
  $P -I t207_d4_ctlcmp.py runs/m_log_dev1/hub12x12_cb_p/run.log $D2/t0c1_pd_log_b4/hub12x12_cb_p/run.log $D2/t0c1_ctl_pd_tauswap8_log/hub12x12_cb_p/run.log 331
  echo "== data.h5 of the dev1 TAUDIAG runs vs the dev2 paired runs (same seeds, production outputs)"
  bash t207_d4_h5cmp.sh runs/m_def_dev1/hub12x12_cb_p/data.h5 $D2/t0c1_pd_def_b4/hub12x12_cb_p/data.h5
  bash t207_d4_h5cmp.sh runs/m_log_dev1/hub12x12_cb_p/data.h5 $D2/t0c1_pd_log_b4/hub12x12_cb_p/data.h5
  echo "== production-like identity (no debug variable), HEAD vs working tree"
  grep h5cmp logs/chain_ident.txt
  grep -h "GPU debug env" runs/id_*/hub12x12_cb_b10_s/run.log | sort | uniq -c
  echo "== GPU minutes (wall of the locked section, all dev4 GPU runs)"
  $P -I - <<'PY'
import glob, json
t = 0
for f in sorted(glob.glob('runs/*/*/provenance.json')):
    p = json.load(open(f)); t += p['wall_s']; print('  %-40s %7.1f s rc %d' % (f.split('/')[1] + '/' + f.split('/')[2], p['wall_s'], p['rc']))
print('  total %.1f s = %.1f GPU-min' % (t, t / 60))
PY
} > logs/report_checks.txt 2>&1
echo report done

#!/bin/bash
# GPU-hour ledger of the 2026-11 production (sacct/eacct are unreliable): node-hours per job from RUN_*/segments.log
# (start .. end, or .. next start / now for a killed or running segment) x #SBATCH --nodes of the run's job.sh, plus a manual
# list of test jobs (ledger_tests.txt: "jobid nodes hours note" per line). 1 node-h = 4 GPU-h (4 Max 1550 cards per node).
# usage: bash ledger.sh [PRODDIR] [TESTLIST]      prints per run, per job, and totals
P=${1:-/hppfs/work/pn93la/di38zuw/prod_2026-11}
T=${2:-$(dirname "$0")/ledger_tests.txt}
now=$(date +%s)
tot=0
for d in "$P"/RUN_*; do
  [ -f "$d/segments.log" ] || continue
  n=$(grep -o '^#SBATCH --nodes=[0-9]*' "$d/job.sh" | head -1 | cut -d= -f2)
  # job id, epoch of every start / end line
  nh=$(grep -E '\| (start|end) ' "$d/segments.log" | while IFS='|' read -r ts job rest; do
         j=${job#*=}; j=${j// /}; e=$(date -d "$ts" +%s); k=$(echo "$rest" | cut -d'|' -f2 | tr -d ' '); echo "$j $k $e"; done |
       awk -v n="$n" -v now="$now" -v run="$(basename "$d")" '
         $2 == "start" && !($1 in s) { s[$1] = $3; ord[++m] = $1 }
         $2 == "end" { e[$1] = $3 }
         END { t = 0
               for (i = 1; i <= m; i++) { j = ord[i]; en = (j in e) ? e[j] : (i < m ? s[ord[i+1]] : now)
                 h = (en - s[j]) / 3600; t += h * n; printf "  %-34s %s %6.2f h x %d = %7.1f node-h\n", run, j, h, n, h * n > "/dev/stderr" }
               printf "%.2f", t }')
  printf '%-36s %8.1f node-h\n' "$(basename "$d")" "$nh"
  tot=$(awk -v a="$tot" -v b="$nh" 'BEGIN { print a + b }')
done
tt=0
if [ -f "$T" ]; then
  tt=$(awk '!/^#/ && NF >= 3 { s += $2 * $3 } END { print s + 0 }' "$T")
fi
printf '%-36s %8.1f node-h\n' "tests (ledger_tests.txt)" "$tt"
awk -v a="$tot" -v b="$tt" 'BEGIN { printf "TOTAL %.0f node-h = %.0f GPU-h (production %.0f, tests %.0f node-h)\n", a + b, 4 * (a + b), a, b }'

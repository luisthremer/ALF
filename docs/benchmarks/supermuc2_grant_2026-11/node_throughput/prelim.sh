#!/bin/bash
# Build the report PDF from the current state of the runs ("Hey give preliminary results").
# usage: ./prelim.sh            -> real data: analyse.py on SuperMUC over all finished runs, figures, tables, PDF
#        ./prelim.sh fake       -> plots/results_fake.csv instead
# Output: report/report_<prelim|fake>.pdf, results.csv (real mode). Runs still in progress are simply not included yet
# (a run directory appears in runs/ only when the run has finished). Needs ssh supermuc2, python with matplotlib +
# scienceplots (PY, default: python3), pdflatex.
set -e
H=$(cd "$(dirname "$0")" && pwd); cd "$H"
PY=${PY:-python3}
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
if [ "$1" = fake ]; then
  MODE=fake; CSV=plots/results_fake.csv; $PY plots/make_fake_data.py > $CSV
  NOTE='\textbf{FAKE DATA.} All numbers and figures in this version are invented, only to check layout and plots.'
  FAKEFLAG=--fake
else
  MODE=prelim; CSV=results.csv
  rsync -a analyse.py supermuc2:$A/jobs/nt/ 2>/dev/null
  timeout 120 ssh supermuc2 "python3 $A/jobs/nt/analyse.py $A/runs" 2>/dev/null | grep -E '^(side|cpu|gpu)' > $CSV
  NRUN=$(($(wc -l < $CSV) - 1))
  QUEUE=$(timeout 60 ssh supermuc2 'squeue -u $USER -h -o "%j %T"' 2>/dev/null | grep -c '^nt_' || true)
  NOTE="\\textbf{Preliminary.} $NRUN finished runs (of 100 planned incl.\\ GPU repeats), $QUEUE jobs still queued or running at $(date '+%Y-%m-%d %H:%M')."
  FAKEFLAG=
fi
rm -rf report/fig; mkdir -p report/fig
$PY plots/plot_all.py $CSV report/fig $FAKEFLAG
$PY plots/make_tables.py $CSV > report/tables.tex
printf '\\newcommand{\\DataState}{%s, %s}\n\\newcommand{\\DataNote}{%s}\n' "$MODE" "$(date '+%Y-%m-%d %H:%M')" "$NOTE" > report/datastate.tex
cd report
pdflatex -interaction=nonstopmode -halt-on-error report.tex > build.log 2>&1 || { tail -30 build.log; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error report.tex > build.log 2>&1
cp report.pdf report_$MODE.pdf
echo "$H/report/report_$MODE.pdf"

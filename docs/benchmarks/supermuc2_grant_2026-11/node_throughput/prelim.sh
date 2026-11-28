#!/bin/bash
# Build the report PDF from the current state of the runs ("Hey give preliminary results").
# usage: ./prelim.sh            -> real data: analyse.py on SuperMUC over all finished runs, figures, tables, PDF
#        ./prelim.sh fake       -> plots/results_fake.csv instead
# Output: report/report_<prelim|fake>.pdf (internal), report/antrag_<prelim|fake>.pdf (proposal version), results.csv (real mode). Runs still in progress are simply not included yet
# (a run directory appears in runs/ only when the run has finished). Needs ssh supermuc2, python with matplotlib +
# scienceplots (PY, default: python3), pdflatex.
set -e
H=$(cd "$(dirname "$0")" && pwd); cd "$H"
PY=${PY:-python3}
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
# SET=betaL: the beta = L series (25 runs, CPU t1 at L=36 out of memory -> 24 counted; t2 required) (runs_betaL/, report_betaL/, results_betaL.csv); default: the beta = 4 series
if [ "$SET" = betaL ]; then RUNS=runs_final_betaL; RDIR=report_betaL; CSVR=results_betaL.csv; NPLAN=34; export NLAY=1 NLAY_L= REQ_TAGS=cpu=t2 BETA_LABEL='\beta = L' CPU_EXP=7 REPORTDIR=report_betaL
else RUNS=runs_final4; RDIR=report; CSVR=results.csv; NPLAN=63; export NLAY=1 NLAY_L= PLANNED_L=16,24,32,36,40 REQ_TAGS= BETA_LABEL='\beta = 4' CPU_EXP=6 REPORTDIR=report; fi
if [ "$1" = fake ]; then
  MODE=fake; CSV=plots/results_fake.csv; $PY plots/make_fake_data.py > $CSV
  NOTE='\textbf{FAKE DATA.} All numbers and figures in this version are invented, only to check layout and plots.'
  FAKEFLAG=--fake
else
  MODE=prelim; CSV=$CSVR
  rsync -a analyse.py supermuc2:$A/jobs/nt/ 2>/dev/null
  timeout 120 ssh supermuc2 "python3 $A/jobs/nt/analyse.py $A/$RUNS" 2>/dev/null | grep -E '^(side|cpu|gpu)' > $CSV
  NRUN=$(($(wc -l < $CSV) - 1))
  QUEUE=$(timeout 60 ssh supermuc2 'squeue -u $USER -h -o "%j %T"' 2>/dev/null | grep -cE '^(nt|nb)_' || true)
  NOTE="\\textbf{Preliminary.} $NRUN finished runs (of $NPLAN planned), $QUEUE jobs still queued or running at $(date '+%Y-%m-%d %H:%M')."
  [ "$NRUN" -ge "$NPLAN" ] && { MODE=final; NOTE="\\textbf{Final.} All $NRUN runs finished ($(date '+%Y-%m-%d %H:%M'))."; }
  FAKEFLAG=
fi
rm -rf $RDIR/fig; mkdir -p $RDIR/fig
# tuned settings (tune.py on SuperMUC) -> tuned.tex
if [ "$1" = fake ]; then printf '\\newcommand{\\TunedTable}{}\n' > $RDIR/tuned.tex
else timeout 60 ssh supermuc2 "cat $A/auto/tuned.json" 2>/dev/null > $RDIR/tuned.json; timeout 60 ssh supermuc2 "cat $A/auto/final_betaL.json" 2>/dev/null > $RDIR/final_betaL.json
  $PY plots/make_tuned.py $RDIR/tuned.json $RDIR/final_betaL.json > $RDIR/tuned.tex; fi
$PY plots/plot_all.py $CSV $RDIR/fig $FAKEFLAG
$PY plots/make_tables.py $CSV > $RDIR/tables.tex
$PY plots/make_antrag.py $CSV > $RDIR/antrag_data.tex
$PY plots/cpu_fit.py $RDIR/antrag_best.dat $CPU_EXP $RDIR/fig/cpu_fit
printf '\\newcommand{\\DataState}{%s, %s}\n\\newcommand{\\DataNote}{%s}\n' "$MODE" "$(date '+%Y-%m-%d %H:%M')" "$NOTE" > $RDIR/datastate.tex
cd $RDIR
pdflatex -interaction=nonstopmode -halt-on-error report.tex > build.log 2>&1 || { tail -30 build.log; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error report.tex > build.log 2>&1
cp report.pdf report_$MODE.pdf
# proposal version (standalone, proposal language); figure first (figures/node_performance.pdf)
(cd figures && cp ../antrag_best.dat ../antrag_fac.dat ../antrag_facx.dat . && pdflatex -interaction=nonstopmode -halt-on-error node_performance.tex > build.log 2>&1 && pdflatex -interaction=nonstopmode -halt-on-error node_performance.tex > build.log 2>&1) || { tail -30 figures/build.log; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error antrag.tex > build_antrag.log 2>&1 || { tail -30 build_antrag.log; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error antrag.tex > build_antrag.log 2>&1
cp antrag.pdf antrag_$MODE.pdf
echo "$H/$RDIR/report_$MODE.pdf"
echo "$H/$RDIR/antrag_$MODE.pdf"

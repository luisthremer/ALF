#!/bin/bash
# L=32 production U-scan: one single-GPU job chain per (U, part): WALKERS ranks of one U on one GPU (MPS), NPARTS parts per U
# (NPARTS x WALKERS walkers per U, seeds unique through setup_pp.py --walker0).
# usage: submit_prod.sh RUNS_ROOT WALKERS NPARTS MPS_PCT NDELAY PARTITION WALLTIME NLINKS [MARGIN_H]   (set DRY=1 to only create the directories)
# Writes RUNS_ROOT/README_jobs.txt (job ids) and RUNS_ROOT/seed_check.txt.
ROOT=$1; W=$2; NP=$3; PCT=$4; ND=$5; PART=$6; WT=$7; NL=$8; MARGIN=${9:-0.3}
HERE=/hnvme/workspace/b133ae29-ALF_GPU/ALF/docs/gpu/scripts/helma
UL="5.0 5.1 5.2 5.25 5.3 5.35 5.4 5.45 5.5 5.6"
mkdir -p $ROOT
for U in $UL; do for p in $(seq 0 $((NP-1))); do
  d=$ROOT/L32_U${U}_p$p
  if [ -e $d ]; then echo "exists: $d (not touched)"; continue; fi
  python3 $HERE/setup_pp.py $d 32 10 $ND 0 $W $U --walker0 $((p*W)) --master helma2026 || exit 1
  if [ -z "$DRY" ]; then $HERE/chain_submit.sh $d 1 $W $PCT $PART $WT $NL $MARGIN | tee -a $ROOT/README_jobs.txt; fi
done; done
cat $ROOT/*/seeds | sort -n | uniq -d | wc -l | xargs echo "duplicate seeds across all runs:" | tee $ROOT/seed_check.txt

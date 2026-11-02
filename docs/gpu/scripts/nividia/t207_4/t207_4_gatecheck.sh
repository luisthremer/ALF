#!/bin/bash
# Ticket 207-4 (development only), nividia: the approved gate (g5 = 4118a1de) on the beta = 10 runs of one path, to check that it accepts
# KEEP 11 and re-processes the 22 kept dumps of every run (all P = 20 positions) bit-identically.  The kept dumps are packed (zstd,
# t207_4_pack.sh); they are unpacked for the gate and packed again.  Done inside ~/alf-gpu/gpu.lock so that no other agent's run writes
# dumps while the unpacked dumps (~0.7 GB per path) are on disk (no GPU use; the lock is held for the CPU-only gate, ~10 min).
#   Gates and bounds: the beta = 4 proposal of 207-3 (inputs, not a decision); controls: var8 and the production control on acc_pTT,
#   eps 1e-12 / 1e-10 / near on the gates named below.  The verdict is expected to show which controls are not killed at beta = 10.
#   usage: t207_4_gatecheck.sh TAG def|log
W=$HOME/alf-gpu/work/t207-4; cd $W || exit 1; export TMPDIR=$W/tmp
TAG=$1; P=$2; C=hub12x12_cb_b10_p
case $P in def) NX=2e-15; NT=2em15;; log) NX=1e-14; NT=1em14;; *) echo bad path; exit 1;; esac
RUNS="m_${P}_dev11 m_${P}_dev12 m_${P}_dev13 c_${P}_dev14 f_${P}_1em12_dev15 f_${P}_1em10_dev16 f_${P}_${NT}_dev18"
exec 9> "$HOME/alf-gpu/gpu.lock"; flock 9; echo "gpu.lock taken $(date -Is)"
fk=$(df -k --output=avail / | tail -1 | tr -d ' ')
if [ "$fk" -lt 5800000 ]; then echo "refusing: / $fk kB free (unpacking ~0.7 GB would approach the 5 GB floor)"; exit 1; fi
for r in $RUNS; do bash t207_4_pack.sh unpack runs/$r/$C | tail -1; done
bash t207_4_gate.sh $TAG $P "dev11 dev12 dev13" "ctl:dev14 flt:1e-12:dev15 flt:1e-10:dev16 flt:$NX:dev18" \
  --gate acc_pTT:pTT:pooled:0.15 --gate pw_emax:emax:pooled:0.15 --gate strat_emax:emax:by_stratum:0.30 --gate q95_emax:emax:pooled:0.60:S_ens_q95 \
  --control var8:variant:var8:acc_pTT --control ctl:run:ctl:acc_pTT --control near:run:flt_$NT:acc_pTT \
  --control f12:run:flt_1em12:pw_emax --control f10:run:flt_1em10:strat_emax
for r in $RUNS; do bash t207_4_pack.sh pack runs/$r/$C | tail -1; done
flock -u 9; exec 9>&-; echo "gpu.lock released $(date -Is)"

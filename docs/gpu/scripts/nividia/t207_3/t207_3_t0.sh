#!/bin/bash
# Ticket 207-3 (development only, dev seeds), nividia: T0.  -> t207-3/logs/t0_<RUN>.txt, runs t207-3/runs/t0<RUN>_<id>/
#   usage: t207_3_t0.sh RUN PART...     PART:
#   selftest  SELFTEST=4, 6 configs x default (nv_def_pp) / LOG (nv_log)
#   faults    perm|prop|prec|blocks|eps:1e-9|eps:1e-11 under SELFTEST=4 (both paths), invalid values
#   off       task 4: tau pass off (Ltau = 0 twin hub6x6_cb_s_lt0) with SELFTEST=4 / TAUTIME=1 / FAULT_TAU / TAUDIAG_* -> must exit non-zero;
#             the same config without them (DEBUG=1 only) exits 0; production refusal unchanged (no DEBUG)
#   ident     normal runs (no debug variable) of the 207-3 build and of base (HEAD before the edits): hub12x12_cb_b10_s (tau pass on)
#             and hub6x6_cb_s_lt0 (tau pass off), default (nv_def) and LOG; data.h5 compared bit for bit
#   off2      (follow-up t207-3b, builds nv_def_pp3 / nv_log3) the added refusals ALF_GPU_TAUSWAP=8 / ALF_GPU_TAU_NS=0 / ALF_GPU_TAUTEST_V=1
#             on the tau-off input (DEBUG=1) -> must exit non-zero with the refusal line; SELFTEST=4 still refused; the old builds
#             nv_def_pp2 / nv_log2 with TAUSWAP=8 on the same input (rc 0 before the fix); the same three variables with the tau pass on
#             (hub6x6_cb_s) -> rc 0 (the refusal fires only when the pass is off)
#   ident3    (follow-up) normal runs (no ALF_GPU_* variable) of nv_def_pp3 / nv_log3 and of the previous nv_def_pp2 / nv_log2 on
#             hub6x6_cb_s (tau pass on) and hub6x6_cb_s_lt0 (off), data.h5 bit for bit; DEBUG=1 without and with eps:1e-6 against the
#             previous agent's runs t03_at_<old>_{0,all} in PREV (default ~/alf-gpu/work/t207-3, read only; same dev1 seeds via SETROOT)
#   epsprod   graded output fault in the sweep: hub6x6_cb, DEBUG=1, eps:1e-6 and eps:1e-10 vs no fault (Precision tau in info)
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
R=$1; shift; out=$W/logs/t0_$R.txt
case "$R" in ""|*/*|*..*) echo bad run label; exit 1;; esac
r() { # id binary config env...
  id=$1; b=$2; c=$3; shift 3
  res=$(TAG=t0${R}_$id bash $W/t207_3_run.sh $b $c "$@")
  echo "== $id: $res" >> $out
  grep -hE "PASSED|FAILED|GPU debug env|GPU fault|refused|invalid|tau self-test|GPU tau pass|branch" $W/runs/t0${R}_$id/$c/run.log | cut -c1-260 | uniq >> $out
}
h5() { echo "h5cmp $1: $(D5W=$W bash $W/t207_d5_h5cmp.sh $2 $3) rc=$?" >> $out; }
E="ALF_GPU_EXIT=1"
for part in "$@"; do case $part in
selftest)
  for c in hub6x6_cb hub6x6_dense hub12x12_cb hub12x12_cb_b10 hub6x6_cb_flux hub6x6_cb_symm; do
    r st4pp_$c nv_def_pp $c $E ALF_GPU_SELFTEST=4; r st4log_$c nv_log $c $E ALF_GPU_SELFTEST=4; done ;;
faults)
  for f in perm prop prec blocks eps:1e-9 eps:1e-11; do
    r fpp_${f/:/_} nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f
    r flog_${f/:/_} nv_log hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done
  for f in xyz eps:abc eps:-1 eps:; do r inv_${f/:/_} nv_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done ;;
off)
  for b in nv_def_pp nv_log; do
    r off_${b}_st4 $b hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4
    for kv in ALF_GPU_TAUTIME=1 ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_FAULT_TAU=perm ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_QREF=3 ALF_GPU_TAUDIAG_FROM=1 \
              ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QVAR=8 ALF_GPU_TAUDIAG_SEED=1 ALF_GPU_TAUDIAG_VAR=8 ALF_GPU_TAUDIAG_AMP=1 ALF_GPU_TAUDIAG_NUL=4; do
      r off_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 $kv; done
    r off_${b}_dbg $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1
    r off_${b}_on_st4 $b hub6x6_cb_s $E ALF_GPU_SELFTEST=4
  done
  r off_prod_ref nv_def hub6x6_cb_s_lt0 ALF_GPU_TAUDIAG=2
  r off_base_st4 base_nv_def hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4          # before 207-3: rc 0 although nothing was tested
  r off_base_fault base_nv_def hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=perm ;;
ident)
  for c in hub12x12_cb_b10_s hub6x6_cb_s_lt0; do
    for b in nv_def base_nv_def nv_log base_nv_log; do r id_${b}_$c $b $c; done
    h5 "def $c" $W/runs/t0${R}_id_nv_def_$c/$c/data.h5 $W/runs/t0${R}_id_base_nv_def_$c/$c/data.h5
    h5 "log $c" $W/runs/t0${R}_id_nv_log_$c/$c/data.h5 $W/runs/t0${R}_id_base_nv_log_$c/$c/data.h5
  done ;;
epsprod)
  for b in nv_def_pp nv_log; do
    r eps_${b}_0 $b hub6x6_cb_s ALF_GPU_DEBUG=1
    for x in 1e-10 1e-6; do r eps_${b}_$x $b hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:$x; done
    for x in 0 1e-10 1e-6; do d=$W/runs/t0${R}_eps_${b}_$x/hub6x6_cb_s; echo "  $b eps $x: $(grep -A1 -i 'Precision tau' $d/info | tr -s ' \n' ' ')" >> $out; done
  done ;;
at)   # localized fault (builds nv_def_pp2 / nv_log2 with ALF_GPU_FAULT_TAU_AT): validation, effect, identity without fault
  r at_inv1 nv_def_pp2 hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU_AT=2:1                          # without eps: refused
  r at_inv2 nv_def_pp2 hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=5:1  # nst > Nstm
  r at_inv3 nv_def_pp2 hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:3  # nfe > N_FL
  r at_inv4 nv_def_pp2 hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2,1  # syntax
  r at_inv5 nv_def_pp2 hub6x6_cb_s ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1                 # production: refused
  r at_off nv_def_pp2 hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1 # tau pass off: refused
  for b in nv_def_pp2 nv_log2; do
    r at_${b}_0 $b hub6x6_cb_s ALF_GPU_DEBUG=1
    r at_${b}_all $b hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6
    r at_${b}_21 $b hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1
    for k in 0 all 21; do d=$W/runs/t0${R}_at_${b}_$k/hub6x6_cb_s; echo "  $b fault $k: $(grep -A1 -i 'Precision tau' $d/info | tr -s ' \n' ' ')" >> $out; done
  done
  h5 "def pp2 vs pp (DEBUG=1, no fault)" $W/runs/t0${R}_at_nv_def_pp2_0/hub6x6_cb_s/data.h5 $W/runs/t02_eps_nv_def_pp_0/hub6x6_cb_s/data.h5
  h5 "log2 vs log (DEBUG=1, no fault)" $W/runs/t0${R}_at_nv_log2_0/hub6x6_cb_s/data.h5 $W/runs/t02_eps_nv_log_0/hub6x6_cb_s/data.h5 ;;
off2)
  for b in nv_def_pp3 nv_log3; do
    r off2_${b}_st4 $b hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4
    for kv in ALF_GPU_TAUSWAP=8 ALF_GPU_TAU_NS=0 ALF_GPU_TAUTEST_V=1; do
      r off2_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 $kv
      r on2_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s ALF_GPU_DEBUG=1 $kv; done
  done
  for b in nv_def_pp2 nv_log2; do r off2_old_${b}_TAUSWAP_8 $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 ALF_GPU_TAUSWAP=8; done ;;
ident3)
  PREV=${PREV:-$HOME/alf-gpu/work/t207-3}
  for pr in nv_def_pp3:nv_def_pp2 nv_log3:nv_log2; do
    n=${pr%%:*}; o=${pr##*:}
    for c in hub6x6_cb_s hub6x6_cb_s_lt0; do
      r id3_${n}_$c $n $c; r id3_${o}_$c $o $c
      h5 "$n vs $o $c (no env)" $W/runs/t0${R}_id3_${n}_$c/$c/data.h5 $W/runs/t0${R}_id3_${o}_$c/$c/data.h5
    done
    r id3_${n}_dbg $n hub6x6_cb_s ALF_GPU_DEBUG=1
    h5 "$n DEBUG=1 vs PREV t03_at_${o}_0" $W/runs/t0${R}_id3_${n}_dbg/hub6x6_cb_s/data.h5 $PREV/runs/t03_at_${o}_0/hub6x6_cb_s/data.h5
    r id3_${n}_eps $n hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6
    h5 "$n DEBUG=1 eps:1e-6 vs PREV t03_at_${o}_all" $W/runs/t0${R}_id3_${n}_eps/hub6x6_cb_s/data.h5 $PREV/runs/t03_at_${o}_all/hub6x6_cb_s/data.h5
  done ;;
esac; done
echo done >> $out

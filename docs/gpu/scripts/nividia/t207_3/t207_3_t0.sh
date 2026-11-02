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
#   at0       (fix round t207-3c, builds nv_def_pp4 / nv_log4; old = nv_def_pp3 / nv_log3, the 9a35bdd6 builds) ALF_GPU_FAULT_TAU_AT=0:<nfe>:
#             bad values refused (0:0, 0:3, 5:2, -1:2, 0:, :2, 0:2 without eps, 0:2 in production); effect per flavour: QREF dumps of every
#             CGR2_2 call (stride 1, passes 46..50) of eps:1e-6 with AT=0:2 against the unfaulted run (same seeds; inputs, B and the
#             flavour-1 outputs must be bit-identical, every flavour-2 output must differ); identity without AT against the 9a35bdd6 binaries
#             (no ALF_GPU_* variable; DEBUG=1 eps:1e-6; both against t207-3b's t04_id3_* runs, same dev1 seeds) and AT=2:1 (old
#             semantics) against t207-3's t03_at_*2_21 (428576dd binaries, same AT code); data.h5 bit for bit; self-test and two
#             fault classes on hub6x6_cb (the hook still works)
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
at0)
  PB=${PREVB:-$HOME/alf-gpu/work/t207-3b}; PA=${PREVA:-$HOME/alf-gpu/work/t207-3}
  F="ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6"
  r at0_inv1 nv_def_pp4 hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:0       # nfe 0
  r at0_inv2 nv_def_pp4 hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:3       # nfe > N_FL_eff
  r at0_inv3 nv_def_pp4 hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=5:2       # nst > Nstm (4)
  r at0_inv4 nv_def_pp4 hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=-1:2      # negative
  r at0_inv5 nv_def_pp4 hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:        # syntax
  r at0_inv6 nv_def_pp4 hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=:2        # syntax
  r at0_inv7 nv_def_pp4 hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU_AT=0:2              # without eps
  r at0_inv8 nv_def_pp4 hub6x6_cb_s ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=0:2  # production run
  Q="ALF_GPU_TAUDIAG_QREF=1 ALF_GPU_TAUDIAG_FROM=46"
  for b in nv_def_pp4 nv_log4; do
    r at0_q_${b}_0 $b hub6x6_cb_s ALF_GPU_DEBUG=1 $Q
    r at0_q_${b}_02 $b hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:2 $Q
    for k in 0 02; do d=$W/runs/t0${R}_at0_q_${b}_$k/hub6x6_cb_s; echo "  $b AT $k: $(grep -A1 -i 'Precision tau' $d/info | tr -s ' \n' ' ')" >> $out; done
    ~/alf-gpu/venv/bin/python3 -I - $W/runs/t0${R}_at0_q_${b}_0/hub6x6_cb_s/qref $W/runs/t0${R}_at0_q_${b}_02/hub6x6_cb_s/qref $b >> $out <<'PY'
import os, struct, sys
import numpy as np
a, b, name = sys.argv[1:4]
fa, fb = sorted(os.listdir(a)), sorted(os.listdir(b))
assert fa == fb and fa, 'dump lists differ or are empty'
st = {}
for f in fa:
    A, B = open(os.path.join(a, f), 'rb').read(), open(os.path.join(b, f), 'rb').read()
    assert len(A) == len(B), f
    h = struct.unpack('<13i', A[:52])
    assert h[0] == 20791 and h[9] == 0 and h[10] == 0 and h[11] == 0 and h == struct.unpack('<13i', B[:52]), (f, h)
    nst, nfe, n, hasb = h[4], h[5], h[7], h[12]
    lo = len(A) - 16 * n * n * (4 + 2 * hasb); hi = lo + 64 * n * n          # the four production outputs G00, GT0, GTT, G0T
    same_in = A[:lo] == B[:lo] and A[hi:] == B[hi:]                          # header, input states, B and Bi
    oa = np.frombuffer(A[lo:hi], dtype='<c16'); ob = np.frombuffer(B[lo:hi], dtype='<c16')
    rel = float(np.max(np.abs(ob - oa)) / np.max(np.abs(oa)))
    s = st.setdefault(nfe, dict(n=0, inputs_identical=0, outputs_identical=0, max_rel=0.0, min_rel=float('inf'), nst=set()))
    s['n'] += 1; s['inputs_identical'] += same_in; s['outputs_identical'] += (A[lo:hi] == B[lo:hi]); s['nst'].add(nst)
    s['max_rel'] = max(s['max_rel'], rel); s['min_rel'] = min(s['min_rel'], rel)
for nfe, s in sorted(st.items()):
    print('  QREF %s AT=0:2 vs no fault, flavour %d: dumps %d (nst %s), inputs+B identical %d, outputs identical %d, output max rel diff %.3e min %.3e'
          % (name, nfe, s['n'], ','.join(map(str, sorted(s['nst']))), s['inputs_identical'], s['outputs_identical'], s['max_rel'], s['min_rel']))
ok = st[1]['outputs_identical'] == st[1]['n'] and st[2]['outputs_identical'] == 0 and all(s['inputs_identical'] == s['n'] for s in st.values())
print('  QREF %s flavour-only check: %s' % (name, 'PASSED (flavour 1 bit-identical, every flavour-2 output perturbed)' if ok else 'FAILED'))
PY
  done
  for pr in nv_def_pp4:nv_def_pp3:nv_def_pp2 nv_log4:nv_log3:nv_log2; do
    n=${pr%%:*}; o=${pr#*:}; o=${o%%:*}; oa=${pr##*:}
    r at0_id_${n} $n hub6x6_cb_s
    h5 "$n vs $o (9a35bdd6) no env" $W/runs/t0${R}_at0_id_${n}/hub6x6_cb_s/data.h5 $PB/runs/t04_id3_${o}_hub6x6_cb_s/hub6x6_cb_s/data.h5
    r at0_eps_${n} $n hub6x6_cb_s $F
    h5 "$n vs $o (9a35bdd6) DEBUG=1 eps:1e-6" $W/runs/t0${R}_at0_eps_${n}/hub6x6_cb_s/data.h5 $PB/runs/t04_id3_${o}_eps/hub6x6_cb_s/data.h5
    r at0_21_${n} $n hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=2:1
    h5 "$n vs $oa (428576dd) DEBUG=1 eps:1e-6 AT=2:1" $W/runs/t0${R}_at0_21_${n}/hub6x6_cb_s/data.h5 $PA/runs/t03_at_${oa}_21/hub6x6_cb_s/data.h5
  done
  for b in nv_def_pp4 nv_log4; do
    r at0_st4_$b $b hub6x6_cb $E ALF_GPU_SELFTEST=4
    for f in perm eps:1e-9; do r at0_f_${b}_${f/:/_} $b hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done
  done ;;
esac; done
echo done >> $out

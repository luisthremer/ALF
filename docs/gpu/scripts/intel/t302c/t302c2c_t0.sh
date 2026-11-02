#!/bin/bash
# Ticket 302c-2 round 3 (development only, dev seeds 'base'), on intel: T0 at HEAD d3e82a36 (final 207-3 b2c71200 merged) on the B580.
# The parts of round 2 (t302c2b_t0.sh, same commands, work dir t302c-2c) plus 'at0' (207-3c t207_3_t0.sh part at0, B580 builds).
#   -> t302c-2c/logs/t0_<RUN>.txt, runs t302c-2c/runs/t0<RUN>_<id>/      usage: t302c2c_t0.sh RUN PART...     PART:
#   selftest  SELFTEST=4, 6 configs x default (in_def_pp) / LOG (in_log)
#   faults    perm|prop|prec|blocks|eps:1e-9|eps:1e-11 under SELFTEST=4 (both paths); invalid values
#   refusal   production runs (no DEBUG) of in_def with each debug variable: must be refused (plus FAULT_TAU_AT=0:2)
#   off       tau pass off (Ltau = 0 twin hub6x6_cb_s_lt0): SELFTEST=4 / TAUTIME / FAULT_TAU(_AT) / TAUDIAG_* / TAUSWAP / TAU_NS /
#             TAUTEST_V must exit non-zero (at d3e82a36 also TAUSWAP, TAU_NS, TAUTEST_V: 207-3 closing round); DEBUG=1 alone exits 0; the
#             tau-on twin with SELFTEST=4 exits 0.  LOG: also TAU_NS and TAUTEST_V (round 2 tried only TAUSWAP there)
#   at0       ALF_GPU_FAULT_TAU_AT=0:<nfe>: bad values refused (0:0, 0:3, 5:2, -1:2, 0:, :2, without eps, production); QREF dumps of every
#             CGR2_2 call (stride 1, passes 46..50) of eps:1e-6 AT=0:2 against the unfaulted run (same seeds): inputs, B and the flavour-1
#             outputs bit-identical, every flavour-2 output perturbed (the dump check of t207_3_t0.sh at0, unchanged); both paths
#   seed      same-seed hub6x6_cb_fd: in_def vs in_none_pp, in_log vs in_none_log (t207_firstdiv.py, <= 1e-9 per bin)
W=$HOME/alf-gpu/work/t302c-2c; R=$1; shift; out=$W/logs/t0_$R.txt
case "$R" in ""|*/*|*..*) echo bad run label; exit 1;; esac
export T302C2_W=$W TMPDIR=$W/tmp
r() { # id binary config env...
  local id=$1 b=$2 c=$3; shift 3
  local res; res=$(SET=base TAG=t0${R}_$id bash $W/t302c2_run.sh $b $c "$@" | tail -1)
  echo "== $id: $res" >> $out
  grep -hE "PASSED|FAILED|GPU debug env|GPU fault|refused|invalid|tau self-test|GPU tau pass|ERROR|error stop|branch" $W/runs/t0${R}_$id/$c/run.log 2>/dev/null | cut -c1-260 | uniq | head -12 >> $out
}
prec() { local d=$W/runs/t0${R}_$1/$2; echo "  $1: $(grep -A1 -i 'Precision tau' $d/info 2>/dev/null | tr -s ' \n' ' ')" >> $out; }
E="ALF_GPU_EXIT=1"
for part in "$@"; do echo "### part $part $(date -Is)" >> $out; case $part in
selftest)
  for c in hub6x6_cb hub6x6_dense hub12x12_cb hub12x12_cb_b10 hub6x6_cb_flux hub6x6_cb_symm; do
    r st4pp_$c in_def_pp $c $E ALF_GPU_SELFTEST=4; r st4log_$c in_log $c $E ALF_GPU_SELFTEST=4; done ;;
faults)
  for f in perm prop prec blocks eps:1e-9 eps:1e-11; do
    r fpp_${f/:/_} in_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f
    r flog_${f/:/_} in_log hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done
  for f in xyz eps:abc eps:-1 eps:; do r inv_${f/:/_} in_def_pp hub6x6_cb $E ALF_GPU_SELFTEST=4 ALF_GPU_FAULT_TAU=$f; done ;;
refusal)
  for kv in ALF_GPU_TAU_NS=0 ALF_GPU_TAUSWAP=8 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=331 ALF_GPU_TAUTIME=1 ALF_GPU_TAUTEST_V=1 \
            ALF_GPU_FAULT_TAU=perm ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_FAULT_TAU=xyz ALF_GPU_TAUDIAG_SEED=1 ALF_GPU_FAULT_TAU_AT=2:1 ALF_GPU_FAULT_TAU_AT=0:2; do
    r ref_$(echo ${kv} | tr "=:" "__") in_def hub6x6_cb_fd $kv; done ;;
off)
  b=in_def_pp
  r off_${b}_st4 $b hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4
  for kv in ALF_GPU_TAUTIME=1 ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_FAULT_TAU=perm ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_QREF=3 ALF_GPU_TAUDIAG_FROM=1 \
            ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QVAR=8 ALF_GPU_TAUDIAG_SEED=1 ALF_GPU_TAUDIAG_VAR=8 ALF_GPU_TAUDIAG_AMP=1 ALF_GPU_TAUDIAG_NUL=4 \
            ALF_GPU_TAUSWAP=8 ALF_GPU_TAU_NS=0 ALF_GPU_TAUTEST_V=1; do
    r off_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 $kv; done
  r off_${b}_at $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1
  r off_${b}_at02 $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=0:2
  r off_${b}_dbg $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1
  r off_${b}_on_st4 $b hub6x6_cb_s $E ALF_GPU_SELFTEST=4
  b=in_log
  r off_${b}_st4 $b hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4
  for kv in ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_TAUDIAG_QREF=3 ALF_GPU_TAUSWAP=8 ALF_GPU_TAU_NS=0 ALF_GPU_TAUTEST_V=1; do
    r off_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 $kv; done
  r off_${b}_dbg $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1
  r off_${b}_on_st4 $b hub6x6_cb_s $E ALF_GPU_SELFTEST=4 ;;
at0)
  F="ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6"
  r at0_inv1 in_def_pp hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:0       # nfe 0
  r at0_inv2 in_def_pp hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:3       # nfe > N_FL_eff
  r at0_inv3 in_def_pp hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=5:2       # nst > Nstm (4)
  r at0_inv4 in_def_pp hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=-1:2      # negative
  r at0_inv5 in_def_pp hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:        # syntax
  r at0_inv6 in_def_pp hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=:2        # syntax
  r at0_inv7 in_def_pp hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU_AT=0:2              # without eps
  r at0_inv8 in_def_pp hub6x6_cb_s ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=0:2  # production run
  Q="ALF_GPU_TAUDIAG_QREF=1 ALF_GPU_TAUDIAG_FROM=46"
  for b in in_def_pp in_log; do
    r at0_q_${b}_0 $b hub6x6_cb_s ALF_GPU_DEBUG=1 $Q
    r at0_q_${b}_02 $b hub6x6_cb_s $F ALF_GPU_FAULT_TAU_AT=0:2 $Q
    for k in 0 02; do d=$W/runs/t0${R}_at0_q_${b}_$k/hub6x6_cb_s; echo "  $b AT $k: $(grep -A1 -i 'Precision tau' $d/info | tr -s ' \n' ' ')" >> $out; done
    ~/alf-gpu/venv/bin/python3 -I - $W/runs/t0${R}_at0_q_${b}_0/hub6x6_cb_s/qref $W/runs/t0${R}_at0_q_${b}_02/hub6x6_cb_s/qref $b >> $out 2>&1 <<'PY'
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
  done ;;
seed)
  r sd_in_def in_def hub6x6_cb_fd; r sd_in_none_pp in_none_pp hub6x6_cb_fd
  r sd_in_log in_log hub6x6_cb_fd; r sd_in_none_log in_none_log hub6x6_cb_fd
  for p in "in_def in_none_pp" "in_log in_none_log"; do set -- $p
    echo "same seed $1 vs $2:" >> $out
    ~/alf-gpu/venv/bin/python3 -I $W/t207_firstdiv.py $W/runs/t0${R}_sd_$1/hub6x6_cb_fd/data.h5 $W/runs/t0${R}_sd_$2/hub6x6_cb_fd/data.h5 2>&1 | head -30 >> $out
  done
  for i in in_def in_none_pp in_log in_none_log; do prec sd_$i hub6x6_cb_fd; done ;;
*) echo "bad part $part" >> $out ;;
esac; done
echo "done $(date -Is)" >> $out

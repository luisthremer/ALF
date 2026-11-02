#!/bin/bash
# Ticket 302c-2 round 2 (development only, dev seeds 'base'), on intel: T0 at HEAD 08559f0b on the B580 (302c-1 review follow-up 1: retest
# after the 207-3 merge), the parts of 302c-1 (t302c_t0a.sh / t0b) and 207-3 (t207_3_t0.sh) with the round-1 run script (load check,
# flock ~/alf-gpu/gpu.lock, B580 env).   -> t302c-2b/logs/t0_<RUN>.txt, runs t302c-2b/runs/t0<RUN>_<id>/
#   usage: t302c2b_t0.sh RUN PART...     PART:
#   selftest  SELFTEST=4, 6 configs x default (in_def_pp) / LOG (in_log)
#   faults    perm|prop|prec|blocks|eps:1e-9|eps:1e-11 under SELFTEST=4 (both paths); invalid values
#   refusal   production runs (no DEBUG) of in_def with each debug variable: must be refused (as 302c-1 t0a 'refusal', plus FAULT_TAU_AT)
#   off       tau pass off (Ltau = 0 twin hub6x6_cb_s_lt0): SELFTEST=4 / TAUTIME / FAULT_TAU(_AT) / TAUDIAG_* / TAUSWAP / TAU_NS must exit
#             non-zero; the same config with DEBUG=1 only exits 0; the tau-on twin with SELFTEST=4 exits 0 (207-3 task 4)
#   epsprod   graded output fault in the sweep: hub6x6_cb_s, DEBUG=1, eps 0 / 1e-10 / 1e-8 / 1e-6, both paths (Precision tau in info)
#   at        localized fault ALF_GPU_FAULT_TAU_AT: invalid forms refused, effect at (2:1) vs all positions, both paths
#   seed      same-seed hub6x6_cb_fd: in_def vs in_none_pp, in_log vs in_none_log (t207_firstdiv.py, <= 1e-9 per bin)
W=$HOME/alf-gpu/work/t302c-2b; R=$1; shift; out=$W/logs/t0_$R.txt
case "$R" in ""|*/*|*..*) echo bad run label; exit 1;; esac
export T302C2_W=$W
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
            ALF_GPU_FAULT_TAU=perm ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_FAULT_TAU=xyz ALF_GPU_TAUDIAG_SEED=1 ALF_GPU_FAULT_TAU_AT=2:1; do
    r ref_$(echo ${kv} | tr "=:" "__") in_def hub6x6_cb_fd $kv; done ;;
off)
  b=in_def_pp
  r off_${b}_st4 $b hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4
  for kv in ALF_GPU_TAUTIME=1 ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_FAULT_TAU=perm ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_QREF=3 ALF_GPU_TAUDIAG_FROM=1 \
            ALF_GPU_TAUDIAG_QHOST=1 ALF_GPU_TAUDIAG_QVAR=8 ALF_GPU_TAUDIAG_SEED=1 ALF_GPU_TAUDIAG_VAR=8 ALF_GPU_TAUDIAG_AMP=1 ALF_GPU_TAUDIAG_NUL=4 \
            ALF_GPU_TAUSWAP=8 ALF_GPU_TAU_NS=0 ALF_GPU_TAUTEST_V=1; do
    r off_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 $kv; done
  r off_${b}_at $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1
  r off_${b}_dbg $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1
  r off_${b}_on_st4 $b hub6x6_cb_s $E ALF_GPU_SELFTEST=4
  b=in_log
  r off_${b}_st4 $b hub6x6_cb_s_lt0 $E ALF_GPU_SELFTEST=4
  for kv in ALF_GPU_FAULT_TAU=eps:1e-9 ALF_GPU_TAUDIAG_QREF=3 ALF_GPU_TAUSWAP=8; do
    r off_${b}_$(echo $kv | tr '=:.' '___') $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1 $kv; done
  r off_${b}_dbg $b hub6x6_cb_s_lt0 ALF_GPU_DEBUG=1
  r off_${b}_on_st4 $b hub6x6_cb_s $E ALF_GPU_SELFTEST=4 ;;
epsprod)
  for b in in_def_pp in_log; do
    r eps_${b}_0 $b hub6x6_cb_s ALF_GPU_DEBUG=1
    for x in 1e-10 1e-8 1e-6; do r eps_${b}_$x $b hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:$x; done
    for x in 0 1e-10 1e-8 1e-6; do prec eps_${b}_$x hub6x6_cb_s; done
  done ;;
at)
  r at_inv1 in_def_pp hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU_AT=2:1                            # without eps: refused
  r at_inv2 in_def_pp hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=5:1  # nst > Nstm
  r at_inv3 in_def_pp hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:3  # nfe > N_FL
  r at_inv4 in_def_pp hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2,1  # syntax
  r at_inv5 in_def_pp hub6x6_cb_s ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1                 # production: refused
  for b in in_def_pp in_log; do
    r at_${b}_21 $b hub6x6_cb_s ALF_GPU_DEBUG=1 ALF_GPU_FAULT_TAU=eps:1e-6 ALF_GPU_FAULT_TAU_AT=2:1
    prec at_${b}_21 hub6x6_cb_s
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

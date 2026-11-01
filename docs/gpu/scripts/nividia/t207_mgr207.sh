#!/bin/bash
# Manager re-run of ticket 207 acceptance on nividia. GPU processes strictly one at a time.
# SUPERSEDED 2026-10-29: historical record of a manager re-run (unchecked exit codes, prep output hidden: vacuous-pass
# sites V12/V14/V16/V20 of the evidence audit). Kept unchanged as evidence; use testsuite/gpu_reference/acceptance/.
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs
L=$W/logs/mgr207.txt; : > $L
say() { echo "$*" | tee -a $L; }
gpu_free() { n=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l); [ "$n" = 0 ] || { say "GPU busy, abort"; exit 9; }; }
cd $W
say "== builds $(date)"
for b in "src_new nv_def nvhpc PGI noMPI HDF5 GPU=nvidia" \
         "src_new nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia" \
         "src_new nv_none_log nvhpc PGI noMPI HDF5 LOG" \
         "src_new nv_none nvhpc PGI noMPI HDF5" \
         "src_head gnu_head gnu GNU noMPI HDF5" \
         "src_new gnu_new gnu GNU noMPI HDF5"; do
  bash t207_build.sh $b 2>&1 | tail -1 | tee -a $L
done
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_build.sh src_new nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia 2>&1 | tail -1 | tee -a $L
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_build.sh src_new nv_none_pp nvhpc PGI noMPI HDF5 2>&1 | tail -1 | tee -a $L
bash t207_prep.sh >/dev/null 2>&1; bash t207_seeds.sh hub12x12_cb_p 1 >> $L
say "== self-tests $(date)"
gpu_free
for n in nv_def_pp nv_log; do
  bash t207_selftest_all.sh $n hub6x6_cb hub6x6_dense hub12x12_cb hub12x12_cb_b10 hub6x6_cb_flux hub6x6_cb_symm
  say "$n: PASSED $(grep -c PASSED logs/selftest_$n.txt) FAILED $(grep -c FAILED logs/selftest_$n.txt)"
  grep -E "^==|worst slice" logs/selftest_$n.txt | cut -c1-220 >> $L
  bash t207_negctl.sh $n > /dev/null
  say "$n faults: $(grep -c 'rc=1' logs/negctl_$n.txt) of 4 rc=1"
done
say "== same seeds, full device path vs CPU oracle $(date)"
for c in hub6x6_cb_fd hub12x12_cb_fd; do
  for p in "nv_def nv_none_pp" "nv_log nv_none_log"; do
    set -- $p; gpu_free
    TAG=m_$1 timeout 900 bash t207_run.sh $1 $c | tee -a $L
    TAG=m_$2 timeout 900 bash t207_run.sh $2 $c | tee -a $L
    say "-- $c $1 vs $2"; . ~/alf-gpu/env/gnu.sh; python3 -I t207_firstdiv.py runs/m_$1/$c/data.h5 runs/m_$2/$c/data.h5 | tail -4 | tee -a $L
  done
done
say "== paired same-trajectory precision test $(date)"
for p in "nv_log hub12x12_cb_p_s1" "nv_def_pp hub12x12_cb_p"; do
  set -- $p; gpu_free
  TAG=pair_$1 timeout 900 bash t207_run.sh $1 $2 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=331 | tee -a $L
  say "-- $1 $2"; grep -i "precision" runs/pair_$1/$2/info | tee -a $L
  python3 -I t207_diagstat.py runs/pair_$1/$2/run.log 331 | tee -a $L
done
say "== GPU=none $(date)"
bash t207_ppcmp.sh | tail -1 | tee -a $L
bash t207_ident.sh 2>&1 | tee -a $L
bash t207_ctest.sh 2>&1 | tail -4 | tee -a $L
say "ALLDONE $(date)"

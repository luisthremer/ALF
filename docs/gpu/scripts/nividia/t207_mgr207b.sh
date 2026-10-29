#!/bin/bash
# Manager follow-up for ticket 207 after the review fixes (fault gating, integer update helper).
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; cd $W
L=$W/logs/mgr207b.txt; : > $L
say() { echo "$*" | tee -a $L; }
[ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" = 0 ] || { say "GPU busy"; exit 9; }
for b in "nv_def_pp" "nv_log"; do :; done
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_build.sh src_new nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia 2>&1 | tail -1 | tee -a $L
bash t207_build.sh src_new nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia 2>&1 | tail -1 | tee -a $L
for n in nv_def_pp nv_log; do
  bash t207_selftest_all.sh $n hub6x6_cb hub6x6_dense hub12x12_cb hub12x12_cb_b10 hub6x6_cb_flux hub6x6_cb_symm
  say "$n: PASSED $(grep -c PASSED logs/selftest_$n.txt) FAILED $(grep -c FAILED logs/selftest_$n.txt)"
  bash t207_negctl.sh $n > /dev/null
  say "$n faults: $(grep -c 'rc=1' logs/negctl_$n.txt) of 4 rc=1"
done
say "-- stray ALF_GPU_FAULT_TAU without the self-test must not change a production run"
TAG=st_plain bash t207_run.sh nv_log hub6x6_cb | tee -a $L
TAG=st_fault bash t207_run.sh nv_log hub6x6_cb ALF_GPU_FAULT_TAU=prop | tee -a $L
bash t207_h5cmp.sh runs/st_plain/hub6x6_cb/data.h5 runs/st_fault/hub6x6_cb/data.h5 | tee -a $L
say "ALLDONE $(date)"

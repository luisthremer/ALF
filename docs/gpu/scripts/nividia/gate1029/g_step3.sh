#!/bin/bash
# Gate 2026-10-29, nividia: step 3 measurements. 3a determinism, 3b MPS trial, 3c CPU jobs alongside a GPU run.
#   usage: g_step3.sh SHA12           (source tree $G/src_SHA12 must exist, shipped with git archive)
# Every GPU run is stamped (g_run.sh). Stops at the first unexpected condition; the MPS daemon is always shut down.
set -o pipefail
S=$1; G=$HOME/alf-gpu/work/gate_1029; cd $G || exit 2; export TMPDIR=$G/tmp
L=$G/logs/step3.txt; : > $L
say() { echo "$*" | tee -a $L; }
die() { say "STOP: $*"; exit 1; }
gpu_idle() {  # nothing on the GPU (ours or anybody else's)
  local n; n=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
  [ "$n" = 0 ] || { nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader | tee -a $L; die "GPU not idle"; }
}
run() { local d=$1; shift; gpu_idle; bash g_run.sh "$d" "$@" || die "run $d rc=$?"; say "ran $d $(grep -E '^(start|end)' $d/provenance.txt | tr '\n' ' ')"; }
same() { python3 -I g_h5same.py "$1" "$2" 2>&1 | tee -a $L; return ${PIPESTATUS[0]}; }
. ~/alf-gpu/env/gnu.sh   # python with h5py for the comparisons
B=$G/bin
say "== builds $(date -Is)"
bash g_build.sh $S nv_def nvhpc PGI noMPI HDF5 GPU=nvidia | tee -a $L || die "build nv_def"
bash g_build.sh $S nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia | tee -a $L || die "build nv_log"
bash g_build.sh $S nv_none_log nvhpc PGI noMPI HDF5 LOG | tee -a $L || die "build nv_none_log"
R=$G/runs
say "== 3a determinism: same binary twice, hub12x12_cb_b10_p, env -i $(date -Is)"
for v in nv_def nv_log; do
  run $R/3a_${v}_1 $B/${S}_$v.out hub12x12_cb_b10_p
  run $R/3a_${v}_2 $B/${S}_$v.out hub12x12_cb_b10_p
  say "-- $v run 1 vs run 2"; same $R/3a_${v}_1 $R/3a_${v}_2; say "   verdict rc=$? (0 = bit-identical)"
done
say "== 3b MPS trial (binary nv_log, hub6x6_cb and hub6x6_cb_fd) $(date -Is)"
M=$HOME/alf-gpu/work/mps_trial; mkdir -p $M/pipe $M/log $M/tmp
MPSENV=(CUDA_MPS_PIPE_DIRECTORY=$M/pipe CUDA_MPS_LOG_DIRECTORY=$M/log)
mps_quit() { env "${MPSENV[@]}" bash -c 'echo quit | nvidia-cuda-mps-control' 2>&1 | tee -a $L; sleep 3
  pgrep -u "$USER" -l -f nvidia-cuda-mps | tee -a $L; say "mps processes left: $(pgrep -u "$USER" -f nvidia-cuda-mps | wc -l)"; }
BL=$B/${S}_nv_log.out
run $R/3b_alone_6x6 $BL hub6x6_cb
run $R/3b_alone_fd $BL hub6x6_cb_fd
# two at once without MPS (time-sliced contexts)
gpu_idle
bash g_run.sh $R/3b_nomps_a $BL hub6x6_cb & p1=$!; bash g_run.sh $R/3b_nomps_b $BL hub6x6_cb & p2=$!
wait $p1; r1=$?; wait $p2; r2=$?; [ $r1 = 0 ] && [ $r2 = 0 ] || die "no-MPS pair rc $r1 $r2"
say "ran 3b_nomps_a/b"
# MPS daemon as this user, own pipe/log dirs, A30 only
gpu_idle
env -i HOME=$HOME USER=$USER PATH=/usr/bin:/bin TMPDIR=$M/tmp CUDA_VISIBLE_DEVICES=0 "${MPSENV[@]}" nvidia-cuda-mps-control -d \
  || die "MPS daemon did not start (rc $?)"
trap mps_quit EXIT
sleep 2; say "mps control: $(pgrep -u "$USER" -l -f nvidia-cuda-mps-control | tr '\n' ' ')"
bash g_run.sh $R/3b_mps_a $BL hub6x6_cb "${MPSENV[@]}" & p1=$!; bash g_run.sh $R/3b_mps_b $BL hub6x6_cb "${MPSENV[@]}" & p2=$!
sleep 4; nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader | tee -a $L
wait $p1; r1=$?; wait $p2; r2=$?; [ $r1 = 0 ] && [ $r2 = 0 ] || die "MPS pair rc $r1 $r2"
say "ran 3b_mps_a/b"
bash g_run.sh $R/3b_mpsfd_a $BL hub6x6_cb_fd "${MPSENV[@]}" & p1=$!; bash g_run.sh $R/3b_mpsfd_b $BL hub6x6_cb_fd "${MPSENV[@]}" & p2=$!
wait $p1; r1=$?; wait $p2; r2=$?; [ $r1 = 0 ] && [ $r2 = 0 ] || die "MPS fd pair rc $r1 $r2"
say "ran 3b_mpsfd_a/b"
# a client that does not know our pipe directory (as another user's job would) while the daemon is up
bash g_run.sh $R/3b_outsider $BL hub6x6_cb || die "outsider client rc $?"
say "ran 3b_outsider (no MPS env, daemon running)"
say "server log tail:"; tail -5 $M/log/control.log 2>/dev/null | tee -a $L
mps_quit; trap - EXIT
gpu_idle
for t in nomps_a nomps_b mps_a mps_b outsider; do say "-- 3b $t vs alone"; same $R/3b_alone_6x6 $R/3b_$t; say "   rc=$?"; done
for t in mpsfd_a mpsfd_b; do say "-- 3b $t vs alone"; same $R/3b_alone_fd $R/3b_$t; say "   rc=$?"; done
say "== 3c CPU oracle runs alongside one GPU run (nv_log hub6x6_cb_fd; CPU nv_none_log hub12x12_cb_fd) $(date -Is)"
BC=$B/${S}_nv_none_log.out
for k in 2 3; do
  gpu_idle; pids=()
  for i in $(seq 1 $k); do bash g_run.sh $R/3c_cpu${k}_$i $BC hub12x12_cb_fd & pids+=($!); done
  sleep 5
  bash g_run.sh $R/3c_gpu_with$k $BL hub6x6_cb_fd || die "3c GPU run rc $?"
  for p in "${pids[@]}"; do wait $p || die "3c CPU run rc $?"; done
  say "-- 3c GPU with $k CPU jobs vs alone"; same $R/3b_alone_fd $R/3c_gpu_with$k; say "   rc=$?"
done
say "-- CPU runs alongside: identical to each other (same binary, config, seeds)?"
same $R/3c_cpu2_1 $R/3c_cpu3_3; say "   rc=$?"
say "== wall times"
for d in $R/*; do s=$(grep '^start' $d/provenance.txt | cut -d' ' -f2); e=$(grep '^end' $d/provenance.txt | cut -d' ' -f2)
  say "$(basename $d) $(( $(date -d $e +%s) - $(date -d $s +%s) )) s rc=$(grep '^exit_code' $d/provenance.txt | cut -d' ' -f2)"; done
gpu_idle
say "ALLDONE $(date -Is)"

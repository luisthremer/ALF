#!/bin/bash
# Manager follow-up for ticket 207 after the audit of 2026-10-29.
#   src_final = b44c3d6e (final), src_pre = bf25ccc3 (before the review fixes), src_wip = 4d7a66e1 (Sonnet final state)
#   (each extracted with git archive of that commit; the commit is in src_*/SOURCE_SHA)
# Every run line carries the source sha and the sha256 of the binary. Raw outputs are kept (no cleanup until the gate).
# Stops at the first failed build, failed run or failed comparison (exit code checked everywhere).
set -o pipefail
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/logs; cd $W || exit 1
L=$W/logs/mgr207c.txt; : > $L
say() { echo "$*" | tee -a $L; }
die() { say "STOP: $*"; exit 1; }
gpu_free() { [ "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" = 0 ] || die "GPU busy"; }
declare -A SRC_OF
build() {  # SRC NAME [ALF_FLAGS_EXT] configure args...
  local src=$1 name=$2; shift 2
  bash t207_build.sh $src $name "$@" > logs/build_$name.out 2>&1 || { tail -5 logs/build_$name.out | tee -a $L; die "build $name"; }
  SRC_OF[$name]=$(cat $src/SOURCE_SHA)
  say "built $name from ${SRC_OF[$name]} sha256 $(sha256sum bin/ALF_$name.out | cut -c1-16)"
}
run() {  # TAG NAME CONFIG [env...]
  local tag=$1 name=$2 c=$3; shift 3; gpu_free
  local out; out=$(TAG=$tag timeout 900 bash t207_run.sh $name $c "$@") || die "run $tag $name $c timed out or failed"
  echo "$out" | grep -q " rc=0 " || die "run: $out"
  say "$out src ${SRC_OF[$name]} bin $(sha256sum bin/ALF_$name.out | cut -c1-16) env [$*]"
}
h5same() { bash t207_h5cmp.sh "$1" "$2" | tee -a $L; local rc=$?; [ $rc = 0 ] && say "  -> bit-identical" || say "  -> DIFFER (rc $rc)"; return 0; }
say "== builds $(date)"
build src_final nv_def nvhpc PGI noMPI HDF5 GPU=nvidia
build src_final nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia
build src_final nv_none_log nvhpc PGI noMPI HDF5 LOG
build src_pre nv_def_pre nvhpc PGI noMPI HDF5 GPU=nvidia
build src_pre nv_log_pre nvhpc PGI noMPI HDF5 LOG GPU=nvidia
build src_wip nv_def_wip nvhpc PGI noMPI HDF5 GPU=nvidia
export ALF_FLAGS_EXT=-DQR_PREPIVOT
build src_final nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia
build src_pre nv_def_pp_pre nvhpc PGI noMPI HDF5 GPU=nvidia
build src_final nv_none_pp nvhpc PGI noMPI HDF5
unset ALF_FLAGS_EXT
bash t207_prep.sh > logs/prep_c.log 2>&1 || die prep
bash t207_seeds.sh hub12x12_cb_p 1 >> $L || die seeds
say "== (1) pre-fix vs final binaries, same seeds, data.h5 bit for bit $(date)"
for p in "nv_def_pre nv_def" "nv_log_pre nv_log" "nv_def_pp_pre nv_def_pp"; do
  set -- $p
  run c_$1 $1 hub6x6_cb_fd; run c_$2 $2 hub6x6_cb_fd
  h5same runs/c_$1/hub6x6_cb_fd/data.h5 runs/c_$2/hub6x6_cb_fd/data.h5
done
say "== (1) default same-seed run with the final binary vs CPU QR_PREPIVOT $(date)"
run c_nv_def nv_def hub12x12_cb_fd; run c_nv_none_pp nv_none_pp hub12x12_cb_fd
. ~/alf-gpu/env/gnu.sh
python3 -I t207_firstdiv.py runs/c_nv_def/hub12x12_cb_fd/data.h5 runs/c_nv_none_pp/hub12x12_cb_fd/data.h5 | tail -5 | tee -a $L || die firstdiv
say "== (1) paired runs with the final binaries $(date)"
for p in "nv_def_pp hub12x12_cb_p" "nv_log hub12x12_cb_p_s1"; do
  set -- $p
  run cpair_$1 $1 $2 ALF_GPU_TAUDIAG=2 ALF_GPU_TAUDIAG_FROM=331
  grep -i "precision" runs/cpair_$1/$2/info | tee -a $L
  python3 -I t207_diagstat.py runs/cpair_$1/$2/run.log 331 | tee -a $L || die "diagstat $1"
done
say "== (5) default-path Green-mean drift: final twice, Sonnet state once, hub12x12_cb_b10_p $(date)"
run c5a nv_def hub12x12_cb_b10_p; run c5b nv_def hub12x12_cb_b10_p; run c5w nv_def_wip hub12x12_cb_b10_p
for t in c5a c5b c5w; do say "$t: $(grep -i 'precision green\|precision tau' runs/$t/hub12x12_cb_b10_p/info | tr -s ' ' | tr '\n' ' ')"; done
say "-- final run A vs final run B"; h5same runs/c5a/hub12x12_cb_b10_p/data.h5 runs/c5b/hub12x12_cb_b10_p/data.h5
say "-- final run A vs Sonnet state"; h5same runs/c5a/hub12x12_cb_b10_p/data.h5 runs/c5w/hub12x12_cb_b10_p/data.h5
say "ALLDONE $(date)"

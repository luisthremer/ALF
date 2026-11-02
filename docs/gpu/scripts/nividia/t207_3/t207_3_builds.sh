#!/bin/bash
# Ticket 207-3 (development only), nividia: builds (sequential, J=4) with t207_d5_build.sh (D5W = t207-3; prebuilt HDF5 in place).
#   gpu:  nv_def_pp (QR_PREPIVOT: the self-test oracle and the QREF runs, as dev5), nv_log (STABLOG), nv_def (production default) from src
#   base: base_nv_def, base_nv_log from base_src (git archive of HEAD before the 207-3 edits: identity of normal runs)
#   fix:  (follow-up t207-3b) nv_def_pp3 (QR_PREPIVOT) and nv_log3 (STABLOG) from src: the refusal list of tau_off_refuse extended
#         (TAUSWAP, TAU_NS, TAUTEST_V); otherwise the source of nv_def_pp2 / nv_log2 (git 428576dd)
#   fix2: (fix round t207-3c) nv_def_pp4 (QR_PREPIVOT) and nv_log4 (STABLOG) from src: comment fixes and ALF_GPU_FAULT_TAU_AT=0:<nfe>
#         (flavour-only fault); otherwise the source of nv_def_pp3 / nv_log3 (git 9a35bdd6)
#   prev: (follow-up) copies of the previous agent's binaries NAMES (default nv_def_pp2 nv_log2) from PREV (default ~/alf-gpu/work/t207-3,
#         read only), sha256 checked against PREV's .sha256 files (for the identity runs of the fix builds); t207-3c:
#         PREV=~/alf-gpu/work/t207-3b NAMES="nv_def_pp3 nv_log3" (the 9a35bdd6 builds: source tree hash 2f62897ab04225ef)
#   q:    the offline quad programs: copies of ~/alf-gpu/work/t207-2/bin/qref_{pp,log}.out (the 207-2 builds of the same
#         t207_d5_qref.F90; sha256 files recomputed here, the source sha256 file copied)
W=${T207_3W:-$HOME/alf-gpu/work/t207-3};   # T207_3W: work dir of a follow-up agent (t207-3b)
cd $W || exit 1
export D5W=$W TMPDIR=$W/tmp J=4
case $1 in
gpu)  ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh src nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh src nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia
      bash t207_d5_build.sh src nv_def nvhpc PGI noMPI HDF5 GPU=nvidia ;;
base) bash t207_d5_build.sh base_src base_nv_def nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh base_src base_nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
fix)  ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh src nv_def_pp3 nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh src nv_log3 nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
fix2) ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh src nv_def_pp4 nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh src nv_log4 nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
prev) PREV=${PREV:-$HOME/alf-gpu/work/t207-3}
      for b in ${NAMES:-nv_def_pp2 nv_log2}; do
        case "$b" in ""|*..*|*/*) echo "refusing name [$b]"; exit 1;; esac
        if [ -e bin/ALF_$b.out ]; then echo "refusing: bin/ALF_$b.out exists"; continue; fi
        want=$(cut -d' ' -f1 $PREV/bin/ALF_$b.out.sha256)
        cp $PREV/bin/ALF_$b.out bin/ALF_$b.out && cp $PREV/bin/ALF_$b.srctree $PREV/bin/ALF_$b.srcstate bin/
        sha256sum bin/ALF_$b.out > bin/ALF_$b.out.sha256
        [ "$(cut -d' ' -f1 bin/ALF_$b.out.sha256)" = "$want" ] && echo "$b copied, sha256 $want" || { echo "$b sha256 MISMATCH"; exit 1; }
      done ;;
q)    for q in qref_pp qref_log; do
        if [ -e bin/$q.out ]; then echo "refusing: bin/$q.out exists"; continue; fi
        cp ~/alf-gpu/work/t207-2/bin/$q.out bin/$q.out && cp ~/alf-gpu/work/t207-2/bin/$q.src.sha256 bin/$q.src.sha256
        sha256sum bin/$q.out > bin/$q.out.sha256
        echo "$q $(cut -c1-16 bin/$q.out.sha256) (207-2: $(cut -c1-16 ~/alf-gpu/work/t207-2/bin/$q.out.sha256))"
      done ;;
*) echo bad; exit 1 ;;
esac
echo "done $1 $(date -Is)"

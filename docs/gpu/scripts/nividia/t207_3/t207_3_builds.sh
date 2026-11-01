#!/bin/bash
# Ticket 207-3 (development only), nividia: builds (sequential, J=4) with t207_d5_build.sh (D5W = t207-3; prebuilt HDF5 in place).
#   gpu:  nv_def_pp (QR_PREPIVOT: the self-test oracle and the QREF runs, as dev5), nv_log (STABLOG), nv_def (production default) from src
#   base: base_nv_def, base_nv_log from base_src (git archive of HEAD before the 207-3 edits: identity of normal runs)
#   q:    the offline quad programs: copies of ~/alf-gpu/work/t207-2/bin/qref_{pp,log}.out (the 207-2 builds of the same
#         t207_d5_qref.F90; sha256 files recomputed here, the source sha256 file copied)
W=$HOME/alf-gpu/work/t207-3; cd $W || exit 1
export D5W=$W TMPDIR=$W/tmp J=4
case $1 in
gpu)  ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh src nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh src nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia
      bash t207_d5_build.sh src nv_def nvhpc PGI noMPI HDF5 GPU=nvidia ;;
base) bash t207_d5_build.sh base_src base_nv_def nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh base_src base_nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
q)    for q in qref_pp qref_log; do
        if [ -e bin/$q.out ]; then echo "refusing: bin/$q.out exists"; continue; fi
        cp ~/alf-gpu/work/t207-2/bin/$q.out bin/$q.out && cp ~/alf-gpu/work/t207-2/bin/$q.src.sha256 bin/$q.src.sha256
        sha256sum bin/$q.out > bin/$q.out.sha256
        echo "$q $(cut -c1-16 bin/$q.out.sha256) (207-2: $(cut -c1-16 ~/alf-gpu/work/t207-2/bin/$q.out.sha256))"
      done ;;
*) echo bad; exit 1 ;;
esac
echo "done $1 $(date -Is)"

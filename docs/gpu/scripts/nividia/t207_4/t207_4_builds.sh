#!/bin/bash
# Ticket 207-4 (development only), nividia: builds (sequential, J=4) with the unchanged g5/t207_d5_build.sh (D5W = t207-4; prebuilt HDF5
# in place), the env of t207_3_builds.sh (nvhpc 26.9 via ~/alf-gpu/env/nvhpc.sh).
#   gpu: nv_def_pp (QR_PREPIVOT: the QREF runs of the default path, as 207-3) and nv_log (STABLOG) from src (git archive 428576dd)
#   q:   the offline binary128 programs: copies of ~/alf-gpu/work/t207-2/bin/qref_{pp,log}.out (the 207-2 builds of t207_d5_qref.F90
#        f8383257..., the same source as g5/ = 4118a1de; sha256 files recomputed here and compared with 207-2's, source sha256 copied)
#   clean: remove the build trees (binaries are kept)
W=$HOME/alf-gpu/work/t207-4; cd $W || exit 1
export D5W=$W TMPDIR=$W/tmp J=4
case $1 in
gpu)  ALF_FLAGS_EXT=-DQR_PREPIVOT bash g5/t207_d5_build.sh src nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia
      bash g5/t207_d5_build.sh src nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
q)    for q in qref_pp qref_log; do
        if [ -e bin/$q.out ]; then echo "refusing: bin/$q.out exists"; continue; fi
        cp ~/alf-gpu/work/t207-2/bin/$q.out bin/$q.out && cp ~/alf-gpu/work/t207-2/bin/$q.src.sha256 bin/$q.src.sha256
        sha256sum bin/$q.out > bin/$q.out.sha256
        a=$(cut -c1-64 bin/$q.out.sha256); b=$(cut -c1-64 ~/alf-gpu/work/t207-2/bin/$q.out.sha256)
        [ "$a" = "$b" ] || { echo "MISMATCH $q: $a vs 207-2 $b"; exit 1; }
        s=$(cut -c1-64 bin/$q.src.sha256); g=$(sha256sum g5/t207_d5_qref.F90 | cut -c1-64)
        [ "$s" = "$g" ] || { echo "MISMATCH $q source: $s vs g5 $g"; exit 1; }
        echo "$q ${a:0:16} (= 207-2), source ${s:0:16} (= g5)"
      done ;;
clean) for t in tree_nv_def_pp tree_nv_log; do [ -d "$W/$t" ] && rm -rf -- "$W/$t" && echo "removed $t"; done ;;
*) echo bad; exit 1 ;;
esac
echo "done $1 $(date -Is)"

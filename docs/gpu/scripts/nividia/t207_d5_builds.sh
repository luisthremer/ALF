#!/bin/bash
# Ticket 207 dev5 (development only), nividia: all builds of dev5, sequentially (J=4).  Phase: gpu | gnu | base | q
#   gpu [SUFFIX]:  nv_def_pp<SUFFIX> (QR_PREPIVOT, default path) and nv_log<SUFFIX> (STABLOG) from src (working tree with the QREF dumps)
#   base: base_nv_def_pp, base_nv_log from base_src (git archive c1cc528e) for the production-identity check
#   gnu:  gnu_pp, gnu_log (GNU noMPI CPU builds; their objects are the CPU oracle of the offline program)
#   q:    the offline quad-reference programs qref_pp, qref_log
cd ~/alf-gpu/work/t207/dev5 || exit 1
export TMPDIR=~/alf-gpu/work/t207/dev5/tmp J=4
case $1 in
gpu)  ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh src nv_def_pp${2:-} nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh src nv_log${2:-} nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
base) ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh base_src base_nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia
      bash t207_d5_build.sh base_src base_nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia ;;
gnu)  ALF_FLAGS_EXT=-DQR_PREPIVOT bash t207_d5_build.sh src gnu_pp gnu GNU noMPI
      bash t207_d5_build.sh src gnu_log gnu GNU noMPI LOG ;;
q)    bash t207_d5_qbuild.sh gnu_pp qref_pp -DQR_PREPIVOT
      bash t207_d5_qbuild.sh gnu_log qref_log -DSTABLOG ;;
*) echo bad; exit 1 ;;
esac
echo "done $1 $(date -Is)"

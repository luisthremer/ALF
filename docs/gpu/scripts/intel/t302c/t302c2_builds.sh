#!/bin/bash
# Ticket 302c-2 (development only), on intel: all builds, sequentially (J = 6 <= 8 cores) -> logs/builds.txt
#   in_def_pp  GPU=intel default path with -DQR_PREPIVOT (the in-run host CGR2_2 of QHOST=1 is then the pre-pivoted algorithm, as nv_def_pp)
#   in_log     GPU=intel STABLOG
#   gnu_pp, gnu_log   GNU noMPI CPU builds (gfortran 15.2 + OpenBLAS); their objects are the CPU oracle ('cpu', null variants) of the
#              offline program
#   qref_pp, qref_log the offline binary128 programs
cd ~/alf-gpu/work/t302c-2 || exit 1
export J=6
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t302c2_build.sh in_def_pp intel IntelLLVM noMPI HDF5 GPU=intel
bash t302c2_build.sh in_log intel IntelLLVM noMPI HDF5 LOG GPU=intel
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t302c2_build.sh gnu_pp gnu GNU noMPI
bash t302c2_build.sh gnu_log gnu GNU noMPI LOG
bash t302c2_qbuild.sh gnu_pp qref_pp -DQR_PREPIVOT
bash t302c2_qbuild.sh gnu_log qref_log -DSTABLOG
echo ALLDONE

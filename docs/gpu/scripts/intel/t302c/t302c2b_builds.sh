#!/bin/bash
# Ticket 302c-2 round 2 (development only), on intel: all builds at HEAD 08559f0b, sequentially (J = 6 <= 8 cores), with the round-1 build
# scripts in the work dir t302c-2b (T302C2_W) -> logs/builds.txt
#   in_def_pp  GPU=intel default path with -DQR_PREPIVOT (QREF runs, self-test oracle; the in-run host CGR2_2 of QHOST=1 is pre-pivoted)
#   in_log     GPU=intel STABLOG
#   in_def     GPU=intel production default (T0: same-seed 6x6 vs the ifx GPU=none oracle, as 302c-1 t0b)
#   in_none_pp, in_none_log   ifx GPU=none oracles (QR_PREPIVOT resp. STABLOG), as 302c-1
#   gnu_pp, gnu_log   GNU noMPI CPU builds (gfortran 15.2 + OpenBLAS); their objects are the CPU oracle of the offline program
#   qref_pp, qref_log the offline binary128 programs (t207_d5_qref.F90 of 4118a1de = that of 5d876291, byte-identical)
cd ~/alf-gpu/work/t302c-2b || exit 1
export J=6 T302C2_W=$HOME/alf-gpu/work/t302c-2b
{ dpkg-query -W "intel-oneapi-compiler-fortran*" "intel-oneapi-mkl*" intel-opencl-icd libze-intel-gpu1 libze1 libigc2 2>&1 | tr "\n" ";"; echo; gfortran --version | head -1; } > logs/dpkg.txt
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t302c2_build.sh in_def_pp intel IntelLLVM noMPI HDF5 GPU=intel
bash t302c2_build.sh in_log intel IntelLLVM noMPI HDF5 LOG GPU=intel
bash t302c2_build.sh in_def intel IntelLLVM noMPI HDF5 GPU=intel
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t302c2_build.sh in_none_pp intel IntelLLVM noMPI HDF5
bash t302c2_build.sh in_none_log intel IntelLLVM noMPI HDF5 LOG
ALF_FLAGS_EXT=-DQR_PREPIVOT bash t302c2_build.sh gnu_pp gnu GNU noMPI
bash t302c2_build.sh gnu_log gnu GNU noMPI LOG
bash t302c2_qbuild.sh gnu_pp qref_pp -DQR_PREPIVOT
bash t302c2_qbuild.sh gnu_log qref_log -DSTABLOG
sha256sum bin/*.out > logs/binaries.sha256.txt
echo ALLDONE

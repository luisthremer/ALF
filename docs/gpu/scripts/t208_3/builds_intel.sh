#!/bin/bash
# intel (J=6 <= 8): GPU builds and GPU=none oracles (ifx) -> ~/alf-gpu/work/t208-3/logs/builds.txt
cd ~/alf-gpu/work/t208-3 || exit 1; export J=6
B=scripts/build.sh
ALF_FLAGS_EXT=-DQR_PREPIVOT bash $B in_def_pp intel IntelLLVM noMPI HDF5 GPU=intel
bash $B in_def  intel IntelLLVM noMPI HDF5 GPU=intel
bash $B in_log  intel IntelLLVM noMPI HDF5 LOG GPU=intel
bash $B inn_def intel IntelLLVM noMPI HDF5
bash $B inn_log intel IntelLLVM noMPI HDF5 LOG
echo BUILDS_DONE

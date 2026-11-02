#!/bin/bash
# nividia (J=4): GPU builds and GPU=none oracles of the same compiler -> ~/alf-gpu/work/t208-2/logs/builds.txt
cd ~/alf-gpu/work/t208-2 || exit 1; export J=4
B=scripts/build.sh
ALF_FLAGS_EXT=-DQR_PREPIVOT bash $B nv_def_pp nvhpc PGI noMPI HDF5 GPU=nvidia
bash $B nv_def  nvhpc PGI noMPI HDF5 GPU=nvidia
bash $B nv_log  nvhpc PGI noMPI HDF5 LOG GPU=nvidia
bash $B nvn_def nvhpc PGI noMPI HDF5
bash $B nvn_log nvhpc PGI noMPI HDF5 LOG
echo BUILDS_DONE

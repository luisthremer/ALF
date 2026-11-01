#!/bin/bash
cd ~/alf-gpu/work/t302c
echo "$(dpkg-query -W "intel-oneapi-compiler-fortran*" "intel-oneapi-mkl*" intel-opencl-icd libze-intel-gpu1 libze1 libigc2 2>&1 | tr "\n" ";")" > logs/dpkg.txt
for v in "" base_; do
  if [ -n "$v" ]; then export SRC=base_src; else unset SRC; fi
  bash t302c_ibuild.sh ${v}in_def IntelLLVM noMPI HDF5 GPU=intel
  bash t302c_ibuild.sh ${v}in_log IntelLLVM noMPI HDF5 LOG GPU=intel
  ALF_FLAGS_EXT=-DQR_PREPIVOT bash t302c_ibuild.sh ${v}in_none_pp IntelLLVM noMPI HDF5
  bash t302c_ibuild.sh ${v}in_none_log IntelLLVM noMPI HDF5 LOG
done
echo ALLDONE

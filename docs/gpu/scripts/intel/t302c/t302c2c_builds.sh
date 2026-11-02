#!/bin/bash
# Ticket 302c-2 round 3 (development only), on intel: all builds at HEAD d3e82a36 (round-2 scripts d0faa7ac + final 207-3 b2c71200),
# sequentially (J = 6 <= 8 cores), with the round-1 build scripts in the work dir t302c-2c (T302C2_W) -> logs/builds.txt
# The same set and flags as round 2 (t302c2b_builds.sh):
#   in_def_pp  GPU=intel default path with -DQR_PREPIVOT (QREF runs, self-test oracle)   in_log  GPU=intel STABLOG
#   in_def     GPU=intel production default (T0 same-seed 6x6 vs in_none_pp)             in_none_pp, in_none_log  ifx GPU=none oracles
#   gnu_pp, gnu_log  GNU noMPI CPU builds (objects = CPU oracle of the offline program)   qref_pp, qref_log  offline binary128 programs
# At the end the sha256 of every binary is compared with round 2's (t302c-2b/bin, read only), for the record (207-3c changed only
# Prog/gpu_tau_mod.F90, which only GPU=intel builds compile).
cd ~/alf-gpu/work/t302c-2c || exit 1
export J=6 T302C2_W=$HOME/alf-gpu/work/t302c-2c TMPDIR=$HOME/alf-gpu/work/t302c-2c/tmp
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
echo "== sha256 vs round 2 (t302c-2b/bin, 08559f0b)"
for f in bin/*.out; do b=$(basename "$f"); o=$HOME/alf-gpu/work/t302c-2b/bin/$b
  [ -e "$o" ] && echo "$b $( [ "$(sha256sum < "$f")" = "$(sha256sum < "$o")" ] && echo same || echo different)"; done
echo ALLDONE

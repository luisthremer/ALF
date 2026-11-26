#!/bin/bash
# Build the gpu-port reference (fork/gpu-port tip, git archive in SRC) on the login node: GPU=intel noMPI, AOT for
# PVC as gpu-port's own Phase 2 builds (spir64 -> spir64_gen + device pvc). usage: build_gp.sh SRC OUT
src=$1; out=$2
here=$(cd "$(dirname "$0")" && pwd)
. "$here/env_smuc2.sh"
rm -rf "$out"; cp -r "$src" "$out"; cd "$out"
source ./configure.sh IntelLLVM noMPI GPU=intel NO-INTERACTIVE > configure.log 2>&1
for v in ALF_FLAGS_MODULES ALF_FLAGS_ANA ALF_FLAGS_PROG ALF_LIB; do
  eval "export $v=\"\${$v//-fopenmp-targets=spir64/-fopenmp-targets=spir64_gen -Xopenmp-target-backend -device -Xopenmp-target-backend pvc}\""
done
export ALF_SYCL_BLK_FLAGS='-fsycl-targets=spir64_gen -Xs "-device pvc"'
make -j16 lib program > make.log 2>&1; echo "gp build rc=$?"; ls -la Prog/ALF.out

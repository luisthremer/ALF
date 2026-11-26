#!/bin/bash
# Build the two binaries of the checks from one source tree (login node): GPU=intel (AOT for PVC) and GPU=none
# with the same ifx (the same-compiler CPU reference). usage: build_smuc2.sh SRC_DIR OUT_DIR
set -eu
src=$1; out=$2
here=$(cd "$(dirname "$0")" && pwd)
. "$here/env_smuc2.sh"
mkdir -p "$out"
for t in intel none; do
  rm -rf "$out/src_$t"; cp -r "$src" "$out/src_$t"
  ( cd "$out/src_$t"
    if [ $t = intel ]; then export ALF_INTEL_AOT=pvc; else unset ALF_INTEL_AOT; fi
    . ./configure.sh IntelLLVM noMPI GPU=$t NO-INTERACTIVE > "$out/configure_$t.log" 2>&1
    make -j16 all > "$out/make_$t.log" 2>&1 )
  cp "$out/src_$t/Prog/ALF.out" "$out/ALF-$t.out"
  echo "built $out/ALF-$t.out"
done

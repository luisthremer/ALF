#!/bin/bash
# Extra GPU=intel builds (AOT for PVC, noMPI) for the gate: usage build_variant.sh SRC OUT NAME [configure args] ; ALF_FLAGS_EXT
# from the environment (e.g. ALF_FLAGS_EXT=-O1 for a layout-shaken build: it follows -O3 on the command line).
src=$1; out=$2; name=$3; shift 3
rm -rf "$out/src_$name"; cp -r "$src" "$out/src_$name"; cd "$out/src_$name"
export ALF_INTEL_AOT=pvc
source ./configure.sh SUPERMUCNG2 noMPI GPU=intel "$@" NO-INTERACTIVE > "$out/configure_$name.log" 2>&1
make -j16 program > "$out/make_$name.log" 2>&1
cp Prog/ALF.out "$out/ALF-$name.out" && echo "built $out/ALF-$name.out" || echo "build $name failed"

#!/bin/bash
# GPU=none / CUDA identity of alf_la_mod.F90 incl. __LINE__: base (git rev, default HEAD) vs the working copy.
# usage (laptop, repo root): docs/gpu/scripts/intel/t302a_ppcheck.sh [BASEREV]   (needs gfortran)
BASE=${1:-HEAD}; D=$(mktemp -d); F=Libraries/Modules/alf_la_mod.F90
git show $BASE:$F > $D/base.F90; cp $F $D/cand.F90
tot=0
for m in "-UNONE" "-DMPI -DHDF5 -DHDF5_ZLIB" "-DHDF5 -DHDF5_ZLIB" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTABLOG" "-DMPI -DHDF5 -DHDF5_ZLIB -DQR_PREPIVOT" "-DMPI -DTEMPERING" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTAB3" "-DMPI -DPARALLEL_PARAMS" "-DMPI -DOBS_LEGACY" "-DALF_LA_CUDA -DALF_GPU" "-DALF_LA_HOST"; do
  for v in base cand; do (cd $D && gfortran -E -cpp -P $m $v.F90 2>$v.err | sed -E 's/"([A-Za-z_0-9]+\.F90)",[0-9]+/"\1",N/' > $v.i); [ -s $D/$v.err ] && { echo "CPP ERR $v [$m]"; head -3 $D/$v.err; }; done
  k=$(diff -B $D/base.i $D/cand.i | wc -l); echo "[$m] differing diff lines: $k"; tot=$((tot+k))
done
echo "TOTAL differing: $tot"; rm -rf $D

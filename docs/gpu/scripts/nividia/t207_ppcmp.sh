#!/bin/bash
# Ticket 207, nividia: GPU=none identity of the changed sources: preprocessed (gfortran -E -cpp -P) main.F90, tau_m_mod.F90,
# control_mod.F90 of the working tree (src_new) vs HEAD (src_head), in identical relative locations (so that __FILE__/__LINE__
# agree), for several macro sets; differences are listed with diff -B (blank lines ignored, as in ticket 206).
# Needs the generated headers (Hamiltonians_interface.h, Hamiltonians_case.h, git.h): taken from tree_nv_none/Prog (built by t207_build.sh).
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp
. ~/alf-gpu/env/gnu.sh
: "${HOME:?}"
for v in head new; do
  rm -rf -- "$W/cmp/$v"; mkdir -p "$W/cmp/$v/Prog"
  src=$W/src_$v
  cp "$W"/tree_nv_none/Prog/*.h "$W/cmp/$v/Prog/" 2>/dev/null
  for f in main.F90 tau_m_mod.F90 control_mod.F90; do cp "$src/Prog/$f" "$W/cmp/$v/Prog/$f"; done
done
tot=0
for m in "-DMPI -DHDF5 -DHDF5_ZLIB" "-DHDF5 -DHDF5_ZLIB" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTABLOG" "-DMPI -DHDF5 -DHDF5_ZLIB -DQR_PREPIVOT" "-DMPI -DTEMPERING" "-DMPI -DHDF5 -DHDF5_ZLIB -DSTAB3"; do
  for f in main tau_m_mod control_mod; do
    ( cd "$W/cmp/head/Prog" && gfortran -E -cpp -P $m $f.F90 > "$W/cmp/head_$f.i" 2> "$W/cmp/head_$f.err" )
    ( cd "$W/cmp/new/Prog"  && gfortran -E -cpp -P $m $f.F90 > "$W/cmp/new_$f.i"  2> "$W/cmp/new_$f.err" )
    n=$(diff -B "$W/cmp/head_$f.i" "$W/cmp/new_$f.i" | wc -l)
    # __LINE__ of runtime-error calls ("file.F90",NNN) shifts when lines are added above them: compare again with the numbers removed
    sed -E 's/"([A-Za-z_0-9]+\.F90)",[0-9]+/"\1",N/' "$W/cmp/head_$f.i" > "$W/cmp/head_$f.n"
    sed -E 's/"([A-Za-z_0-9]+\.F90)",[0-9]+/"\1",N/' "$W/cmp/new_$f.i" > "$W/cmp/new_$f.n"
    n2=$(diff -B "$W/cmp/head_$f.n" "$W/cmp/new_$f.n" | wc -l)
    echo "$f [$m]: $(wc -l < "$W/cmp/new_$f.i") lines, diff -B lines: $n, after removing the __LINE__ numbers: $n2, preprocessor errors: $(cat "$W/cmp/new_$f.err" | wc -l)"
    tot=$((tot+n2))
  done
done
echo "total differing diff lines (after removing the __LINE__ numbers): $tot"

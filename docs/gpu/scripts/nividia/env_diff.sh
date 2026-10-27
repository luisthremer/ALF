#!/bin/bash
# Ticket 201, acceptance 1: GPU=none identity of configure.sh.
# usage: env_diff.sh TREE OUTDIR
#   TREE must contain the NEW configure.sh and HDF5/<compiler> builds (reused, never rebuilt here).
#   The unmodified script is taken from 'git show <ORIG_REV>:configure.sh' (ORIG_REV, default 54ef3173) -- pass a
#   ready file via ORIG_CONFIGURE=/path to override. It is copied into TREE as configure_orig.sh (it uses $PWD files).
# For every MACHINE x MODE x STAB x HDF5 x Devel it dumps `env | grep '^ALF_' | sort` after sourcing
#   orig / new (no GPU=) / new GPU=none, and diffs the three. Combos where the ORIGINAL fails are listed as SKIP.
TREE=$1; OUT=$2; mkdir -p "$OUT"
[ -n "$TMPDIR" ] || { echo "set TMPDIR"; exit 1; }
cd "$TREE" || exit 1
if [ -n "$ORIG_CONFIGURE" ]; then cp "$ORIG_CONFIGURE" configure_orig.sh; else
  git -C "$TREE" show ${ORIG_REV:-54ef3173}:configure.sh > configure_orig.sh; fi
dump() { # dump FILE CONFIGURE ARGS...  (fresh shell state each time via subshell)
  local f=$1 c=$2; shift 2
  ( . "./$c" "$@" NO-INTERACTIVE >"$f.log" 2>&1 || exit 1; env | grep '^ALF_' | sort > "$f" ) ; return $?
}
printf "%-10s %-6s %-5s %-5s %-5s %s\n" MACHINE MODE STAB HDF5 DEVEL RESULT | tee "$OUT/table.txt"
for M in GNU IntelLLVM PGI; do
 ( case $M in GNU) . ~/alf-gpu/env/gnu.sh;; IntelLLVM) . ~/alf-gpu/env/intel.sh;; PGI) . ~/alf-gpu/env/nvhpc.sh;; esac
  for MODE in noMPI MPI; do for STAB in default LOG; do for H in nohdf5 HDF5; do for D in nodevel Devel; do
    A=(); [ $STAB = LOG ] && A+=(LOG); [ $H = HDF5 ] && A+=(HDF5); [ $D = Devel ] && A+=(Devel)
    n="$OUT/${M}_${MODE}_${STAB}_${H}_${D}"
    if ! dump "$n.orig" configure_orig.sh $M $MODE "${A[@]}"; then r="SKIP (unmodified configure fails)"
    else
      dump "$n.new" configure.sh $M $MODE "${A[@]}" ; dump "$n.none" configure.sh $M $MODE "${A[@]}" GPU=none
      if [ -s "$n.orig" ] && cmp -s "$n.orig" "$n.new" && cmp -s "$n.orig" "$n.none"; then r="IDENTICAL ($(wc -l < $n.orig) vars)"; else r="DIFFERENT"; fi
    fi
    printf "%-10s %-6s %-5s %-5s %-5s %s\n" $M $MODE $STAB $H $D "$r" | tee -a "$OUT/table.txt"
  done; done; done; done )
done

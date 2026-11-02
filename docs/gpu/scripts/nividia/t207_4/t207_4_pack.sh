#!/bin/bash
# Ticket 207-4 (development only), nividia: lossless packing of the KEPT dumps of processed runs (disk: KEEP 11 keeps 22 dumps = 131 MB per
# run; 16+ runs would exceed the 1.5 GB per-ticket evidence limit and / was at 6-7 GB free).  zstd -3 (measured ~32 % of the size).
#   pack   RUNDIR...  every qref/q*.bin -> q*.bin.zst; the .bin is removed only after zstd -dc of the .zst gives the sha256 that
#                     qref_manifest.txt lists for it (the processing record); nothing else is touched
#   unpack RUNDIR...  the reverse (before t207_4_gate.sh re-processes the kept dumps); the .zst is removed only after the restored .bin
#                     has the manifest sha256
# Only directories under ~/alf-gpu/work/t207-4/runs/ that hold a passing qref_proc.json are accepted.
W=$HOME/alf-gpu/work/t207-4; export TMPDIR=$W/tmp
op=$1; shift
for R in "$@"; do
  R=$(cd "$R" 2>/dev/null && pwd) || { echo "no dir"; exit 1; }
  case "$R" in $W/runs/*) ;; *) echo "refusing: $R"; exit 1;; esac
  grep -q '"check_rc": 0' $R/qref_proc.json 2>/dev/null || { echo "refusing: $R not processed"; exit 1; }
  cd $R || exit 1; n=0
  case $op in
  pack)
    for f in qref/q*.bin; do [ -e "$f" ] || continue
      want=$(awk -v f="$f" '$2 == f {print $1}' qref_manifest.txt); [ -n "$want" ] || { echo "FAIL: $f not in manifest"; exit 2; }
      zstd -q -3 -T1 -f "$f" -o "$f.zst" || exit 2
      got=$(zstd -dc "$f.zst" | sha256sum | cut -d' ' -f1)
      [ "$got" = "$want" ] && rm -f -- "$f" && n=$((n+1)) || { echo "FAIL: $f.zst does not restore the manifest sha256"; rm -f -- "$f.zst"; exit 2; }
    done ;;
  unpack)
    for z in qref/q*.bin.zst; do [ -e "$z" ] || continue; f=${z%.zst}
      want=$(awk -v f="$f" '$2 == f {print $1}' qref_manifest.txt)
      zstd -q -d -f "$z" -o "$f" || exit 2
      got=$(sha256sum "$f" | cut -d' ' -f1)
      [ "$got" = "$want" ] && rm -f -- "$z" && n=$((n+1)) || { echo "FAIL: $f restored with sha256 $got, manifest $want"; rm -f -- "$f"; exit 2; }
    done ;;
  *) echo "usage: t207_4_pack.sh pack|unpack RUNDIR..."; exit 1 ;;
  esac
  echo "$op $R: $n dumps, $(du -sh qref | cut -f1)"
done

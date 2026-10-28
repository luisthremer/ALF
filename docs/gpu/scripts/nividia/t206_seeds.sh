#!/bin/bash
# Ticket 206, nividia: make K copies of a prepared run directory (runs/base/BASE) with rotated seeds (independent streams):
#   runs/base/BASE_s1 ... BASE_sK ; seeds rotated by 37*k lines (pyALF's only_prep writes the same seeds file for every case).
#   usage: t206_seeds.sh BASE K
W=$HOME/alf-gpu/work/t206
B=$1; K=$2
case "$B" in ""|*..*|*/*) echo "unsafe name"; exit 1;; esac
: "${HOME:?}"
for k in $(seq 1 "$K"); do
  D="$W/runs/base/${B}_s$k"
  rm -rf -- "$D"; mkdir -p "$D"; cp -r "$W/runs/base/$B/." "$D/"
  n=$(wc -l < "$D/seeds"); r=$(( (37*k) % n ))
  { tail -n +$((r+1)) "$D/seeds"; head -n $r "$D/seeds"; } > "$D/seeds.new" && mv "$D/seeds.new" "$D/seeds"
done
echo "made ${B}_s1..s$K"

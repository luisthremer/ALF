#!/bin/bash
# Phase 0 CPU reference batch. Each run < 15 min; 16 MPI ranks unless noted.
cd ~/alf-gpu
echo "== Intel test_vs_ed"; timeout 1000 ./run_ed.sh intel-DEFAULT hubbard_finT_chain 2>&1 | grep -vE "^Prepare directory|^\./|No orbital"
ref() { # tag outdir config
  t=${1%%-*}; . ~/alf-gpu/env/$t.sh; cd ~/alf-gpu/gpu_reference
  rm -rf $2/$3; python prep.py $2 $3 --alf-dir ~/alf-gpu/src/ALF-$t >/dev/null 2>&1
  timeout 900 bash run.sh ~/alf-gpu/bin/ALF-$1.out 16 $2/$3
  python summarize.py $2/$3 2>&1 | tail -1
}
for c in hub6x6_cb hub6x6_dense hub12x12_cb hub12x12_dense hub20x20_cb hub20x20_dense; do
  echo "== ref gnu-DEFAULT $c"; ref gnu-DEFAULT ~/alf-gpu/runs/ref_gnu-DEFAULT $c
done
for tag in intel-DEFAULT gnu-LOG; do
  echo "== check $tag hub12x12_cb"; ref $tag ~/alf-gpu/runs/ref_$tag hub12x12_cb
  python ~/alf-gpu/gpu_reference/compare.py ~/alf-gpu/runs/ref_gnu-DEFAULT/hub12x12_cb/reference.json ~/alf-gpu/runs/ref_$tag/hub12x12_cb/reference.json
done
echo BATCH_DONE

#!/bin/bash
# Ticket 207-3 (development only), on the LAPTOP from the worktree root: rsync the working tree to nividia ~/alf-gpu/work/t207-3/src/
# (source state HEAD + sha256 of the uncommitted diff in src/.d5_srcstate, copied into every build), optionally base_src = git archive
# of HEAD (the code before the 207-3 edits), the scripts of this ticket and the QREF pipeline scripts t207_d5_* (used unchanged).
#   usage: [T207_3D=t207-3b] t207_3_sync.sh [base|pipe]      (T207_3D: work dir name under ~/alf-gpu/work, default t207-3)
#   pipe (fix round t207-3c): g3/ = the t207_d5_* scripts of gpu-port-207-2 305488ec (reader of the runs processed by it, spec version 2)
#        and g5/ = those of 4118a1de (the approved pipeline, spec version 3), both by git archive; refuses if the directory exists
set -e
D=${T207_3D:-t207-3}; case "$D" in t207-3|t207-3[a-z]) ;; *) echo "bad T207_3D"; exit 1;; esac
ssh nividia "mkdir -p ~/alf-gpu/work/$D/{src,base_src,tmp,logs,bin,runs}"
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD | shasum -a 256 | cut -c1-64)"; git status --short; } > .d5_srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/$D/src/
rm -f .d5_srcstate
rsync -a docs/gpu/scripts/nividia/t207_3/ docs/gpu/scripts/nividia/t207_d5_* nividia:alf-gpu/work/$D/
if [ "$1" = base ]; then
  ssh nividia "rm -rf -- ~/alf-gpu/work/$D/base_src && mkdir -p ~/alf-gpu/work/$D/base_src"
  git archive HEAD | ssh nividia "tar -x -C ~/alf-gpu/work/$D/base_src"
  echo "head $(git rev-parse HEAD) (git archive, no diff)" | ssh nividia "cat > ~/alf-gpu/work/$D/base_src/.d5_srcstate"
fi
if [ "$1" = pipe ]; then
  for gv in g3:305488ec g5:4118a1de; do
    g=${gv%%:*}; c=${gv##*:}
    if ssh nividia "test -e ~/alf-gpu/work/$D/$g"; then echo "refusing: $D/$g exists"; continue; fi
    git archive $c docs/gpu/scripts/nividia | ssh nividia "mkdir -p ~/alf-gpu/work/$D/$g && tar -x --strip-components=4 -C ~/alf-gpu/work/$D/$g --wildcards 'docs/gpu/scripts/nividia/t207_d5_*' && cd ~/alf-gpu/work/$D/$g && echo $g $c qparse \$(sha256sum t207_d5_qparse.py | cut -c1-16) proc \$(sha256sum t207_d5_proc.sh | cut -c1-16)"
  done
fi
echo synced

#!/bin/bash
# Ticket 207-1 (development only), on the LAPTOP from the checkout root (worktree t207r): rsync the working tree to nividia
# ~/alf-gpu/work/t207-1/src/ (and base_src = git archive 0df7217c = gpu-port before the merge), the scripts of this ticket, and
# changed_files.txt (git diff --name-status 0df7217c of the working tree).  The source state (HEAD, sha256 of the uncommitted diff)
# goes to src/.srcstate.     usage: t207_1_sync.sh SCRIPTDIR [base]
set -e
SD=$1
ssh nividia 'mkdir -p ~/alf-gpu/work/t207-1/{src,base_src,tmp,logs,bin,runs}'
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD | shasum -a 256 | cut -c1-64)"; git status --short; } > .srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation \
      --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
      --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ nividia:alf-gpu/work/t207-1/src/
rm -f .srcstate
git diff --name-status 0df7217c | ssh nividia 'cat > ~/alf-gpu/work/t207-1/changed_files.txt'
rsync -a $SD/t207_1_* $SD/t207_d5_qref.F90 docs/gpu/scripts/nividia/t207_firstdiv.py docs/gpu/scripts/nividia/t207_diagstat.py nividia:alf-gpu/work/t207-1/
for a in "${@:2}"; do
  if [ "$a" = base ]; then
    ssh nividia 'rm -rf -- ~/alf-gpu/work/t207-1/base_src && mkdir -p ~/alf-gpu/work/t207-1/base_src'
    git archive 0df7217c | ssh nividia 'tar -x -C ~/alf-gpu/work/t207-1/base_src'
  fi
done
echo synced

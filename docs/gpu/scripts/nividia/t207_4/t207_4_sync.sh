#!/bin/bash
# Ticket 207-4 (development only), on the LAPTOP from the worktree root: put the β = 10 calibration inputs on nividia ~/alf-gpu/work/t207-4/.
#   src/  = git archive of SRC (default 428576dd = gpu-port-207-3: sweep output fault, localized fault), no working-tree diff;
#           src/.d5_srcstate names the commit (copied into every build record)
#   g5/   = the approved QREF processing / gate tools t207_d5_* of gpu-port-207-2 4118a1de (git archive, unchanged; g5/SOURCE.txt)
#   ./    = the scripts of this ticket (docs/gpu/scripts/nividia/t207_4/)
#   usage: t207_4_sync.sh [src|g5|scripts]...   (default: scripts)
set -e
SRC=${SRC:-428576dd}; GATE=4118a1de
R=alf-gpu/work/t207-4
ssh nividia "mkdir -p ~/$R/{src,g5,tmp,logs,bin,runs,gate}"
for what in "${@:-scripts}"; do case $what in
src)
  ssh nividia "test -z \"\$(ls -A ~/$R/src)\"" || { echo "refusing: $R/src not empty"; exit 1; }
  git archive "$SRC" | ssh nividia "tar -x -C ~/$R/src"
  echo "head $(git rev-parse "$SRC") (git archive, no diff)" | ssh nividia "cat > ~/$R/src/.d5_srcstate" ;;
g5)
  ssh nividia "test -z \"\$(ls -A ~/$R/g5)\"" || { echo "refusing: $R/g5 not empty"; exit 1; }
  git archive "$GATE" docs/gpu/scripts/nividia | ssh nividia "tar -x -C ~/$R/g5 --strip-components=4 --wildcards 'docs/gpu/scripts/nividia/t207_d5_*'"
  echo "git archive $(git rev-parse $GATE) docs/gpu/scripts/nividia/t207_d5_* (gpu-port-207-2, approved gate)" | ssh nividia "cat > ~/$R/g5/SOURCE.txt" ;;
scripts)
  rsync -a docs/gpu/scripts/nividia/t207_4/ nividia:$R/ ;;
*) echo "bad item $what"; exit 1 ;;
esac; done
echo synced "${@:-scripts}"

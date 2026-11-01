#!/bin/bash
# Ticket 302c-2 (development only), on the LAPTOP: copy the candidate source and the QREF pipeline to intel:~/alf-gpu/work/t302c-2/.
#   src/            the worktree (branch gpu-port-302c-2; code = 302c-1's 0e656fc8), with .srcstate (HEAD + sha256 of the diff)
#   t207_d5_{proc.sh,qparse.py,agg.py,qref.F90}   the reviewed 207-2 pipeline, taken unchanged from commit 5d876291 (branch gpu-port-207-2),
#                   not from this branch (the copies in docs/gpu/scripts/nividia/ of this branch are older); pipeline.sha256 lists them
#   t302c2_*        the scripts of this directory
# usage: t302c2_sync.sh [PIPELINE_COMMIT]   (default 5d876291)
set -e
WT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)
PC=${1:-5d876291}
cd "$WT"
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD -- Prog Libraries | shasum -a 256 | cut -c1-64)"; git status --short -- Prog Libraries; } > .srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation --exclude doxygen \
  --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' \
  --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ intel:alf-gpu/work/t302c-2/src/
rm -f .srcstate
T=$(mktemp -d "${TMPDIR:-/tmp}/t302c2_sync.XXXXXX")
git archive -o "$T/p.tar" "$PC" docs/gpu/scripts/nividia/t207_d5_proc.sh docs/gpu/scripts/nividia/t207_d5_qparse.py \
  docs/gpu/scripts/nividia/t207_d5_agg.py docs/gpu/scripts/nividia/t207_d5_qref.F90
tar -xf "$T/p.tar" -C "$T"
( cd "$T/docs/gpu/scripts/nividia" && shasum -a 256 t207_d5_* > pipeline.sha256 && echo "commit $(git -C "$WT" rev-parse "$PC")" >> pipeline.sha256 )
rsync -a "$T/docs/gpu/scripts/nividia/" intel:alf-gpu/work/t302c-2/
rsync -a "$WT/docs/gpu/scripts/intel/t302c/" intel:alf-gpu/work/t302c-2/
rm -rf -- "$T"
echo "synced (pipeline from $PC)"

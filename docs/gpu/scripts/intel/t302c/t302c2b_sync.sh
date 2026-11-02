#!/bin/bash
# Ticket 302c-2 round 2 (development only), on the LAPTOP: copy the candidate source and the QREF pipeline to intel:~/alf-gpu/work/t302c-2b/.
#   src/            the worktree (branch gpu-port-302c-2, HEAD 08559f0b = 302c-1 device tau + 207-3 sweep output fault), with .srcstate
#   t207_d5_{proc.sh,qparse.py,agg.py,qref.F90,gate.py}   the approved 207-2 pipeline, unchanged from commit 4118a1de (gate fix round 5,
#                   branch gpu-port-207-2: proc.sh sources gnu.sh only if it exists, which intel lacks); pipeline.sha256 lists them
#   t207_3_cal.py   207-3's calibration table script (unchanged, from this branch), the generator of the NVIDIA beta = 4 numbers
#   t207_firstdiv.py  same-seed first-difference tool (unchanged, docs/gpu/scripts/nividia/ of this branch)
#   t302c2_*, t302c2b_*  the scripts of this directory (round 1 scripts take T302C2_W=~/alf-gpu/work/t302c-2b)
# Round 1's ~/alf-gpu/work/t302c-2/ is not touched.   usage: t302c2b_sync.sh [PIPELINE_COMMIT]   (default 4118a1de)
set -e
WT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)
PC=${1:-4118a1de}
cd "$WT"
ssh intel 'mkdir -p ~/alf-gpu/work/t302c-2b/{src,tmp,logs,bin,runs,tab}'
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD -- Prog Libraries | shasum -a 256 | cut -c1-64)"; git status --short -- Prog Libraries; } > .srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation --exclude doxygen \
  --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' \
  --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ intel:alf-gpu/work/t302c-2b/src/
rm -f .srcstate
T=$(mktemp -d "${TMPDIR:-/tmp}/t302c2b_sync.XXXXXX")
git archive -o "$T/p.tar" "$PC" docs/gpu/scripts/nividia/t207_d5_proc.sh docs/gpu/scripts/nividia/t207_d5_qparse.py \
  docs/gpu/scripts/nividia/t207_d5_agg.py docs/gpu/scripts/nividia/t207_d5_qref.F90 docs/gpu/scripts/nividia/t207_d5_gate.py
tar -xf "$T/p.tar" -C "$T"
( cd "$T/docs/gpu/scripts/nividia" && shasum -a 256 t207_d5_* > pipeline.sha256 && echo "commit $(git -C "$WT" rev-parse "$PC")" >> pipeline.sha256 )
rsync -a "$T/docs/gpu/scripts/nividia/" intel:alf-gpu/work/t302c-2b/
rsync -a "$WT/docs/gpu/scripts/nividia/t207_3/t207_3_cal.py" "$WT/docs/gpu/scripts/nividia/t207_firstdiv.py" intel:alf-gpu/work/t302c-2b/
rsync -a "$WT/docs/gpu/scripts/intel/t302c/" intel:alf-gpu/work/t302c-2b/
rm -rf -- "$T"
echo "synced (pipeline from $PC)"

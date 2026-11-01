#!/bin/bash
# Ticket 302c-1 (development only), on the LAPTOP: rsync the worktree (candidate src) to intel:~/alf-gpu/work/t302c/src/.
cd "${REPO:?set REPO to the checkout}" || exit 1
{ echo "head $(git rev-parse HEAD)"; echo "diff_sha256 $(git diff HEAD | shasum -a 256 | cut -c1-64)"; git status --short; } > .srcstate
rsync -a --delete --exclude .git --exclude .alf_gpu_stamp --exclude HDF5 --exclude gpu_bench --exclude Documentation --exclude doxygen --exclude Mathematica --exclude .DS_Store --exclude '*.o' --exclude '*.mod' --exclude '*.a' --exclude '*.so' --exclude '*.out' --exclude .deps.mk --exclude tests --exclude '*.pdf' ./ intel:alf-gpu/work/t302c/src/
rm -f .srcstate
echo synced

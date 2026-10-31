#!/bin/bash
# T0 development only (ticket 302b): ship the WORKING TREE (tracked files, current content) and the base 71af3e07 to a
# dev candidate root on intel, mimicking acc_ship.sh's layout so that the harness scripts work on it:
#   intel:~/alf-gpu/work/t302b/dev/cand/dev/<TAG>/{src,base_src}{,.tar}, .alf_candidate, tools/ (acc_build/acc_ctest with the
#   candidate root moved into dev/).  Not a candidate in the sense of VERIFICATION.md (no commit): T0 evidence only.
#   usage (laptop, repo root): docs/gpu/scripts/intel/t302b_sync.sh TAG [HOST]  (HOST intel or nividia)
set -e
TAG=$1; [ -n "$TAG" ] || { echo usage; exit 2; }
HOST=${2:-intel}; BASE=71af3e076322fe5115451de1b73482702a6b17e0
W=$(git rev-parse --show-toplevel); cd "$W"
T=$(mktemp -d "${TMPDIR:-/tmp}/t302b.XXXXXX")
git ls-files > "$T/files"
COPYFILE_DISABLE=1 tar -cf "$T/src.tar" -T "$T/files"
git archive --format=tar "$BASE" > "$T/base_src.tar"
SH=$(shasum -a 256 "$T/src.tar" | cut -d' ' -f1); BH=$(shasum -a 256 "$T/base_src.tar" | cut -d' ' -f1)
{ echo "sha devtree$SH"; echo "tree worktree-of-$(git rev-parse --short HEAD)"; echo "tar_sha256 $SH"; echo "ticket 302b"
  echo "base_sha $BASE"; echo "base_tar_sha256 $BH"; echo "who dev"; echo "shipped $(date -u +%FT%TZ)"; } > "$T/.alf_candidate"
R=alf-gpu/work/t302b/dev/cand/dev/$TAG
ssh $HOST "mkdir -p ~/alf-gpu/work/t302b/dev/tmp ~/alf-gpu/work/t302b/dev/cand/dev && mkdir ~/$R" 
scp -q "$T/src.tar" "$T/base_src.tar" "$T/.alf_candidate" "$HOST:$R/"
ssh $HOST "cd ~/$R && mkdir src base_src tools && tar -xf src.tar -C src && tar -xf base_src.tar -C base_src && 
  for f in acc_build.sh acc_ctest.sh; do sed -e 's#C=\$HOME/alf-gpu/work/t\$TICKET/cand/\$WHO/\$S12#C=\$HOME/alf-gpu/work/t302b/dev/cand/\$WHO/\$S12#' -e 's#TMPDIR=\$HOME/alf-gpu/work/t\$TICKET/tmp#TMPDIR=\$HOME/alf-gpu/work/t302b/dev/tmp#' src/testsuite/gpu_reference/acceptance/\$f > tools/\$f; done;
  grep -c 'work/t302b/dev' tools/acc_build.sh tools/acc_ctest.sh"
echo "shipped devtree$SH as $TAG"

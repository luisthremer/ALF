#!/bin/bash
# Acceptance harness, LAPTOP side: ship a frozen candidate commit to a test machine with `git archive`.
#   usage: acc_ship.sh HOST TICKET CANDIDATE_SHA [BASE_SHA]
#   Result on HOST: ~/alf-gpu/work/t<TICKET>/cand/<sha12>/src/ (read-only) and .alf_candidate (sha, tree, tar sha256),
#   and, if BASE_SHA is given, cand/<sha12>/base_src/ for the GPU=none comparison.
# Never ships a working tree: only committed content. Exit code 0 only if every hash check passes.
set -euo pipefail
HOST=$1; TICKET=$2; SHA=$(git rev-parse --verify "$3^{commit}"); BASE=${4:-}
S12=${SHA:0:12}; TREE=$(git rev-parse "$SHA^{tree}")
case "$TICKET" in *[!0-9a-z]*|"") echo "bad ticket id"; exit 2;; esac
TMP=$(mktemp -d "${TMPDIR:-/tmp}/accship.XXXXXX"); trap 'rm -rf -- "$TMP"' EXIT
ship() {  # SHA DESTNAME
  local sha=$1 dst=$2 tar=$TMP/$2.tar h
  git archive --format=tar "$sha" > "$tar"
  h=$(shasum -a 256 "$tar" | cut -d' ' -f1)
  scp -q "$tar" "$HOST:alf-gpu/work/t$TICKET/cand/$S12/$dst.tar"
  ssh "$HOST" "bash -c 'set -e; cd ~/alf-gpu/work/t$TICKET/cand/$S12; [ ! -e $dst ] || { echo \"refusing: $dst exists\"; exit 3; };
     h=\$(sha256sum $dst.tar | cut -d\" \" -f1); [ \"\$h\" = \"$h\" ] || { echo \"TAR HASH MISMATCH\"; exit 4; };
     mkdir $dst; tar -xf $dst.tar -C $dst; chmod -R a-w $dst; echo shipped $dst'"
  echo "$h"
}
ssh "$HOST" "mkdir -p ~/alf-gpu/work/t$TICKET/cand/$S12"
H=$(ship "$SHA" src | tail -1)
{ echo "sha $SHA"; echo "tree $TREE"; echo "tar_sha256 $H"; echo "ticket $TICKET"; echo "shipped $(date -u +%FT%TZ)"
  if [ -n "$BASE" ]; then B=$(git rev-parse --verify "$BASE^{commit}"); echo "base_sha $B"; echo "base_tar_sha256 $(ship "$B" base_src | tail -1)"; fi
} > "$TMP/.alf_candidate"
scp -q "$TMP/.alf_candidate" "$HOST:alf-gpu/work/t$TICKET/cand/$S12/.alf_candidate"
ssh "$HOST" "cat ~/alf-gpu/work/t$TICKET/cand/$S12/.alf_candidate"

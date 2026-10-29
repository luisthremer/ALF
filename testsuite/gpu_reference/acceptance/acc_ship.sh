#!/bin/bash
# Acceptance harness, LAPTOP side: ship a frozen candidate commit to a test machine with `git archive`.
#   usage: acc_ship.sh HOST TICKET CANDIDATE_SHA BASE_SHA PREREG_SHA
#     PREREG_SHA: the commit in which testsuite/gpu_reference/acceptance/<TICKET>.json was pre-registered (recorded in the
#     ticket and DECISIONS.md). The candidate must carry that file, the harness (acceptance/), compare2.py and configs*.py
#     unchanged since PREREG_SHA: the tools that judge a candidate are not the candidate's to change.
#   Result on HOST: ~/alf-gpu/work/t<TICKET>/cand/<sha12>/{src,base_src} (read-only), {src,base_src}.tar,
#   changed_files.txt (git diff --name-status --no-renames base..candidate) and .alf_candidate (read-only):
#   sha, tree, tar_sha256, base_sha, base_tar_sha256, prereg_sha, acceptance_sha256.
# Every hash is checked on the machine after the copy. Exit code 0 only if every step and check succeeds.
set -euo pipefail
export LC_ALL=C
[ $# = 5 ] || { echo "usage: acc_ship.sh HOST TICKET CANDIDATE_SHA BASE_SHA PREREG_SHA"; exit 2; }
HOST=$1; TICKET=$2
case "$TICKET" in *[!0-9a-z]*|"") echo "bad ticket id"; exit 2;; esac
SHA=$(git rev-parse --verify "$3^{commit}"); BASE=$(git rev-parse --verify "$4^{commit}"); PRE=$(git rev-parse --verify "$5^{commit}")
S12=${SHA:0:12}; TREE=$(git rev-parse "$SHA^{tree}")
ACCF=testsuite/gpu_reference/acceptance/$TICKET.json
git cat-file -e "$PRE:$ACCF" || { echo "$ACCF does not exist in the pre-registration commit $PRE"; exit 3; }
if ! git diff --quiet "$PRE" "$SHA" -- "$ACCF" testsuite/gpu_reference/acceptance testsuite/gpu_reference/compare2.py \
      testsuite/gpu_reference/configs.py testsuite/gpu_reference/configs_stat.py; then
  echo "REFUSED: the acceptance file or the judging tools differ from the pre-registration commit $PRE:"
  git diff --stat "$PRE" "$SHA" -- "$ACCF" testsuite/gpu_reference/acceptance testsuite/gpu_reference/compare2.py \
      testsuite/gpu_reference/configs.py testsuite/gpu_reference/configs_stat.py
  exit 3
fi
ACCH=$(git show "$PRE:$ACCF" | shasum -a 256 | cut -d' ' -f1)
TMP=$(mktemp -d "${TMPDIR:-/tmp}/accship.XXXXXX"); trap 'rm -rf -- "$TMP"' EXIT
R=alf-gpu/work/t$TICKET/cand/$S12
ssh "$HOST" "mkdir -p ~/alf-gpu/work/t$TICKET/cand && mkdir ~/$R" || { echo "refusing: $R exists on $HOST"; exit 3; }
ship() {  # SHA NAME -> writes $TMP/NAME.sha256; returns nonzero on any failure
  local sha=$1 dst=$2 tar=$TMP/$2.tar h
  git archive --format=tar "$sha" > "$tar" || return 1
  h=$(shasum -a 256 "$tar" | cut -d' ' -f1) || return 1
  scp -q "$tar" "$HOST:$R/$dst.tar" || return 1
  ssh "$HOST" "bash -c 'set -e; cd ~/$R; h=\$(sha256sum $dst.tar | cut -d\" \" -f1); [ \"\$h\" = \"$h\" ] || { echo \"TAR HASH MISMATCH $dst\"; exit 4; };
     mkdir $dst; tar -xf $dst.tar -C $dst; chmod -R a-w $dst $dst.tar'" || return 1
  echo "$h" > "$TMP/$dst.sha256"
}
ship "$SHA" src || { echo "shipping the candidate failed"; exit 4; }
ship "$BASE" base_src || { echo "shipping the base failed"; exit 4; }
git diff --name-status --no-renames "$BASE" "$SHA" > "$TMP/changed_files.txt"
{ echo "sha $SHA"; echo "tree $TREE"; echo "tar_sha256 $(cat "$TMP/src.sha256")"; echo "ticket $TICKET"
  echo "base_sha $BASE"; echo "base_tar_sha256 $(cat "$TMP/base_src.sha256")"
  echo "prereg_sha $PRE"; echo "acceptance_sha256 $ACCH"; echo "shipped $(date -u +%FT%TZ)"; } > "$TMP/.alf_candidate"
scp -q "$TMP/.alf_candidate" "$TMP/changed_files.txt" "$HOST:$R/"
cat > "$TMP/postcheck.sh" <<'PC'
set -e; cd "$1"
for k in src base_src; do
  key=tar_sha256; [ $k = base_src ] && key=base_tar_sha256
  w=$(awk -v k=$key '$1==k {print $2}' .alf_candidate); h=$(sha256sum $k.tar | cut -d' ' -f1)
  [ -n "$w" ] && [ "$h" = "$w" ] || { echo "POST-CHECK MISMATCH $k"; exit 4; }
done
chmod a-w .alf_candidate changed_files.txt
cat .alf_candidate
PC
ssh "$HOST" "bash -s ~/$R" < "$TMP/postcheck.sh"

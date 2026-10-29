#!/bin/bash
# Acceptance harness, LAPTOP side: ship a frozen candidate commit to a test machine with `git archive`.
#   usage: acc_ship.sh HOST TICKET WHO CANDIDATE_SHA BASE_SHA PREREG_SHA
#     WHO: mgr | agent - each party gets its own candidate root, so T2 never reuses T1's binaries.
#     BASE_SHA: the accepted commit the ticket starts from; its harness is the approved one.
#     PREREG_SHA: the commit in which acceptance/<TICKET>*.json was pre-registered (recorded in the ticket and
#     DECISIONS.md). Checks: PREREG is a strict ancestor of the candidate; between BASE and PREREG the judging tools
#     (acceptance/, compare2.py, configs*.py) change only by adding <TICKET>*.json; between PREREG and the candidate they
#     do not change at all, and existing test files under testsuite/ are not modified or deleted unless the acceptance
#     file lists them in "allowed_test_changes".
#   Result on HOST: ~/alf-gpu/work/t<TICKET>/cand/<WHO>/<sha12>/{src,base_src} (read-only), {src,base_src}.tar,
#   changed_files.txt (git diff --name-status --no-renames base..candidate) and .alf_candidate (read-only):
#   sha, tree, tar_sha256, base_sha, base_tar_sha256, prereg_sha, acceptance_sha256.
# Every hash is checked on the machine after the copy. Exit code 0 only if every step and check succeeds.
set -euo pipefail
export LC_ALL=C
[ $# = 6 ] || { echo "usage: acc_ship.sh HOST TICKET WHO CANDIDATE_SHA BASE_SHA PREREG_SHA"; exit 2; }
HOST=$1; TICKET=$2; WHO=$3
case "$TICKET" in *[!0-9a-z]*|"") echo "bad ticket id"; exit 2;; esac
case "$WHO" in mgr|agent) ;; *) echo "WHO must be mgr or agent"; exit 2;; esac
SHA=$(git rev-parse --verify "$4^{commit}"); BASE=$(git rev-parse --verify "$5^{commit}"); PRE=$(git rev-parse --verify "$6^{commit}")
S12=${SHA:0:12}; TREE=$(git rev-parse "$SHA^{tree}")
ACCF=testsuite/gpu_reference/acceptance/$TICKET.json
TOOLS="testsuite/gpu_reference/acceptance testsuite/gpu_reference/compare2.py testsuite/gpu_reference/configs.py testsuite/gpu_reference/configs_stat.py"
git cat-file -e "$PRE:$ACCF" || { echo "$ACCF does not exist in the pre-registration commit $PRE"; exit 3; }
[ "$PRE" != "$SHA" ] && git merge-base --is-ancestor "$PRE" "$SHA" || { echo "REFUSED: PREREG must be a strict ancestor of the candidate"; exit 3; }
git merge-base --is-ancestor "$BASE" "$PRE" || { echo "REFUSED: BASE must be an ancestor of PREREG"; exit 3; }
ABASE=$(git show "$PRE:$ACCF" | python3 -c "import json,sys; print(json.load(sys.stdin).get('base_sha', ''))")
case "$ABASE" in ????????????*) ;; *) echo "REFUSED: the acceptance file's base_sha must have >= 12 hex digits (got '$ABASE')"; exit 3;; esac
[ "${BASE#$ABASE}" != "$BASE" ] || { echo "REFUSED: BASE $BASE is not the acceptance file's base_sha $ABASE"; exit 3; }
git fetch -q fork gpu-port 2>/dev/null || true
git merge-base --is-ancestor "$BASE" fork/gpu-port || { echo "REFUSED: BASE is not on fork/gpu-port (not an accepted state)"; exit 3; }
other=$(git diff --name-status --no-renames "$BASE" "$PRE" -- $TOOLS | grep -v -E "^A[[:space:]]testsuite/gpu_reference/acceptance/$TICKET(\.[a-z0-9]+)*\.json$" || true)
[ -z "$other" ] || { echo "REFUSED: between BASE and PREREG the judging tools changed beyond adding $TICKET*.json:"; echo "$other"; exit 3; }
allowed=$(git show "$PRE:$ACCF" | python3 -c "import json,sys; print(' '.join(json.load(sys.stdin).get('allowed_test_changes', [])))")
tch=$(git diff --name-status --no-renames "$PRE" "$SHA" -- testsuite | awk '$1=="M"||$1=="D" {print $2}' | grep -v -E '^testsuite/gpu_reference/acceptance/' || true)
for f in $tch; do case " $allowed " in *" $f "*) ;; *) echo "REFUSED: existing test file changed after pre-registration: $f"; exit 3;; esac; done
if ! git diff --quiet "$PRE" "$SHA" -- $TOOLS; then
  echo "REFUSED: the acceptance file or the judging tools differ from the pre-registration commit $PRE:"
  git diff --stat "$PRE" "$SHA" -- $TOOLS
  exit 3
fi
ACCH=$(git show "$PRE:$ACCF" | shasum -a 256 | cut -d' ' -f1)
PARTS=$(git ls-tree --name-only "$PRE" testsuite/gpu_reference/acceptance/ | grep -E "/$TICKET(\.[a-z0-9]+)?\.json$" || true)
TMP=$(mktemp -d "${TMPDIR:-/tmp}/accship.XXXXXX"); trap 'rm -rf -- "$TMP"' EXIT
R=alf-gpu/work/t$TICKET/cand/$WHO/$S12
ssh "$HOST" "mkdir -p ~/alf-gpu/work/t$TICKET/cand/$WHO && mkdir ~/$R" || { echo "refusing: $R exists on $HOST"; exit 3; }
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
  echo "prereg_sha $PRE"; echo "acceptance_sha256 $ACCH"; echo "who $WHO"
  for f in $PARTS; do echo "acceptance_part $(basename $f) $(git show "$PRE:$f" | shasum -a 256 | cut -d' ' -f1)"; done
  echo "shipped $(date -u +%FT%TZ)"; } > "$TMP/.alf_candidate"
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

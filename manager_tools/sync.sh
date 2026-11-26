#!/bin/sh
# Copy manager_tools/ from one branch to the other main lines, so the folder is the same on all of them.
#   ./manager_tools/sync.sh [SOURCE]      SOURCE: branch whose manager_tools/ is copied (default: current branch)
# Targets: Code-modularity, Code-modularity-CPU, Code-modularity-GPU (local branches; they must equal fork/<name>).
# Only manager_tools/ changes; one commit per target with [skip ci]; the working tree is not touched.
# If Code-modularity-GPU equalled Code-modularity-CPU before, it is fast-forwarded to the new CPU commit.
# Then pushes the changed branches to the remote "fork" (never origin).
set -eu
cd "$(git rev-parse --show-toplevel)"
src="${1:-$(git symbolic-ref --short HEAD)}"
srctree="$(git rev-parse "$src:manager_tools")"
git fetch -q fork
cpu_old="$(git rev-parse Code-modularity-CPU)"
gpu_old="$(git rev-parse Code-modularity-GPU)"
changed=""
for b in Code-modularity Code-modularity-CPU Code-modularity-GPU; do
  [ "$(git rev-parse "$b")" = "$(git rev-parse "fork/$b")" ] || { echo "$b differs from fork/$b: pull or push first" >&2; exit 1; }
  if [ "$b" = Code-modularity-GPU ] && [ "$gpu_old" = "$cpu_old" ] && [ "$(git rev-parse Code-modularity-CPU)" != "$cpu_old" ]; then
    git update-ref refs/heads/Code-modularity-GPU "$(git rev-parse Code-modularity-CPU)" "$gpu_old"
    changed="$changed $b"; echo "$b: fast-forwarded to Code-modularity-CPU"; continue
  fi
  old="$(git rev-parse "$b")"
  [ "$(git rev-parse "$b:manager_tools" 2>/dev/null || true)" = "$srctree" ] && { echo "$b: already the same"; continue; }
  idx="$(mktemp)"
  GIT_INDEX_FILE="$idx" git read-tree "$b"
  GIT_INDEX_FILE="$idx" git rm -r -q -f --cached --ignore-unmatch manager_tools
  GIT_INDEX_FILE="$idx" git read-tree --prefix=manager_tools/ "$srctree"
  tree="$(GIT_INDEX_FILE="$idx" git write-tree)"
  rm -f "$idx"
  c="$(git commit-tree "$tree" -p "$old" -m "manager_tools: sync from $src [skip ci]")"
  git update-ref "refs/heads/$b" "$c" "$old"
  changed="$changed $b"; echo "$b: $(git rev-parse --short "$c")"
done
# a checked-out target keeps its old files in the working tree until reset; say so instead of touching it
cur="$(git symbolic-ref --short HEAD 2>/dev/null || true)"
case " $changed " in *" $cur "*) echo "note: $cur is checked out here; run 'git reset --hard $cur' only if you have no local changes" ;; esac
[ -n "$changed" ] && git push fork $changed

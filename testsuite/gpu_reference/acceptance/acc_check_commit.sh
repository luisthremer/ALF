#!/bin/bash
# Acceptance harness, LAPTOP side: before committing, check that the commit equals the accepted candidate except docs/**.
#   usage: acc_check_commit.sh ACCEPTED_SHA [COMMIT=HEAD]
# Exit 0 only if every changed path between the two lies under docs/.
set -euo pipefail
A=$(git rev-parse --verify "$1^{commit}"); B=$(git rev-parse --verify "${2:-HEAD}^{commit}")
# --no-renames: a file moved into docs/ shows up at its old path too; docs/gpu/scripts/ is executable tooling
bad=$(git diff --name-only --no-renames "$A" "$B" | grep -v '^docs/' ; git diff --name-only --no-renames "$A" "$B" | grep '^docs/gpu/scripts/' || true)
bad=$(echo "$bad" | sed '/^$/d')
if [ -n "$bad" ]; then
  echo "NOT the accepted candidate: non-docs changes since $A:"; echo "$bad"; exit 1
fi
echo "OK: $B equals accepted $A up to docs/**"

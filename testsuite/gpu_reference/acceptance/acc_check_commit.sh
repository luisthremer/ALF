#!/bin/bash
# Acceptance harness, LAPTOP side: before committing, check that the commit equals the accepted candidate except docs/**.
#   usage: acc_check_commit.sh ACCEPTED_SHA [COMMIT=HEAD]
# Exit 0 only if every changed path between the two lies under docs/.
set -euo pipefail
A=$(git rev-parse --verify "$1^{commit}"); B=$(git rev-parse --verify "${2:-HEAD}^{commit}")
bad=$(git diff --name-only "$A" "$B" | grep -v '^docs/' || true)
if [ -n "$bad" ]; then
  echo "NOT the accepted candidate: non-docs changes since $A:"; echo "$bad"; exit 1
fi
echo "OK: $B equals accepted $A up to docs/**"

#!/bin/bash
# rsync the worktree to ~/alf-gpu/work/t302a/dev/src on a machine: t302a_sync.sh intel|nividia  (run on the laptop)
set -e
cd "$(dirname "$0")/../../../.."
ssh $1 'mkdir -p ~/alf-gpu/work/t302a/dev/src ~/alf-gpu/work/t302a/dev/tmp'
rsync -a --delete --exclude .git --exclude '*.swp' ./ $1:alf-gpu/work/t302a/dev/src/

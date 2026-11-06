#!/bin/bash
# sync Prog from the worktree and rebuild on nividia
cd /Users/luis/Documents/_WORK/Coding/CLAUDE/ALF/.claude/worktrees/gpu-performance-bottlenecks-bc7e27 || exit 1
rsync -a Prog/*.F90 Prog/Makefile nividia:alf-gpu/work/dlsub/src/Prog/ || exit 1
ssh nividia bash -s <<'R'
cd ~/alf-gpu/work/dlsub/src && export TMPDIR=~/alf-gpu/work/dlsub/tmp && . ~/alf-gpu/env/nvhpc.sh && . ./configure.sh PGI noMPI GPU=nvidia NO-INTERACTIVE >/dev/null 2>&1
cd Prog && make ALF.out > ../../make.log 2>&1; echo build rc=$?; grep -i "NVFORTRAN-S\|Error" ../../make.log | head
R

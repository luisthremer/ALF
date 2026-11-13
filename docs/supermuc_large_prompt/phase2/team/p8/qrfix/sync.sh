#!/bin/bash
# push p8/qrfix (scripts, test source, header copies) and, if given, the git-archive tree of HEAD to $SCRATCH/large_prompt/team/p8/qrfix
# usage: sync.sh [path-of-extracted-git-archive]
D=$(cd $(dirname $0) && pwd)
F='grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"'
timeout 60 ssh supermuc2 'mkdir -p $SCRATCH/large_prompt/team/p8/qrfix' >/dev/null 2>&1
timeout 120 rsync -a $D/ supermuc2:'$SCRATCH/large_prompt/team/p8/qrfix/' 2>&1 | eval $F
[ -n "$1" ] && timeout 300 rsync -a --delete $1/ supermuc2:'$SCRATCH/large_prompt/team/p8/qrfix/src_head/' 2>&1 | eval $F
echo synced

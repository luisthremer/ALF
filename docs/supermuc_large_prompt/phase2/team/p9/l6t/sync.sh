#!/bin/bash
# push p9/l6t scripts to $SCRATCH/large_prompt/team/p9/l6t (src_head was pushed once from `git archive HEAD` without docs/, see NOTES.md)
D=$(cd $(dirname $0) && pwd)
timeout 60 ssh supermuc2 'mkdir -p $SCRATCH/large_prompt/team/p9/l6t/jobs' >/dev/null 2>&1
timeout 120 rsync -a $D/jobs/ supermuc2:'$SCRATCH/large_prompt/team/p9/l6t/jobs/' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"
echo synced

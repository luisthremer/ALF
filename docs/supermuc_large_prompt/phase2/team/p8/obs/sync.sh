#!/bin/bash
# push p8/obs scripts + bench to the cluster
R=/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_GPU_SUBMATRIX_POINT_DEFECT/ALF_wt_lp_team_combined; D=$R/docs/supermuc_large_prompt/phase2/team/p8/obs
F='grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"'
timeout 60 ssh supermuc2 'mkdir -p $SCRATCH/large_prompt/team/p8/obs/{bench,jobs}' >/dev/null 2>&1
timeout 120 rsync -a $D/bench/ supermuc2:'$SCRATCH/large_prompt/team/p8/obs/bench/' 2>&1 | eval $F
timeout 120 rsync -a $D/jobs/ supermuc2:'$SCRATCH/large_prompt/team/p8/obs/jobs/' 2>&1 | eval $F
echo synced

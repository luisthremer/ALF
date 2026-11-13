#!/bin/bash
# push the p8/blk job scripts and (optionally) the source tree to $SCRATCH/large_prompt/team/p8/blk
R=/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_GPU_SUBMATRIX_POINT_DEFECT/ALF_wt_lp_team_combined
timeout 120 rsync -a $R/docs/supermuc_large_prompt/phase2/team/p8/blk/jobs/ supermuc2:'$SCRATCH/large_prompt/team/p8/blk/jobs/' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack|INFO"
echo synced

#!/bin/bash
# push the Phase 6 prototype sources and scripts to $SCRATCH/large_prompt/team/p6 (src_alf = copy of p5/src_alf + the three prototype files)
R=/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_GPU_SUBMATRIX_POINT_DEFECT/ALF_wt_lp_team_combined
F='grep -v -E module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack'
timeout 60 ssh supermuc2 '[ -d $SCRATCH/large_prompt/team/p6/src_alf ] || cp -a $SCRATCH/large_prompt/team/p5/src_alf $SCRATCH/large_prompt/team/p6/src_alf' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"
timeout 120 rsync -a $R/Prog/gpu_wrap_mod.F90 $R/Prog/gpu_update_mod.F90 $R/Prog/main.F90 supermuc2:'$SCRATCH/large_prompt/team/p6/src_alf/Prog/' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"
timeout 120 rsync -a $R/docs/gpu/sessions/supermuc_large_prompt/phase2/team/p6/*.sh $R/docs/gpu/sessions/supermuc_large_prompt/phase2/team/p6/*.py $R/docs/gpu/sessions/supermuc_large_prompt/phase2/team/p6/*.F90 supermuc2:'$SCRATCH/large_prompt/team/p6/' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"
echo synced

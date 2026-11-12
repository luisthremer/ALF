#!/bin/bash
# push the sources of Phase 5 to $SCRATCH/large_prompt/team/p5 (src/ = Libraries/Modules files, jobs/ = scripts)
R=/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_GPU_SUBMATRIX_POINT_DEFECT/ALF_wt_lp_team_combined
timeout 60 ssh supermuc2 'mkdir -p $SCRATCH/large_prompt/team/p5/{src,jobs,queue}' >/dev/null 2>&1
timeout 120 rsync -a $R/Libraries/Modules/alf_la_sycl_devqr.hpp $R/Libraries/Modules/alf_la_sycl_devqr.cpp $R/Libraries/Modules/alf_la_sycl_qr2.hpp $R/Libraries/Modules/alf_la_sycl_shim.cpp supermuc2:'$SCRATCH/large_prompt/team/p5/src/' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"
timeout 120 rsync -a $R/docs/gpu/sessions/supermuc_large_prompt/phase2/team/p5/*.sh $R/docs/gpu/sessions/supermuc_large_prompt/phase2/team/p5/*.cpp supermuc2:'$SCRATCH/large_prompt/team/p5/jobs/' 2>&1 | grep -v -E "module|Unload|WARNING|ERROR|EAR|SLURM|^$|spack"
echo synced

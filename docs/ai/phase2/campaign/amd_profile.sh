#!/bin/bash
# wait until the MI210 is idle (other sessions use it), then rocprofv3 of the Hubbard L = 24 / 32 runs
idle=0; while [ $idle -lt 3 ]; do u=$(rocm-smi --showuse 2>/dev/null | grep -m1 "GPU use" | grep -o "[0-9]*$"); if [ "${u:-100}" -lt 5 ]; then idle=$((idle+1)); else idle=0; fi; sleep 20; done
. ~/p2_amd/env_run.sh
cd ~/p2_amd/camp && python3 docs/ai/phase2/campaign/profile.py rocprof $PWD/Prog/ALF.out ~/p2_amd/campaign/prof_work ~/p2_amd/campaign/out mi210 > ~/p2_amd/campaign/prof.out 2>&1
echo PROF-DONE > ~/p2_amd/campaign/PROF_DONE

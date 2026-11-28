#!/bin/bash
# MI210 (GPU=amd) job set, one at a time: tol, then bench.  usage: amd_run.sh   (run on host amd, build in ~/p2_amd/camp)
. ~/p2_amd/env_run.sh
cd ~/p2_amd/camp; C=docs/ai/phase2/campaign/campaign.py; W=~/p2_amd/campaign; mkdir -p $W/out
python3 $C tol --exe $PWD/Prog/ALF.out --machine mi210 --work $W/tol_work --out $W/out --models hub-mz,hub-su2 --jobs 6 --nsweep 2 --timeout 7200
python3 $C tol --exe $PWD/Prog/ALF.out --machine mi210-refused --work $W/tol_work_other --out $W/out --models kondo,tV,sp --sizes 8 --jobs 3 --nsweep 2 --timeout 3600
python3 $C bench --exe $PWD/Prog/ALF.out --machine mi210 --work $W/bench_work --out $W/out --models hub-mz --timeout 1800
python3 $C bench --exe $PWD/Prog/ALF.out --machine mi210-refused --work $W/bench_work_other --out $W/out --models kondo,tV,sp --sizes 8 --timeout 1800
echo ALL-DONE > $W/DONE

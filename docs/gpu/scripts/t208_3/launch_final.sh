#!/bin/bash
# 208-3: one launch for both machines.  Run from the repository root (branch gpu-port-208-3, committed):  bash docs/gpu/scripts/t208_3/launch_final.sh
# Syncs HEAD to ~/alf-gpu/work/t208-3f on nividia and intel and starts final.sh detached on both.  Result: logs/SUMMARY.txt there.
set -e
for m in nividia intel; do
  P=$([ $m = nividia ] && echo nv || echo in)
  ssh $m 'W=$HOME/alf-gpu/work/t208-3f; rm -rf $W/src; mkdir -p $W/src $W/scripts $W/logs $W/tmp $W/bin $W/runs $W/ref'
  git archive HEAD | ssh $m 'tar -x -C $HOME/alf-gpu/work/t208-3f/src'
  for f in build builds_nvidia builds_intel t0 obs ctest; do
    sed 's/t208-3/t208-3f/g' docs/gpu/scripts/t208_3/$f.sh | ssh $m "cat > \$HOME/alf-gpu/work/t208-3f/scripts/$f.sh"
  done
  cat docs/gpu/scripts/t208_3/final.sh | ssh $m "cat > \$HOME/alf-gpu/work/t208-3f/scripts/final.sh"     # already uses t208-3f
  scp -q docs/gpu/scripts/t208_3/t207_firstdiv.py ${m}:alf-gpu/work/t208-3f/scripts/
  if [ $m = nividia ]; then
    ssh $m 'cd $HOME/alf-gpu/work/t208-3f; ln -sfn $HOME/alf-gpu/runs/ref_gnu-DEFAULT ref/ref_gnu-DEFAULT; ln -sfn $HOME/alf-gpu/runs/ref_gnu-LOG ref/ref_gnu-LOG; ln -sfn $HOME/alf-gpu/work/t207-final/src_up src_up'
  else
    ssh $m 'cd $HOME/alf-gpu/work/t208-3f; ln -sfn $HOME/alf-gpu/work/t207-final/ref/ref_gnu-DEFAULT ref/ref_gnu-DEFAULT; ln -sfn $HOME/alf-gpu/work/t207-final/ref/ref_gnu-LOG ref/ref_gnu-LOG'
  fi
  ssh $m "cd \$HOME/alf-gpu/work/t208-3f && setsid nohup bash scripts/final.sh $P > logs/final.out 2>&1 < /dev/null & echo started on $m"
done

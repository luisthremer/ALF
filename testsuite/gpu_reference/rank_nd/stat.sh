#!/bin/bash
# Ticket 401-B (development seeds only: the committed pyALF default seeds rotated by --seed-offset): statistical runs.
#   usage: stat.sh GPUBIN CPUBIN      (e.g. stat.sh nv cpu)    4 independent runs per Symm, NBin 60 x NSweep 20, Ltau 1
# GPU runs: seed offsets 0..3 (gpu lock on nividia); CPU runs: seed offsets 10..13.  Dirs runs/st_<bin>/pvs<Symm>_k<i>.
W=$HOME/alf-gpu/work/t401b; R=$W/src/testsuite/gpu_reference/rank_nd; G=$1; C=$2
if [ -d /opt/nvidia ]; then . ~/alf-gpu/env/gnu.sh; LK=1; else . ~/alf-gpu/env/intel.sh; LK=0; fi
. ~/alf-gpu/venv/bin/activate 2>/dev/null; export TMPDIR=$W/tmp; cd $W || exit 1
for s in 1 0; do for k in 0 1 2 3; do
  S=$([ $s = 1 ] && echo True || echo False)
  for b in $G $C; do
    off=$k; [ $b = $C ] && off=$((k+10))
    n=pvs${s}_k${k}_$b
    [ -d runs/base/$n ] || python3 $R/mkrun.py src runs/base $n Symm=$S Ltau=1 NSweep=20 NBin=60 --seed-offset $off > /dev/null 2>&1
    if [ $b = $G ]; then GPULOCK=$LK TAG=st_$b $R/run.sh $b $n; else GPULOCK=0 TAG=st_$b $R/run.sh $b $n; fi
  done
done; done

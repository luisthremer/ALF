#!/bin/bash
# Ticket 207, nividia: run directories for the statistical comparison (compare2.py), config from configs_stat.py, with
# DIFFERENT explicit seeds for the GPU run (runs/stat/<cfg>_gpu) and the GPU=none run (runs/stat/<cfg>_cpu).
#   usage: t207_stat_prep.sh CONFIG NSWEEP NBIN
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
C=$1; NS=$2; NB=$3
cd $W/src_new/testsuite/gpu_reference || exit 1
: "${HOME:?}"
rm -rf -- "$W/runs/stat"; mkdir -p "$W/runs/stat"
python3 prep.py $W/runs/stat $C --nsweep $NS --nbin $NB --alf-dir $W/src_new > /dev/null
cp -r "$W/runs/stat/$C" "$W/runs/stat/${C}_gpu"; mv "$W/runs/stat/$C" "$W/runs/stat/${C}_cpu"
python3 -I - "$W/runs/stat/${C}_gpu/seeds" "$W/runs/stat/${C}_cpu/seeds" <<'PY'
import sys
# explicit, different seed lists (one integer per line, as ALF's seeds file)
a = [1357911, 2468013, 3141592, 2718281, 1618033, 1414213, 1732050, 2236067, 2645751, 3162277]
b = [9876543, 8765432, 7654321, 6543210, 5432109, 4321098, 3210987, 2109876, 1098765, 9988776]
for f, v in zip(sys.argv[1:], (a, b)):
    open(f, "w").write("".join(f"{x:12d}\n" for x in v))
PY
echo prepared "$C" NSweep=$NS NBin=$NB

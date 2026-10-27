#!/bin/bash
# Ticket 202, run on nividia after `t202_build.sh cuda` (or any tool): run the alf_la timing table.
#   usage: t202_timing.sh cuda|pgi|gnu ...    (uses ~/alf-gpu/work/t202/tree_<tool>/testsuite/tests/alf_la.tests/t_timing)
T=$1; W=$HOME/alf-gpu/work/t202; export TMPDIR=$W/tmp
case $T in gnu) . ~/alf-gpu/env/gnu.sh;; intel) . ~/alf-gpu/env/intel.sh;; *) . ~/alf-gpu/env/nvhpc.sh;; esac
[ $T = cuda ] && { nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv; export OMP_TARGET_OFFLOAD=MANDATORY; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 840 $W/tree_$T/testsuite/tests/alf_la.tests/t_timing

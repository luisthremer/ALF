#!/bin/bash
# mpiexec wrapper of the "gpu" run of run_cmp_gg.py: ALFSH_BLKPACK=$PACK_N (default 1) (+ ALFSH_BLKPACK_PLANT=$PLANT_N), everything else as the "cpu" (reference) run
export ALFSH_BLKPACK=${PACK_N:-1}
[ -n "$PLANT_N" ] && export ALFSH_BLKPACK_PLANT=$PLANT_N
exec "$@"

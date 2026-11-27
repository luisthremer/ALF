#!/bin/bash
# Per-rank pinning wrapper: mpiexec -n NR ./pin.sh EXE   (Intel MPI pinning is switched off, I_MPI_PIN=0)
# Node: 2 NUMA domains x 56 cores (0-55 | 56-111), hyperthreads +112. GPUs 0,1 (tiles 0-3) on NUMA 0, GPUs 2,3 (tiles 4-7) on NUMA 1.
# MODE=cpu: rank r gets T consecutive physical cores inside one NUMA domain (R ranks split evenly over the 2 sockets).
# MODE=gpu: rank r -> tile r / W (ZE_AFFINITY_MASK, flat hierarchy); bound to that tile's 14 nearby cores and their
#           hyperthreads (the GPU runtime's helper threads need room; pinning a rank to one cpu made it ~6x slower, jobs 560229-34).
r=${MPI_LOCALRANKID:?MPI_LOCALRANKID not set}
if [ "$MODE" = cpu ]; then
  rps=$(( NR / 2 )); s=$(( r / rps )); i=$(( r % rps ))
  c0=$(( s * 56 + i * T )); c1=$(( c0 + T - 1 ))
  exec numactl --physcpubind=$c0-$c1 --membind=$s "$@"
else
  tile=$(( r / W )); s=$(( tile / 4 ))
  c0=$(( s * 56 + (tile % 4) * 14 )); c1=$(( c0 + 13 ))
  export ZE_AFFINITY_MASK=$tile
  exec numactl --physcpubind=$c0-$c1,$(( c0 + 112 ))-$(( c1 + 112 )) --membind=$s "$@"
fi

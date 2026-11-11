#!/bin/bash
# one rank per GPU card (tile 0 of the card): flat hierarchy, devices are numbered card*2+tile
export ZE_AFFINITY_MASK=$(( ${MPI_LOCALRANKID:-$PMI_RANK} * 2 ))
exec "$@"

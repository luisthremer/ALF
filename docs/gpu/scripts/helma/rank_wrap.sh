#!/bin/bash
# per-rank wrapper: bind rank to one GPU (ALF picks the device as rank mod #devices, so each rank must see exactly one GPU)
# env (set by run_pp.sh): RPG ranks per GPU, MPS_PCT (0 = no MPS), MPS_DIR, ALF_BIN
r=${OMPI_COMM_WORLD_RANK:?}
g=$(( r / RPG ))
export CUDA_VISIBLE_DEVICES=$g
if [ "$MPS_PCT" != 0 ]; then
  export CUDA_MPS_PIPE_DIRECTORY=$MPS_DIR/pipe$g CUDA_MPS_LOG_DIRECTORY=$MPS_DIR/log$g CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=$MPS_PCT
fi
exec "$ALF_BIN"

#!/bin/bash
# per-rank wrapper for the prototype timing runs: WBCFG=mixed (one configuration per rank, below) or uni (WB_B, WB_M for every rank; WB_B empty = default)
r=${PMI_RANK:-0}
export ALF_GPU_WB_TIME=1 ALF_GPU_DEBUG=1   # ALF refuses ALF_GPU_* variables in production runs unless ALF_GPU_DEBUG=1 (allow-list, no numerical effect)
case "$WBCFG" in
  mixed) case $r in
    0) ;;                                       # default, unbatched (baseline)
    1) export ALF_GPU_WB=1  ALF_GPU_WB_MODE=1;;
    2) export ALF_GPU_WB=2  ALF_GPU_WB_MODE=1;;
    3) export ALF_GPU_WB=4  ALF_GPU_WB_MODE=1;;
    4) export ALF_GPU_WB=8  ALF_GPU_WB_MODE=1;;
    5) export ALF_GPU_WB=4  ALF_GPU_WB_MODE=2;;
    6) export ALF_GPU_WB=4  ALF_GPU_WB_MODE=3;; # separate gemms (control)
    7) export ALF_GPU_WB=16 ALF_GPU_WB_MODE=1;;
  esac;;
  uni) [ -n "$WB_B" ] && export ALF_GPU_WB=$WB_B ALF_GPU_WB_MODE=${WB_M:-1};;
esac
exec "$@"

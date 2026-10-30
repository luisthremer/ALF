#!/bin/bash
# run a built binary on the B580 with the mandatory environment: t302a_run.sh BINARY [args]
ARGS=("$@")
. $HOME/alf-gpu/env/intel.sh 2>/dev/null
export ONEAPI_DEVICE_SELECTOR=level_zero:0 OMP_TARGET_OFFLOAD=MANDATORY FORT_BUFFERED=false
exec "${ARGS[@]}"

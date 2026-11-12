#!/bin/bash
# preload the PTI tracer into the application process only (not into mpiexec/hydra)
PTI=/lrz/sys/intel/oneapi_2025.2.1/pti/latest
export LD_PRELOAD=$SCRATCH/large_prompt/team/jobs/libptitrace.so:$PTI/lib/libpti_view.so
exec "$@"

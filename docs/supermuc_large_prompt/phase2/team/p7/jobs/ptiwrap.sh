#!/bin/bash
# preload the PTI kernel tracer (Phase 4/team ptitrace.cpp, built in team/jobs) into the application only
PTI=${PTIDIR:-/lrz/sys/intel/oneapi_2025.1.0/pti/latest}
export LD_PRELOAD=$SCRATCH/large_prompt/team/jobs/libptitrace.so:$PTI/lib/libpti_view.so
exec "$@"

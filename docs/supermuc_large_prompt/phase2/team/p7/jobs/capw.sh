#!/bin/bash
# per-rank wrapper: own stdout/stderr file per rank (PLUGIN_PROFILE tables, ALFP7/ALFSH/ALFDEVQR_PROF lines); optional LD_PRELOAD (P7PRE) for the application only
[ -n "$P7PRE" ] && export LD_PRELOAD=$P7PRE
exec "$@" > rank${PMI_RANK:-${PMIX_RANK:-0}}.log 2>&1

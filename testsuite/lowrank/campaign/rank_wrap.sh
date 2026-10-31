#!/bin/bash
# rank_wrap.sh MODE OUTDIR CMD ...   per-rank wrapper, started by the launcher once per MPI rank (bash 3.2, no python).
#   MODE alf    record rank provenance, then run CMD under /usr/bin/time (max RSS per rank: OUTDIR/rss_<rank>.txt)
#   MODE probe  record rank provenance, then run CMD as is (the MKL_VERBOSE probe)
# fakeintel: LD_PRELOAD is set HERE only (never in the job environment), so it reaches exactly the rank processes.
# Writes OUTDIR/rank_<MODE>_<rank>.prov: rank, host, pid, affinity (taskset -cp), ld_preload, fakeintel_maps.
# Exit 90: fakeintel on but the library is missing; 91: library not mapped (/proc/self/maps check).
mode="$1"
out="$2"
shift 2
. "$(cd "$(dirname "$0")" && pwd)/provenance.sh"

# per-rank variables of the MPI library first: under mpiexec (Intel MPI / MPICH hydra) SLURM_PROCID is the proxy's
# (or the batch step's) value, identical for all ranks of a node; with direct srun launch it is the rank.
rank="${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-0}}}}"
f="$out/rank_${mode}_${rank}.prov"
prov_kv "$f" rank "$rank"
prov_kv "$f" host "$(hostname)"
prov_kv "$f" pid "$$"
if command -v taskset >/dev/null 2>&1; then
    prov_kv "$f" affinity "$(taskset -cp $$ 2>&1 | sed 's/.*: *//')"
else
    prov_kv "$f" affinity "unavailable(no taskset)"
fi

if [ "${FAKEINTEL:-off}" = on ]; then
    if [ -z "$FAKEINTEL_LIB" ] || [ ! -r "$FAKEINTEL_LIB" ]; then
        prov_kv "$f" fakeintel_maps "library missing: '$FAKEINTEL_LIB'"
        exit 90
    fi
    export LD_PRELOAD="$FAKEINTEL_LIB${LD_PRELOAD:+:$LD_PRELOAD}"
    maps="${PROV_PROC_ROOT}/proc/self/maps"
    if [ -r "$maps" ]; then
        # grep inherits LD_PRELOAD, so its own /proc/self/maps shows whether the dynamic loader honours the preload
        n="$(grep -F -c "$(basename "$FAKEINTEL_LIB")" "$maps")"
        if [ "${n:-0}" -gt 0 ]; then
            prov_kv "$f" fakeintel_maps "loaded ($n mappings in /proc/self/maps)"
        else
            prov_kv "$f" fakeintel_maps "NOT loaded (0 mappings in /proc/self/maps)"
            exit 91
        fi
    else
        prov_kv "$f" fakeintel_maps "unavailable(no /proc/self/maps)"
    fi
else
    prov_kv "$f" fakeintel_maps "n/a (fakeintel off)"
fi
prov_kv "$f" ld_preload "${LD_PRELOAD:-unset}"

if [ "$mode" = alf ]; then
    case "${TIME_STYLE:-none}" in
        gnu) [ -x /usr/bin/time ] && exec /usr/bin/time -v -o "$out/rss_$rank.txt" "$@" ;;
        bsd) /usr/bin/time -l "$@" 2> "$out/rss_$rank.txt"; exit $? ;;
    esac
fi
exec "$@"

#!/bin/bash
# Re-judge the pinning of finished runs that were marked INVALID only because the I_MPI_DEBUG pin table could not be
# parsed (Intel MPI 2019 cuts "Pin cpu" lists after ~30 entries; fixed in prov_pin_check, which now falls back to the
# OMP_AFFINITY lines). Only runs whose status starts with "INVALID pinning: BAD unparsable I_MPI_DEBUG pin line" are
# touched; every other check of the run (launcher rc, info file, scheme, placement) passed when that status was set,
# because the job template checks those first. For each such run the latest attempt's pin_map.txt / omp_affinity.txt
# are judged again with the current provenance.sh; on OK the run becomes OK, and run.prov keeps the original verdict
# (pin_check_original, status_original) next to the new one (pin_check, status, pin_rejudged). Otherwise nothing changes.
# Runs on the cluster (GNU sed -i). Usage: testsuite/lowrank/campaign/rejudge_pinning.sh CAMPAIGN_DIR [--dry-run]
HERE="$(cd "$(dirname "$0")" && pwd)"
C="${1:?campaign directory}"; DRY="${2:-}"
. "$HERE/provenance.sh"
psha="$(prov_sha256 "$HERE/provenance.sh")"
nok=0; nbad=0
for R in "$C"/runs/*/; do
    R="${R%/}"; st="$(cat "$R/status" 2>/dev/null)"
    case "$st" in "INVALID pinning: BAD unparsable I_MPI_DEBUG pin line"*) ;; *) continue ;; esac
    A="$(ls -d "$R"/a[0-9]* 2>/dev/null | sort -V | tail -1)"
    lay="$(basename "$R" | sed -n 's/.*_\([0-9][0-9]*x[0-9][0-9]*\)_r[0-9].*/\1/p')"
    [ -n "$A" ] && [ -n "$lay" ] || { echo "SKIP $(basename "$R"): no attempt or layout"; continue; }
    ranks="${lay%x*}"; threads="${lay#*x}"
    pin="$(prov_pin_check "$A/pin_map.txt" "$ranks" "$threads" 1 "$A/omp_affinity.txt")"; rc=$?
    if [ "$rc" -ne 0 ]; then echo "STILL_BAD $(basename "$R"): $pin"; nbad=$((nbad + 1)); continue; fi
    echo "OK $(basename "$R"): $pin"; nok=$((nok + 1))
    [ "$DRY" = --dry-run ] && continue
    f="$A/run.prov"
    sed -i -e 's/^pin_check=/pin_check_original=/' -e 's/^status=/status_original=/' "$f"
    prov_kv "$f" pin_check "$pin"
    prov_kv "$f" pin_rejudged "$(date '+%Y-%m-%dT%H:%M:%S') rejudge_pinning.sh, provenance.sh sha256 $psha"
    prov_kv "$f" status OK
    echo OK > "$R/status"
    echo "$(date '+%Y-%m-%dT%H:%M:%S') rejudge_pinning.sh OK ($pin)" >> "$R/status.history"
done
echo "rejudge_pinning: $nok OK, $nbad still bad${DRY:+ ($DRY)}"

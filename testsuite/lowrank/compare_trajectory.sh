#!/usr/bin/env bash
# Trajectory comparison of two Sequential_scheme values.
#
# Usage:  testsuite/lowrank/compare_trajectory.sh [--rebuild] SET SCHEME_A SCHEME_B [KEY=VALUE ...]
#   SET                 one of the sets of testsuite/lowrank/sets.sh: SETS (start, explicit, doped, density,
#                       continuous, tv, kondo, pv_adiab, z2gt, tv_flux, tv_proj, z2, spin_peierls) or TRAJ_ONLY_SETS
#                       (jz); only its common keys are used (the scheme is set explicitly).
#   SCHEME_A, SCHEME_B  standard | delayed | submatrix (any case)
#   KEY=VALUE ...       extra parameters for both runs (e.g. N_delay=8, SEED=42)
#   --rebuild           force a rebuild of the trace binary.
#   ALF_TRACE_BIN=path  (environment, for testing the harness) use that binary instead of building one.
#
# Two runs with identical seeds and parameters (NSweep=3 NBin=2, OMP_NUM_THREADS=1) differ only in
# Sequential_scheme and must make the same accept/reject decisions. The runs use one shared trace binary
# builds/ALF_trace.out (ALF built with -DALF_TRACE_UPDATES; it appends every Metropolis decision to
# trace_updates.txt) and live in runs/traj/<SET>/<scheme> (<scheme>_2 for the second of two identical schemes).
# trace_diff.py walks both traces: identical decisions -> compare_runs.py requires identical confout_0 and bins
# within 1e-8; the first differing toggle is classified ROUNDOFF (both runs within 1e-10 of the Metropolis
# threshold: prints "DIVERGED-ROUNDOFF at proposal k (margin ...)", skips the confout/bin comparison, exit 0)
# or BUG (exit 1). Non-standard schemes must also report "Sequential scheme: <NAME>" in info.
# The normal binary is restored on any exit (trap) if the trace build has touched the tree.
# Exit: 0 pass, 1 mismatch or BUG, 2 harness failure (build, run, usage).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
MAKE_PARAMS="$HERE/make_params.py"
TRACE_BIN="$ALFDIR/builds/ALF_trace.out"
TRACE_BIN_OVERRIDE="${ALF_TRACE_BIN:-}"   # testing hook: use this binary, no build/staleness check
OUTROOT="$ALFDIR/runs/traj"
SHORT="NSweep=3 NBin=2"

REBUILD=0
if [ "${1:-}" = "--rebuild" ]; then REBUILD=1; shift; fi
if [ $# -lt 3 ]; then sed -n '2,22p' "${BASH_SOURCE[0]}" >&2; exit 2; fi
SET="$1"; SCHEME_A="$2"; SCHEME_B="$3"; shift 3
EXTRA=("$@")

source "$HERE/sets.sh"
common=""; found=0
for entry in "${SETS[@]}" "${TRAJ_ONLY_SETS[@]}"; do
    IFS='|' read -r name c _ <<<"$entry"
    if [ "$name" = "$SET" ]; then common="$c"; found=1; fi
done
[ $found = 1 ] || { echo "unknown set '$SET'" >&2; exit 2; }

mkdir -p "$ALFDIR/builds" "$ALFDIR/runs"

build_in_tree() {  # $1 = extra flag or empty; leaves Prog/ALF.out built from scratch with that flag
    ( set +eu; cd "$ALFDIR"
      if [ -n "$1" ]; then export ALF_FLAGS_EXT="$1"; fi
      source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1
      make -C Prog clean >/dev/null 2>&1 && make program -j8 ) >"$ALFDIR/runs/build_trace.log" 2>&1 \
        || { echo "build failed, see $ALFDIR/runs/build_trace.log" >&2; return 2; }
}

# While the tree holds traced objects (between the two builds) TRACE_TREE_DIRTY=1. Any exit (failure,
# interrupt) then rebuilds the normal binary, so later incremental builds never reuse traced objects.
TRACE_TREE_DIRTY=0
restore_normal_build() {
    local rc=$?
    trap - EXIT INT TERM
    if [ "$TRACE_TREE_DIRTY" = 1 ]; then
        echo "restoring the normal build ..." >&2
        build_in_tree "" || rc=2
        TRACE_TREE_DIRTY=0
    fi
    exit $rc
}
trap restore_normal_build EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

trace_binary_stale() {
    [ -x "$TRACE_BIN" ] || return 0
    # Generated files (touched by every clean+rebuild or commit) are not sources.
    stale=$(find "$ALFDIR/Prog" "$ALFDIR/Libraries" -type f \( -name '*.F90' -o -name '*.h' -o -name Makefile \) \
        -newer "$TRACE_BIN" ! -name 'git.h' ! -name 'git_status.h' ! -name 'Hamiltonians_case.h' \
        ! -name 'Hamiltonians_interface.h' ! -name '*_read_write_parameters.F90' | head -1)
    [ -n "$stale" ]
}

if [ -n "$TRACE_BIN_OVERRIDE" ]; then
    TRACE_BIN="$TRACE_BIN_OVERRIDE"; echo "using $TRACE_BIN (ALF_TRACE_BIN)"
elif [ $REBUILD = 1 ] || trace_binary_stale; then
    echo "building trace binary (-DALF_TRACE_UPDATES) ..."
    TRACE_TREE_DIRTY=1
    build_in_tree "-DALF_TRACE_UPDATES" || exit 2
    cp "$ALFDIR/Prog/ALF.out" "$TRACE_BIN"
    if [ -n "${TRAJ_TEST_FAIL_AFTER_TRACE:-}" ]; then echo "simulated failure after the trace build" >&2; exit 1; fi
    echo "rebuilding the normal binary ..."
    build_in_tree "" || exit 2
    TRACE_TREE_DIRTY=0
else
    echo "reusing $TRACE_BIN"
fi

export OMP_NUM_THREADS=1
dirs=()
for scheme in "$SCHEME_A" "$SCHEME_B"; do
    dir="$OUTROOT/$SET/$scheme"
    [ -d "$dir" ] && [ "${#dirs[@]}" = 1 ] && [ "$dir" = "${dirs[0]}" ] && dir="${dir}_2"
    /bin/rm -rf "$dir"
    keys="$SHORT $common ${EXTRA[*]:-}"
    read -r -a keyarr <<<"$keys"
    python3 -I "$MAKE_PARAMS" "$dir" "${keyarr[@]}" "Sequential_scheme=\"$scheme\""
    (cd "$dir" && "$TRACE_BIN" >stdout.log 2>&1) \
        || { echo "run failed ($scheme), see $dir/stdout.log" >&2; exit 2; }
    dirs+=("$dir")
done
DIR_A="${dirs[0]}"; DIR_B="${dirs[1]}"

status=0
td_out=$(python3 -I "$HERE/trace_diff.py" "$DIR_A/trace_updates.txt" "$DIR_B/trace_updates.txt") || td_rc=$?
td_rc=${td_rc:-0}
echo "$td_out"
result=$(printf '%s\n' "$td_out" | sed -n 's/^RESULT: //p' | tail -1)
case "$result" in
    ROUNDOFF*)
        read -r _ k margin <<<"$result"
        echo "DIVERGED-ROUNDOFF at proposal $k (margin $margin)"
        python3 -I "$HERE/compare_runs.py" --info-only "$DIR_A" "$DIR_B" "$SCHEME_A" "$SCHEME_B" || status=1 ;;
    IDENTICAL*)
        python3 -I "$HERE/compare_runs.py" "$DIR_A" "$DIR_B" "$SCHEME_A" "$SCHEME_B" || status=1 ;;
    BUG*)
        echo "TRAJECTORY BUG: $result"; status=1 ;;
    *)
        echo "trace_diff failed: $result" >&2; exit 2 ;;
esac
if [ $status = 0 ]; then echo "$SET $SCHEME_A vs $SCHEME_B: PASS"; else echo "$SET $SCHEME_A vs $SCHEME_B: FAIL"; fi
exit $status

#!/usr/bin/env bash
# One benchmark job: Mz Hubbard, square L x L, U=4, Beta=2, Dtau=0.1, Ltau=0, one measurement slice, NBin=2.
#
# Usage:  testsuite/lowrank/bench.sh SCHEME L N_delay [THREADS]
#   SCHEME   standard | delayed | submatrix     (Sequential_scheme)
#   L        linear lattice size (N = L*L sites)
#   N_delay  max pending accepted vertices (ignored by 'standard' but recorded)
#   THREADS  OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS (default 1)
# Environment:
#   BENCH_BIN     ALF binary (default runs/bench/ALF.out, a frozen copy so rebuilds do not change it mid-scan)
#   NSWEEP        overrides the per-L sweep count (calibration / re-runs; a warning is printed to stderr)
#   TAG           optional suffix of the run directory name (e.g. TAG=_rep2)
#   BENCH_COMMIT  optional expected commit; a warning is printed if the binary reports another one
#   SEED          first seed (default: the one in Start/seeds); identical for all schemes => identical trajectories
# Run directory: runs/bench/<SCHEME>_L<L>_nd<N_delay>_t<THREADS>_ns<NSweep><TAG>
# Output (one line, stdout):  scheme L N N_delay threads sec_per_sweep acceptance nsweep commit date pgreen_mean pgreen_max pphase_max
#   (the last three are the "Precision Green Mean, Max" and "Precision Phase, Max" lines of info; stability diagnostics)
#   sec_per_sweep = "CPU Time" in info / (NSweep*NBin). "CPU Time" is WALL time (system_clock) of the sweep
#   phase only: the timer is restarted after the initial Green function computation (Prog/main.F90 calls
#   Control_init again before the sweeps, see Prog/control_mod.F90), so setup is excluded.
#   acceptance = the "Acceptance" line of info.
#   commit = "This executable represents commit ..." line of info, date = day of "START TIME" in info.
#   A warning goes to stderr if the commit differs from $BENCH_COMMIT or from the previous bench.sh run
#   (recorded in runs/bench/.last_commit); bench_report.py refuses to mix commits in one table.
#
# NSweep rule (same for all three schemes at a given L): the table below targets ~50 s for the standard scheme
# (1 thread), from the calibration runs (standard: 6.9 ms/sweep at L=6 ... 2.46 s/sweep at L=20, ~ N^2.4).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
[ $# -ge 3 ] || { sed -n 2,25p "${BASH_SOURCE[0]}" >&2; exit 2; }
SCHEME="$1"; L="$2"; ND="$3"; THREADS="${4:-1}"
case "$SCHEME" in
    standard|delayed|submatrix) ;;
    *) echo "unknown scheme '$SCHEME' (standard|delayed|submatrix)" >&2; exit 2 ;;
esac

# NSweep per L (NBin=2 => 2*NSweep sweeps in total). Calibrated 2026-10-31 on M1 Pro, 1 thread, standard scheme.
nsweep_for() {
    case "$1" in
        6)  echo 3600 ;;   # standard ~7 ms/sweep
        8)  echo 1000 ;;
        10) echo 350 ;;    # 71 ms/sweep
        12) echo 140 ;;
        14) echo 66 ;;     # 0.38 s/sweep
        16) echo 35 ;;
        18) echo 14 ;;     # added in Task 15 (standard ~1.8 s/sweep expected)
        20) echo 12 ;;     # 2.46 s/sweep
        *) echo "no NSweep calibrated for L=$1 (set NSWEEP=...)" >&2; return 1 ;;
    esac
}
NSW="${NSWEEP:-$(nsweep_for "$L")}"
NBIN=2
BIN="${BENCH_BIN:-$ALFDIR/runs/bench/ALF.out}"
[ -x "$BIN" ] || { echo "binary $BIN not found" >&2; exit 2; }

[ -z "${NSWEEP:-}" ] || echo "bench.sh: NSWEEP=$NSWEEP overrides the calibrated per-L value" >&2
DIR="$ALFDIR/runs/bench/${SCHEME}_L${L}_nd${ND}_t${THREADS}_ns${NSW}${TAG:-}"
/bin/rm -rf "$DIR"
mkdir -p "$DIR"
KEYS=(L1="$L" L2="$L" ham_U=4.d0 Beta=2.d0 Dtau=0.1d0 Ltau=0 LOBS_ST=1 LOBS_EN=1
      NSweep="$NSW" NBin="$NBIN" Sequential_scheme="\"$SCHEME\"" N_delay="$ND")
[ -z "${SEED:-}" ] || KEYS+=(SEED="$SEED")
python3 -I "$HERE/make_params.py" "$DIR" "${KEYS[@]}"

export OMP_NUM_THREADS="$THREADS" VECLIB_MAXIMUM_THREADS="$THREADS"
cd "$DIR"
"$BIN" >stdout.log 2>&1 || { echo "bench $SCHEME L=$L: run failed, see $DIR/stdout.log" >&2; exit 1; }
if [ "$SCHEME" != standard ]; then
    want="Sequential scheme: $(printf '%s' "$SCHEME" | tr '[:lower:]' '[:upper:]')"
    grep -q -F "$want" info || { echo "bench $SCHEME L=$L: expected line '$want' not in $DIR/info (scheme not active?)" >&2; exit 1; }
fi
cpu="$(awk '/^ *CPU Time *:/ {print $4}' info)"
acc="$(awk '/^ *Acceptance *:/ {print $3}' info)"
[ -n "$cpu" ] && [ -n "$acc" ] || { echo "bench $SCHEME L=$L: could not parse info in $DIR" >&2; exit 1; }
pg="$(awk '/^ *Precision Green *Mean, Max *:/ {print $6, $7}' info)"
pp="$(awk '/^ *Precision Phase, Max *:/ {print $5}' info)"
[ -n "$pg" ] && [ -n "$pp" ] || { echo "bench $SCHEME L=$L: no Precision lines in info of $DIR" >&2; exit 1; }
commit="$(sed -n 's/^ *This executable represents commit \([^ ]*\) .*/\1/p' info)"
day="$(sed -n 's/^ *START TIME: \([0-9-]*\)T.*/\1/p' info)"
[ -n "$commit" ] && [ -n "$day" ] || { echo "bench $SCHEME L=$L: no commit/START TIME line in $DIR/info" >&2; exit 1; }
if [ -n "${BENCH_COMMIT:-}" ] && [ "$commit" != "$BENCH_COMMIT" ]; then
    echo "bench.sh WARNING: binary reports commit $commit, expected $BENCH_COMMIT" >&2
fi
LAST="$ALFDIR/runs/bench/.last_commit"
if [ -f "$LAST" ] && [ "$(cat "$LAST")" != "$commit" ]; then
    echo "bench.sh WARNING: commit $commit differs from the previous run's commit $(cat "$LAST")" >&2
fi
echo "$commit" > "$LAST"
awk -v s="$SCHEME" -v l="$L" -v nd="$ND" -v t="$THREADS" -v cpu="$cpu" -v acc="$acc" -v ns="$NSW" -v nb="$NBIN" \
    -v c="$commit" -v d="$day" -v pg="$pg" -v pp="$pp" \
    'BEGIN { printf "%s %d %d %d %d %.6g %.6f %d %s %s %s %s\n", s, l, l*l, nd, t, cpu/(ns*nb), acc, ns, c, d, pg, pp }'

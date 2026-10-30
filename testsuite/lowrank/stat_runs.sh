#!/usr/bin/env bash
# One statistical job: parameters from STAT_SETS (sets.sh), run ALF.out, analyse with ana.out.
#
# Usage:  testsuite/lowrank/stat_runs.sh SET SCHEME SEED [KEY=VALUE ...]
#   SET     half | doped | density            (sets.sh, STAT_SETS)
#   SCHEME  standard | delayed | submatrix    (Sequential_scheme)
#   SEED    integer, first seed of the seeds file (independent chains need different seeds)
#   KEY=VALUE ...  extra parameters, applied last (override NSweep=50 NBin=40 Ltau=0 n_skip=5 too)
#   STAT_TAG=<s>   env: suffix of the run directory (e.g. for a modified-parameter run)
# Run directory: runs/stat/<SET>/<SCHEME>_s<SEED><STAT_TAG>. Start it with run_in_background for long jobs.
# Compare two directories with stat_compare.py.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
[ $# -ge 3 ] || { sed -n 2,11p "${BASH_SOURCE[0]}" >&2; exit 2; }
SET="$1"; SCHEME="$2"; SEED="$3"; shift 3
EXTRA=("$@")

source "$HERE/sets.sh"
keys=""
for entry in "${STAT_SETS[@]}"; do
    IFS='|' read -r name k <<<"$entry"
    [ "$name" = "$SET" ] && keys="$k"
done
[ -n "$keys" ] || { echo "unknown set '$SET'" >&2; exit 2; }

DIR="$ALFDIR/runs/stat/$SET/${SCHEME}_s${SEED}${STAT_TAG:-}"
/bin/rm -rf "$DIR"
mkdir -p "$DIR"
read -r -a keyarr <<<"NSweep=50 NBin=40 Ltau=0 n_skip=5 N_rebin=1 $keys Sequential_scheme=\"$SCHEME\" SEED=$SEED ${EXTRA[*]:-}"
python3 -I "$HERE/make_params.py" "$DIR" "${keyarr[@]}"

export OMP_NUM_THREADS=1
cd "$DIR"
"$ALFDIR/Prog/ALF.out" >stdout.log 2>&1 || { echo "$SET/$SCHEME: run failed, see $DIR/stdout.log" >&2; exit 1; }
for obs in Ener_scal Pot_scal Kin_scal Part_scal SpinZ_eq Green_eq; do
    "$ALFDIR/Analysis/ana.out" "$obs" >"ana_$obs.log" 2>&1 || { echo "ana.out $obs failed in $DIR" >&2; exit 1; }
done
echo "done $DIR"

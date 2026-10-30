#!/usr/bin/env bash
# Bit-identity regression: branch ALF.out vs the unmodified baseline commit.
#
# Usage:  testsuite/lowrank/regress_vs_baseline.sh [KEY=VALUE ...]
#   KEY=VALUE ...      extra parameters applied to every set, baseline and branch.
#   BRANCH_EXTRA="K=V ..."   extra parameters applied to the branch side only
#                      (keys the baseline cannot parse, or a deliberate difference).
# For each parameter set it runs the baseline and the branch binary from
# identical parameters/seeds in runs/regress/<set>/{base,branch} and cmp's every
# output file except info* (commit and time stamps) and the inputs parameters, seeds.
# Prints "<set> IDENTICAL|DIFFERENT"; exit 0 iff all sets are identical.
set -euo pipefail

BASELINE_COMMIT=f81b5cbc
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
BASEDIR="$ALFDIR/../ALF_baseline"
MAKE_PARAMS="$HERE/make_params.py"
LOGDIR="$ALFDIR/runs"
OUTROOT="$LOGDIR/regress"
SHORT="NSweep=5 NBin=2"
BRANCH_EXTRA="${BRANCH_EXTRA:-}"
COMMON_EXTRA=("$@")

build() {  # build ALF.out in the tree $1, log to $2
    # configure.sh is not safe under set -eu, so relax them inside this subshell
    (set +eu; cd "$1" && source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 \
        && make program -j8 >"$2" 2>&1) \
        || { echo "build failed in $1, see $2" >&2; exit 2; }
}

# Parameter sets: "name|common keys|branch-only keys"
SETS=(
  "start||"
  # explicit: the new keys at their defaults must not change anything (branch-only: baseline cannot parse them)
  "explicit||Sequential_scheme=\"standard\" N_delay=32"
  "doped|ham_chem=-0.5 L1=4 L2=4 Beta=4|"
  "density|Mz=.F.|"
  "continuous|Continuous=.T.|"
  # tV: the Start default ham_V=4 aborts at baseline (unstable, positive V gives an imaginary
  # coupling), so use the parameters of testsuite/test_branch_parameters.json (ham_V=-1, N_SUN=1).
  "tv|ham_name=\"tV\" Model=\"tV\" ham_V=-1.d0 N_SUN=1 L1=4 L2=4 Beta=4|"
)

mkdir -p "$LOGDIR"

# Baseline worktree (detached, sources never modified).
if [ ! -d "$BASEDIR" ]; then
    git -C "$ALFDIR" worktree add --detach "$BASEDIR" "$BASELINE_COMMIT" >&2
fi
base_head="$(git -C "$BASEDIR" rev-parse HEAD)"
case "$base_head" in
    "$BASELINE_COMMIT"*) ;;
    *) echo "$BASEDIR is at $base_head, expected $BASELINE_COMMIT" >&2; exit 2 ;;
esac
[ -x "$BASEDIR/Prog/ALF.out" ] || build "$BASEDIR" "$LOGDIR/build_baseline.log"
build "$ALFDIR" "$LOGDIR/build_branch.log"

export OMP_NUM_THREADS=1
/bin/rm -rf "$OUTROOT"
status=0
for entry in "${SETS[@]}"; do
    IFS='|' read -r name common branch_only <<<"$entry"
    for side in base branch; do
        dir="$OUTROOT/$name/$side"
        keys="$SHORT $common"
        [ ${#COMMON_EXTRA[@]} -gt 0 ] && keys="$keys ${COMMON_EXTRA[*]}"
        exe="$BASEDIR/Prog/ALF.out"
        # Each side starts from its own tree's Start template: the branch template carries keys
        # (Sequential_scheme, N_delay) that the baseline binary cannot parse.
        template="$BASEDIR/Scripts_and_Parameters_files/Start"
        if [ "$side" = branch ]; then
            template="$ALFDIR/Scripts_and_Parameters_files/Start"
            keys="$keys $branch_only $BRANCH_EXTRA"
            exe="$ALFDIR/Prog/ALF.out"
        fi
        mkdir -p "$dir"
        read -r -a keyarr <<<"$keys"   # plain word split: quotes in values stay literal
        python3 -I "$MAKE_PARAMS" --template "$template" "$dir" "${keyarr[@]}"
        (cd "$dir" && "$exe" >stdout.log 2>&1) \
            || { echo "$name: $side run failed, see $dir/stdout.log" >&2; exit 2; }
        /bin/rm -f "$dir/stdout.log"
    done
    # Sanity: both runs must have produced outputs.
    for side in base branch; do
        d="$OUTROOT/$name/$side"
        if [ ! -s "$d/confout_0" ] || ! compgen -G "$d/*_scal" >/dev/null; then
            echo "$name: $side run produced no outputs (need confout_0 and *_scal) in $d" >&2
            exit 2
        fi
    done
    same=1
    # Compare outputs only: skip info* (commit/time stamps) and the inputs parameters, seeds
    # (they may differ textually when branch-only keys are used).
    files=$( (cd "$OUTROOT/$name/base" && ls; cd "$OUTROOT/$name/branch" && ls) | sort -u \
             | { grep -v -x -e 'parameters' -e 'seeds' -e 'info.*' || true; } )
    [ -n "$files" ] || { echo "$name: no output files to compare" >&2; exit 2; }
    for f in $files; do
        if ! cmp -s "$OUTROOT/$name/base/$f" "$OUTROOT/$name/branch/$f"; then
            same=0
            echo "  differs: $name/$f" >&2
        fi
    done
    if [ $same = 1 ]; then echo "$name IDENTICAL"; else echo "$name DIFFERENT"; status=1; fi
done
exit $status

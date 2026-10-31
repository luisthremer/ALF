#!/usr/bin/env bash
# Builds the timing-split binary (ALF with -DALF_TIMERS) in a separate copy of the tree.
#
# Usage:  testsuite/lowrank/timers_build.sh [--rebuild] [--mpi | --tempering]
#   --rebuild    rebuild even if the binary is newer than every source
#   --mpi        MPI build (configure.sh GNU MPI)       -> builds/ALF_timers_mpi.out
#   --tempering  TEMPERING build (configure.sh GNU Tempering) -> builds/ALF_timers_tempering.out
#   default      noMPI build (configure.sh GNU noMPI)   -> builds/ALF_timers.out
# Environment: MAKE_JOBS=n  parallel make jobs (default 2).
#
# The copy builds/timers[_mpi|_tempering]/src holds the working tree's files as they are now (tracked and untracked,
# not ignored: `git ls-files -co --exclude-standard`, i.e. uncommitted changes included), so Prog/*.o of
# the main tree are never touched (a flagged build there would poison test 36 and later incremental builds).
# The copy lies inside the repository, so Git_config.sh records the commit (with -dirty if the tree has
# uncommitted changes) and info carries the "This executable represents commit" line as usual.
# The binary is reused while it is newer than every *.F90, *.h, Makefile under Prog/ and Libraries/ of the
# main tree (generated files ignored, as compare_trajectory.sh) and the build flavour is the same.
# Prints the path of the binary on the last line. Exit 0 ok, 2 build failure or usage.
# A run of the binary writes timers.txt next to info (format: Prog/Timers_mod.F90, README.md).
# No `set -u` (configure.sh uses unset variables).
set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
REBUILD=0; MODE=noMPI
for a in "$@"; do
    case "$a" in
        --rebuild) REBUILD=1 ;;
        --mpi) MODE=MPI ;;
        --tempering) MODE=Tempering ;;
        *) sed -n '2,20p' "${BASH_SOURCE[0]}" >&2; exit 2 ;;
    esac
done
if [ $MODE = MPI ]; then
    COPY="$ALFDIR/builds/timers_mpi"; BIN="$ALFDIR/builds/ALF_timers_mpi.out"
elif [ $MODE = Tempering ]; then
    COPY="$ALFDIR/builds/timers_tempering"; BIN="$ALFDIR/builds/ALF_timers_tempering.out"
else
    COPY="$ALFDIR/builds/timers"; BIN="$ALFDIR/builds/ALF_timers.out"
fi
mkdir -p "$ALFDIR/builds" "$ALFDIR/runs"
LOG="$ALFDIR/runs/build_timers_$MODE.log"

stale() {
    [ -x "$BIN" ] || return 0
    s=$(find "$ALFDIR/Prog" "$ALFDIR/Libraries" -type f \( -name '*.F90' -o -name '*.h' -o -name Makefile \) \
        -newer "$BIN" ! -name 'git.h' ! -name 'git_status.h' ! -name 'Hamiltonians_case.h' \
        ! -name 'Hamiltonians_interface.h' ! -name '*_read_write_parameters.F90' | head -1)
    [ -n "$s" ]
}

if [ $REBUILD = 0 ] && ! stale; then
    echo "reusing $BIN"
    echo "$BIN"
    exit 0
fi

echo "building the timing-split binary (-DALF_TIMERS, $MODE) in $COPY ..."
/bin/rm -rf "$COPY"
mkdir -p "$COPY/src"
# Working-tree copy: tracked + untracked (not ignored) files that exist
( cd "$ALFDIR" && git ls-files -z -co --exclude-standard ) \
  | while IFS= read -r -d '' f; do [ -e "$ALFDIR/$f" ] && printf '%s\0' "$f"; done \
  | ( cd "$ALFDIR" && tar --null -T - -cf - ) | tar -xf - -C "$COPY/src"
(
    cd "$COPY/src" || exit 1
    export ALF_FLAGS_EXT="-DALF_TIMERS"
    # shellcheck disable=SC1091
    source configure.sh GNU "$MODE" NO-INTERACTIVE >/dev/null 2>&1
    make program -j"${MAKE_JOBS:-2}"
) >"$LOG" 2>&1 || { echo "build failed, see $LOG" >&2; exit 2; }
[ -x "$COPY/src/Prog/ALF.out" ] || { echo "build produced no binary, see $LOG" >&2; exit 2; }
grep -q -- "-DALF_TIMERS" "$LOG" || { echo "flag -DALF_TIMERS not in the compile lines, see $LOG" >&2; exit 2; }
cp "$COPY/src/Prog/ALF.out" "$BIN"
echo "$BIN"

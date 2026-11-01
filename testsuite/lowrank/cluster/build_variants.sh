#!/usr/bin/env bash
# Build every ALF variant the cluster validation (validate.sh) needs, on a login node (compiling only, no runs).
#
# Usage:  ALF_MACHINE=SuperMUC-NG ALF_ENV_SH=<workdir>/env_alf.sh testsuite/lowrank/cluster/build_variants.sh [VARIANT ...]
#   VARIANT  any of: normal trace ctest mpi tempering hdf5 debug timers timers_mpi  (default: all, in this order)
#   normal      Prog/ALF.out + Analysis/ana.out (configure.sh <machine> noMPI; used by ctest, regress, the trajectory harness)
#   trace       builds/ALF_trace.out (-DALF_TRACE_UPDATES, built in-tree as compare_trajectory.sh does; the normal build is restored)
#   mpi tempering hdf5 debug   archived copies of HEAD in builds/cluster/<variant>/src (configure.sh <machine> MPI | Tempering |
#               noMPI HDF5 | noMPI with ALF_FLAGS_EXT="-DALF_LOWRANK_DEBUG $ALF_DEBUG_EXTRA")
#   ctest       the CTest tree testsuite/tests (cmake + make against the normal build's Prog/*.o; compute nodes may lack make/cmake);
#               ALF_TEST_LDFLAGS = extra link flags of the tests (ifort: -qopt-matmul, the unit-test harness uses matmul)
#   timers timers_mpi   testsuite/lowrank/timers_build.sh [--mpi]
# Environment: ALF_MACHINE (default GNU), ALF_ENV_SH (module setup, sourced before configure.sh), MAKE_J (default 4),
#   ALF_DEBUG_EXTRA (extra flags of the debug build, e.g. -qopt-matmul- for ifort: the debug shadow uses matmul).
# Logs: runs/cluster_build/<variant>.log; one line "<variant> OK|FAILED <binary>" per variant on stdout; exit 1 if any failed.
# No `set -u` (configure.sh uses unset variables).
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../../.." && pwd)"
LOG="$ALFDIR/runs/cluster_build"
MAKE_J="${MAKE_J:-4}"
mkdir -p "$LOG" "$ALFDIR/builds/cluster"
cfg() {  # cfg FLAGS MODE...: module setup + configure.sh in the current shell; FLAGS -> ALF_FLAGS_EXT (exported before sourcing)
    local flags="$1"; shift
    { [ -z "${ALF_ENV_SH:-}" ] || source "$ALF_ENV_SH" || true; }
    if [ -n "$flags" ]; then export ALF_FLAGS_EXT="$flags"; else unset ALF_FLAGS_EXT; fi
    source configure.sh "${ALF_MACHINE:-GNU}" "$@" NO-INTERACTIVE
}
in_tree() {  # in_tree FLAGS: clean + build Prog/ALF.out and Analysis with ALF_FLAGS_EXT=FLAGS (empty = normal)
    ( cd "$ALFDIR" && cfg "$1" noMPI && make -C Prog clean && make all -j"$MAKE_J" )
}
copy_build() {  # copy_build VARIANT FLAGS MODE...: archived HEAD copy, built there
    local v="$1" flags="$2"; shift 2
    local d="$ALFDIR/builds/cluster/$v"
    /bin/rm -rf "$d"; mkdir -p "$d/src"
    git -C "$ALFDIR" archive HEAD | tar -x -C "$d/src"
    ( cd "$d/src" && cfg "$flags" "$@" && make all -j"$MAKE_J" )
    git -C "$ALFDIR" rev-parse HEAD >"$d/commit"
}
status=0
report() { if [ -s "$2" ]; then echo "$1 OK $2"; else echo "$1 FAILED (log $LOG/$1.log)"; status=1; fi; }
VARIANTS="${*:-normal trace ctest mpi tempering hdf5 debug timers timers_mpi}"
for v in $VARIANTS; do
    case $v in
        normal)     in_tree "" >"$LOG/$v.log" 2>&1; report $v "$ALFDIR/Prog/ALF.out" ;;
        trace)      mkdir -p "$ALFDIR/builds"; /bin/rm -f "$ALFDIR/builds/ALF_trace.out"
                    in_tree "-DALF_TRACE_UPDATES" >"$LOG/$v.log" 2>&1 && cp "$ALFDIR/Prog/ALF.out" "$ALFDIR/builds/ALF_trace.out"
                    in_tree "" >>"$LOG/$v.log" 2>&1
                    report $v "$ALFDIR/builds/ALF_trace.out" ;;
        ctest)      /bin/rm -f "$ALFDIR/testsuite/tests/ctest_built.ok"
                    ( cd "$ALFDIR" && cfg "" noMPI && { module load cmake >/dev/null 2>&1 || true; } && mkdir -p testsuite/tests \
                      && cd testsuite/tests && cmake .. -DCMAKE_EXE_LINKER_FLAGS="${ALF_TEST_LDFLAGS:-}" && make -j"$MAKE_J" \
                      && date >ctest_built.ok ) >"$LOG/$v.log" 2>&1
                    report $v "$ALFDIR/testsuite/tests/ctest_built.ok" ;;
        mpi)        copy_build $v "" MPI >"$LOG/$v.log" 2>&1; report $v "$ALFDIR/builds/cluster/$v/src/Prog/ALF.out" ;;
        tempering)  copy_build $v "" Tempering >"$LOG/$v.log" 2>&1; report $v "$ALFDIR/builds/cluster/$v/src/Prog/ALF.out" ;;
        hdf5)       copy_build $v "" noMPI HDF5 >"$LOG/$v.log" 2>&1; report $v "$ALFDIR/builds/cluster/$v/src/Prog/ALF.out" ;;
        debug)      copy_build $v "-DALF_LOWRANK_DEBUG ${ALF_DEBUG_EXTRA:-}" noMPI >"$LOG/$v.log" 2>&1
                    report $v "$ALFDIR/builds/cluster/$v/src/Prog/ALF.out" ;;
        timers)     "$ALFDIR/testsuite/lowrank/timers_build.sh" >"$LOG/$v.log" 2>&1; report $v "$ALFDIR/builds/ALF_timers.out" ;;
        timers_mpi) "$ALFDIR/testsuite/lowrank/timers_build.sh" --mpi >"$LOG/$v.log" 2>&1; report $v "$ALFDIR/builds/ALF_timers_mpi.out" ;;
        *)          echo "unknown variant $v" >&2; exit 2 ;;
    esac
done
exit $status

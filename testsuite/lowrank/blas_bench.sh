#!/usr/bin/env bash
# BLAS ceiling microbenchmark of one layout: builds blas_bench.F90 with ALF's compiler/flags/BLAS (configure.sh), runs one
# process per rank concurrently, writes ONE log with the provenance block, the per-rank result lines and the summary.
#
# Usage:  testsuite/lowrank/blas_bench.sh [--smoke] [--build-only] [KEY=VALUE ...]
#   LAYOUT=RxT       R ranks (concurrent processes) x T threads each (default 1x1); MKL_NUM_THREADS = OMP_NUM_THREADS =
#                    VECLIB_MAXIMUM_THREADS = T, MKL_DYNAMIC=FALSE (the profile's profile_setup may add pinning variables)
#   N=400,1600,3600  matrix sizes;  KS=16,64,256  inner dimensions of ZGEMM N x N x k
#   KERNELS=...      subset of zgemm_k,zgemm_n,zgeru,zgemv,zgeqp3;  REPS=7  WARMUP=2  MINTIME=0.02  CACHE_MB=256
#   TIMER_REPS=2000000  push/pop pairs per timer-cost sample (0 = skip)    (all of these are passed to blas_bench.out,
#                    whose header documents them: testsuite/lowrank/blas_bench.F90)
#   TAG=name         part of the log name (default <host>_<layout>_<date>);   LOG=file  the log (default runs/blas_bench/blas_bench_<TAG>.txt)
#   --smoke          fast self-check: N=64 REPS=2 WARMUP=1 KS=16,64 CACHE_MB=8 TIMER_REPS=20000 MINTIME=0.002 (< 30 s)
#   --build-only     build builds/blas_bench/blas_bench.out and stop
# Environment:
#   BLAS_BENCH_MACHINE  configure.sh machine name (default GNU; e.g. SUPERMUC-NG, HELMA...), always built "noMPI" (no MPI is
#                       used; the ranks are independent processes synchronised by a start-time file on a shared file system)
#   BLAS_BENCH_BIN      use this prebuilt blas_bench.out instead of building (provenance then lacks the build flags)
#   BLAS_BENCH_PROFILE  machine profile of testsuite/lowrank/campaign/profiles/<name>.sh (default local): TIME_STYLE,
#                       MPI_VERSION_CMD, EAR_SETTING, profile_setup (pinning/OpenMP variables), profile_launch
#   BLAS_BENCH_LAUNCH   how the R processes are started (default "local": R background processes on this machine, rank in
#                       PMI_RANK). "profile": the profile's profile_launch (srun ... on the clusters). Anything else is a
#                       command prefix started ONCE for the whole layout, @RANKS@ and @THREADS@ replaced, e.g.
#                       BLAS_BENCH_LAUNCH='srun -n @RANKS@ --cpus-per-task=@THREADS@ --cpu-bind=cores'
#                       BLAS_BENCH_LAUNCH='mpiexec -n @RANKS@'     (Intel MPI / MPICH; rank from PMI_RANK / PMIX_RANK /
#                       OMPI_COMM_WORLD_RANK / SLURM_PROCID, as campaign/rank_wrap.sh)
#   BLAS_BENCH_SYNC_DELAY  seconds between the launch and the common start time (default 3; use 20-30 under srun/mpiexec on
#                       many ranks); a process released late shows up as lag_s in its SYNC line and in the log (WARN)
#   BLAS_BENCH_MAXCORES refuse a layout with R x T above this many cores (a guard for shared machines)
#   FAKEINTEL=on|off, FAKEINTEL_LIB  as in the campaign (LD_PRELOAD only inside campaign/rank_wrap.sh)
# Log (one file, plain text, every line "KEY ..."): provenance "PROV key=value" lines (campaign/provenance.sh prov_job:
#   node, CPU, clock state, MPI, environment; plus build, binary and run facts, and the MKL_VERBOSE probe), then per rank
#   META, SYNC, TIMERCOST (ns per start/stop of the T1 timers, and of the clocks), RESULT (rank threads kernel N k reps inner
#   t_median_s t_min_s t_max_s gflops gflops_best gbs), then SUMMARY layout kernel N k nranks t_median_s gflops_med gflops_sum
#   gbs_med gbs_sum (median over the ranks; sum = node throughput) and TIMERCOST_SUMMARY layout what nranks t_median_ns.
# Exit 0 iff every rank finished and wrote at least one RESULT line. No `set -u` (configure.sh uses unset variables).
# Files: builds/blas_bench/ (copy of Prog/Timers_mod.F90 compiled with -DALF_TIMERS, blas_bench.out, build_info.txt),
#   runs/blas_bench/<TAG>.d (per rank output and provenance side files). Nothing is ever written into Prog/.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
CDIR="$HERE/campaign"
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"

# ---- per-rank mode (started by the launcher once per rank): rank from the MPI/Slurm variables, as campaign/rank_wrap.sh ----
if [ "$1" = "--rank-exec" ]; then
    outdir="$2"; startfile="$3"; threads="$4"; bin="$5"; shift 5
    rank="${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-0}}}}"
    "$CDIR/rank_wrap.sh" blas "$outdir" "$bin" "$@" RANK="$rank" THREADS="$threads" START="$startfile" \
        > "$outdir/rank_$rank.txt" 2> "$outdir/rank_$rank.err"
    exit $?
fi

SMOKE=0; BUILD_ONLY=0
LAYOUT=1x1; TAG=""; LOG=""
PASS=""      # KEY=VALUE arguments passed to blas_bench.out
for a in "$@"; do
    case "$a" in
        --smoke) SMOKE=1 ;;
        --build-only) BUILD_ONLY=1 ;;
        --help|-h) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
        LAYOUT=*) LAYOUT="${a#LAYOUT=}" ;;
        TAG=*) TAG="${a#TAG=}" ;;
        LOG=*) LOG="${a#LOG=}" ;;
        N=*|KS=*|KERNELS=*|REPS=*|WARMUP=*|MINTIME=*|CACHE_MB=*|TIMER_REPS=*) PASS="$PASS $a" ;;
        *) echo "blas_bench.sh: unknown argument $a" >&2; sed -n '2,12p' "${BASH_SOURCE[0]}" >&2; exit 2 ;;
    esac
done
if [ "$SMOKE" = 1 ]; then
    # smoke defaults first; explicit arguments (later on the command line of blas_bench.out) override them
    PASS="N=64 REPS=2 WARMUP=1 KS=16,64 CACHE_MB=8 TIMER_REPS=20000 MINTIME=0.002$PASS"
fi
case "$LAYOUT" in
    [1-9]*x[1-9]*) RANKS="${LAYOUT%%x*}"; THREADS="${LAYOUT##*x}" ;;
    *) echo "blas_bench.sh: LAYOUT must be RxT, got '$LAYOUT'" >&2; exit 2 ;;
esac
if [ -n "$BLAS_BENCH_MAXCORES" ] && [ $((RANKS * THREADS)) -gt "$BLAS_BENCH_MAXCORES" ]; then
    echo "blas_bench.sh: layout $LAYOUT needs $((RANKS * THREADS)) cores > BLAS_BENCH_MAXCORES=$BLAS_BENCH_MAXCORES" >&2; exit 2
fi

BDIR="$ALFDIR/builds/blas_bench"
BIN="${BLAS_BENCH_BIN:-$BDIR/blas_bench.out}"
BINFO="$BDIR/build_info.txt"

# ---- build (own directory, never Prog/) ----------------------------------------------------------------------------------
build() {
    mkdir -p "$BDIR"
    cp "$ALFDIR/Prog/Timers_mod.F90" "$BDIR/Timers_mod.F90"
    cp "$HERE/blas_bench.F90" "$BDIR/blas_bench.F90"
    # configure.sh works relative to $PWD and must be sourced in this shell (it exports ALF_FC, ALF_FLAGS_PROG, LIB_BLAS_LAPACK)
    cd "$ALFDIR" || return 1
    source configure.sh "${BLAS_BENCH_MACHINE:-GNU}" noMPI NO-INTERACTIVE >/dev/null 2>&1
    if [ -z "$ALF_FC" ] || [ -z "$ALF_FLAGS_PROG" ]; then echo "blas_bench.sh: configure.sh did not set ALF_FC / ALF_FLAGS_PROG" >&2; return 1; fi
    cd "$BDIR" || return 1
    /bin/rm -f Timers_mod.o timers_mod.mod blas_bench.o blas_bench.out
    local flags="$ALF_FLAGS_PROG -DALF_TIMERS"
    echo "build: $ALF_FC $flags" > "$BDIR/build.log"
    { $ALF_FC $flags -c Timers_mod.F90 \
        && $ALF_FC $flags -c blas_bench.F90 \
        && $ALF_FC $flags -o blas_bench.out Timers_mod.o blas_bench.o $LIB_BLAS_LAPACK; } >> "$BDIR/build.log" 2>&1 || {
        echo "blas_bench.sh: build failed, see $BDIR/build.log" >&2; tail -20 "$BDIR/build.log" >&2; return 1; }
    {
        echo "machine=${BLAS_BENCH_MACHINE:-GNU}"
        echo "compiler=$ALF_FC"
        echo "compiler_version=$($ALF_FC --version 2>&1 | head -1)"
        echo "flags=$flags"
        echo "blas_link=$LIB_BLAS_LAPACK"
        echo "timers_mod_sha256=$(. "$CDIR/provenance.sh"; prov_sha256 "$BDIR/Timers_mod.F90")"
        echo "blas_bench_f90_sha256=$(. "$CDIR/provenance.sh"; prov_sha256 "$BDIR/blas_bench.F90")"
        echo "built=$(date '+%Y-%m-%dT%H:%M:%S%z')"
    } > "$BINFO"
}
stale() {
    [ -x "$BIN" ] || return 0
    [ -n "$BLAS_BENCH_BIN" ] && return 1
    [ -r "$BINFO" ] || return 0
    cmp -s "$ALFDIR/Prog/Timers_mod.F90" "$BDIR/Timers_mod.F90" || return 0
    cmp -s "$HERE/blas_bench.F90" "$BDIR/blas_bench.F90" || return 0
    [ "${BLAS_BENCH_MACHINE:-GNU}" = "$(sed -n 's/^machine=//p' "$BINFO")" ] || return 0
    return 1
}
if stale; then build || exit 2; fi
[ -x "$BIN" ] || { echo "blas_bench.sh: binary $BIN not found" >&2; exit 2; }
if [ "$BUILD_ONLY" = 1 ]; then echo "$BIN"; exit 0; fi

. "$CDIR/provenance.sh"
PROFILE_NAME="${BLAS_BENCH_PROFILE:-local}"
. "$CDIR/profiles/$PROFILE_NAME.sh" || { echo "blas_bench.sh: no profile $PROFILE_NAME" >&2; exit 2; }
HOST="$(hostname -s 2>/dev/null || hostname)"
STAMP="$(date '+%Y%m%d-%H%M%S')"
[ -n "$TAG" ] || TAG="${HOST}_${LAYOUT}_${STAMP}"
RUNS="$ALFDIR/runs/blas_bench"
OUTDIR="$RUNS/$TAG.d"
[ -n "$LOG" ] || LOG="$RUNS/blas_bench_$TAG.txt"
if [ -e "$OUTDIR" ]; then echo "blas_bench.sh: $OUTDIR exists (choose another TAG)" >&2; exit 2; fi
mkdir -p "$OUTDIR" "$(dirname "$LOG")" || exit 2

# ---- environment of every rank ----------------------------------------------------------------------------------------------
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MKL_DYNAMIC=FALSE
unset MKL_VERBOSE LD_PRELOAD
FAKEINTEL="${FAKEINTEL:-$FAKEINTEL_DEFAULT}"
[ -n "$FAKEINTEL" ] || FAKEINTEL=off
export FAKEINTEL FAKEINTEL_LIB TIME_STYLE PROV_SYSFS_ROOT PROV_PROC_ROOT
if type profile_setup >/dev/null 2>&1; then profile_setup || { echo "blas_bench.sh: profile_setup failed" >&2; exit 2; }; fi

# ---- provenance block -----------------------------------------------------------------------------------------------------------
CAMPAIGN_ID="blas_bench"; JOB_NAME="blas_bench_$TAG"
PF="$OUTDIR/job.prov"
prov_job "$OUTDIR" || { echo "blas_bench.sh: prov_job failed" >&2; exit 2; }
prov_kv "$PF" layout "$LAYOUT"
prov_kv "$PF" ranks "$RANKS"
prov_kv "$PF" threads "$THREADS"
prov_kv "$PF" blas_bench_args "${PASS# }"
prov_kv "$PF" launcher "${BLAS_BENCH_LAUNCH:-local}"
prov_kv "$PF" sync_delay_s "${BLAS_BENCH_SYNC_DELAY:-3}"
prov_kv "$PF" fakeintel "$FAKEINTEL"
prov_kv "$PF" blas_bench_sh_sha256 "$(prov_sha256 "$SELF")"
prov_kv "$PF" blas_bench_f90_sha256 "$(prov_sha256 "$HERE/blas_bench.F90")"
prov_kv "$PF" binary "$BIN"
prov_kv "$PF" binary_sha256 "$(prov_sha256 "$BIN")"
if [ -r "$BINFO" ] && [ -z "$BLAS_BENCH_BIN" ]; then
    while IFS= read -r line; do prov_kv "$PF" "build_${line%%=*}" "${line#*=}"; done < "$BINFO"
else
    prov_kv "$PF" build_info "unavailable(prebuilt binary BLAS_BENCH_BIN, no build_info.txt)"
fi
if command -v git >/dev/null 2>&1 && [ -n "$(git -C "$ALFDIR" rev-parse --show-toplevel 2>/dev/null)" ]; then
    prov_kv "$PF" git_commit "$(git -C "$ALFDIR" rev-parse HEAD)"
    prov_kv "$PF" git_dirty "$(git -C "$ALFDIR" status --porcelain --untracked-files=no | wc -l | tr -d ' ') modified tracked files"
else
    prov_kv "$PF" git_commit "unavailable(not a git tree)"
fi
# what this machine can tell beyond prov_job (which reads lscpu/numactl/free, not on macOS)
prov_kv "$PF" uname "$(uname -a)"
if command -v sysctl >/dev/null 2>&1 && sysctl -n machdep.cpu.brand_string >/dev/null 2>&1; then
    prov_kv "$PF" cpu_model_sysctl "$(sysctl -n machdep.cpu.brand_string)"
    prov_kv "$PF" cpus_sysctl "$(sysctl -n hw.ncpu) logical, $(sysctl -n hw.physicalcpu) physical"
    prov_kv "$PF" mem_bytes_sysctl "$(sysctl -n hw.memsize)"
    prov_kv "$PF" perf_cores_sysctl "$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo 'unavailable(no perflevel0)')"
fi
if command -v otool >/dev/null 2>&1; then blibs="$(otool -L "$BIN" 2>&1 | sed '1d' | awk '{print $1}')"
elif command -v ldd >/dev/null 2>&1; then blibs="$(ldd "$BIN" 2>&1 | awk '{print $1 "->" $3}')"
else blibs="unavailable(no otool/ldd)"; fi
prov_kv "$PF" binary_libs "$blibs"
# other load on the machine (the numbers are only as good as the machine was quiet)
prov_kv "$PF" load_start "$(uptime | sed 's/.*load averages*: *//')"
prov_kv "$PF" busy_procs_start "$(ps -Ao pcpu=,comm= 2>/dev/null | sort -rn | head -6 | awk '{n = split($2, p, "/"); printf "%s%%%s ", $1, p[n]}')"
# MKL probe (as campaign/job_template.sh): the first MKL_VERBOSE version line of the --probe mode
mline="$(MKL_VERBOSE=1 "$BIN" --probe 2>&1 | grep -E -m1 'MKL_VERBOSE.*(oneMKL|Intel\(R\) MKL|Math Kernel Library)')"
if [ -n "$mline" ]; then
    prov_mkl_parse "$mline"
else
    MKL_VERSION="n/a (no MKL_VERBOSE line: BLAS is not MKL; see binary_libs)"; MKL_ISA="$MKL_VERSION"; MKL_THREADING="$MKL_VERSION"
fi
prov_kv "$PF" mkl_probe_line "${mline:-none}"
prov_kv "$PF" mkl_version "$MKL_VERSION"
prov_kv "$PF" mkl_isa "$MKL_ISA"
prov_kv "$PF" mkl_threading "$MKL_THREADING"
if [ -n "$mline" ]; then
    prov_mkl_threads_check "$THREADS" "$MKL_THREADING" || { echo "blas_bench.sh: THREADS=$THREADS but MKL runs sequential" >&2; exit 2; }
    prov_mkl_check "$FAKEINTEL" "$MKL_ISA" || { echo "blas_bench.sh: fakeintel=on but MKL took the $MKL_ISA path" >&2; exit 2; }
fi

# ---- launch: R processes released at a common wall-clock time ---------------------------------------------------------------------
DELAY="${BLAS_BENCH_SYNC_DELAY:-3}"
START=$(( $(date +%s) + DELAY ))
echo "$START" > "$OUTDIR/start_time.txt"
prov_kv "$PF" start_epoch_s "$START"
t_launch="$(prov_now)"
LAUNCH="${BLAS_BENCH_LAUNCH:-local}"
RC=0
RANKCMD="$SELF --rank-exec $OUTDIR $OUTDIR/start_time.txt $THREADS $BIN $PASS"
case "$LAUNCH" in
    local)
        pids=""
        r=0
        while [ $r -lt "$RANKS" ]; do
            PMI_RANK=$r $RANKCMD < /dev/null &
            pids="$pids $!"
            r=$((r + 1))
        done
        for p in $pids; do wait "$p" || RC=1; done ;;
    profile)
        profile_launch "$RANKS" "$THREADS" $RANKCMD < /dev/null || RC=1 ;;
    *)
        pre="${LAUNCH//@RANKS@/$RANKS}"; pre="${pre//@THREADS@/$THREADS}"
        $pre $RANKCMD < /dev/null || RC=1 ;;
esac
t_end="$(prov_now)"
prov_kv "$PF" run_wall_s "$(awk -v a="$t_launch" -v b="$t_end" 'BEGIN {printf "%.3f", b - a}')"
prov_kv "$PF" launch_rc "$RC"
prov_kv "$PF" load_end "$(uptime | sed 's/.*load averages*: *//')"
prov_kv "$PF" busy_procs_end "$(ps -Ao pcpu=,comm= 2>/dev/null | sort -rn | head -6 | awk '{n = split($2, p, "/"); printf "%s%%%s ", $1, p[n]}')"
prov_kv "$PF" rank_hosts "$(cat "$OUTDIR"/rank_blas_*.prov 2>/dev/null | sed -n 's/^host=//p' | sort -u | tr '\n' ' ')"
prov_kv "$PF" rank_affinity_preinit "$(for p in "$OUTDIR"/rank_blas_*.prov; do [ -r "$p" ] && printf 'r%s:%s ' "$(sed -n 's/^rank=//p' "$p")" "$(sed -n 's/^affinity=//p' "$p")"; done)"

# ---- the log -----------------------------------------------------------------------------------------------------------------------
{
    echo "# blas_bench.sh log (testsuite/lowrank/blas_bench.sh): layout $LAYOUT, tag $TAG"
    echo "# PROV key=value: provenance (campaign/provenance.sh); RESULT/TIMERCOST/SYNC/META per rank; SUMMARY lines at the end"
    sed 's/^/PROV /' "$PF"
    echo "PROV_END"
    r=0
    while [ $r -lt "$RANKS" ]; do
        f="$OUTDIR/rank_$r.txt"
        echo "# ---- rank $r"
        if [ -s "$f" ]; then cat "$f"; else echo "# rank $r: no output"; fi
        if [ -s "$OUTDIR/rank_$r.err" ]; then sed "s/^/# STDERR rank $r: /" "$OUTDIR/rank_$r.err"; fi
        r=$((r + 1))
    done
    cat "$OUTDIR"/rank_*.txt 2>/dev/null | awk -v layout="$LAYOUT" '
        function med(arr, n,   i, j, t, y) {
            for (i = 1; i <= n; i++) y[i] = arr[i]
            for (i = 2; i <= n; i++) { t = y[i]; j = i - 1; while (j >= 1 && y[j] > t) { y[j+1] = y[j]; j-- } y[j+1] = t }
            return (n % 2) ? y[(n+1)/2] : 0.5 * (y[n/2] + y[n/2+1])
        }
        $1 == "RESULT" {
            key = $4 " " $5 " " $6
            if (!(key in seen)) { seen[key] = 1; order[++nkeys] = key }
            n[key]++; i = n[key]
            tm[key, i] = $9; gf[key, i] = $12; gb[key, i] = $14
            gfs[key] += $12; gbs[key] += $14
        }
        $1 == "TIMERCOST" {
            key = $4
            if (!(key in tseen)) { tseen[key] = 1; torder[++ntk] = key }
            tn[key]++; tv[key, tn[key]] = $7
        }
        $1 == "SYNC" { if ($5 + 0 > maxlag) maxlag = $5 + 0; nsync++ }
        END {
            for (q = 1; q <= nkeys; q++) {
                key = order[q]; m = n[key]
                for (i = 1; i <= m; i++) { a[i] = tm[key, i]; b[i] = gf[key, i]; c[i] = gb[key, i] }
                tmed = med(a, m); gmed = med(b, m); cmed = med(c, m)
                printf "SUMMARY %s %s %d %.6E %.3f %.3f %.3f %.3f\n", layout, key, m, tmed, gmed, gfs[key], cmed, gbs[key]
            }
            for (q = 1; q <= ntk; q++) {
                key = torder[q]; m = tn[key]
                for (i = 1; i <= m; i++) a[i] = tv[key, i]
                printf "TIMERCOST_SUMMARY %s %s %d %.2f\n", layout, key, m, med(a, m)
            }
            printf "SYNC_SUMMARY %s ranks_synced %d max_lag_s %.3f\n", layout, nsync, maxlag
            if (maxlag > 0.5) print "# WARN: a rank was released more than 0.5 s after the common start time (late launch): the ranks did not run concurrently; raise BLAS_BENCH_SYNC_DELAY"
        }'
    nres="$(cat "$OUTDIR"/rank_*.txt 2>/dev/null | grep -c '^RESULT ')"
    echo "# end: launch rc=$RC, $nres RESULT lines from $RANKS ranks"
} > "$LOG"
echo "$LOG"
grep -E '^(SUMMARY|TIMERCOST_SUMMARY|SYNC_SUMMARY|# WARN|# end)' "$LOG"
[ "$RC" -eq 0 ] && [ "$(grep -c '^RESULT ' "$LOG")" -ge "$RANKS" ] || exit 1
exit 0

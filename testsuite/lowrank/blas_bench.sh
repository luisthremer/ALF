#!/usr/bin/env bash
# BLAS ceiling microbenchmark of one layout: builds blas_bench.F90 with ALF's compiler/flags/BLAS (configure.sh), runs one
# process per rank concurrently (all ranks measure the same point at the same time), writes ONE log with the provenance
# block, the per-rank result lines and the summary.
#
# Usage:  testsuite/lowrank/blas_bench.sh [--smoke] [--build-only] [--plan] [KEY=VALUE ...]
#   LAYOUT=RxT       R ranks (concurrent processes) x T threads each (default 1x1); MKL_NUM_THREADS = OMP_NUM_THREADS =
#                    VECLIB_MAXIMUM_THREADS = T, MKL_DYNAMIC=FALSE (the profile's profile_setup may add pinning variables)
#   N=400,1600,3600  matrix sizes;  KS=16,64,256  inner dimensions of ZGEMM N x N x k
#   KERNELS=...      subset of zgemm_k,zgemm_n,zgeru,zgemv,zgeqp3;  REPS=7  WARMUP=2  MINTIME=0.02
#   CACHE_MB=256     floor of the rotating-buffer size per rank (see "Cache"); TIMER_REPS=2000000  push/pop pairs per timer-cost
#                    sample (0 = skip)    (all of these are passed to the binary, whose header documents them:
#                    testsuite/lowrank/blas_bench.F90)
#   TAG=name         part of the log name (default <host>_<layout>_<date>);   LOG=file  the log (default runs/blas_bench/blas_bench_<TAG>.txt)
#   --smoke          fast self-check: N=64 REPS=2 WARMUP=1 KS=16,64 CACHE_MB=8 TIMER_REPS=20000 MINTIME=0.002 (< 30 s)
#   --build-only     build builds/blas_bench/blas_bench.out and blas_bench_mpi.out and stop
#   --plan           print the sizing (cache, memory, SKIP lines) and stop
# Binaries (builds/blas_bench/): blas_bench.out (noMPI; launcher "local") and blas_bench_mpi.out (-DBLAS_BENCH_MPI, built with
#   the MPI wrapper of `configure.sh <machine> MPI`: MPI_Init, so Intel MPI / srun pin the ranks as for ALF, MPI_Barrier before every
#   point). The launchers `profile` and any srun/mpiexec prefix REQUIRE the MPI binary (error otherwise).
# Synchronisation: before EVERY point (timer cost, each kernel x N x k) all ranks wait for each other: MPI_Barrier (MPI binary) or
#   ready files runs/blas_bench/<TAG>.d/ready/ready_<point>_<rank> (local launcher; a timeout of BLAS_BENCH_READY_TIMEOUT s = 900
#   FAILS the run). The first point is released at a common wall-clock time (start-time file; a rank released > 0.5 s late FAILS).
#   Each RESULT line carries the release time; SUMMARY has the spread (skew_s; > 1 s FAILS).
# Pinning: local launcher on Linux sets, for rank r, OMP_PLACES="{r*T}:T" (cores r*T .. r*T+T-1) and OMP_PROC_BIND=close; on macOS no
#   pinning exists (logged as pinning=none). With MPI launchers the MPI library / srun pins. OMP_DISPLAY_AFFINITY=TRUE: the runtime's
#   affinity lines (stderr of each rank) are kept in the log as "# STDERR rank r: ...".
# Memory: per rank peak(N) = 16 B x max( 3 N^2 [zgemm_n], N^2 + 2 N kmax [zgemm_k], nbuf N^2 [zgeru/zgemv], 2 N^2 + 64 N [zgeqp3] )
#   + 200 MB, nbuf = max(1, ceil(cache/(16 N^2))), over the selected kernels. N is dropped when R x peak(N) > 0.8 x available memory
#   (MemAvailable of /proc/meminfo, or hw.memsize on the Mac; BLAS_BENCH_MEM_AVAIL_MB overrides) with the log line
#   "SKIP N=<n> reason=memory est_gb=<R x peak> avail_gb=<available>" (never silently).
# Cache: per-rank rotating-buffer size = max(CACHE_MB, 4 x node L3 total / R) MB (L3 total: sysfs cache/index*/ with level 3, summed
#   over the distinct shared_cpu_list groups, else lscpu, else sysctl hw.l3cachesize, else 0), BLAS_BENCH_CACHE_MB overrides the whole
#   expression; capped at 25 % of 0.8 x available memory / R (logged). The value used is logged (cache_mb_used).
# Environment:
#   BLAS_BENCH_MACHINE  configure.sh machine name (default GNU; e.g. SUPERMUC-NG, HELMA...)
#   BLAS_BENCH_BIN, BLAS_BENCH_BIN_MPI   use these prebuilt binaries instead of building (provenance then lacks the build flags)
#   BLAS_BENCH_PROFILE  machine profile of testsuite/lowrank/campaign/profiles/<name>.sh (default local): TIME_STYLE,
#                       MPI_VERSION_CMD, EAR_SETTING, profile_setup (pinning/OpenMP variables), profile_launch
#   BLAS_BENCH_LAUNCH   how the R processes are started (default "local": R background processes on this machine, rank in
#                       PMI_RANK, noMPI binary). "profile": the profile's profile_launch (srun ... on the clusters). Anything else
#                       is a command prefix started ONCE for the whole layout, @RANKS@ and @THREADS@ replaced, e.g.
#                       BLAS_BENCH_LAUNCH='srun -n @RANKS@ --cpus-per-task=@THREADS@ --cpu-bind=cores'
#                       BLAS_BENCH_LAUNCH='mpiexec -n @RANKS@'     (rank from PMI_RANK / PMIX_RANK / OMPI_COMM_WORLD_RANK /
#                       SLURM_PROCID, as campaign/rank_wrap.sh; on the Mac Homebrew mpiexec needs FI_PROVIDER=tcp)
#   BLAS_BENCH_SYNC_DELAY  seconds between the launch and the common first start time (default 3 for local, 30 otherwise)
#   BLAS_BENCH_READY_TIMEOUT, BLAS_BENCH_MAXCORES (refuse R x T above this), BLAS_BENCH_REQUIRE_MKL (default 0 for the local
#                       profile, else 1: abort when the probe shows no MKL_VERBOSE line)
#   FAKEINTEL=on|off, FAKEINTEL_LIB  as in the campaign (LD_PRELOAD only inside campaign/rank_wrap.sh, also for the MKL probe)
# Log (one file, plain text, every line "KEY ..."): provenance "PROV key=value" lines (campaign/provenance.sh prov_job: node, CPU,
#   clock state, MPI, environment; plus build, binary, sizing and run facts, and the MKL_VERBOSE probe run like job_template.sh
#   through the launcher and rank_wrap.sh), SKIP lines, then per rank META, SYNC, TIMERCOST (ns per start/stop of the T1 timers and
#   of the clocks), RESULT (rank threads kernel N k reps inner t_median_s t_min_s t_max_s gflops gflops_best gbs t_release nbuf),
#   then SUMMARY layout kernel N k nranks t_median_s gflops_med gflops_sum gbs_med gbs_sum skew_s (median over the ranks; sum = node
#   throughput) and TIMERCOST_SUMMARY layout what nranks t_median_ns.
# Exit 0 iff every rank 0..R-1 wrote exactly one output with the full expected RESULT count (kernels x N x k minus SKIPs), R SYNC
#   lines, R distinct rank records, no late release and no large skew; otherwise non-zero with a message and "# FAIL" lines.
# No `set -u` (configure.sh uses unset variables). The rank command is built as a bash array (paths with spaces are fine).
# Files: builds/blas_bench/ (copies of Prog/Timers_mod.F90 and blas_bench.F90 compiled with -DALF_TIMERS; mpi/ for the MPI objects,
#   blas_bench.out, blas_bench_mpi.out, build_info*.txt), runs/blas_bench/<TAG>.d (per rank output, provenance side files).
#   Nothing is ever written into Prog/.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
CDIR="$HERE/campaign"
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"

# ---- per-rank mode (started by the launcher once per rank): rank from the MPI/Slurm variables, as campaign/rank_wrap.sh ----
if [ "$1" = "--rank-exec" ]; then
    outdir="$2"; startfile="$3"; threads="$4"; nranks="$5"; bin="$6"; shift 6
    rank="${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-0}}}}"
    if [ "$BLAS_BENCH_PIN" = linux ]; then   # local launcher on Linux: rank r on cores r*T .. r*T+T-1
        export OMP_PLACES="{$((rank * threads))}:$threads" OMP_PROC_BIND=close
    fi
    # output names unique even if the launcher gives every process the same rank (detected by the collector)
    "$CDIR/rank_wrap.sh" blas "$outdir" "$bin" "$@" RANK="$rank" NRANKS="$nranks" THREADS="$threads" START="$startfile" \
        > "$outdir/out_${rank}_p$$.txt" 2> "$outdir/out_${rank}_p$$.err"
    rc=$?
    [ "$rc" -eq 0 ] || { [ -d "$outdir/ready" ] && : > "$outdir/ready/ABORT"; }   # the other ranks stop waiting at the barrier
    exit $rc
fi

SMOKE=0; BUILD_ONLY=0; PLAN_ONLY=0
LAYOUT=1x1; TAG=""; LOG=""
V_N=400,1600,3600; V_KS=16,64,256; V_KERNELS=zgemm_k,zgemm_n,zgeru,zgemv,zgeqp3
V_REPS=""; V_WARMUP=""; V_MINTIME=""; V_CACHE_MB=256; V_TIMER_REPS=""
for a in "$@"; do
    case "$a" in
        --smoke) SMOKE=1 ;;
        --build-only) BUILD_ONLY=1 ;;
        --plan) PLAN_ONLY=1 ;;
        --help|-h) sed -n '2,62p' "${BASH_SOURCE[0]}"; exit 0 ;;
        LAYOUT=*) LAYOUT="${a#LAYOUT=}" ;;
        TAG=*) TAG="${a#TAG=}" ;;
        LOG=*) LOG="${a#LOG=}" ;;
        N=*) V_N="${a#N=}" ;;
        KS=*) V_KS="${a#KS=}" ;;
        KERNELS=*) V_KERNELS="${a#KERNELS=}" ;;
        REPS=*) V_REPS="${a#REPS=}" ;;
        WARMUP=*) V_WARMUP="${a#WARMUP=}" ;;
        MINTIME=*) V_MINTIME="${a#MINTIME=}" ;;
        CACHE_MB=*) V_CACHE_MB="${a#CACHE_MB=}" ;;
        TIMER_REPS=*) V_TIMER_REPS="${a#TIMER_REPS=}" ;;
        *) echo "blas_bench.sh: unknown argument $a" >&2; sed -n '2,12p' "${BASH_SOURCE[0]}" >&2; exit 2 ;;
    esac
done
if [ "$SMOKE" = 1 ]; then   # smoke defaults only where the command line did not say otherwise
    case " $* " in *" N="*) ;; *) V_N=64 ;; esac
    case " $* " in *" KS="*) ;; *) V_KS=16,64 ;; esac
    [ -n "$V_REPS" ] || V_REPS=2
    [ -n "$V_WARMUP" ] || V_WARMUP=1
    [ -n "$V_MINTIME" ] || V_MINTIME=0.002
    [ -n "$V_TIMER_REPS" ] || V_TIMER_REPS=20000
    case " $* " in *" CACHE_MB="*) ;; *) V_CACHE_MB=8 ;; esac
fi
case "$LAYOUT" in
    *[!0-9x]*|x*|*x|*x*x*|"") echo "blas_bench.sh: LAYOUT must be RxT (digits), got '$LAYOUT'" >&2; exit 2 ;;
    *x*) RANKS="${LAYOUT%%x*}"; THREADS="${LAYOUT##*x}" ;;
    *) echo "blas_bench.sh: LAYOUT must be RxT, got '$LAYOUT'" >&2; exit 2 ;;
esac
if [ "$RANKS" -lt 1 ] || [ "$THREADS" -lt 1 ]; then echo "blas_bench.sh: LAYOUT needs R >= 1 and T >= 1" >&2; exit 2; fi
if [ -n "$BLAS_BENCH_MAXCORES" ] && [ $((RANKS * THREADS)) -gt "$BLAS_BENCH_MAXCORES" ]; then
    echo "blas_bench.sh: layout $LAYOUT needs $((RANKS * THREADS)) cores > BLAS_BENCH_MAXCORES=$BLAS_BENCH_MAXCORES" >&2; exit 2
fi
LAUNCH="${BLAS_BENCH_LAUNCH:-local}"

BDIR="$ALFDIR/builds/blas_bench"
BIN="${BLAS_BENCH_BIN:-$BDIR/blas_bench.out}"
BIN_MPI="${BLAS_BENCH_BIN_MPI:-$BDIR/blas_bench_mpi.out}"

# ---- build (own directory, never Prog/); flavour noMPI or mpi, each in a subshell (configure.sh exports a lot) ---------------
build() {   # build noMPI|mpi
    local flavour="$1" dir="$BDIR" cfg=noMPI extra="" out="$BDIR/blas_bench.out" info="$BDIR/build_info.txt"
    if [ "$flavour" = mpi ]; then
        dir="$BDIR/mpi"; cfg=MPI; extra="-DBLAS_BENCH_MPI"; out="$BDIR/blas_bench_mpi.out"; info="$BDIR/build_info_mpi.txt"
    fi
    mkdir -p "$dir"
    (
        cp "$ALFDIR/Prog/Timers_mod.F90" "$dir/Timers_mod.F90"
        cp "$HERE/blas_bench.F90" "$dir/blas_bench.F90"
        # configure.sh works relative to $PWD and must be sourced in this shell (exports ALF_FC, ALF_FLAGS_PROG, LIB_BLAS_LAPACK)
        cd "$ALFDIR" || exit 1
        source configure.sh "${BLAS_BENCH_MACHINE:-GNU}" $cfg NO-INTERACTIVE >/dev/null 2>&1
        if [ -z "$ALF_FC" ] || [ -z "$ALF_FLAGS_PROG" ]; then echo "blas_bench.sh: configure.sh ($cfg) did not set ALF_FC / ALF_FLAGS_PROG" >&2; exit 1; fi
        cd "$dir" || exit 1
        /bin/rm -f Timers_mod.o timers_mod.mod blas_bench.o "$out"
        flags="$ALF_FLAGS_PROG -DALF_TIMERS $extra"
        echo "build ($flavour): $ALF_FC $flags" > "$dir/build.log"
        { $ALF_FC $flags -c Timers_mod.F90 \
            && $ALF_FC $flags -c blas_bench.F90 \
            && $ALF_FC $flags -o "$out" Timers_mod.o blas_bench.o $LIB_BLAS_LAPACK; } >> "$dir/build.log" 2>&1 || {
            echo "blas_bench.sh: $flavour build failed, see $dir/build.log" >&2; tail -20 "$dir/build.log" >&2; exit 1; }
        {
            echo "flavour=$flavour"
            echo "machine=${BLAS_BENCH_MACHINE:-GNU}"
            echo "compiler=$ALF_FC"
            echo "compiler_version=$($ALF_FC --version 2>&1 | head -1)"
            echo "flags=$flags"
            echo "blas_link=$LIB_BLAS_LAPACK"
            echo "timers_mod_sha256=$(. "$CDIR/provenance.sh"; prov_sha256 "$dir/Timers_mod.F90")"
            echo "blas_bench_f90_sha256=$(. "$CDIR/provenance.sh"; prov_sha256 "$dir/blas_bench.F90")"
            echo "built=$(date '+%Y-%m-%dT%H:%M:%S%z')"
        } > "$info"
    )
}
stale() {   # stale noMPI|mpi
    local dir="$BDIR" b="$BIN" info="$BDIR/build_info.txt" pre=""
    if [ "$1" = mpi ]; then dir="$BDIR/mpi"; b="$BIN_MPI"; info="$BDIR/build_info_mpi.txt"; fi
    [ -x "$b" ] || return 0
    if [ "$1" = noMPI ] && [ -n "$BLAS_BENCH_BIN" ]; then return 1; fi
    if [ "$1" = mpi ] && [ -n "$BLAS_BENCH_BIN_MPI" ]; then return 1; fi
    [ -r "$info" ] || return 0
    cmp -s "$ALFDIR/Prog/Timers_mod.F90" "$dir/Timers_mod.F90" || return 0
    cmp -s "$HERE/blas_bench.F90" "$dir/blas_bench.F90" || return 0
    [ "${BLAS_BENCH_MACHINE:-GNU}" = "$(sed -n 's/^machine=//p' "$info")" ] || return 0
    return 1
}
NEED_MPI=0
[ "$LAUNCH" != local ] && NEED_MPI=1
[ "$BUILD_ONLY" = 1 ] && NEED_MPI=1
if [ "$PLAN_ONLY" != 1 ]; then
    if stale noMPI; then build noMPI || exit 2; fi
    [ -x "$BIN" ] || { echo "blas_bench.sh: binary $BIN not found" >&2; exit 2; }
    if [ "$NEED_MPI" = 1 ]; then
        if stale mpi; then build mpi || exit 2; fi
        [ -x "$BIN_MPI" ] || { echo "blas_bench.sh: MPI binary $BIN_MPI not found (the launcher '$LAUNCH' needs it)" >&2; exit 2; }
    fi
    if [ "$BUILD_ONLY" = 1 ]; then echo "$BIN"; echo "$BIN_MPI"; exit 0; fi
fi

. "$CDIR/provenance.sh"
PROFILE_NAME="${BLAS_BENCH_PROFILE:-local}"
. "$CDIR/profiles/$PROFILE_NAME.sh" || { echo "blas_bench.sh: no profile $PROFILE_NAME" >&2; exit 2; }

# ---- sizing: memory guard and rotating-buffer size -----------------------------------------------------------------------------
mem_avail_mb() {
    local v
    if [ -n "$BLAS_BENCH_MEM_AVAIL_MB" ]; then echo "$BLAS_BENCH_MEM_AVAIL_MB"; return; fi
    if [ -r "${PROV_PROC_ROOT}/proc/meminfo" ]; then
        v="$(awk '/^MemAvailable:/ {printf "%d", $2 / 1024}' "${PROV_PROC_ROOT}/proc/meminfo")"
        [ -n "$v" ] && { echo "$v"; return; }
    fi
    if command -v sysctl >/dev/null 2>&1 && v="$(sysctl -n hw.memsize 2>/dev/null)" && [ -n "$v" ]; then
        echo $((v / 1048576)); return
    fi
    echo ""
}
# node L3 total in MB (0 if unknown)
l3_total_mb() {
    local base="${PROV_SYSFS_ROOT}/sys/devices/system/cpu" d lvl v tot
    if [ -n "$BLAS_BENCH_L3_MB" ]; then echo "$BLAS_BENCH_L3_MB"; return; fi
    tot="$(for d in "$base"/cpu[0-9]*/cache/index*; do
            [ -r "$d/level" ] && [ "$(head -1 "$d/level")" = 3 ] && [ -r "$d/shared_cpu_list" ] && [ -r "$d/size" ] || continue
            echo "$(head -1 "$d/shared_cpu_list") $(head -1 "$d/size")"
        done 2>/dev/null | sort -u | awk '{v = $2; m = 1; if (v ~ /[Kk]$/) m = 1 / 1024; else if (v ~ /[Gg]$/) m = 1024; sub(/[KkMmGg]$/, "", v); s += v * m} END {if (NR) printf "%d", s}')"
    if [ -n "$tot" ]; then echo "$tot"; return; fi
    if command -v lscpu >/dev/null 2>&1; then
        # util-linux >= 2.34: "L3 cache: 64 MiB (2 instances)" is the total; older: per instance in K, times Socket(s)
        tot="$(lscpu 2>/dev/null | awk '
            /^L3 cache:/ {v = $3; u = $4; if (u ~ /^K/) m = 1 / 1024; else if (u ~ /^G/) m = 1024; else m = 1; tot = v * m; inst = ($0 ~ /instance/)}
            /^Socket\(s\):/ {sock = $2}
            END {if (tot > 0) {if (!inst && sock > 0) tot *= sock; printf "%d", tot}}')"
        if [ -n "$tot" ]; then echo "$tot"; return; fi
    fi
    if command -v sysctl >/dev/null 2>&1 && v="$(sysctl -n hw.l3cachesize 2>/dev/null)" && [ -n "$v" ]; then echo $((v / 1048576)); return; fi
    echo 0
}
L3_MB="$(l3_total_mb)"
MEM_MB="$(mem_avail_mb)"
if [ -n "$BLAS_BENCH_CACHE_MB" ]; then
    CACHE_USED="$BLAS_BENCH_CACHE_MB"; CACHE_NOTE="BLAS_BENCH_CACHE_MB"
else
    CACHE_USED="$(awk -v f="$V_CACHE_MB" -v l3="$L3_MB" -v r="$RANKS" 'BEGIN {c = 4 * l3 / r; if (f > c) c = f; printf "%d", c + 0.999}')"
    CACHE_NOTE="max(CACHE_MB=$V_CACHE_MB, 4 x L3 $L3_MB MB / $RANKS ranks)"
fi
if [ -n "$MEM_MB" ]; then
    cap="$(awk -v m="$MEM_MB" -v r="$RANKS" 'BEGIN {printf "%d", 0.25 * 0.8 * m / r}')"
    if [ "$CACHE_USED" -gt "$cap" ]; then CACHE_NOTE="$CACHE_NOTE, capped from $CACHE_USED MB at 25% of 0.8 x memory / ranks"; CACHE_USED="$cap"; fi
fi
[ "$CACHE_USED" -ge 1 ] || CACHE_USED=1
KMAX="$(echo "$V_KS" | tr ',' '\n' | sort -n | tail -1)"
sel() { case ",$V_KERNELS," in *",$1,"*) return 0 ;; *) return 1 ;; esac; }
SKIPS=""; NLIST=""
for n in $(echo "$V_N" | tr ',' ' '); do
    peak="$(awk -v n="$n" -v kmax="$KMAX" -v cache="$CACHE_USED" -v zk="$(sel zgemm_k && echo 1 || echo 0)" \
        -v zn="$(sel zgemm_n && echo 1 || echo 0)" -v b2="$( (sel zgeru || sel zgemv) && echo 1 || echo 0)" \
        -v qp="$(sel zgeqp3 && echo 1 || echo 0)" 'BEGIN {
            m = 16 * n * n; nbuf = int(cache * 1048576 / m); if (nbuf * m < cache * 1048576) nbuf++; if (nbuf < 1) nbuf = 1
            p = 0
            if (zk && 16 * (n * n + 2 * n * kmax) > p) p = 16 * (n * n + 2 * n * kmax)
            if (zn && 3 * m > p) p = 3 * m
            if (b2 && nbuf * m > p) p = nbuf * m
            if (qp && 2 * m + 16 * 64 * n > p) p = 2 * m + 16 * 64 * n
            printf "%.1f", p / 1048576 + 200 }')"
    if [ -n "$MEM_MB" ] && awk -v p="$peak" -v r="$RANKS" -v m="$MEM_MB" 'BEGIN {exit !(r * p > 0.8 * m)}'; then
        SKIPS="$SKIPS
SKIP N=$n reason=memory est_gb=$(awk -v p="$peak" -v r="$RANKS" 'BEGIN {printf "%.1f", r * p / 1024}') avail_gb=$(awk -v m="$MEM_MB" 'BEGIN {printf "%.1f", m / 1024}')"
    else
        NLIST="$NLIST,$n"
    fi
done
NLIST="${NLIST#,}"
if [ "$PLAN_ONLY" = 1 ]; then
    echo "layout=$LAYOUT l3_total_mb=$L3_MB mem_avail_mb=${MEM_MB:-unknown} cache_mb_used=$CACHE_USED ($CACHE_NOTE) n_run=${NLIST:-none}"
    [ -n "$SKIPS" ] && echo "${SKIPS#
}"
    exit 0
fi
if [ -z "$NLIST" ]; then echo "blas_bench.sh: no N fits into the memory:${SKIPS}" >&2; exit 2; fi
# arguments of the binary
FARGS=("N=$NLIST" "KS=$V_KS" "KERNELS=$V_KERNELS" "CACHE_MB=$CACHE_USED")
[ -n "$V_REPS" ] && FARGS+=("REPS=$V_REPS")
[ -n "$V_WARMUP" ] && FARGS+=("WARMUP=$V_WARMUP")
[ -n "$V_MINTIME" ] && FARGS+=("MINTIME=$V_MINTIME")
[ -n "$V_TIMER_REPS" ] && FARGS+=("TIMER_REPS=$V_TIMER_REPS")
# expected RESULT lines per rank
NNN="$(echo "$NLIST" | tr ',' '\n' | wc -l | tr -d ' ')"
NKK="$(echo "$V_KS" | tr ',' '\n' | wc -l | tr -d ' ')"
PER_N=0
sel zgemm_k && PER_N=$((PER_N + NKK))
for k in zgemm_n zgeru zgemv zgeqp3; do sel $k && PER_N=$((PER_N + 1)); done
EXPECT=$((NNN * PER_N))
HAS_TIMERCOST=1; [ "$V_TIMER_REPS" = 0 ] && HAS_TIMERCOST=0

HOST="$(hostname -s 2>/dev/null || hostname)"
STAMP="$(date '+%Y%m%d-%H%M%S')"
[ -n "$TAG" ] || TAG="${HOST}_${LAYOUT}_${STAMP}"
RUNS="$ALFDIR/runs/blas_bench"
OUTDIR="$RUNS/$TAG.d"
[ -n "$LOG" ] || LOG="$RUNS/blas_bench_$TAG.txt"
if [ -e "$OUTDIR" ]; then echo "blas_bench.sh: $OUTDIR exists (choose another TAG)" >&2; exit 2; fi
mkdir -p "$OUTDIR" "$OUTDIR/ready" "$(dirname "$LOG")" || exit 2

# ---- environment of every rank ----------------------------------------------------------------------------------------------
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MKL_DYNAMIC=FALSE
export OMP_DISPLAY_AFFINITY=TRUE OMP_AFFINITY_FORMAT='OMP_AFFINITY host=%H pid=%P thread=%n affinity=%A'
unset MKL_VERBOSE LD_PRELOAD
FAKEINTEL="${FAKEINTEL:-$FAKEINTEL_DEFAULT}"
[ -n "$FAKEINTEL" ] || FAKEINTEL=off
export FAKEINTEL FAKEINTEL_LIB TIME_STYLE PROV_SYSFS_ROOT PROV_PROC_ROOT
if type profile_setup >/dev/null 2>&1; then profile_setup || { echo "blas_bench.sh: profile_setup failed" >&2; exit 2; }; fi
BLAS_BENCH_PIN=none
PIN_NOTE="by the launcher ($LAUNCH): MPI library / srun binding, see pin_env"
if [ "$LAUNCH" = local ]; then
    if [ "$(uname -s)" = Linux ]; then
        BLAS_BENCH_PIN=linux
        PIN_NOTE="explicit per rank: OMP_PLACES={r*T}:T, OMP_PROC_BIND=close (local launcher, Linux)"
    else
        PIN_NOTE="none: $(uname -s) has no thread/core pinning that the runtime honours (local launcher)"
    fi
fi
export BLAS_BENCH_PIN

# ---- provenance block -----------------------------------------------------------------------------------------------------------
CAMPAIGN_ID="blas_bench"; JOB_NAME="blas_bench_$TAG"
PF="$OUTDIR/job.prov"
prov_job "$OUTDIR" || { echo "blas_bench.sh: prov_job failed" >&2; exit 2; }
if [ "$LAUNCH" = local ]; then RBIN="$BIN"; else RBIN="$BIN_MPI"; fi
prov_kv "$PF" layout "$LAYOUT"
prov_kv "$PF" ranks "$RANKS"
prov_kv "$PF" threads "$THREADS"
prov_kv "$PF" blas_bench_args "${FARGS[*]}"
prov_kv "$PF" launcher "$LAUNCH"
prov_kv "$PF" mpi_binary "$([ "$LAUNCH" = local ] && echo 'no (noMPI binary, ready-file barrier)' || echo 'yes (MPI_Init, MPI_Barrier before every point)')"
prov_kv "$PF" pinning "$PIN_NOTE"
prov_kv "$PF" cache_mb_used "$CACHE_USED ($CACHE_NOTE)"
prov_kv "$PF" l3_total_mb "$L3_MB"
prov_kv "$PF" mem_avail_mb "${MEM_MB:-unavailable(no /proc/meminfo or sysctl)}"
prov_kv "$PF" expected_results_per_rank "$EXPECT"
prov_kv "$PF" fakeintel "$FAKEINTEL"
prov_kv "$PF" blas_bench_sh_sha256 "$(prov_sha256 "$SELF")"
prov_kv "$PF" blas_bench_f90_sha256 "$(prov_sha256 "$HERE/blas_bench.F90")"
prov_kv "$PF" binary "$RBIN"
prov_kv "$PF" binary_sha256 "$(prov_sha256 "$RBIN")"
if [ "$LAUNCH" = local ]; then BINFO="$BDIR/build_info.txt"; [ -z "$BLAS_BENCH_BIN" ] || BINFO=""
else BINFO="$BDIR/build_info_mpi.txt"; [ -z "$BLAS_BENCH_BIN_MPI" ] || BINFO=""; fi
if [ -n "$BINFO" ] && [ -r "$BINFO" ]; then
    while IFS= read -r line; do prov_kv "$PF" "build_${line%%=*}" "${line#*=}"; done < "$BINFO"
else
    prov_kv "$PF" build_info "unavailable(prebuilt binary BLAS_BENCH_BIN*, no build_info)"
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
if command -v otool >/dev/null 2>&1; then blibs="$(otool -L "$RBIN" 2>&1 | sed '1d' | awk '{print $1}')"
elif command -v ldd >/dev/null 2>&1; then blibs="$(ldd "$RBIN" 2>&1 | awk '{print $1 "->" $3}')"
else blibs="unavailable(no otool/ldd)"; fi
prov_kv "$PF" binary_libs "$blibs"
# other load on the machine (the numbers are only as good as the machine was quiet)
prov_kv "$PF" load_start "$(uptime | sed 's/.*load averages*: *//')"
prov_kv "$PF" busy_procs_start "$(ps -Ao pcpu=,comm= 2>/dev/null | sort -rn | head -6 | awk '{n = split($2, p, "/"); printf "%s%%%s ", $1, p[n]}')"

# ---- MKL probe, exactly as campaign/job_template.sh: through the launcher and rank_wrap.sh (fakeintel applied) -------------------
P="$OUTDIR/mkl_probe"
mkdir -p "$P"
case "$LAUNCH" in
    local)   ( cd "$P" && MKL_VERBOSE=1 PMI_RANK=0 "$CDIR/rank_wrap.sh" probe "$P" "$BIN" --probe ) > "$P/out.txt" 2>&1 < /dev/null ;;
    profile) ( cd "$P" && MKL_VERBOSE=1 profile_launch "$RANKS" "$THREADS" "$CDIR/rank_wrap.sh" probe "$P" "$BIN_MPI" --probe ) > "$P/out.txt" 2>&1 < /dev/null ;;
    *)       pre="${LAUNCH//@RANKS@/$RANKS}"; pre="${pre//@THREADS@/$THREADS}"
             ( cd "$P" && MKL_VERBOSE=1 $pre "$CDIR/rank_wrap.sh" probe "$P" "$BIN_MPI" --probe ) > "$P/out.txt" 2>&1 < /dev/null ;;
esac
prc=$?
mline="$(grep -E -m1 'MKL_VERBOSE.*(oneMKL|Intel\(R\) MKL|Math Kernel Library)' "$P/out.txt")"
if [ -n "$mline" ]; then
    prov_mkl_parse "$mline"
else
    MKL_VERSION="n/a (no MKL_VERBOSE line, probe rc=$prc: BLAS is not MKL or the probe failed; see binary_libs)"
    MKL_ISA="$MKL_VERSION"; MKL_THREADING="$MKL_VERSION"
fi
prov_kv "$PF" mkl_probe_rc "$prc"
prov_kv "$PF" mkl_probe_line "${mline:-none}"
prov_kv "$PF" mkl_version "$MKL_VERSION"
prov_kv "$PF" mkl_isa "$MKL_ISA"
prov_kv "$PF" mkl_threading "$MKL_THREADING"
prov_kv "$PF" fakeintel_maps_probe "$(cat "$P"/rank_probe_*.prov 2>/dev/null | sed -n 's/^fakeintel_maps=//p' | sort -u | tr '\n' ' ')"
REQ_MKL="${BLAS_BENCH_REQUIRE_MKL:-$([ "$PROFILE_NAME" = local ] && echo 0 || echo 1)}"
if [ "$prc" -ne 0 ] && [ "$REQ_MKL" = 1 ]; then echo "blas_bench.sh: the probe through the launcher failed (rc=$prc), see $P/out.txt" >&2; exit 2; fi
if [ -z "$mline" ] && [ "$REQ_MKL" = 1 ]; then echo "blas_bench.sh: the probe gave no MKL_VERBOSE version line, see $P/out.txt" >&2; exit 2; fi
if [ -n "$mline" ]; then
    prov_mkl_threads_check "$THREADS" "$MKL_THREADING" || { echo "blas_bench.sh: THREADS=$THREADS but MKL runs sequential" >&2; exit 2; }
    prov_mkl_check "$FAKEINTEL" "$MKL_ISA" || { echo "blas_bench.sh: fakeintel=on but MKL took the $MKL_ISA path" >&2; exit 2; }
fi

# ---- launch: R processes released at a common wall-clock time, then a barrier before every point ---------------------------------
if [ "$LAUNCH" = local ]; then DEFAULT_DELAY=3; else DEFAULT_DELAY=30; fi
DELAY="${BLAS_BENCH_SYNC_DELAY:-$DEFAULT_DELAY}"
prov_kv "$PF" sync_delay_s "$DELAY"
START=$(( $(date +%s) + DELAY ))
echo "$START" > "$OUTDIR/start_time.txt"
prov_kv "$PF" start_epoch_s "$START"
t_launch="$(prov_now)"
RC=0
EXTRA=()
[ "$LAUNCH" = local ] && EXTRA=("READYDIR=$OUTDIR/ready" "READY_TIMEOUT=${BLAS_BENCH_READY_TIMEOUT:-900}")
RANKCMD=("$SELF" --rank-exec "$OUTDIR" "$OUTDIR/start_time.txt" "$THREADS" "$RANKS" "$RBIN" "${FARGS[@]}" ${EXTRA[@]+"${EXTRA[@]}"})
case "$LAUNCH" in
    local)
        pids=""
        r=0
        while [ $r -lt "$RANKS" ]; do
            PMI_RANK=$r "${RANKCMD[@]}" < /dev/null &
            pids="$pids $!"
            r=$((r + 1))
        done
        for p in $pids; do wait "$p" || RC=1; done ;;
    profile)
        profile_launch "$RANKS" "$THREADS" "${RANKCMD[@]}" < /dev/null || RC=1 ;;
    *)
        pre="${LAUNCH//@RANKS@/$RANKS}"; pre="${pre//@THREADS@/$THREADS}"
        $pre "${RANKCMD[@]}" < /dev/null || RC=1 ;;
esac
t_end="$(prov_now)"
prov_kv "$PF" run_wall_s "$(awk -v a="$t_launch" -v b="$t_end" 'BEGIN {printf "%.3f", b - a}')"
prov_kv "$PF" launch_rc "$RC"
prov_kv "$PF" load_end "$(uptime | sed 's/.*load averages*: *//')"
prov_kv "$PF" busy_procs_end "$(ps -Ao pcpu=,comm= 2>/dev/null | sort -rn | head -6 | awk '{n = split($2, p, "/"); printf "%s%%%s ", $1, p[n]}')"
prov_kv "$PF" rank_hosts "$(cat "$OUTDIR"/rank_blas_*.prov 2>/dev/null | sed -n 's/^host=//p' | sort -u | tr '\n' ' ')"
prov_kv "$PF" rank_affinity_preinit "$(for p in "$OUTDIR"/rank_blas_*.prov; do [ -r "$p" ] && printf 'r%s:%s ' "$(sed -n 's/^rank=//p' "$p")" "$(sed -n 's/^affinity=//p' "$p")"; done)"

# ---- checks: every rank present exactly once with the full result count ----------------------------------------------------------
FAILS=""
fail() { FAILS="$FAILS
# FAIL $1"; echo "blas_bench.sh: FAIL: $1" >&2; }
[ "$RC" -eq 0 ] || fail "the launcher returned non-zero (a rank failed)"
r=0
while [ $r -lt "$RANKS" ]; do
    nf="$(ls "$OUTDIR"/out_${r}_p*.txt 2>/dev/null | wc -l | tr -d ' ')"
    if [ "$nf" -ne 1 ]; then
        fail "rank $r has $nf output files (expected 1; the launcher gave several processes the same rank, or none ran)"
    else
        f="$(ls "$OUTDIR"/out_${r}_p*.txt)"
        nr="$(grep -c '^RESULT ' "$f")"
        [ "$nr" -eq "$EXPECT" ] || fail "rank $r wrote $nr RESULT lines, expected $EXPECT"
        ns="$(grep -c '^SYNC ' "$f")"
        [ "$ns" -eq 1 ] || fail "rank $r wrote $ns SYNC lines, expected 1"
        if [ "$HAS_TIMERCOST" = 1 ]; then
            nt="$(grep -c '^TIMERCOST ' "$f")"
            [ "$nt" -eq 4 ] || fail "rank $r wrote $nt TIMERCOST lines, expected 4"
        fi
        awk -v r="$r" '$1 == "RESULT" && $2 != r {bad = 1} END {exit bad}' "$f" || fail "rank $r output contains lines of another rank"
    fi
    r=$((r + 1))
done
nprov="$(ls "$OUTDIR"/rank_blas_*.prov 2>/dev/null | wc -l | tr -d ' ')"
[ "$nprov" -eq "$RANKS" ] || fail "$nprov distinct rank records (rank_blas_*.prov), expected $RANKS"
MAXLAG="$(cat "$OUTDIR"/out_*_p*.txt 2>/dev/null | awk '$1 == "SYNC" {if ($5 + 0 > m) m = $5 + 0} END {printf "%.3f", m}')"
awk -v m="$MAXLAG" 'BEGIN {exit !(m > 0.5)}' && fail "a rank was released $MAXLAG s after the common start time (late launch): the ranks did not start together; raise BLAS_BENCH_SYNC_DELAY"
MAXSKEW="$(cat "$OUTDIR"/out_*_p*.txt 2>/dev/null | awk '$1 == "RESULT" {k = $4 " " $5 " " $6; t = $15 + 0; if (!(k in mn) || t < mn[k]) mn[k] = t; if (!(k in mx) || t > mx[k]) mx[k] = t}
    END {for (k in mn) if (mx[k] - mn[k] > m) m = mx[k] - mn[k]; printf "%.3f", m}')"
awk -v m="$MAXSKEW" 'BEGIN {exit !(m > 1.0)}' && fail "release skew of a point is $MAXSKEW s (> 1 s): the ranks did not run the points concurrently"

# ---- the log -----------------------------------------------------------------------------------------------------------------------
{
    echo "# blas_bench.sh log (testsuite/lowrank/blas_bench.sh): layout $LAYOUT, tag $TAG"
    echo "# PROV key=value: provenance (campaign/provenance.sh); RESULT/TIMERCOST/SYNC/META per rank; SUMMARY lines at the end"
    sed 's/^/PROV /' "$PF"
    echo "PROV_END"
    [ -n "$SKIPS" ] && echo "${SKIPS#
}"
    r=0
    while [ $r -lt "$RANKS" ]; do
        echo "# ---- rank $r"
        for f in "$OUTDIR"/out_${r}_p*.txt; do
            if [ -s "$f" ]; then cat "$f"; else echo "# rank $r: no output"; fi
            if [ -s "${f%.txt}.err" ]; then sed "s/^/# STDERR rank $r: /" "${f%.txt}.err"; fi
        done
        r=$((r + 1))
    done
    cat "$OUTDIR"/out_*_p*.txt 2>/dev/null | awk -v layout="$LAYOUT" '
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
            t = $15 + 0
            if (!(key in tmin) || t < tmin[key]) tmin[key] = t
            if (!(key in tmax) || t > tmax[key]) tmax[key] = t
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
                printf "SUMMARY %s %s %d %.6E %.3f %.3f %.3f %.3f %.3f\n", layout, key, m, tmed, gmed, gfs[key], cmed, gbs[key], tmax[key] - tmin[key]
            }
            for (q = 1; q <= ntk; q++) {
                key = torder[q]; m = tn[key]
                for (i = 1; i <= m; i++) a[i] = tv[key, i]
                printf "TIMERCOST_SUMMARY %s %s %d %.2f\n", layout, key, m, med(a, m)
            }
            printf "SYNC_SUMMARY %s ranks_synced %d max_lag_s %.3f\n", layout, nsync, maxlag
        }'
    [ -n "$FAILS" ] && echo "${FAILS#
}"
    nres="$(cat "$OUTDIR"/out_*_p*.txt 2>/dev/null | grep -c '^RESULT ')"
    echo "# end: launch rc=$RC, $nres RESULT lines from $RANKS ranks (expected $((EXPECT * RANKS))), max release skew $MAXSKEW s"
} > "$LOG"
echo "$LOG"
grep -E '^(SKIP|SUMMARY|TIMERCOST_SUMMARY|SYNC_SUMMARY|# FAIL|# end)' "$LOG"
[ -z "$FAILS" ] || exit 1
exit 0

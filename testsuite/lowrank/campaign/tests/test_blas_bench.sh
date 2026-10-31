#!/bin/bash
# Smoke test of the BLAS ceiling microbenchmark (testsuite/lowrank/blas_bench.{F90,sh}), bash 3.2, < 60 s on the Mac.
# Builds builds/blas_bench/ (noMPI and MPI binary) if needed, then checks: --probe ("probe ok", exit 0, < 1 s with sub-second
# timing); the sizing (memory guard with its SKIP line, L3 from a sysfs fixture and from an lscpu shim, cache override); a smoke
# run of layout 2x1 with the noMPI binary (ready-file barrier) and with the MPI binary (mpiexec -n 2, FI_PROVIDER=tcp on the
# Mac): the provenance block has no empty field and the fields the machine can provide are filled, every kernel x N x k x rank
# has a RESULT line with positive times (min <= median <= max) and positive GFLOP/s and GB/s, the TIMERCOST lines are positive
# (wtime < push_pop), the SUMMARY lines (median over ranks, node sum, skew) and the SYNC lines are there; ZGEMM N=400 (1 thread) is
# > 1 GFLOP/s; a layout whose launcher starts one process for two ranks fails (non-zero exit, FAIL lines). Exit 0 iff all pass.
# Usage: testsuite/lowrank/campaign/tests/test_blas_bench.sh [SCRATCH_TAG]      (logs in runs/blas_bench/<tag>*)
HERE="$(cd "$(dirname "$0")" && pwd)"
ALFDIR="$(cd "$HERE/../../../.." && pwd)"
SH="$ALFDIR/testsuite/lowrank/blas_bench.sh"
TAG="${1:-smoke_test}"
RUNS="$ALFDIR/runs/blas_bench"
fails=0
check() { if [ "$2" = "$3" ]; then echo "PASS $1 ($2)"; else echo "FAIL $1: got '$2', want '$3'"; fails=$((fails + 1)); fi; }
ok() { if [ "$2" = 1 ]; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }
filled() { case "$1" in ""|unavailable*) echo 0 ;; *) echo 1 ;; esac; }
now() { python3 -I -c 'import time; print("%.3f" % time.time())'; }
# MPICH of Homebrew hangs at the first collective over the default ofi provider on the Mac; tcp works
if [ "$(uname -s)" = Darwin ] && [ -z "$FI_PROVIDER" ]; then export FI_PROVIDER=tcp; fi

"$SH" --build-only > /dev/null || { echo "FAIL build"; exit 1; }
BIN="$ALFDIR/builds/blas_bench/blas_bench.out"
BINMPI="$ALFDIR/builds/blas_bench/blas_bench_mpi.out"

# --probe, timed with sub-second resolution
t0=$(now); out="$("$BIN" --probe)"; rc=$?; t1=$(now)
check "probe output" "$out" "probe ok"
check "probe rc" "$rc" 0
ok "probe < 1 s ($(awk -v a="$t0" -v b="$t1" 'BEGIN {printf "%.2f", b - a}') s)" "$(awk -v a="$t0" -v b="$t1" 'BEGIN {print (b - a < 1) ? 1 : 0}')"
check "probe MPI binary" "$("$BINMPI" --probe)" "probe ok"
"$BIN" KERNELS=nope > /dev/null 2>&1; check "unknown kernel rc" "$?" 2
for bad in 2x4abc x2 2x 2x2x2 0x1 1x0 ""; do
    "$SH" --plan LAYOUT="$bad" > /dev/null 2>&1; check "LAYOUT '$bad' rejected" "$?" 2
done

# sizing: memory guard (SKIP line), L3 from a sysfs fixture and from an lscpu shim, overrides
export BLAS_BENCH_MEM_AVAIL_MB=98304
plan="$("$SH" --plan LAYOUT=48x1 N=400,8000)"
echo "$plan" | sed 's/^/  plan: /'
ok "memory guard drops N=8000 with a SKIP line" "$(echo "$plan" | grep -q '^SKIP N=8000 reason=memory est_gb=' && echo 1 || echo 0)"
ok "memory guard keeps N=400" "$(echo "$plan" | grep -q 'n_run=400 *$' && echo 1 || echo 0)"
check "N=8000 fits one rank" "$("$SH" --plan LAYOUT=1x1 N=8000 | grep -c '^SKIP')" 0
FX="$RUNS/$TAG.sysfs"; /bin/rm -rf "$FX"
for c in 0 1 2 3; do
    d="$FX/sys/devices/system/cpu/cpu$c/cache"
    mkdir -p "$d/index2" "$d/index3"
    echo 2 > "$d/index2/level"; echo 3 > "$d/index3/level"
    echo 65536K > "$d/index3/size"; echo 1048576K > "$d/index2/size"
    if [ $c -lt 2 ]; then echo 0-1 > "$d/index3/shared_cpu_list"; else echo 2-3 > "$d/index3/shared_cpu_list"; fi
    echo $c > "$d/index2/shared_cpu_list"
done
p="$(LOCAL_SYSFS_ROOT="$FX" "$SH" --plan LAYOUT=1x1)"
check "L3 total from sysfs (2 instances x 64 MB)" "$(echo "$p" | sed -n 's/.*l3_total_mb=\([0-9]*\).*/\1/p')" 128
check "cache = 4 x L3 / R (1 rank)" "$(echo "$p" | sed -n 's/.*cache_mb_used=\([0-9]*\).*/\1/p')" 512
p="$(LOCAL_SYSFS_ROOT="$FX" "$SH" --plan LAYOUT=4x1)"
check "cache floor 256 wins (4 ranks: 128)" "$(echo "$p" | sed -n 's/.*cache_mb_used=\([0-9]*\).*/\1/p')" 256
p="$(LOCAL_SYSFS_ROOT="$FX" BLAS_BENCH_CACHE_MB=64 "$SH" --plan LAYOUT=1x1)"
check "BLAS_BENCH_CACHE_MB overrides" "$(echo "$p" | sed -n 's/.*cache_mb_used=\([0-9]*\).*/\1/p')" 64
SH2="$RUNS/$TAG.shim"; /bin/rm -rf "$SH2"; mkdir -p "$SH2" "$RUNS/$TAG.emptysys"
printf '#!/bin/sh\necho "L3 cache:                       384 MiB (2 instances)"\necho "Socket(s):                       2"\n' > "$SH2/lscpu"; chmod +x "$SH2/lscpu"
p="$(PATH="$SH2:$PATH" LOCAL_SYSFS_ROOT="$RUNS/$TAG.emptysys" "$SH" --plan LAYOUT=2x1)"
check "L3 total from lscpu shim" "$(echo "$p" | sed -n 's/.*l3_total_mb=\([0-9]*\).*/\1/p')" 384
check "cache = 4 x 384 / 2" "$(echo "$p" | sed -n 's/.*cache_mb_used=\([0-9]*\).*/\1/p')" 768
unset BLAS_BENCH_MEM_AVAIL_MB

# one smoke run, checked in full: $1 = tag suffix, rest = environment words for env(1)
check_smoke_log() {   # check_smoke_log NAME LOG
    local name="$1" LOG="$2" P bad k v kern rk w
    ok "$name: log exists" "$([ -s "$LOG" ] && echo 1 || echo 0)"
    awk '/^PROV_END/ {exit} /^PROV / {print substr($0, 6)}' "$LOG" > "$LOG.prov"
    P="$LOG.prov"
    check "$name: provenance lines without '='" "$(grep -vc '=' "$P")" 0
    check "$name: provenance empty values" "$(grep -c '=$' "$P")" 0
    check "$name: provenance duplicate keys" "$(sed 's/=.*//' "$P" | sort | uniq -d | wc -l | tr -d ' ')" 0
    for k in hostname profile job_start layout ranks threads blas_bench_args binary_sha256 build_compiler_version build_flags \
             build_blas_link git_commit binary_libs omp_num_threads mkl_num_threads load_start mkl_version mkl_threading \
             pinning cache_mb_used mem_avail_mb launcher mpi_binary; do
        v="$(sed -n "s/^$k=//p" "$P")"
        ok "$name: provenance $k filled" "$(filled "$v")"
    done
    if [ "$(uname -s)" = Darwin ]; then
        for k in cpu_model_sysctl cpus_sysctl mem_bytes_sysctl; do
            ok "$name: provenance $k filled (macOS)" "$(filled "$(sed -n "s/^$k=//p" "$P")")"
        done
    fi
    check "$name: provenance ranks" "$(sed -n 's/^ranks=//p' "$P")" 2
    check "$name: omp = mkl threads" "$(sed -n 's/^omp_num_threads=//p' "$P")" "$(sed -n 's/^mkl_num_threads=//p' "$P")"
    check "$name: RESULT lines" "$(grep -c '^RESULT ' "$LOG" | tr -d ' ')" 12
    bad="$(grep '^RESULT ' "$LOG" | awk '
        NF != 16 {print "fields: " $0; next}
        !($9 > 0 && $10 > 0 && $11 > 0 && $12 > 0 && $13 > 0 && $14 > 0 && $15 > 1e9 && $16 >= 1) {print "nonpositive: " $0; next}
        ($10 > $9 || $9 > $11) {print "order min<=med<=max: " $0; next}
        $7 < 2 {print "reps: " $0}')"
    check "$name: RESULT schema and positivity" "$bad" ""
    for kern in zgemm_k zgemm_n zgeru zgemv zgeqp3; do
        for rk in 0 1; do
            ok "$name: RESULT rank $rk $kern present" "$([ "$(grep -c "^RESULT $rk 1 $kern 64 " "$LOG")" -ge 1 ] && echo 1 || echo 0)"
        done
    done
    for w in push_pop switch system_clock wtime; do
        check "$name: TIMERCOST $w lines" "$(grep -c "^TIMERCOST [01] 1 $w " "$LOG" | tr -d ' ')" 2
    done
    check "$name: TIMERCOST positive" "$(grep '^TIMERCOST ' "$LOG" | awk '!($7 > 0 && $8 > 0 && $9 > 0) {print}')" ""
    check "$name: wtime (clock alone) < push_pop pair" "$(awk '$1 == "TIMERCOST_SUMMARY" {v[$3] = $5} END {print (v["wtime"] < v["push_pop"]) ? "yes" : "no"}' "$LOG")" yes
    check "$name: SUMMARY lines" "$(grep -c '^SUMMARY ' "$LOG" | tr -d ' ')" 6
    check "$name: SUMMARY ranks" "$(grep '^SUMMARY ' "$LOG" | awk '$6 != 2 {print}')" ""
    check "$name: SUMMARY positive and skew < 1 s" "$(grep '^SUMMARY ' "$LOG" | awk '!($7 > 0 && $8 > 0 && $9 > 0 && $10 > 0 && $11 > 0 && $12 >= 0 && $12 < 1) {print}')" ""
    check "$name: TIMERCOST_SUMMARY lines" "$(grep -c '^TIMERCOST_SUMMARY ' "$LOG" | tr -d ' ')" 4
    check "$name: SYNC lines" "$(grep -c '^SYNC ' "$LOG" | tr -d ' ')" 2
    check "$name: SYNC lag in [0, 0.5)" "$(grep '^SYNC ' "$LOG" | awk '!($5 >= 0 && $5 < 0.5) {print}')" ""
    check "$name: no FAIL lines" "$(grep -c '^# FAIL' "$LOG" | tr -d ' ')" 0
}

# noMPI binary, local launcher, two concurrent ranks (ready-file barrier before every point)
/bin/rm -rf "$RUNS/$TAG.d" "$RUNS/blas_bench_$TAG.txt"
"$SH" --smoke LAYOUT=2x1 TAG="$TAG" > /dev/null; check "noMPI smoke rc" "$?" 0
check_smoke_log noMPI "$RUNS/blas_bench_$TAG.txt"
check "noMPI: ready files of every point (2 ranks x 7 points)" "$(ls "$RUNS/$TAG.d/ready" | grep -c '^ready_')" 14

# MPI binary under mpiexec -n 2: MPI_Init, MPI_Barrier before every point
if command -v mpiexec >/dev/null 2>&1; then
    /bin/rm -rf "$RUNS/${TAG}_mpi.d" "$RUNS/blas_bench_${TAG}_mpi.txt"
    BLAS_BENCH_LAUNCH='mpiexec -n @RANKS@' BLAS_BENCH_SYNC_DELAY=5 "$SH" --smoke LAYOUT=2x1 TAG="${TAG}_mpi" > /dev/null
    check "MPI smoke rc" "$?" 0
    check_smoke_log MPI "$RUNS/blas_bench_${TAG}_mpi.txt"
    check "MPI: META says mpi=yes" "$(grep -c '^META mpi=yes' "$RUNS/blas_bench_${TAG}_mpi.txt" | tr -d ' ')" 2
else
    echo "FAIL mpiexec not found"; fails=$((fails + 1))
fi

# a launcher that starts ONE process for two ranks must fail (missing rank), with FAIL lines in the log
/bin/rm -rf "$RUNS/${TAG}_one.d" "$RUNS/blas_bench_${TAG}_one.txt"
BLAS_BENCH_LAUNCH='env' BLAS_BENCH_SYNC_DELAY=1 "$SH" --smoke LAYOUT=2x1 TAG="${TAG}_one" > /dev/null 2>&1
ok "one process for two ranks -> non-zero exit" "$([ $? -ne 0 ] && echo 1 || echo 0)"
ok "... and FAIL lines in the log" "$([ "$(grep -c '^# FAIL' "$RUNS/blas_bench_${TAG}_one.txt" 2>/dev/null)" -ge 1 ] && echo 1 || echo 0)"

# ZGEMM N=400, 1 thread: plausible rate
/bin/rm -rf "$RUNS/${TAG}_400.d" "$RUNS/blas_bench_${TAG}_400.txt"
"$SH" LAYOUT=1x1 TAG="${TAG}_400" N=400 KERNELS=zgemm_n REPS=3 WARMUP=1 TIMER_REPS=0 > /dev/null
check "ZGEMM N=400 run rc" "$?" 0
G="$(grep '^RESULT 0 1 zgemm_n 400 400 ' "$RUNS/blas_bench_${TAG}_400.txt" | awk '{print $12}')"
ok "ZGEMM N=400 > 1 GFLOP/s (got $G)" "$(awk -v g="$G" 'BEGIN {print (g > 1) ? 1 : 0}')"
echo "failures: $fails"
[ "$fails" -eq 0 ]

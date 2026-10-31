#!/bin/bash
# Smoke test of the BLAS ceiling microbenchmark (testsuite/lowrank/blas_bench.{F90,sh}), bash 3.2, < 30 s on the Mac.
# Builds builds/blas_bench/ if needed, then checks: --probe ("probe ok", exit 0, < 1 s); a smoke run of layout 2x1 (two
# concurrent ranks): the provenance block has no empty field and the fields the Mac can provide are filled, every
# kernel x N x k x rank has a RESULT line with positive times (min <= median <= max) and positive GFLOP/s and GB/s, the
# four TIMERCOST lines are positive, the SUMMARY lines (median over ranks, node sum) and the SYNC lines (ranks released
# together) are there; ZGEMM N=400 (1 thread) is > 1 GFLOP/s. Exit 0 iff all pass.
# Usage: testsuite/lowrank/campaign/tests/test_blas_bench.sh [SCRATCH_TAG]      (logs in runs/blas_bench/<tag>*)
HERE="$(cd "$(dirname "$0")" && pwd)"
ALFDIR="$(cd "$HERE/../../../.." && pwd)"
SH="$ALFDIR/testsuite/lowrank/blas_bench.sh"
TAG="${1:-smoke_test}"
fails=0
check() { if [ "$2" = "$3" ]; then echo "PASS $1 ($2)"; else echo "FAIL $1: got '$2', want '$3'"; fails=$((fails + 1)); fi; }
filled() { case "$1" in ""|unavailable*) echo 0 ;; *) echo 1 ;; esac; }
ok() { if [ "$2" = 1 ]; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }

"$SH" --build-only > /dev/null || { echo "FAIL build"; exit 1; }
BIN="$ALFDIR/builds/blas_bench/blas_bench.out"

# --probe
t0=$(date +%s)
out="$("$BIN" --probe)"; rc=$?
t1=$(date +%s)
check "probe output" "$out" "probe ok"
check "probe rc" "$rc" 0
ok "probe < 1 s" "$([ $((t1 - t0)) -le 1 ] && echo 1 || echo 0)"

# smoke run, two concurrent ranks
/bin/rm -rf "$ALFDIR/runs/blas_bench/$TAG.d" "$ALFDIR/runs/blas_bench/blas_bench_$TAG.txt"
"$SH" --smoke LAYOUT=2x1 TAG="$TAG" > /dev/null; rc=$?
check "blas_bench.sh smoke rc" "$rc" 0
LOG="$ALFDIR/runs/blas_bench/blas_bench_$TAG.txt"
ok "log exists" "$([ -s "$LOG" ] && echo 1 || echo 0)"

# provenance: key=value, no empty value, no duplicate key, mandatory fields present and filled
awk '/^PROV_END/ {exit} /^PROV / {print substr($0, 6)}' "$LOG" > "$ALFDIR/runs/blas_bench/$TAG.prov"
P="$ALFDIR/runs/blas_bench/$TAG.prov"
check "provenance lines without '='" "$(grep -vc '=' "$P")" 0
check "provenance empty values" "$(grep -c '=$' "$P")" 0
check "provenance duplicate keys" "$(sed 's/=.*//' "$P" | sort | uniq -d | wc -l | tr -d ' ')" 0
for k in hostname profile job_start layout ranks threads blas_bench_args binary_sha256 build_compiler_version build_flags \
         build_blas_link git_commit binary_libs omp_num_threads mkl_num_threads load_start mkl_version mkl_threading; do
    v="$(sed -n "s/^$k=//p" "$P")"
    ok "provenance $k filled" "$(filled "$v")"
done
if [ "$(uname -s)" = Darwin ]; then
    for k in cpu_model_sysctl cpus_sysctl mem_bytes_sysctl; do
        ok "provenance $k filled (macOS)" "$([ -n "$(sed -n "s/^$k=//p" "$P")" ] && echo 1 || echo 0)"
    done
fi
check "provenance ranks" "$(sed -n 's/^ranks=//p' "$P")" 2
check "provenance omp = mkl threads" "$(sed -n 's/^omp_num_threads=//p' "$P")" "$(sed -n 's/^mkl_num_threads=//p' "$P")"

# RESULT lines: 2 ranks x (zgemm_k 16, zgemm_k 64, zgemm_n, zgeru, zgemv, zgeqp3) at N=64
check "RESULT lines" "$(grep -c '^RESULT ' "$LOG" | tr -d ' ')" 12
bad="$(grep '^RESULT ' "$LOG" | awk '
    NF != 14 {print "fields: " $0; next}
    !($9 > 0 && $10 > 0 && $11 > 0 && $12 > 0 && $13 > 0 && $14 > 0) {print "nonpositive: " $0; next}
    ($10 > $9 || $9 > $11) {print "order min<=med<=max: " $0; next}
    $7 < 2 {print "reps: " $0}')"
check "RESULT schema and positivity" "$bad" ""
for kern in zgemm_k zgemm_n zgeru zgemv zgeqp3; do
    for rk in 0 1; do
        n="$(grep -c "^RESULT $rk 1 $kern 64 " "$LOG")"
        ok "RESULT rank $rk $kern present" "$([ "$n" -ge 1 ] && echo 1 || echo 0)"
    done
done
# TIMERCOST: push_pop, switch, system_clock, wtime for each rank, positive ns
for w in push_pop switch system_clock wtime; do
    check "TIMERCOST $w lines" "$(grep -c "^TIMERCOST [01] 1 $w " "$LOG" | tr -d ' ')" 2
done
check "TIMERCOST positive" "$(grep '^TIMERCOST ' "$LOG" | awk '!($7 > 0 && $8 > 0 && $9 > 0) {print}')" ""
# SUMMARY: one per kernel/N/k, 2 ranks, node sum = about 2 x median
check "SUMMARY lines" "$(grep -c '^SUMMARY ' "$LOG" | tr -d ' ')" 6
check "SUMMARY ranks" "$(grep '^SUMMARY ' "$LOG" | awk '$6 != 2 {print}')" ""
check "SUMMARY positive" "$(grep '^SUMMARY ' "$LOG" | awk '!($7 > 0 && $8 > 0 && $9 > 0 && $10 > 0 && $11 > 0) {print}')" ""
check "TIMERCOST_SUMMARY lines" "$(grep -c '^TIMERCOST_SUMMARY ' "$LOG" | tr -d ' ')" 4
# SYNC: both ranks released together (lag >= 0, < 0.5 s)
check "SYNC lines" "$(grep -c '^SYNC ' "$LOG" | tr -d ' ')" 2
check "SYNC lag in [0, 0.5)" "$(grep '^SYNC ' "$LOG" | awk '!($5 >= 0 && $5 < 0.5) {print}')" ""
check "no WARN" "$(grep -c '^# WARN' "$LOG" | tr -d ' ')" 0

# ZGEMM N=400, 1 thread: plausible rate
/bin/rm -rf "$ALFDIR/runs/blas_bench/${TAG}_400.d" "$ALFDIR/runs/blas_bench/blas_bench_${TAG}_400.txt"
"$SH" LAYOUT=1x1 TAG="${TAG}_400" N=400 KERNELS=zgemm_n REPS=3 WARMUP=1 TIMER_REPS=0 > /dev/null
check "ZGEMM N=400 run rc" "$?" 0
G="$(grep '^RESULT 0 1 zgemm_n 400 400 ' "$ALFDIR/runs/blas_bench/blas_bench_${TAG}_400.txt" | awk '{print $12}')"
ok "ZGEMM N=400 > 1 GFLOP/s (got $G)" "$(awk -v g="$G" 'BEGIN {print (g > 1) ? 1 : 0}')"
# a bad kernel name / size must fail loudly
"$BIN" KERNELS=nope > /dev/null 2>&1; check "unknown kernel rc" "$?" 2
echo "failures: $fails"
[ "$fails" -eq 0 ]

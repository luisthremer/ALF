#!/bin/bash
# End-to-end local test of the campaign generator with the 'local' profile and the PATH shims (no cluster).
# Usage: testsuite/lowrank/campaign/tests/run_local_test.sh MPI_BINARY NOMPI_BINARY
#   frozen binaries with <binary>.build_info.txt next to them (build_info.sh), e.g.
#   runs/campaign_test/bin/ALF_mpi.out (2-rank layouts need MPI) and runs/campaign_test/bin/ALF_nompi.out.
# Output under runs/campaign_test/ (git-excluded); every test uses fresh campaign ids (time stamp), nothing is deleted.
# Checks: placeholder refusal (supermuc, helma); layout and memory refusal; generator idempotency; mini campaign
# (Hubbard L=8, standard + delayed N_delay 8, layouts 1x1 and 2x2, 2 repetitions) end to end with collect --strict
# (every provenance field filled from shims/fixtures, none unknown/unavailable); an injected FAILED run is excluded
# and listed; a second submit_all.sh skips finished runs and re-runs only the failed one; MAX_ATTEMPTS gives up;
# the queue limit stops submission; the SIGTERM trap marks the interrupted run FAILED rc=143 and the rest NOT_RUN;
# fakeintel=on aborts the job when MKL reports AVX2 and passes with AVX-512; collect refuses a missing or duplicate
# provenance key; report refuses to mix binaries/commits. Core budget: one job at a time, at most 2 MPI ranks.
HERE="$(cd "$(dirname "$0")" && pwd)"
CDIR="$(dirname "$HERE")"
ALF="$(cd "$CDIR/../../.." && pwd)"
MPI_BIN="$(cd "$(dirname "${1:?MPI binary}")" && pwd)/$(basename "$1")"
NOMPI_BIN="$(cd "$(dirname "${2:?noMPI binary}")" && pwd)/$(basename "$2")"
STAMP="$(date +%Y%m%d%H%M%S)"
T="$ALF/runs/campaign_test"
ROOT="$T/campaigns"
export PATH="$HERE/shims:$PATH" SHIM_STATE="$T/shim_state_$STAMP" LOCAL_SYSFS_ROOT="$HERE/fixtures/sysfs"
export SHIM_MKL_FIXTURE=supermuc
mkdir -p "$ROOT" "$SHIM_STATE"
MK="python3 -I $CDIR/make_campaign.py"
COLLECT="python3 -I $CDIR/collect.py"
REPORT="python3 -I $CDIR/report.py"
fails=0
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; fails=$((fails + 1)); }
expect() { if [ "$2" = "$3" ]; then pass "$1 ($2)"; else fail "$1: got '$2', want '$3'"; fi; }
st() { awk '{print $1; exit}' "$1/runs/$2/status" 2>/dev/null; }
natt() { ls -d "$1/runs/$2"/a[0-9]* 2>/dev/null | wc -l | tr -d ' '; }

echo "== unit checks"
python3 -I "$HERE/test_collect.py" | tail -1 | grep -q 'ALL PASS' && pass "test_collect.py" || fail "test_collect.py"
"$HERE/test_provenance.sh" "$T/unit_$STAMP" | tail -1 | grep -q 'ALL PASS' && pass "test_provenance.sh" || fail "test_provenance.sh"

echo "== generator refusals"
for p in supermuc helma; do
    out="$(CAMPAIGN_TEST_PROFILE=$p $MK "$HERE/campaign_placeholder.py" --root "$ROOT" 2>&1)"; rc=$?
    expect "placeholder refusal $p rc" "$rc" 1
    printf '%s\n' "$out" | grep -q "still contains placeholders" && pass "placeholder message $p ($(printf '%s\n' "$out" | grep -c 'line '))" \
        || fail "placeholder message $p: $out"
done
echo "$out" | head -4
export CAMPAIGN_TEST_BIN="$MPI_BIN"
TMPDESC="$T/refuse_$STAMP.py"
sed 's/"partial_node": True/"partial_node": False/' "$HERE/campaign_test.py" > "$TMPDESC"
out="$(CAMPAIGN_TEST_ID=t3refuse_$STAMP $MK "$TMPDESC" --dry-run 2>&1)"; rc=$?
expect "layout refusal (1x1 on a 4-core node without partial_node) rc" "$rc" 1; echo "  $out"
sed 's/"partial_node": True/"partial_node": True, "mem_gb_per_rank": 9/' "$HERE/campaign_test.py" > "$TMPDESC"
out="$(CAMPAIGN_TEST_ID=t3refuse_$STAMP $MK "$TMPDESC" --dry-run 2>&1)"; rc=$?
expect "memory refusal (2 ranks x 9 GB > 90% of 16 GB) rc" "$rc" 1; echo "  $out"

echo "== main campaign (MPI binary): generate, idempotency"
export CAMPAIGN_TEST_ID="t3main_$STAMP"
C="$ROOT/$CAMPAIGN_TEST_ID"
$MK "$HERE/campaign_test.py" --root "$ROOT" || fail "generate main"
out="$($MK "$HERE/campaign_test.py" --root "$ROOT" 2>&1)"; rc=$?
expect "regenerate same description rc" "$rc" 0
printf '%s\n' "$out" | grep -q unchanged && pass "regenerate: unchanged" || fail "regenerate: $out"
sed 's/"nbin": 3/"nbin": 4/' "$HERE/campaign_test.py" > "$TMPDESC"
out="$($MK "$TMPDESC" --root "$ROOT" 2>&1)"; rc=$?
expect "regenerate changed description under the same id rc" "$rc" 1
expect "jobs" "$(wc -l < "$C/jobs.list" | tr -d ' ')" 4

echo "== queue limit"
for i in 1 2 3 4; do echo "c_${CAMPAIGN_TEST_ID}_fake$i"; done > "$SHIM_STATE/queue"
out="$("$C/submit_all.sh" 2>&1)"
printf '%s\n' "$out" | grep -q "queue limit reached" && pass "queue limit stops submission" || fail "queue limit: $out"
expect "no job submitted while the queue is full" "$(ls "$C/jobs"/*.submissions 2>/dev/null | wc -l | tr -d ' ')" 0
mv "$SHIM_STATE/queue" "$SHIM_STATE/queue.used"

echo "== main campaign: first submission with an injected failure (OOM-like rc 137 in hub_L8_2x2_r2_delayed_nd8)"
FAILRUN=hub_L8_2x2_r2_delayed_nd8
SHIM_SRUN_FAIL_MATCH="/$FAILRUN/" "$C/submit_all.sh"
for r in $(cut -d' ' -f2- "$C/jobs.list"); do
    if [ "$r" = "$FAILRUN" ]; then expect "status $r" "$(st "$C" "$r")" FAILED; else expect "status $r" "$(st "$C" "$r")" OK; fi
done
cat "$C/runs/$FAILRUN/status"
$COLLECT "$C" --strict > "$T/collect1_$STAMP.txt" 2>&1; rc=$?
cat "$T/collect1_$STAMP.txt"
expect "collect --strict rc (every field filled)" "$rc" 0
grep -q "EXCLUDED $FAILRUN: FAILED rc=137" "$T/collect1_$STAMP.txt" && pass "failed run excluded and listed" || fail "failed run not listed"
expect "collected OK rows" "$(python3 -I -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["rows"]))' "$C/collected.json")" 7
n_unknown="$(python3 -I -c '
import json, sys
rows = json.load(open(sys.argv[1]))["rows"]
bad = [(r["run_id"], k, v) for r in rows for k, v in r.items() if isinstance(v, str) and ("unknown" in v.lower() or v.startswith("unavailable")) and k != "job.build_git_commit" and k != "job.build_git_dirty"]
for b in bad: print(b, file=sys.stderr)
print(len(bad))' "$C/collected.json")"
expect "row values unknown/unavailable (except build git of the archived MPI tree)" "$n_unknown" 0
$REPORT "$C/collected.json" --out "$C/report.md" && pass "report" || fail "report"

echo "== second submission: finished runs skipped, failed run re-run (still failing), then MAX_ATTEMPTS gives up"
before="$(for r in $(cut -d' ' -f2- "$C/jobs.list"); do printf '%s:%s ' "$r" "$(natt "$C" "$r")"; done)"
SHIM_SRUN_FAIL_MATCH="/$FAILRUN/" "$C/submit_all.sh" | tee "$T/submit2_$STAMP.txt"
expect "second submit: jobs submitted" "$(grep -c '^submit ' "$T/submit2_$STAMP.txt")" 1
expect "second submit: jobs done" "$(grep -c '^done ' "$T/submit2_$STAMP.txt")" 3
for r in $(cut -d' ' -f2- "$C/jobs.list"); do
    want=1; [ "$r" = "$FAILRUN" ] && want=2
    expect "attempts $r" "$(natt "$C" "$r")" "$want"
done
grep -q "skip hub_L8_2x2_r2_standard (status OK)" "$C"/logs/c_${CAMPAIGN_TEST_ID}_hub_L8_2x2_r2.*.out \
    && pass "job skips the OK run of the same job" || fail "OK run not skipped in the job log"
"$C/submit_all.sh" | tee "$T/submit3_$STAMP.txt"
grep -q "^GAVE UP" "$T/submit3_$STAMP.txt" && pass "MAX_ATTEMPTS: gave up" || fail "no GAVE UP"
expect "third submit: nothing submitted" "$(grep -c '^submit ' "$T/submit3_$STAMP.txt")" 0
$COLLECT "$C" --strict | tee "$T/collect2_$STAMP.txt"
grep -q "EXCLUDED $FAILRUN: FAILED rc=137 .*attempts 2" "$T/collect2_$STAMP.txt" && pass "failed run listed with 2 attempts" || fail "collect2"
$REPORT "$C/collected.json" --out "$C/report.md" && pass "report after resubmission" || fail "report 2"

echo "== collect refuses a missing / duplicate provenance key (on a copy)"
CP="$T/copy_$STAMP"; cp -R "$C" "$CP"
p="$CP/runs/hub_L8_1x1_r1_standard/a1/run.prov"
grep -v '^max_rss_kb=' "$p" > "$p.new" && mv "$p.new" "$p"
$COLLECT "$CP" > "$T/collect_missing_$STAMP.txt" 2>&1; rc=$?
expect "missing key rc" "$rc" 1; cat "$T/collect_missing_$STAMP.txt"
echo "max_rss_kb=1" >> "$p"; echo "max_rss_kb=2" >> "$p"
$COLLECT "$CP" > "$T/collect_dup_$STAMP.txt" 2>&1; rc=$?
expect "duplicate key rc" "$rc" 1; cat "$T/collect_dup_$STAMP.txt"

echo "== SIGTERM trap (noMPI binary, 1x1, 2 repetitions; time limit simulated in hub_L8_1x1_r1_standard)"
export CAMPAIGN_TEST_ID="t3nompi_$STAMP" CAMPAIGN_TEST_BIN="$NOMPI_BIN" CAMPAIGN_TEST_LAYOUTS="1x1"
C2="$ROOT/$CAMPAIGN_TEST_ID"
$MK "$HERE/campaign_test.py" --root "$ROOT" || fail "generate nompi"
SHIM_SRUN_TERM_MATCH="/hub_L8_1x1_r1_standard/" "$C2/submit_all.sh"
cat "$C2/runs/hub_L8_1x1_r1_standard/status" "$C2/runs/hub_L8_1x1_r1_delayed_nd8/status"
grep -q '^FAILED rc=143 .*SIGTERM' "$C2/runs/hub_L8_1x1_r1_standard/status" && pass "trap: interrupted run FAILED rc=143" || fail "trap status"
expect "trap: later run of the job" "$(st "$C2" hub_L8_1x1_r1_delayed_nd8)" NOT_RUN
expect "other job unaffected" "$(st "$C2" hub_L8_1x1_r2_standard)$(st "$C2" hub_L8_1x1_r2_delayed_nd8)" OKOK
$COLLECT "$C2" --strict | tee "$T/collect_nompi_$STAMP.txt"
"$C2/submit_all.sh"
$COLLECT "$C2" --strict && pass "nompi collect after resubmission" || fail "nompi collect"
$REPORT "$C2/collected.json" --out "$C2/report.md" && pass "report nompi" || fail "report nompi"
$REPORT "$C/collected.json" "$C2/collected.json" > "$T/report_mixed_$STAMP.txt" 2>&1; rc=$?
expect "report refuses to mix binaries/commits rc" "$rc" 2; cat "$T/report_mixed_$STAMP.txt"

echo "== fakeintel=on (noMPI binary, standard only)"
lib="$T/fake_$STAMP/libfakeintel.so"; mkdir -p "$(dirname "$lib")"; : > "$lib"
export CAMPAIGN_TEST_FAKEINTEL=on CAMPAIGN_TEST_SCHEMES=standard CAMPAIGN_TEST_REPS=1 LOCAL_FAKEINTEL_LIB="$lib"
export LOCAL_PROC_ROOT="$HERE/fixtures"
export CAMPAIGN_TEST_ID="t3fi_$STAMP"
C3="$ROOT/$CAMPAIGN_TEST_ID"
$MK "$HERE/campaign_test.py" --root "$ROOT" || fail "generate fakeintel"
SHIM_MKL_FIXTURE=helma_nofakeintel "$C3/submit_all.sh"
grep -h ABORT "$C3"/logs/*.out
grep -q "ABORT: fakeintel=on but MKL took the AVX2 path" "$C3"/logs/*.out && pass "fakeintel on + AVX2: job aborts" || fail "fakeintel abort"
expect "fakeintel abort: run not started" "$(st "$C3" hub_L8_1x1_r1_standard)" NOT_RUN
SHIM_MKL_FIXTURE=helma_fakeintel "$C3/submit_all.sh"
expect "fakeintel on + AVX-512: run" "$(st "$C3" hub_L8_1x1_r1_standard)" OK
$COLLECT "$C3" --strict && pass "fakeintel collect --strict" || fail "fakeintel collect"
grep -E '^(fakeintel|mkl_isa|fakeintel_lib_sha256|fakeintel_maps_probe)=' "$C3"/jobs/*/job_*/job.prov | tail -4
grep -E '^(fakeintel_maps|ld_preload)=' "$C3"/runs/hub_L8_1x1_r1_standard/a1/run.prov

echo "run_local_test: $([ $fails = 0 ] && echo ALL PASS || echo FAILED) ($fails failures); campaigns: $C $C2 $C3"
[ $fails = 0 ]

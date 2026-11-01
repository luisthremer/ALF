#!/bin/bash
# Local test of packed multi-node jobs (description key pack_nodes) with the 'local' profile and the PATH shims (no cluster).
# Usage: testsuite/lowrank/campaign/tests/run_pack_test.sh MPI_BINARY NOMPI_BINARY
#   (the frozen binaries of run_local_test.sh; the packed campaign runs the noMPI binary, layout 1x1: at most 3 processes)
# Output under runs/campaign_test/ (git-excluded); fresh campaign ids (time stamp), nothing is deleted.
# Checks: pack_nodes = 1 given explicitly gives the same job scripts / runs as omitting it (one node, --ntasks); refusals
# (pack_nodes < 1, a single run longer than the partition limit); pack_nodes above the partition's node limit is cut with a
# warning; --dry-run prints the packing; the first-fit-decreasing slot assignment (4 slots, jobs of 3 + 1 nodes, run order
# inside a slot rotated by the repetition) against an independent re-implementation; the generated job script (nodes,
# --ntasks-per-node, --exclusive, host list, background slots, wait, per-slot status); end to end: slots run concurrently,
# runs of a slot one after the other, slot i on host i, a failing run does not stop the others, a slot that cannot start
# (MKL probe fails on its host) marks only its own runs NOT_RUN, per-slot status lines, collect.py --strict and report.py
# read the packed campaign, a resubmission re-runs only the pending runs; SIGTERM (time limit) marks the running runs of
# every slot FAILED rc=143, the rest NOT_RUN, and leaves no process behind; the host check aborts all slots when the
# launcher does not start them on their host.
HERE="$(cd "$(dirname "$0")" && pwd)"
CDIR="$(dirname "$HERE")"
ALF="$(cd "$CDIR/../../.." && pwd)"
MPI_BIN="$(cd "$(dirname "${1:?MPI binary}")" && pwd)/$(basename "$1")"
NOMPI_BIN="$(cd "$(dirname "${2:?noMPI binary}")" && pwd)/$(basename "$2")"
STAMP="$(date +%Y%m%d%H%M%S)"
T="$ALF/runs/campaign_test"
ROOT="$T/campaigns"
export PATH="$HERE/shims:$PATH" SHIM_STATE="$T/shim_state_pack_$STAMP" LOCAL_SYSFS_ROOT="$HERE/fixtures/sysfs"
export SHIM_MKL_FIXTURE=supermuc
mkdir -p "$ROOT" "$SHIM_STATE"
MK="python3 -I $CDIR/make_campaign.py"
COLLECT="python3 -I $CDIR/collect.py"
REPORT="python3 -I $CDIR/report.py"
DESC="$HERE/campaign_pack_test.py"
fails=0
pass() { echo "PASS $*"; }
fail() { echo "FAIL $*"; fails=$((fails + 1)); }
expect() { if [ "$2" = "$3" ]; then pass "$1 ($2)"; else fail "$1: got '$2', want '$3'"; fi; }
st() { awk '{print $1; exit}' "$1/runs/$2/status" 2>/dev/null; }
has() { if printf '%s\n' "$2" | grep -q -F -- "$3"; then pass "$1"; else fail "$1: no '$3' in: $2"; fi; }
hasf() { if grep -q -F -- "$3" "$2"; then pass "$1"; else fail "$1: no '$3' in $2"; fi; }

echo "== static checks (bash 3.2 syntax, no bash-4 constructs)"
for f in job_template_packed.sh node_prov.sh provenance.sh profiles/local.sh profiles/supermuc.sh profiles/helma.sh; do
    /bin/bash -n "$CDIR/$f" && pass "bash -n $f" || fail "bash -n $f"
    if grep -n -E 'declare -A|mapfile|readarray|\$\{[A-Za-z_]+(,,|\^\^)|&>>|\|&' "$CDIR/$f" > /dev/null; then fail "bash-4 construct in $f"; fi
done
/bin/bash --version | head -1

echo "== pack_nodes = 1 (explicit) == omitted"
export CAMPAIGN_TEST_BIN="$MPI_BIN" CAMPAIGN_TEST_ID="t4one_$STAMP"
sed 's/"seed": 790789,/"seed": 790789, "pack_nodes": 1,/' "$HERE/campaign_test.py" > "$T/pack1_desc_$STAMP.py"
$MK "$HERE/campaign_test.py" --root "$ROOT/omit_$STAMP" > /dev/null || fail "generate without pack_nodes"
$MK "$T/pack1_desc_$STAMP.py" --root "$ROOT/one_$STAMP" > /dev/null || fail "generate with pack_nodes = 1"
A="$ROOT/omit_$STAMP/$CAMPAIGN_TEST_ID"; B="$ROOT/one_$STAMP/$CAMPAIGN_TEST_ID"
diff -r "$A/jobs" "$B/jobs" > /dev/null && pass "job scripts identical" || fail "job scripts differ"
diff -r "$A/runs" "$B/runs" > /dev/null && pass "run inputs identical" || fail "run inputs differ"
cmp -s "$A/jobs.list" "$B/jobs.list" && cmp -s "$A/campaign.env" "$B/campaign.env" && pass "jobs.list, campaign.env identical" || fail "jobs.list / campaign.env differ"
ls "$B/scripts" | tr '\n' ' ' | grep -q node_prov && fail "node_prov.sh copied for pack_nodes = 1" || pass "no node_prov.sh for pack_nodes = 1"
hasf "one-node job: --nodes=1" "$B/jobs/c_${CAMPAIGN_TEST_ID}_hub_L8_1x1_r1.sbatch" "#SBATCH --nodes=1"
hasf "one-node job: --ntasks" "$B/jobs/c_${CAMPAIGN_TEST_ID}_hub_L8_1x1_r1.sbatch" "#SBATCH --ntasks=1"
grep -q 'slot' "$B/jobs/c_${CAMPAIGN_TEST_ID}_hub_L8_1x1_r1.sbatch" && fail "slot code in a one-node job" || pass "no slot code in a one-node job"
grep -q '"pack_nodes"' "$B/campaign.json" && fail "pack_nodes in campaign.json for pack_nodes = 1" || pass "campaign.json has no pack_nodes entry for 1"

echo "== refusals and the node limit"
export CAMPAIGN_TEST_BIN="$NOMPI_BIN" CAMPAIGN_TEST_ID="t4ref_$STAMP"
out="$(CAMPAIGN_TEST_PACK=0 $MK "$DESC" --dry-run 2>&1)"; rc=$?
expect "pack_nodes = 0 refused rc" "$rc" 1; echo "  $out"
out="$(CAMPAIGN_TEST_SETUP8=3500 $MK "$DESC" --dry-run 2>&1)"; rc=$?
expect "a single run longer than the partition limit refused rc" "$rc" 1; echo "  $out"
has "single-run refusal message" "$out" "a single run"
out="$(CAMPAIGN_TEST_PACK=9 $MK "$DESC" --dry-run 2>&1)"; rc=$?
expect "pack_nodes above the partition node limit (4): accepted rc" "$rc" 0
has "node-limit warning" "$out" "WARNING: pack_nodes = 9 exceeds the 4 nodes of partition local"
expect "cut at the node limit: one job" "$(printf '%s\n' "$out" | grep -c '^job ')" 1
has "cut at the node limit: 4 nodes" "$out" "nodes=4"

echo "== dry-run of the packed campaign"
out="$($MK "$DESC" --dry-run 2>&1)"; rc=$?
expect "dry-run rc" "$rc" 0
printf '%s\n' "$out" | head -8
has "dry-run summary" "$out" "2 jobs, 12 runs"
has "dry-run summary: packed" "$out" "packed: pack_nodes = 3, at most 3 slots per job, 4 nodes in total"
expect "dry-run: job lines" "$(printf '%s\n' "$out" | grep -c '^job ')" 2
expect "dry-run: slot lines" "$(printf '%s\n' "$out" | grep -c '^  slot ')" 4
has "dry-run: first job has 3 nodes" "$out" "pack_1x1_01: layout 1x1, nodes=3"
has "dry-run: second job has 1 node" "$out" "pack_1x1_02: layout 1x1, nodes=1"
has "dry-run: slot wall" "$out" "estimated wall  32.1 min (54 min with safety and margin; limit 60)"
expect "dry-run generated nothing" "$(ls "$ROOT" | grep -c t4ref_)" 0

echo "== generate the packed campaign"
export CAMPAIGN_TEST_ID="t4main_$STAMP"
C="$ROOT/$CAMPAIGN_TEST_ID"
$MK "$DESC" --root "$ROOT" > /dev/null || fail "generate packed"
out="$($MK "$DESC" --root "$ROOT" 2>&1)"; rc=$?
expect "regenerate same description rc" "$rc" 0; has "regenerate: unchanged" "$out" unchanged
J1="$C/jobs/c_${CAMPAIGN_TEST_ID}_pack_1x1_01.sbatch"; J2="$C/jobs/c_${CAMPAIGN_TEST_ID}_pack_1x1_02.sbatch"
expect "jobs.list" "$(wc -l < "$C/jobs.list" | tr -d ' ')" 2
for k in "--nodes=3" "--ntasks-per-node=1" "--cpus-per-task=1" "--exclusive" "--no-requeue" "--time=0:54:00"; do hasf "job 1 header $k" "$J1" "#SBATCH $k"; done
hasf "job 2: --nodes=1" "$J2" "#SBATCH --nodes=1"
hasf "job script: host list" "$J1" 'scontrol show hostnames "$SLURM_JOB_NODELIST"'
hasf "job script: slots in the background" "$J1" 'run_slot "$i" "$host" &'
hasf "job script: waits" "$J1" "wait_slots"
hasf "job script: per-slot status" "$J1" 'slot $i host=$(host_of "$i") rc=$rc runs=$n ok=$ok'
hasf "job script: launch on the slot's host" "$J1" 'profile_launch_on "$HOST" "$RANKS" "$THREADS"'
hasf "job script: run table with the slot as first field" "$J1" "0 hub_L8_1x1_r1_standard standard n/a"
expect "node_prov.sh copied" "$(ls "$C/scripts" | grep -c '^node_prov.sh$')" 1
python3 -I - "$C" <<'EOF' && pass "slot assignment = independent first-fit-decreasing; order inside the slots" || fail "slot assignment (see above)"
import json, math, sys
c = json.load(open(sys.argv[1] + "/campaign.json"))
est = {6: 300 + 20 * 3 * 0.03, 8: 900 + 20 * 3 * 0.03}
def wall(s): return int(math.ceil((s + 120) * 1.5 / 60.0)) + 5
runs = {r["run_id"]: r for r in c["runs"]}
slots = []
for rid in sorted(runs, key=lambda r: (-est[runs[r]["L"]], r)):
    e = est[runs[rid]["L"]]
    for s in slots:
        if wall(sum(est[runs[x]["L"]] for x in s) + e) <= 60:
            s.append(rid); break
    else:
        slots.append([rid])
jobs = [slots[i:i + 3] for i in range(0, len(slots), 3)]
assert len(slots) == 4 and [len(j) for j in jobs] == [3, 1], ([len(j) for j in jobs])
assert len(c["jobs"]) == 2 and c["pack_nodes"] == 3, c["jobs"]
for j, exp in zip(c["jobs"], jobs):
    assert j["nodes"] == len(exp), j["nodes"]
    assert [sorted(s["runs"]) for s in j["slots"]] == [sorted(s) for s in exp], (j["slots"], exp)
    assert j["wall_min"] <= 60 and all(s["wall_min"] <= 60 for s in j["slots"])
want = [["hub_L8_1x1_r1_standard", "hub_L8_1x1_r1_delayed_nd8"], ["hub_L8_1x1_r1_submatrix_nd8", "hub_L8_1x1_r2_delayed_nd8"],
        ["hub_L8_1x1_r2_submatrix_nd8", "hub_L8_1x1_r2_standard"]]
assert [s["runs"] for s in c["jobs"][0]["slots"]] == want, c["jobs"][0]["slots"]
assert [s["runs"] for s in c["jobs"][1]["slots"]] == [[
    "hub_L6_1x1_r1_standard", "hub_L6_1x1_r1_delayed_nd8", "hub_L6_1x1_r1_submatrix_nd8",
    "hub_L6_1x1_r2_delayed_nd8", "hub_L6_1x1_r2_submatrix_nd8", "hub_L6_1x1_r2_standard"]]
ids = [r for j in c["jobs"] for r in j["runs"]]
assert sorted(ids) == sorted(runs) and len(set(ids)) == 12            # every run in exactly one slot
assert all(runs[r]["job"] == j["job"] and runs[r]["slot"] == i for j in c["jobs"] for i, s in enumerate(j["slots"]) for r in s["runs"])
EOF

echo "== first submission: a failing run in slot 0, an MKL probe failure on the host of slot 1 (shimnode2)"
FAILRUN=hub_L8_1x1_r1_standard
SHIM_SRUN_FAIL_MATCH="/$FAILRUN/" SHIM_PROBE_FAIL_HOST=shimnode2 "$C/submit_all.sh" 2>&1 | grep -v WARNING
JD="$(ls -d "$C"/jobs/c_${CAMPAIGN_TEST_ID}_pack_1x1_01/job_* | head -1)"
cat "$JD/job.status"
expect "slot 0 failing run" "$(st "$C" hub_L8_1x1_r1_standard)" FAILED
expect "slot 0: the next run still runs" "$(st "$C" hub_L8_1x1_r1_delayed_nd8)" OK
expect "slot 1 (probe failed): first run" "$(st "$C" hub_L8_1x1_r1_submatrix_nd8)" NOT_RUN
expect "slot 1 (probe failed): second run" "$(st "$C" hub_L8_1x1_r2_delayed_nd8)" NOT_RUN
expect "slot 2 unaffected: first run" "$(st "$C" hub_L8_1x1_r2_submatrix_nd8)" OK
expect "slot 2 unaffected: second run" "$(st "$C" hub_L8_1x1_r2_standard)" OK
for r in $(cut -d' ' -f2- "$C/jobs.list" | tail -1); do expect "job 2 (one slot) $r" "$(st "$C" "$r")" OK; done
hasf "slot 0 status line" "$JD/job.status" "slot 0 host=shimnode1 rc=1 runs=2 ok=1 failed=1 invalid=0 not_run=0"
hasf "slot 1 status line" "$JD/job.status" "slot 1 host=shimnode2 rc=3 runs=2 ok=0 failed=0 invalid=0 not_run=2"
hasf "slot 2 status line" "$JD/job.status" "slot 2 host=shimnode3 rc=0 runs=2 ok=2 failed=0 invalid=0 not_run=0"
hasf "job rc" "$JD/job.status" "rc=1"
hasf "slot 1 abort reason" "$JD/slot_1/abort_reason.txt" "MKL probe gave no MKL_VERBOSE version line"
grep -q "NOT_RUN (job .* slot 1 ended before this run" "$C/runs/hub_L8_1x1_r2_delayed_nd8/status" && pass "NOT_RUN names the slot" || fail "NOT_RUN text"
cat "$C/runs/$FAILRUN/status"
python3 -I - "$C" "$SHIM_STATE/launch_hosts" <<'EOF' && pass "slot i ran on host i (launcher host per run); slots concurrent, runs of a slot sequential; node provenance per slot" || fail "hosts / concurrency (see above)"
import calendar, json, os, subprocess, sys, time
c = json.load(open(sys.argv[1] + "/campaign.json"))
runs = {r["run_id"]: r for r in c["runs"]}
jobs = {j["job"]: j for j in c["jobs"]}
real = subprocess.run(["hostname", "-s"], capture_output=True, text=True).stdout.strip()
launch = {}
for l in open(sys.argv[2]):
    cwd, host = l.split()
    launch.setdefault(cwd, []).append(host)
def prov(p):
    return dict(l.rstrip("\n").split("=", 1) for l in open(p) if "=" in l)
iv = {}
for rid, r in runs.items():
    a = os.path.join(sys.argv[1], "runs", rid, "a1")
    if not os.path.isdir(a):
        continue
    want = "shimnode%d" % (r["slot"] + 1) if jobs[r["job"]]["nodes"] > 1 else real      # sbatch shim: 1 node = the real host name
    assert launch[a] == [want], (rid, launch.get(a), r["slot"])    # the launcher got host i for slot i
    p = prov(a + "/run.prov")
    assert p["slot"] == str(r["slot"]) and p["slot_host"] == want
    t0 = calendar.timegm(time.strptime(p["run_start"][:19], "%Y-%m-%dT%H:%M:%S"))      # whole seconds, local time = same shift for all
    iv[rid] = (r["job"], r["slot"], t0, t0 + float(p["wall_s"]))
assert len(iv) == 10, len(iv)                                  # 12 runs minus the 2 never started
job1 = {k: v for k, v in iv.items() if v[0].endswith("_01")}
overl = [(a, b) for a in job1 for b in job1 if a < b and job1[a][1] != job1[b][1] and job1[a][2] < job1[b][3] and job1[b][2] < job1[a][3]]
assert overl, "no overlap between runs of different slots: not concurrent"
for a in iv:
    for b in iv:
        if a < b and iv[a][:2] == iv[b][:2]:                   # same slot: sequential (start times have 1 s resolution)
            lo, hi = (iv[a], iv[b]) if iv[a][2] <= iv[b][2] else (iv[b], iv[a])
            assert hi[2] >= lo[3] - 1.5, (a, b, lo, hi)
print("  %d overlapping pairs of runs of different slots in job 1" % len(overl))
# node provenance: one node.prov per started slot (host requested = host i), the slot job.prov = job part + node part
jd = [d for d in os.listdir(sys.argv[1] + "/jobs/" + c["jobs"][0]["job"])][0]
base = sys.argv[1] + "/jobs/" + c["jobs"][0]["job"] + "/" + jd
for i in (0, 2):
    n = prov(base + "/slot_%d/node.prov" % i)
    assert n["slot_host_requested"] == "shimnode%d" % (i + 1) and n["cpu_model"], n
    j = prov(base + "/slot_%d/job.prov" % i)
    assert j["slot_host_requested"] == n["slot_host_requested"] and j["job_script_sha256"] == prov(base + "/job.prov")["job_script_sha256"]
    assert "cpu_model" not in prov(base + "/job.prov"), "node-level key in the job-level file"
EOF
$COLLECT "$C" --strict > "$T/collect_pack1_$STAMP.txt" 2>&1; rc=$?
cat "$T/collect_pack1_$STAMP.txt" | grep -v 'WARNING  .*build_info git commit'
expect "collect --strict rc (packed campaign)" "$rc" 0
has "collect: 9 OK rows, 1 excluded, 2 missing" "$(cat "$T/collect_pack1_$STAMP.txt")" "9 OK rows, 1 excluded, 2 missing"
has "failed run excluded" "$(cat "$T/collect_pack1_$STAMP.txt")" "EXCLUDED $FAILRUN: FAILED rc=137"
$REPORT "$C/collected.json" --out "$C/report.md" > /dev/null && pass "report.py on the packed campaign" || fail "report.py"
expect "collected rows carry the slot" "$(python3 -I -c 'import json,sys; print(sum("run.slot_host" in r for r in json.load(open(sys.argv[1]))["rows"]))' "$C/collected.json")" 9

echo "== second submission: only the pending runs are executed (failed run again, slot 1 now runs, slot 2 skips)"
"$C/submit_all.sh" 2>&1 | grep -v WARNING
for r in hub_L8_1x1_r1_standard hub_L8_1x1_r1_delayed_nd8 hub_L8_1x1_r1_submatrix_nd8 hub_L8_1x1_r2_delayed_nd8 hub_L8_1x1_r2_submatrix_nd8 hub_L8_1x1_r2_standard; do
    expect "status $r" "$(st "$C" "$r")" OK
done
expect "attempts of the earlier FAILED run" "$(ls -d "$C"/runs/$FAILRUN/a[0-9]* | wc -l | tr -d ' ')" 2
expect "attempts of a run that was OK" "$(ls -d "$C"/runs/hub_L8_1x1_r2_standard/a[0-9]* | wc -l | tr -d ' ')" 1
grep -q "skip hub_L8_1x1_r2_standard (status OK)" "$C"/logs/c_${CAMPAIGN_TEST_ID}_pack_1x1_01.*.out && pass "slot 2 skips its OK runs" || fail "no skip line"
$COLLECT "$C" --strict > "$T/collect_pack2_$STAMP.txt" 2>&1; rc=$?
expect "collect --strict rc after the resubmission" "$rc" 0
has "collect: all 12 runs" "$(cat "$T/collect_pack2_$STAMP.txt")" "12 OK rows, 0 excluded, 0 missing"
$REPORT "$C/collected.json" --out "$C/report.md" > /dev/null && pass "report.py after the resubmission" || fail "report.py 2"

echo "== SIGTERM (time limit) during the first runs of all slots (simulated in hub_L8_1x1_r1_standard, slot 0)"
export CAMPAIGN_TEST_ID="t4term_$STAMP"
C2="$ROOT/$CAMPAIGN_TEST_ID"
$MK "$DESC" --root "$ROOT" > /dev/null || fail "generate term"
SHIM_SRUN_TERM_MATCH="/hub_L8_1x1_r1_standard/" "$C2/submit_all.sh" 2>&1 | grep -v WARNING
JT="$(ls -d "$C2"/jobs/c_${CAMPAIGN_TEST_ID}_pack_1x1_01/job_* | head -1)"
cat "$JT/job.status"
grep -q '^FAILED rc=143 .*SIGTERM' "$C2/runs/hub_L8_1x1_r1_standard/status" && pass "interrupted run FAILED rc=143 (SIGTERM)" || fail "term status"
expect "slot 0: the run after the interrupted one" "$(st "$C2" hub_L8_1x1_r1_delayed_nd8)" NOT_RUN
bad=0
for r in hub_L8_1x1_r1_standard hub_L8_1x1_r1_delayed_nd8 hub_L8_1x1_r1_submatrix_nd8 hub_L8_1x1_r2_delayed_nd8 hub_L8_1x1_r2_submatrix_nd8 hub_L8_1x1_r2_standard; do
    case "$(st "$C2" $r)" in FAILED|NOT_RUN) ;; *) bad=$((bad + 1)); echo "  unexpected status of $r: $(cat "$C2/runs/$r/status")" ;; esac
done
expect "after the TERM every run of job 1 is FAILED or NOT_RUN (none OK, none RUNNING)" "$bad" 0
hasf "job.status: rc=143 SIGTERM" "$JT/job.status" "rc=143 SIGTERM"
hasf "slot status line with the signal" "$JT/job.status" "slot 0 host=shimnode1 rc=143"
sleep 1
expect "no process of the campaign left running" "$(pgrep -f "$C2/runs" | wc -l | tr -d ' ')" 0
expect "job 2 (other L, other job) unaffected" "$(st "$C2" hub_L6_1x1_r1_standard)" OK

echo "== host check: the launcher does not start the slots on their hosts (names are fake here) -> every slot aborts"
export CAMPAIGN_TEST_ID="t4host_$STAMP" CAMPAIGN_TEST_REPS=1
C3="$ROOT/$CAMPAIGN_TEST_ID"
$MK "$DESC" --root "$ROOT" > /dev/null || fail "generate host"
LOCAL_NODE_HOST_CHECK=1 "$C3/submit_all.sh" 2>&1 | grep -v WARNING
JH="$(ls -d "$C3"/jobs/c_${CAMPAIGN_TEST_ID}_pack_1x1_01/job_* | head -1)"
cat "$JH/job.status"
grep -h "ABORT" "$JH"/slot_0/abort_reason.txt
hasf "host check message" "$JH/slot_0/abort_reason.txt" "instead of 'shimnode1'"
n_notrun=0; for r in $(cut -d' ' -f2- "$C3/jobs.list"); do [ "$(st "$C3" "$r")" = NOT_RUN ] && n_notrun=$((n_notrun + 1)); done
expect "all runs NOT_RUN (nothing was started on a wrong host)" "$n_notrun" "$(cut -d' ' -f2- "$C3/jobs.list" | wc -w | tr -d ' ')"
expect "no attempt directory created" "$(ls -d "$C3"/runs/*/a[0-9]* 2>/dev/null | wc -l | tr -d ' ')" 0

echo "run_pack_test: $([ $fails = 0 ] && echo ALL PASS || echo FAILED) ($fails failures); campaigns: $C $C2 $C3"
[ $fails = 0 ]

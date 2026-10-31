#!/bin/bash
# submit_all.sh [--dry-run] [--wait SECONDS]   (copied into every generated campaign directory; bash 3.2, no python)
# Submits the campaign's jobs (jobs.list: "<job_name> <run_id> ..."), idempotently:
#   - a job whose runs all have status OK is skipped;
#   - a job that is pending/running (squeue, by name) is skipped;
#   - a run with MAX_ATTEMPTS attempts (runs/<id>/a<k>) that is still not OK is given up (listed, never resubmitted);
#   - a job submitted 2 x MAX_ATTEMPTS times is given up as well (job-level aborts before any run, e.g. the fakeintel
#     ISA check, create no attempt directory);
#   - at most MAX_JOBS_IN_QUEUE jobs of this campaign are queued at any time: without --wait the script stops when the
#     limit is reached (run it again later); with --wait it polls squeue every SECONDS until everything is submitted.
# Nothing is ever deleted. Submissions are appended to jobs/<job_name>.submissions.
cd "$(dirname "$0")" || exit 2
. ./campaign.env || exit 2
DRY=0; WAIT=0
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY=1 ;;
        --wait) WAIT="$2"; shift ;;
        *) echo "usage: $0 [--dry-run] [--wait SECONDS]" >&2; exit 2 ;;
    esac
    shift
done

queued_names() { squeue -h -u "${USER:-$(id -un)}" -o '%j' 2>/dev/null | grep "^${JOB_PREFIX}_" ; }

n_sub=0; n_done=0; n_queued=0; n_gaveup=0
while read -r job runs; do
    [ -n "$job" ] || continue
    need=""; gave=""
    for r in $runs; do
        st="$(awk '{print $1; exit}' "runs/$r/status" 2>/dev/null)"
        [ "$st" = OK ] && continue
        a="$(ls -d "runs/$r"/a[0-9]* 2>/dev/null | wc -l | tr -d ' ')"
        if [ "$a" -ge "$MAX_ATTEMPTS" ]; then gave="$gave $r($st)"; else need="$need $r"; fi
    done
    if [ -n "$gave" ]; then echo "GAVE UP  $job:$gave (MAX_ATTEMPTS=$MAX_ATTEMPTS)"; n_gaveup=$((n_gaveup + 1)); fi
    if [ -z "$need" ]; then
        [ -n "$gave" ] || { echo "done     $job"; n_done=$((n_done + 1)); }
        continue
    fi
    nsubm="$(grep -c ' rc=0 ' "jobs/$job.submissions" 2>/dev/null)"
    if [ "${nsubm:-0}" -ge $((2 * MAX_ATTEMPTS)) ]; then
        echo "GAVE UP  $job: submitted $nsubm times (job-level aborts? see logs/$job.*.out)"; n_gaveup=$((n_gaveup + 1)); continue
    fi
    if queued_names | grep -q -x "$job"; then echo "queued   $job"; n_queued=$((n_queued + 1)); continue; fi
    while :; do
        inq="$(queued_names | wc -l | tr -d ' ')"
        [ "$inq" -lt "$MAX_JOBS_IN_QUEUE" ] && break
        if [ "$WAIT" -gt 0 ] 2>/dev/null; then sleep "$WAIT"
        else
            echo "queue limit reached ($inq >= MAX_JOBS_IN_QUEUE=$MAX_JOBS_IN_QUEUE): stopping; run $0 again later"
            echo "summary: submitted $n_sub, done $n_done, queued $n_queued, gave up $n_gaveup (stopped at $job)"
            exit 0
        fi
    done
    if [ "$DRY" = 1 ]; then echo "would submit $job (runs:$need)"; n_sub=$((n_sub + 1)); continue; fi
    out="$(sbatch "jobs/$job.sbatch" < /dev/null 2>&1)"; rc=$?
    echo "$(date '+%Y-%m-%dT%H:%M:%S') rc=$rc $out" >> "jobs/$job.submissions"
    if [ $rc -ne 0 ]; then echo "sbatch failed for $job: $out" >&2; exit 1; fi
    echo "submit   $job: $(printf '%s\n' "$out" | tail -1) (runs:$need)"
    n_sub=$((n_sub + 1))
done < jobs.list
echo "summary: submitted $n_sub, done $n_done, queued $n_queued, gave up $n_gaveup"

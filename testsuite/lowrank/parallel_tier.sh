#!/usr/bin/env bash
# Local parallel tier (tier 6): MPI rank count, OpenMP threads, parallel tempering, PARALLEL_PARAMS.
#
# Usage:  testsuite/lowrank/parallel_tier.sh [--rebuild] [--no-replicate] [PHASE ...] [KEY=VALUE ...]
#   PHASE  mpi | threads | tempering | pp     (default: all four, in this order)
#   --rebuild       rebuild the variant trees even if the sources did not change
#   --no-replicate  do not run the pre-registered replicate after an array-only failure
#   KEY=VALUE       extra ALF parameters applied last to every run (e.g. NSweep=10 NBin=20 for a smoke test)
#   env PT_SCHEMES="standard delayed", PT_SETS="half", PT_NP="1 4"  restrict the schemes / MPI-phase sets / MPI-phase rank counts
#   env PT_TAG=<s>  suffix of the output tree runs/parallel<s> (smoke tests, so the real results are kept)
#
# Builds (once per source state, in a COPY of the committed tree: git archive HEAD; the main tree is never touched):
#   builds/<variant>/src  for variant noMPI | MPI | TEMPERING | PARALLEL_PARAMS
#   (`source configure.sh GNU <mode> NO-INTERACTIVE && make all -j4`; log runs/parallel/build_<variant>.log).
# Runs live in runs/parallel/<phase>/<case>/..., references in runs/parallel/ref/<name>_<tag>, logs there too.
# Prints one summary line per check (PASS / FAIL / FAIL-ARRAY with the decisive numbers); the lines are also written to
# runs/parallel/summary_<timestamp>.txt. Exit 0 iff every check passed (after the replicate protocol, below).
#
# Checks (all comparisons use stat_compare.py unchanged, via parallel_compare.py):
#   mpi        sets `half` (6x6 U=4 Beta=5, the "start" set) and `tv` (STAT_SETS), each scheme, on 1, 2 and 4 MPI ranks
#              (MPI build) against ONE independent noMPI standard chain (different seed family). NSweep/NBin are the
#              stat_runs.sh values for every rank count: ALF with MPI pools the per-rank bins (every bin is the average
#              over the ranks of that bin), so the number of bins is the same and the N-rank run has N times the statistics.
#   threads    noMPI, OMP_NUM_THREADS=4 vs 1 (independent seed), `half`, each scheme.
#   tempering  TEMPERING build, mpirun -n 2, &VAR_TEMP mpi_per_parameter_set=1 (2 replicas), `half` with ham_U = 3 and 4
#              (Temp_0, Temp_1), each scheme; each replica vs an independent noMPI standard run of the same parameters
#              (the ham_U=3 reference has its own seed, reference + 500). N_Tempering_frequency=2 (the template's 10 gives
#              only ~16 accepted exchanges per run at ham_U 3 vs 4; 2 gives ~100), N_exchange_steps=6.
#   pp         PARALLEL_PARAMS build, same layout (QMC keys read from Temp_<g>/parameters), same references. No exchange steps
#              in this build: the replicas are independent chains (replica 1 is bit-identical to the noMPI run with the same seed).
#   For every non-standard run Temp_<g>/info (info) must contain the line "Sequential scheme: <NAME>" exactly once.
# Seeds: rank r of an N-rank run reads line N-r of `seeds` (Set_Random_number_Generator), so the script writes N distinct
#   seeds per run: base + 7919*i, i = 1..16. Seed family bases: reference 1000, candidate 2000 (+100 delayed, +200 submatrix: a
#   different chain per scheme), replicate candidate 3000 (+100/+200), replicate reference 4000 (references and candidates never share a seed).
# Pre-registered protocol (testsuite/lowrank/README.md): a scalar |z| >= 3 is a failure with no replicate. If ONLY an array
#   rule fails (FAIL-ARRAY), ONE replicate of the whole case is run with fresh seeds (candidate 3000, reference 4000); the
#   case passes iff the replicate passes everything (tempering/pp: the replicate concerns the replica(s) whose round-1
#   comparison failed; only those are judged in round 2). Expected false-fail rate of the array rule: ~11 % per comparison.
# Core budget: at most 4 cores busy (MPI up to 4 ranks, threads 4, tempering/PP 2 ranks).
# No `set -u` (configure.sh uses unset variables).

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALFDIR="$(cd "$HERE/../.." && pwd)"
PAR="$ALFDIR/runs/parallel${PT_TAG:-}"
BUILDS="$ALFDIR/builds${PT_TAG:-}"
RETRY_G=""
SCHEMES=(${PT_SCHEMES:-standard delayed submatrix}); PSETS=(${PT_SETS:-half tv}); NPS=(${PT_NP:-1 2 4})
REBUILD=0; REPLICATE=1
PHASES=(); EXTRA=()
for a in "$@"; do
    case "$a" in
        --rebuild) REBUILD=1 ;;
        --no-replicate) REPLICATE=0 ;;
        mpi|threads|tempering|pp) PHASES+=("$a") ;;
        *=*) EXTRA+=("$a") ;;
        *) echo "unknown argument '$a'" >&2; sed -n 2,9p "${BASH_SOURCE[0]}" >&2; exit 2 ;;
    esac
done
[ ${#PHASES[@]} -gt 0 ] || PHASES=(mpi threads tempering pp)

source "$HERE/sets.sh"
mkdir -p "$PAR"
SUMMARY="$PAR/summary_$(date +%Y%m%d_%H%M%S).txt"
: >"$SUMMARY"
NFAIL=0; NPASS=0

say() { echo "$*" | tee -a "$SUMMARY"; }

# ---------------------------------------------------------------- builds
variant_mode() { case "$1" in noMPI) echo noMPI ;; MPI) echo MPI ;; TEMPERING) echo Tempering ;; PARALLEL_PARAMS) echo PARALLEL_PARAMS ;; esac; }
src_stamp() { git -C "$ALFDIR" rev-parse HEAD:Prog HEAD:Libraries HEAD:Analysis HEAD:Makefile HEAD:configure.sh | tr '\n' ' '; }

build_variant() {   # build_variant VARIANT: copy of the committed tree, built once per source state
    local v="$1" dir="$BUILDS/$1" mode stamp
    mode="$(variant_mode "$v")"
    stamp="$(src_stamp)"
    if [ $REBUILD = 0 ] && [ -x "$dir/src/Prog/ALF.out" ] && [ "$(cat "$dir/stamp" 2>/dev/null)" = "$stamp" ]; then
        return 0
    fi
    say "build $v (configure.sh GNU $mode, copy of HEAD $(git -C "$ALFDIR" rev-parse --short HEAD)) ..."
    /bin/rm -rf "$dir"; mkdir -p "$dir/src"
    git -C "$ALFDIR" archive HEAD | tar -x -C "$dir/src"
    (
        cd "$dir/src" || exit 1
        # shellcheck disable=SC1091
        source configure.sh GNU "$mode" NO-INTERACTIVE >/dev/null 2>&1 && make all -j4
    ) >"$PAR/build_$v.log" 2>&1
    if [ ! -x "$dir/src/Prog/ALF.out" ] || [ ! -x "$dir/src/Analysis/ana.out" ]; then
        say "ERROR: build of $v failed, see $PAR/build_$v.log"; exit 2
    fi
    echo "$stamp" >"$dir/stamp"
}
bin_of() { echo "$BUILDS/$1/src/Prog/ALF.out"; }
ANA() { echo "$BUILDS/noMPI/src/Analysis/ana.out"; }

# ---------------------------------------------------------------- runs
BASEKEYS="NSweep=50 NBin=40 Ltau=0 n_skip=5 N_rebin=1"
set_keys() { local k=""; for e in "${STAT_SETS[@]}"; do IFS='|' read -r n kk <<<"$e"; [ "$n" = "$1" ] && k="$kk"; done; echo "$k"; }

write_seeds() {     # write_seeds DIR BASE: 16 distinct seeds (the template seeds file otherwise)
    awk -v b="$2" 'NR<=16{printf "%12d\n", b+7919*NR; next}{print}' "$ALFDIR/Scripts_and_Parameters_files/Start/seeds" >"$1/seeds"
}

mkparams() {        # mkparams DIR SEEDBASE KEY=VALUE ...   (parameters from the Start template, seeds overwritten)
    local dir="$1" sb="$2"; shift 2
    /bin/rm -rf "$dir"; mkdir -p "$dir"
    python3 -I "$HERE/make_params.py" "$dir" "$@" || return 1
    write_seeds "$dir" "$sb"
}

analyse() {         # analyse DIR  (same observable list as stat_runs.sh)
    local dir="$1" obs list=""
    ( cd "$dir" || exit 1
      for obs in Ener_scal Pot_scal Kin_scal Part_scal SpinZ_eq Green_eq; do
          if [ -e "$obs" ]; then list="$list $obs"
          elif [ -e Ener_scal ]; then echo "standard observable $obs missing in $dir" >&2; exit 1; fi
      done
      if [ ! -e Ener_scal ]; then for f in *_scal; do case " $list " in *" $f "*) ;; *) list="$list $f" ;; esac; done; fi
      for obs in $list; do "$(ANA)" "$obs" >"ana_$obs.log" 2>&1 || { echo "ana.out $obs failed in $dir" >&2; exit 1; }; done )
}

scheme_line_once() {  # scheme_line_once INFO SCHEME -> 0 iff "Sequential scheme: NAME" occurs exactly once
    local want n
    want="$(printf '%s' "$2" | tr '[:lower:]' '[:upper:]')"
    n=$(grep -c -E "^[[:space:]]*Sequential scheme: $want(,|[[:space:]]*$)" "$1" 2>/dev/null)
    [ "$n" = 1 ]
}

# run_job DIR NPROCS VARIANT [OMP]: run ALF in DIR (noMPI binary directly; MPI-type builds under mpirun -n NPROCS).
run_job() {
    local dir="$1" np="$2" v="$3" omp="${4:-1}" t0 t1 rc
    t0=$(date +%s)
    if [ "$v" = noMPI ]; then
        ( cd "$dir" && OMP_NUM_THREADS="$omp" "$(bin_of "$v")" >stdout.log 2>&1 ); rc=$?
    else
        ( cd "$dir" && OMP_NUM_THREADS="$omp" mpirun -n "$np" "$(bin_of "$v")" >stdout.log 2>&1 ); rc=$?
    fi
    t1=$(date +%s)
    echo $((t1 - t0)) >"$dir/wall_seconds"
    return $rc
}

# ensure_ref NAME SEEDBASE KEY=VALUE ...: independent noMPI standard chain (cached; the cache key is the build stamp + keys)
ensure_ref() {
    local name="$1" sb="$2"; shift 2
    local dir="$PAR/ref/$name" sig
    sig="$(src_stamp)|$*|$sb|${EXTRA[*]}"
    if [ -e "$dir/ref_done" ] && [ "$(cat "$dir/ref_done")" = "$sig" ]; then return 0; fi
    build_variant noMPI
    mkparams "$dir" "$sb" $BASEKEYS "$@" Sequential_scheme=\"standard\" "${EXTRA[@]}" || return 1
    run_job "$dir" 1 noMPI 1 || { say "ERROR: reference run $name failed, see $dir/stdout.log"; return 1; }
    analyse "$dir" || return 1
    echo "$sig" >"$dir/ref_done"
}

# check LABEL DIR_REF DIR_CAND -> prints/records the summary line; returns 0 pass, 1 hard fail, 3 array-only fail, 2 error
check() {
    local label="$1" a="$2" b="$3" out rc
    out="$(python3 -I "$HERE/parallel_compare.py" "$a" "$b" "$b/stat_compare.log" "$label")"; rc=$?
    say "$out"
    return $rc
}

worst() { # worst RC...: hard fail (1/2) > array-only (3) > pass
    local r=0 x
    for x in "$@"; do
        case "$x" in 1|2) echo 1; return ;; 3) r=3 ;; esac
    done
    echo $r
}

# run_case PHASE-LABEL FUNCTION ARGS...: FUNCTION ROUND ARGS... runs one case (ROUND 1 or 2 = replicate) and returns the worst rc
run_case() {
    local label="$1" fn="$2"; shift 2
    local rc
    "$fn" 1 "$@"; rc=$?
    if [ $rc = 3 ] && [ $REPLICATE = 1 ]; then
        say "replicate (array-only failure, fresh seeds) for $label:"
        "$fn" 2 "$@"; rc=$?
        [ $rc = 0 ] && say "REPLICATE-PASS $label (round 1 failed an array rule only; the replicate passes everything)"
    fi
    if [ $rc = 0 ]; then NPASS=$((NPASS + 1)); else NFAIL=$((NFAIL + 1)); say "CASE-FAIL $label (rc=$rc)"; fi
}

round_seeds() {  # round_seeds ROUND SCHEME: sets CSEED RSEED RTAG; the candidate seed family also depends on the scheme
    local off=0
    case "$2" in delayed) off=100 ;; submatrix) off=200 ;; esac
    if [ "$1" = 1 ]; then CSEED=$((2000 + off)); RSEED=1000; RTAG=r1; else CSEED=$((3000 + off)); RSEED=4000; RTAG=r2; fi
}

# ---------------------------------------------------------------- phase mpi
case_mpi() {   # ROUND SET SCHEME NRANKS
    local round="$1" set="$2" scheme="$3" n="$4" keys dir lab
    round_seeds "$round" "$scheme"
    keys="$(set_keys "$set")"
    ensure_ref "${set}_${RTAG}" "$RSEED" $keys || return 2
    dir="$PAR/mpi/$set/${scheme}_np${n}_${RTAG}"
    lab="mpi/$set/$scheme/np$n"
    mkparams "$dir" "$CSEED" $BASEKEYS $keys Sequential_scheme=\"$scheme\" "${EXTRA[@]}" || return 2
    run_job "$dir" "$n" MPI 1 || { say "ERROR $lab: run failed, see $dir/stdout.log"; return 2; }
    analyse "$dir" || { say "ERROR $lab: analysis failed"; return 2; }
    local rc=0
    if [ "$scheme" != standard ]; then
        if scheme_line_once "$dir/info" "$scheme"; then :; else say "FAIL $lab: info lacks exactly one 'Sequential scheme: $(echo $scheme | tr a-z A-Z)' line"; return 1; fi
    fi
    check "$lab" "$PAR/ref/${set}_${RTAG}" "$dir"
}

phase_mpi() {
    build_variant noMPI; build_variant MPI
    local set scheme n
    for set in "${PSETS[@]}"; do
        for scheme in "${SCHEMES[@]}"; do
            for n in "${NPS[@]}"; do
                run_case "mpi/$set/$scheme/np$n" case_mpi "$set" "$scheme" "$n"
            done
        done
    done
}

# ---------------------------------------------------------------- phase threads
case_threads() {  # ROUND SCHEME
    local round="$1" scheme="$2" keys d1 d4 lab rc
    round_seeds "$round" "$scheme"
    keys="$(set_keys half)"
    d1="$PAR/threads/${scheme}_omp1_${RTAG}"; d4="$PAR/threads/${scheme}_omp4_${RTAG}"
    lab="threads/half/$scheme/omp4-vs-omp1"
    mkparams "$d1" "$RSEED" $BASEKEYS $keys Sequential_scheme=\"$scheme\" "${EXTRA[@]}" || return 2
    mkparams "$d4" "$CSEED" $BASEKEYS $keys Sequential_scheme=\"$scheme\" "${EXTRA[@]}" || return 2
    run_job "$d1" 1 noMPI 1 || { say "ERROR $lab: omp1 run failed, see $d1/stdout.log"; return 2; }
    run_job "$d4" 1 noMPI 4 || { say "ERROR $lab: omp4 run failed, see $d4/stdout.log"; return 2; }
    analyse "$d1" && analyse "$d4" || { say "ERROR $lab: analysis failed"; return 2; }
    if [ "$scheme" != standard ]; then
        scheme_line_once "$d4/info" "$scheme" || { say "FAIL $lab: info lacks exactly one scheme line"; return 1; }
    fi
    check "$lab" "$d1" "$d4"
}

phase_threads() {
    build_variant noMPI
    local scheme
    for scheme in "${SCHEMES[@]}"; do run_case "threads/half/$scheme" case_threads "$scheme"; done
}

# ---------------------------------------------------------------- phases tempering / pp
# Layout (both builds): top-level parameters (&VAR_TEMP with mpi_per_parameter_set=1, plus the QMC keys: TEMPERING reads
# &VAR_QMC and &VAR_HAM_NAME there; PARALLEL_PARAMS reads them from Temp_<g>/parameters), Temp_<g>/parameters with
# ham_U = U_g (and the QMC keys, read by PARALLEL_PARAMS), one top-level seeds file (one seed per MPI rank).
case_temper() {   # ROUND VARIANT SCHEME
    local round="$1" v="$2" scheme="$3" keys dir lab g U rc=0 rcs=() tag refname rseed judge
    round_seeds "$round" "$scheme"
    keys="$(set_keys half)"
    tag=tempering; [ "$v" = PARALLEL_PARAMS ] && tag=pp
    dir="$PAR/$tag/${scheme}_${RTAG}"
    lab="$tag/half/$scheme"
    # Pre-registered replicate = of the comparison that failed: in round 2 only the replicas whose round-1 comparison
    # failed an array rule are judged (the other replica ran too, but its round-1 verdict stands).
    judge="0 1"; [ "$round" = 2 ] && judge="$RETRY_G"
    local common="$BASEKEYS $keys Sequential_scheme=\"$scheme\" N_Tempering_frequency=2 ${EXTRA[*]}"
    local tmpl="$ALFDIR/Scripts_and_Parameters_files/Start_Tempering"
    /bin/rm -rf "$dir"; mkdir -p "$dir"
    # shellcheck disable=SC2086
    python3 -I "$HERE/make_params.py" --template "$tmpl" "$dir" $common ham_U=4.d0 mpi_per_parameter_set=1 || return 2
    write_seeds "$dir" "$CSEED"
    for g in 0 1; do
        U=$((3 + g))
        # shellcheck disable=SC2086
        python3 -I "$HERE/make_params.py" --template "$tmpl" "$dir/Temp_$g" $common ham_U=$U.d0 mpi_per_parameter_set=1 || return 2
        /bin/rm -f "$dir/Temp_$g/seeds"
    done
    # references: ham_U = 4 is the MPI-phase reference (seed family RSEED); ham_U = 3 gets its own seed (RSEED + 500), so
    # the two reference chains are not the same random stream at nearby parameters (their fluctuations would be correlated)
    ensure_ref "half_${RTAG}" "$RSEED" $keys || return 2
    ensure_ref "half_U3_${RTAG}" "$((RSEED + 500))" $keys ham_U=3.d0 || return 2
    run_job "$dir" 2 "$v" 1 || { say "ERROR $lab: run failed, see $dir/stdout.log"; return 2; }
    [ "$round" = 1 ] && RETRY_G=""
    for g in 0 1; do
        U=$((3 + g))
        analyse "$dir/Temp_$g" || { say "ERROR $lab: analysis of Temp_$g failed"; return 2; }
        case " $judge " in *" $g "*) ;; *) continue ;; esac
        if [ "$scheme" != standard ]; then
            if ! scheme_line_once "$dir/Temp_$g/info" "$scheme"; then
                say "FAIL $lab/Temp_$g: info lacks exactly one 'Sequential scheme: $(echo $scheme | tr a-z A-Z)' line"; rcs+=(1); continue
            fi
        fi
        refname="half_U${U}_${RTAG}"; [ $U = 4 ] && refname="half_${RTAG}"
        check "$lab/Temp_$g(U=$U)" "$PAR/ref/$refname" "$dir/Temp_$g"; rc=$?
        rcs+=($rc)
        [ "$round" = 1 ] && [ $rc = 3 ] && RETRY_G="$RETRY_G $g"
    done
    return "$(worst "${rcs[@]}")"
}

phase_tp() {  # phase_tp VARIANT
    local v="$1" scheme
    build_variant noMPI; build_variant "$v"
    for scheme in "${SCHEMES[@]}"; do
        run_case "$v/half/$scheme" case_temper "$v" "$scheme"
    done
}

# ---------------------------------------------------------------- main
say "parallel tier on HEAD $(git -C "$ALFDIR" rev-parse --short HEAD), phases: ${PHASES[*]}, extra: ${EXTRA[*]:-none}, replicate=$REPLICATE"
T0=$(date +%s)
for ph in "${PHASES[@]}"; do
    t0=$(date +%s)
    case "$ph" in
        mpi) phase_mpi ;;
        threads) phase_threads ;;
        tempering) phase_tp TEMPERING ;;
        pp) phase_tp PARALLEL_PARAMS ;;
    esac
    say "phase $ph: $(( $(date +%s) - t0 )) s"
done
say "parallel tier: $NPASS cases passed, $NFAIL failed, $(( $(date +%s) - T0 )) s total; summary $SUMMARY"
[ $NFAIL = 0 ]

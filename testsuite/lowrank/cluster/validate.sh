#!/usr/bin/env bash
# Cluster validation (Phase 7), all checks concurrently on ONE node (e.g. inside one Slurm job; binaries from build_variants.sh and
# `parallel_tier.sh --build-only`, built beforehand on the login node; nothing is compiled here except the ctest tree if missing).
#
# Usage:  testsuite/lowrank/cluster/validate.sh [GROUP ...]      (default: all groups)
#   core     ctest; regress_vs_baseline.sh (13 sets byte-identical to the baseline commit)
#   traj     compare_trajectory.sh (trace binary builds/ALF_trace.out): sets start tv kondo jz pv_adiab tv_proj tv_flux z2 z2gt
#            spin_peierls, each: standard vs delayed, standard vs submatrix, and both again with N_delay=-1 (z2gt: N_delay=-1 with
#            Global_tau); continuous with HMC (Leapfrog_steps=10, dt 0.05: acceptance ~1) + sequential low-rank, both schemes;
#            the HMC acceptance of each run is reported
#   extras   HDF5 attributes (h5dump -A, delayed vs standard); debug shadow runs (tv kondo jz z2 spin_peierls x delayed/submatrix);
#            timers build (timers.txt); MPI trajectory identity (4 ranks: delayed, submatrix, both with Lowrank_cross_slice=.T., vs
#            standard; confout_0 and bins, compare_runs.py)
#   ptier    parallel_tier.sh phases concurrently, each in its own output tree: mpi (half and tv separately), threads, tempering, pp
# Environment: ALF_MACHINE, ALF_ENV_SH (module setup), HDF5_MODULE (module providing h5dump, default hdf5/1.10.7-intel21),
#   VAL_JOBS (concurrent trajectory sets, default 11; e.g. 4 on a laptop),
#   MPIRUN (default "mpirun"; single node: I_MPI_HYDRA_BOOTSTRAP=fork, I_MPI_PIN=0 so concurrent runs do not pin to the same cores).
# Output: runs/cluster_val/<group>/...; one line "CHECK <name> PASS|FAIL <detail>" per check in runs/cluster_val/summary.txt and on
#   stdout; exit 0 iff every check passed. No `set -u` (configure.sh uses unset variables).
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LR="$(cd "$HERE/.." && pwd)"
ALFDIR="$(cd "$LR/../.." && pwd)"
OUT="$ALFDIR/runs/cluster_val"
SUM="$OUT/summary.txt"
mkdir -p "$OUT"
{ [ -z "${ALF_ENV_SH:-}" ] || source "$ALF_ENV_SH" || true; }
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP="${I_MPI_HYDRA_BOOTSTRAP:-fork}" I_MPI_PIN="${I_MPI_PIN:-0}"
MPIRUN="${MPIRUN:-mpirun}"
check() { echo "CHECK $1 $2 ${3:-}" | tee -a "$SUM"; }
GROUPS_RUN="${*:-core traj extras ptier}"

g_core() {  # regress first (its incremental `make program` must not race the ctest link of Prog/*.o), then ctest
    local d="$OUT/core" rc n; mkdir -p "$d"
    "$LR/regress_vs_baseline.sh" >"$d/regress.log" 2>&1; rc=$?
    n=$(grep -c " IDENTICAL" "$d/regress.log")
    if [ $rc = 0 ]; then check regress PASS "$n IDENTICAL"; else check regress FAIL "rc=$rc $(grep -v IDENTICAL "$d/regress.log" | tail -3 | tr '\n' ' ')"; fi
    ( cd "$ALFDIR" && source configure.sh "${ALF_MACHINE:-GNU}" noMPI NO-INTERACTIVE >/dev/null 2>&1
      module load cmake >/dev/null 2>&1
      mkdir -p testsuite/tests && cd testsuite/tests && cmake .. >"$d/cmake.log" 2>&1 && make -j8 >"$d/make_tests.log" 2>&1
      ctest -j8 >"$d/ctest.log" 2>&1 )
    if grep -q "100% tests passed" "$d/ctest.log"; then check ctest PASS "$(grep 'tests passed' "$d/ctest.log")"
    else check ctest FAIL "$(grep -E 'tests passed|Failed|\*\*\*' "$d/ctest.log" | head -5 | tr '\n' ' ')"; fi
}

traj_one() {  # traj_one SET: all comparisons of one set, serially (they share runs/traj/<SET>)
    local s="$1" d="$OUT/traj" a b extra rc
    for spec in "standard delayed" "standard submatrix" "standard delayed N_delay=-1" "standard submatrix N_delay=-1"; do
        set -- $spec; a=$1; b=$2; shift 2; extra="$*"
        [ "$s" = continuous ] && extra="${extra:+$extra }HMC=.T. Leapfrog_steps=10 Delta_t_Langevin_HMC=0.05"   # Leapfrog_steps=0 (template) accepts nothing
        ALF_TRACE_BIN="$ALFDIR/builds/ALF_trace.out" "$LR/compare_trajectory.sh" "$s" "$a" "$b" $extra \
            >"$d/${s}_${b}${extra// /_}.log" 2>&1
        rc=$?
        local hmc=""; [ "$s" = continuous ] && hmc=" HMC acceptance $(awk '/Acceptance_HMC/{print $3}' "$ALFDIR/runs/traj/$s/$b/info")"
        if [ $rc = 0 ]; then check "traj:$s:$a-$b${extra:+:$extra}" PASS "$(grep -E 'IDENTICAL|ROUNDOFF' "$d/${s}_${b}${extra// /_}.log" | tail -1)$hmc"
        else check "traj:$s:$a-$b${extra:+:$extra}" FAIL "rc=$rc $(tail -2 "$d/${s}_${b}${extra// /_}.log" | tr '\n' ' ')"; fi
    done
}
g_traj() {
    mkdir -p "$OUT/traj"
    local pids=() s
    for s in start tv kondo jz pv_adiab tv_proj tv_flux z2 z2gt spin_peierls continuous; do
        traj_one "$s" & pids+=($!)
        if [ ${#pids[@]} -ge "${VAL_JOBS:-11}" ]; then wait "${pids[@]}"; pids=(); fi
    done
    [ ${#pids[@]} = 0 ] || wait "${pids[@]}"
}

mkrun() {  # mkrun DIR KEY=VALUE...: run directory from the Start template
    /bin/rm -rf "$1"; python3 -I "$LR/make_params.py" "$@" >/dev/null
}
ex_hdf5() {
    local d="$OUT/extras/hdf5" bin="$ALFDIR/builds/cluster/hdf5/src/Prog/ALF.out" s ok=1
    module load "${HDF5_MODULE:-hdf5/1.10.7-intel21}" >/dev/null 2>&1
    for s in standard delayed submatrix; do
        mkrun "$d/$s" Sequential_scheme=\"$s\" NSweep=10 NBin=2
        ( cd "$d/$s" && "$bin" >stdout.log 2>&1 && h5dump -A -g / data.h5 >attrs.txt 2>&1 )
    done
    grep -q "Sequential_scheme" "$d/standard/attrs.txt" && ok=0
    for s in delayed submatrix; do
        for a in Sequential_scheme N_delay_requested N_delay_used; do grep -q "ATTRIBUTE \"$a\"" "$d/$s/attrs.txt" || ok=0; done
    done
    if [ $ok = 1 ]; then check hdf5_attributes PASS "delayed/submatrix: Sequential_scheme N_delay_requested N_delay_used; standard: none"
    else check hdf5_attributes FAIL "see $d/*/attrs.txt"; fi
}
ex_debug() {
    local d="$OUT/extras/debug" bin="$ALFDIR/builds/cluster/debug/src/Prog/ALF.out" s sc common pids=() name c
    source "$LR/sets.sh"
    for s in tv kondo jz z2 spin_peierls; do
        for entry in "${SETS[@]}" "${TRAJ_ONLY_SETS[@]}"; do IFS='|' read -r name c _ <<<"$entry"; [ "$name" = "$s" ] && common="$c"; done
        for sc in delayed submatrix; do
            ( mkrun "$d/${s}_$sc" $common Sequential_scheme=\"$sc\" NSweep=3 NBin=2
              cd "$d/${s}_$sc" && "$bin" >stdout.log 2>&1; echo $? >rc ) & pids+=($!)
        done
    done
    wait "${pids[@]}"
    for s in tv kondo jz z2 spin_peierls; do for sc in delayed submatrix; do
        local r="$d/${s}_$sc"
        if [ "$(cat "$r/rc")" = 0 ] && [ -f "$r/lowrank_debug_summary.txt" ]; then
            check "debug:$s:$sc" PASS "$(grep -E 'checks|max relative' "$r/lowrank_debug_summary.txt" | tr -s ' ' | tr '\n' ' ')"
        else check "debug:$s:$sc" FAIL "rc=$(cat "$r/rc") $(tail -3 "$r/stdout.log" | tr '\n' ' ')"; fi
    done; done
}
ex_timers() {
    local d="$OUT/extras/timers/delayed_L12"
    mkrun "$d" Sequential_scheme=\"delayed\" L1=12 L2=12 Beta=2.d0 Ltau=0 NSweep=5 NBin=2
    ( cd "$d" && "$ALFDIR/builds/ALF_timers.out" >stdout.log 2>&1 )
    if [ -s "$d/timers.txt" ]; then check timers PASS "$(wc -l <"$d/timers.txt") lines in timers.txt"
    else check timers FAIL "no timers.txt ($(tail -2 "$d/stdout.log" | tr '\n' ' '))"; fi
}
ex_mpi_identity() {
    local d="$OUT/extras/mpi4" bin="$ALFDIR/builds/cluster/mpi/src/Prog/ALF.out" v rc
    for v in "standard|" "delayed|" "submatrix|" "delayed_xs|Lowrank_cross_slice=.T." "submatrix_xs|Lowrank_cross_slice=.T."; do
        IFS='|' read -r name kv <<<"$v"
        mkrun "$d/$name" Sequential_scheme=\"${name%_xs}\" NSweep=10 NBin=2 $kv
        ( cd "$d/$name" && $MPIRUN -n 4 "$bin" >stdout.log 2>&1 )
    done
    for name in delayed submatrix delayed_xs submatrix_xs; do
        python3 -I "$LR/compare_runs.py" "$d/standard" "$d/$name" standard "${name%_xs}" >"$d/cmp_$name.log" 2>&1; rc=$?
        if [ $rc = 0 ]; then check "mpi4_identity:$name" PASS "$(tail -1 "$d/cmp_$name.log")"
        else check "mpi4_identity:$name" FAIL "$(tail -2 "$d/cmp_$name.log" | tr '\n' ' ')"; fi
    done
}
g_extras() {
    mkdir -p "$OUT/extras"
    ex_hdf5 & local p1=$!
    ex_debug & local p2=$!
    ex_timers & local p3=$!
    ex_mpi_identity & local p4=$!
    wait $p1 $p2 $p3 $p4
}

g_ptier() {  # exit code of parallel_tier.sh = verdict (0 iff every case passed, after the pre-registered replicate protocol)
    local d="$OUT/ptier" pids=() t
    mkdir -p "$d"
    ( PT_SETS=half PT_TAG=_cl_mpi_half "$LR/parallel_tier.sh" mpi >"$d/mpi_half.log" 2>&1; echo $? >"$d/mpi_half.rc" ) & pids+=($!)
    ( PT_SETS=tv   PT_TAG=_cl_mpi_tv   "$LR/parallel_tier.sh" mpi >"$d/mpi_tv.log" 2>&1; echo $? >"$d/mpi_tv.rc" ) & pids+=($!)
    for t in threads tempering pp; do
        ( PT_TAG=_cl_$t "$LR/parallel_tier.sh" $t >"$d/$t.log" 2>&1; echo $? >"$d/$t.rc" ) & pids+=($!)
    done
    wait "${pids[@]}"
    for t in mpi_half mpi_tv threads tempering pp; do
        if [ "$(cat "$d/$t.rc")" = 0 ]; then check "ptier:$t" PASS "$(grep 'parallel tier:' "$d/$t.log" | tail -1)"
        else check "ptier:$t" FAIL "rc=$(cat "$d/$t.rc") $(grep -E 'FAIL|ERROR|parallel tier:' "$d/$t.log" | head -3 | tr '\n' ' ')"; fi
    done
}

echo "# validate.sh $(date -u +%FT%TZ) host $(hostname) commit $(git -C "$ALFDIR" rev-parse --short HEAD) groups: $GROUPS_RUN" | tee -a "$SUM"
pids=()
for g in $GROUPS_RUN; do "g_$g" & pids+=($!); done
wait "${pids[@]}"
nfail=$(grep -c "^CHECK .* FAIL" "$SUM")
echo "# done $(date -u +%FT%TZ): $(grep -c '^CHECK .* PASS' "$SUM") PASS, $nfail FAIL" | tee -a "$SUM"
[ "$nfail" = 0 ]

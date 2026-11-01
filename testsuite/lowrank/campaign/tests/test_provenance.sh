#!/bin/bash
# Unit checks of provenance.sh / rank_wrap.sh (bash 3.2): MKL_VERBOSE fixtures (SuperMUC AVX-512, Helma with and
# without fakeintel), the fakeintel ISA rule, the per-rank preload and /proc/self/maps check. Exit 0 iff all pass.
# Usage: testsuite/lowrank/campaign/tests/test_provenance.sh SCRATCH_DIR
HERE="$(cd "$(dirname "$0")" && pwd)"
CDIR="$(dirname "$HERE")"
TMP="${1:?scratch directory}"
mkdir -p "$TMP"
. "$CDIR/provenance.sh"
fails=0
check() { if [ "$2" = "$3" ]; then echo "PASS $1 ($2)"; else echo "FAIL $1: got '$2', want '$3'"; fails=$((fails + 1)); fi; }

for fx in supermuc:AVX-512:"oneMKL 2025.0 Update 3" helma_nofakeintel:AVX2:"oneMKL 2024.0 Update 2" \
          helma_fakeintel:AVX-512:"oneMKL 2024.0 Update 2"; do
    name="${fx%%:*}"; rest="${fx#*:}"; isa="${rest%%:*}"; ver="${rest#*:}"
    prov_mkl_parse "$(head -1 "$HERE/fixtures/mkl_verbose_$name.txt")"
    check "mkl version $name" "$MKL_VERSION" "$ver"
    check "mkl isa $name" "$MKL_ISA" "$isa"
    check "mkl threading $name" "$MKL_THREADING" "intel_thread"
done
prov_mkl_parse "MKL_VERBOSE Intel(R) MKL 2020.0 Update 4 Product build 20200917 for Intel(R) 64 architecture Intel(R) Advanced Vector Extensions 2 (Intel(R) AVX2) enabled processors, Lnx 2.50GHz lp64 sequential"
check "mkl legacy version" "$MKL_VERSION" "Intel(R) MKL 2020.0 Update 4"
check "mkl legacy threading" "$MKL_THREADING" "sequential"
prov_mkl_parse "MKL_VERBOSE something unexpected"
check "mkl unparsable -> unavailable" "${MKL_ISA%%(*}" "unavailable"

r=0; prov_mkl_check on AVX-512 || r=1;  check "fakeintel on + AVX-512 passes" "$r" 0
r=0; prov_mkl_check on AVX2 || r=1;     check "fakeintel on + AVX2 fails" "$r" 1
r=0; prov_mkl_check off AVX2 || r=1;    check "fakeintel off + AVX2 passes" "$r" 0

# rank wrapper: preload only in the wrapper, maps check against the fixture /proc
lib="$TMP/libfakeintel.so"; : > "$lib"
d="$TMP/rw1"; mkdir -p "$d"
FAKEINTEL=on FAKEINTEL_LIB="$lib" PROV_PROC_ROOT="$HERE/fixtures" PMI_RANK=3 \
    "$CDIR/rank_wrap.sh" probe "$d" /bin/sh -c 'echo "LD_PRELOAD=$LD_PRELOAD" > '"$d/child_env.txt"
check "wrapper rc (fakeintel on, mapped)" "$?" 0
check "wrapper rank from PMI_RANK" "$(sed -n 's/^rank=//p' "$d/rank_probe_3.prov")" 3
check "wrapper maps check" "$(sed -n 's/^fakeintel_maps=//p' "$d/rank_probe_3.prov")" "loaded (2 mappings in /proc/self/maps)"
check "child sees LD_PRELOAD" "$(cat "$d/child_env.txt")" "LD_PRELOAD=$lib"
d="$TMP/rw2"; mkdir -p "$d"
lib2="$TMP/libotherpreload.so"; : > "$lib2"
FAKEINTEL=on FAKEINTEL_LIB="$lib2" PROV_PROC_ROOT="$HERE/fixtures" "$CDIR/rank_wrap.sh" probe "$d" /usr/bin/true
check "wrapper rc (library not mapped)" "$?" 91
d="$TMP/rw3"; mkdir -p "$d"
FAKEINTEL=on FAKEINTEL_LIB="$TMP/missing.so" "$CDIR/rank_wrap.sh" probe "$d" /usr/bin/true
check "wrapper rc (library missing)" "$?" 90
d="$TMP/rw4"; mkdir -p "$d"
FAKEINTEL=off "$CDIR/rank_wrap.sh" probe "$d" /bin/sh -c 'echo "LD_PRELOAD=${LD_PRELOAD:-unset}" > '"$d/child_env.txt"
check "fakeintel off: no preload" "$(cat "$d/child_env.txt")" "LD_PRELOAD=unset"
d="$TMP/rw5"; mkdir -p "$d"; : > "$d/x.prov"
prov_kv "$d/x.prov" a 1; prov_kv "$d/x.prov" a 2 2>/dev/null
check "prov_kv refuses a duplicate key" "$?" 1
prov_kv "$d/x.prov" b ""
check "empty value -> unavailable(empty)" "$(sed -n 's/^b=//p' "$d/x.prov")" "unavailable(empty)"
# pinning verification (Intel MPI I_MPI_DEBUG table and srun cpu-bind masks; good and bad cases)
pc() { prov_pin_check "$HERE/fixtures/pin_$1.txt" "$2" "$3" "${4:-1}" >/dev/null; echo $?; }
pc_rc() { prov_pin_check "$@" >/dev/null; echo $?; }
check "pin intel 4x12 good" "$(pc intel_good 4 12)" 0
check "pin intel overlapping domains" "$(pc intel_overlap 4 12)" 1
check "pin intel wrong rank count" "$(pc intel_good 8 6)" 1
check "pin srun 8x48 good (384-bit masks)" "$(pc srun_good 8 48)" 0
check "pin srun overlapping masks" "$(pc srun_overlap 2 2)" 1
check "pin srun 4 cpus per rank for THREADS=2, no SMT" "$(pc srun_wrongsize 2 2 1)" 1
check "pin srun 4 cpus per rank for THREADS=2, 2 hw threads per core" "$(pc srun_wrongsize 2 2 2)" 0
check "pin message" "$(prov_pin_check "$HERE/fixtures/pin_srun_overlap.txt" 2 2)" "BAD overlapping domains (srun cpu-bind masks): h01c01:cpu1(r0,r1)"
: > "$TMP/empty_pin.txt"
check "pin: no record is bad" "$(prov_pin_check "$TMP/empty_pin.txt" 1 1)" "BAD no parsable I_MPI_DEBUG pin lines or cpu-bind masks"
check "pin intel 2 nodes, same cpu sets per node: valid" "$(pc intel_2node_good 4 24)" 0
check "pin intel same-node overlap" "$(pc intel_samenode_overlap 2 24)" 1
check "pin srun 2 nodes, same mask per node: valid" "$(pc srun_2node_good 2 2)" 0
check "pin intel lists/ranges without braces" "$(pc intel_nobrace_good 4 4)" 0
check "pin intel unparsable line quoted" "$(prov_pin_check "$HERE/fixtures/pin_intel_unparsable.txt" 2 4)" \
    "BAD unparsable I_MPI_DEBUG pin line: [0] MPI startup(): 0       201234   i01r01c01s01  cpu0-cpu3"
check "pin both sources: Intel table judged, srun masks informational" "$(prov_pin_check "$HERE/fixtures/pin_both.txt" 2 2)" \
    "OK I_MPI_DEBUG pin table (srun masks informational): 2 ranks x 2 cpus, disjoint"
F="$HERE/fixtures"
check "pin intel truncated list, OMP affinity 1x48 decides" "$(prov_pin_check "$F/pin_intel_truncated.txt" 1 48 1 "$F/omp_1x48_good.txt")" \
    "OK OpenMP thread affinity (I_MPI_DEBUG pin list truncated): 1 ranks x 48 cpus, disjoint"
check "pin intel truncated list, OMP affinity 2x24 decides" "$(prov_pin_check "$F/pin_intel_truncated_2r.txt" 2 24 1 "$F/omp_2x24_good.txt")" \
    "OK OpenMP thread affinity (I_MPI_DEBUG pin list truncated): 2 ranks x 24 cpus, disjoint"
check "pin intel truncated list, unpinned OMP threads" "$(pc_rc "$F/pin_intel_truncated_2r.txt" 2 24 1 "$F/omp_2x24_unpinned.txt")" 1
check "pin intel truncated list, no OMP record" "$(prov_pin_check "$F/pin_intel_truncated.txt" 1 48 1)" \
    "BAD truncated I_MPI_DEBUG pin list and no OMP_AFFINITY record: [0] MPI startup(): 0       972908   i01r01c03s02  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,"
check "pin intel complete list ignores OMP record" "$(prov_pin_check "$F/pin_both.txt" 2 2 1 "$F/omp_2x24_unpinned.txt")" \
    "OK I_MPI_DEBUG pin table (srun masks informational): 2 ranks x 2 cpus, disjoint"
# threading layer
r=0; prov_mkl_threads_check 2 sequential || r=1;   check "THREADS=2 + sequential MKL fails" "$r" 1
r=0; prov_mkl_threads_check 1 sequential || r=1;   check "THREADS=1 + sequential MKL passes" "$r" 0
r=0; prov_mkl_threads_check 4 intel_thread || r=1; check "THREADS=4 + intel_thread passes" "$r" 0
echo "test_provenance: $([ $fails = 0 ] && echo ALL PASS || echo FAILED) ($fails failures)"
[ $fails = 0 ]

# provenance.sh: provenance helpers, sourced by the generated job scripts and by rank_wrap.sh (bash 3.2, no python).
#
# Every field is written as one line KEY=VALUE (prov_kv refuses a duplicate key). A source that this machine cannot
# provide is written as "unavailable(<reason>)", never left empty; collect.py --strict rejects such values, so the local
# test (shims + fixtures) proves that every field has a real source. Multi-line dumps go to separate files next to the
# .prov file. PROV_SYSFS_ROOT / PROV_PROC_ROOT (default empty = the real /sys and /proc) let the local test read fixtures.

prov_sha256() {
    if [ ! -r "$1" ]; then echo "unavailable(no file $1)"; return 0; fi
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | awk '{print $1}'
    else echo "unavailable(no sha256sum/shasum)"; fi
}

# epoch seconds, with nanoseconds where date supports %N (GNU date on Linux), else whole seconds (BSD date)
prov_now() {
    local t
    t="$(date +%s.%N)"
    case "$t" in *N*|*.) t="$(date +%s)" ;; esac
    printf '%s' "$t"
}

# one line, tabs/newlines -> ' ; ', trimmed; empty -> unavailable(empty)
prov_flat() {
    local v
    v="$(printf '%s' "$1" | tr '\t' ' ' | tr -s ' ' | awk 'NF {gsub(/^ +| +$/, ""); out = (out == "" ? $0 : out " ; " $0)} END {print out}')"
    [ -n "$v" ] || v="unavailable(empty)"
    printf '%s' "$v"
}

# prov_kv FILE KEY VALUE
prov_kv() {
    if [ -f "$1" ] && grep -q "^$2=" "$1"; then
        echo "provenance: duplicate key '$2' in $1" >&2
        return 1
    fi
    printf '%s=%s\n' "$2" "$(prov_flat "$3")" >> "$1"
}

# prov_cmd OUTFILE CMD ...: run CMD, keep stdout+stderr in OUTFILE; return its status (127 if not found)
prov_cmd() {
    local out="$1"
    shift
    if ! command -v "$1" >/dev/null 2>&1 && ! type "$1" >/dev/null 2>&1; then
        echo "command not found: $1" > "$out"
        return 127
    fi
    "$@" > "$out" 2>&1
}

# prov_field FILE SED_EXPR: first match of a sed -n expression in FILE, or unavailable(...)
prov_field() {
    local v=""
    [ -r "$1" ] && v="$(sed -n "$2" "$1" | head -1)"
    [ -n "$v" ] || v="unavailable(not in $(basename "$1"))"
    printf '%s' "$v"
}

prov_sysfs() {
    local f="${PROV_SYSFS_ROOT}$1"
    if [ -r "$f" ]; then head -1 "$f"; else printf 'unavailable(no %s)' "$1"; fi
}

# prov_freq FILE TAG: clock state (governor, driver, current frequency min/mean/max over all CPUs, max frequency, boost)
prov_freq() {
    local f="$1" t="$2" base="${PROV_SYSFS_ROOT}/sys/devices/system/cpu" stats boost
    prov_kv "$f" "governor_$t" "$(prov_sysfs /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)"
    prov_kv "$f" "freq_driver_$t" "$(prov_sysfs /sys/devices/system/cpu/cpu0/cpufreq/scaling_driver)"
    stats="$(cat "$base"/cpu[0-9]*/cpufreq/scaling_cur_freq 2>/dev/null | awk '
        {n++; s += $1; if (n == 1 || $1 < mn) mn = $1; if (n == 1 || $1 > mx) mx = $1}
        END {if (n) printf "min=%d mean=%d max=%d kHz over %d cpus", mn, s / n, mx, n}')"
    [ -n "$stats" ] || stats="unavailable(no scaling_cur_freq)"
    prov_kv "$f" "cur_freq_$t" "$stats"
    prov_kv "$f" "max_freq_$t" "$(prov_sysfs /sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq) kHz (scaling_max_freq)"
    if [ -r "$base/cpufreq/boost" ]; then boost="boost=$(head -1 "$base/cpufreq/boost") (cpufreq/boost)"
    elif [ -r "$base/intel_pstate/no_turbo" ]; then boost="no_turbo=$(head -1 "$base/intel_pstate/no_turbo") (intel_pstate)"
    else boost="unavailable(no cpufreq/boost or intel_pstate/no_turbo)"; fi
    prov_kv "$f" "boost_$t" "$boost"
}

# prov_mkl_parse LINE: sets MKL_VERSION, MKL_ISA, MKL_THREADING from the first MKL_VERBOSE line (the version line)
prov_mkl_parse() {
    local l="$1"
    MKL_VERSION="$(printf '%s' "$l" | grep -o -E '(oneMKL|Intel\(R\) MKL|Intel\(R\) oneAPI Math Kernel Library) ([Vv]ersion )?[0-9][0-9.]*( Update [0-9]+)?' | head -1)"
    case "$l" in
        *AVX-512*|*AVX512*) MKL_ISA="AVX-512" ;;
        *"AVX2"*)           MKL_ISA="AVX2" ;;
        *"(Intel(R) AVX)"*) MKL_ISA="AVX" ;;
        *SSE4*)             MKL_ISA="SSE4" ;;
        *)                  MKL_ISA="" ;;
    esac
    MKL_THREADING="$(printf '%s' "$l" | grep -o -E '(intel_thread|gnu_thread|tbb_thread|pgi_thread|sequential)' | tail -1)"
    [ -n "$MKL_VERSION" ] || MKL_VERSION="unavailable(no version in MKL_VERBOSE line)"
    [ -n "$MKL_ISA" ] || MKL_ISA="unavailable(no ISA in MKL_VERBOSE line)"
    [ -n "$MKL_THREADING" ] || MKL_THREADING="unavailable(no threading layer in MKL_VERBOSE line)"
}

# prov_mkl_check FAKEINTEL ISA: non-zero if fakeintel is on and MKL did not take the AVX-512 path
prov_mkl_check() {
    [ "$1" = on ] || return 0
    [ "$2" = "AVX-512" ]
}

# prov_rss FILE STYLE: max resident set size in kB from /usr/bin/time output (gnu: -v, bsd: -l), empty if absent
prov_rss() {
    [ -r "$1" ] || return 0
    case "$2" in
        gnu) sed -n 's/^[[:space:]]*Maximum resident set size (kbytes):[[:space:]]*\([0-9]*\).*/\1/p' "$1" | head -1 ;;
        bsd) awk '/maximum resident set size/ {printf "%d\n", $1 / 1024; exit}' "$1" ;;
    esac
}

# prov_job DIR: job-level provenance (node, CPU, memory, clock, Slurm, MPI, modules, environment) into DIR/job.prov
prov_job() {
    local d="$1" f="$1/job.prov" v
    prov_kv "$f" campaign_id "$CAMPAIGN_ID" || return 1
    prov_kv "$f" job_name "$JOB_NAME"
    prov_kv "$f" slurm_job_id "${SLURM_JOB_ID:-unavailable(no SLURM_JOB_ID)}"
    prov_kv "$f" cluster "${SLURM_CLUSTER_NAME:-unavailable(no SLURM_CLUSTER_NAME)}"
    prov_kv "$f" node "${SLURMD_NODENAME:-unavailable(no SLURMD_NODENAME)}"
    prov_kv "$f" hostname "$(hostname)"
    prov_kv "$f" profile "$PROFILE_NAME"
    prov_kv "$f" job_start "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    # Slurm view of the job
    prov_cmd "$d/scontrol.txt" scontrol show job "${SLURM_JOB_ID:-0}"
    prov_kv "$f" node_list "$(prov_field "$d/scontrol.txt" 's/.*[[:space:]]NodeList=\([^[:space:]]*\).*/\1/p')"
    prov_kv "$f" partition "$(prov_field "$d/scontrol.txt" 's/.*[[:space:]]Partition=\([^[:space:]]*\).*/\1/p')"
    v="$(prov_field "$d/scontrol.txt" 's/.*OverSubscribe=\([^[:space:]]*\).*/\1/p')"
    case "$v" in NO) v="yes (OverSubscribe=NO)" ;; unavailable*) ;; *) v="no (OverSubscribe=$v)" ;; esac
    prov_kv "$f" exclusive "$v"
    prov_kv "$f" ear "$EAR_SETTING"
    # CPU, NUMA, memory
    prov_cmd "$d/lscpu.txt" lscpu
    prov_kv "$f" cpu_model "$(prov_field "$d/lscpu.txt" 's/^Model name:[[:space:]]*//p')"
    prov_kv "$f" cpus "$(prov_field "$d/lscpu.txt" 's/^CPU(s):[[:space:]]*//p')"
    prov_kv "$f" sockets "$(prov_field "$d/lscpu.txt" 's/^Socket(s):[[:space:]]*//p')"
    prov_kv "$f" cores_per_socket "$(prov_field "$d/lscpu.txt" 's/^Core(s) per socket:[[:space:]]*//p')"
    prov_kv "$f" threads_per_core "$(prov_field "$d/lscpu.txt" 's/^Thread(s) per core:[[:space:]]*//p')"
    prov_kv "$f" lscpu_flags_avx512f "$( [ -r "$d/lscpu.txt" ] && { grep -q '^Flags:.*avx512f' "$d/lscpu.txt" && echo yes || echo no; } )"
    prov_cmd "$d/numactl.txt" numactl -H
    prov_kv "$f" numa_nodes "$(prov_field "$d/numactl.txt" 's/^available:[[:space:]]*\([0-9]*\) nodes.*/\1/p')"
    prov_cmd "$d/free.txt" free -m
    prov_kv "$f" mem_total_mib "$(prov_field "$d/free.txt" 's/^Mem:[[:space:]]*\([0-9]*\).*/\1/p')"
    prov_freq "$f" job_start
    # MPI, modules, environment of the runs
    prov_cmd "$d/mpi_version.txt" $MPI_VERSION_CMD
    v="$(grep -v '^[[:space:]]*$' "$d/mpi_version.txt" | head -1)"
    case "$v" in *[Vv]ersion*) ;; *) v="$v $(grep -m1 -i 'version' "$d/mpi_version.txt")" ;; esac
    prov_kv "$f" mpi_version "$v"
    prov_kv "$f" i_mpi_root "${I_MPI_ROOT:-unset}"
    if type module >/dev/null 2>&1; then module list > "$d/modules.txt" 2>&1; v="$(cat "$d/modules.txt")"
    else v="unavailable(no module command)"; echo "$v" > "$d/modules.txt"; fi
    prov_kv "$f" modules_job "$v"
    env | LC_ALL=C sort > "$d/env.txt"
    prov_kv "$f" pin_env "$(env | LC_ALL=C sort | grep -E '^(I_MPI_PIN|I_MPI_DEBUG|KMP_|OMP_PLACES|OMP_PROC_BIND|SLURM_CPU_BIND)' | tr '\n' ' ')"
    prov_kv "$f" omp_num_threads "${OMP_NUM_THREADS:-unset}"
    prov_kv "$f" mkl_num_threads "${MKL_NUM_THREADS:-unset}"
    prov_kv "$f" mkl_dynamic "${MKL_DYNAMIC:-unset}"
    prov_kv "$f" omp_display_affinity "${OMP_DISPLAY_AFFINITY:-unset}"
    prov_kv "$f" mkl_verbose_in_runs "${MKL_VERBOSE:-unset}"
    prov_kv "$f" prov_sysfs_root "${PROV_SYSFS_ROOT:-/ (real)}"
    prov_kv "$f" prov_proc_root "${PROV_PROC_ROOT:-/ (real)}"
}

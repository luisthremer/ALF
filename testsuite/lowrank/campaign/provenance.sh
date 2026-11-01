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

# The job provenance is built from three parts. prov_job (one node per job) writes all of them into DIR/job.prov;
# a packed job (several nodes, job_template_packed.sh) writes the job-level parts once (prov_job_common, prov_job_env)
# and the node-level part (prov_node_hw, called by node_prov.sh on each host) once per host.

# prov_job_common DIR FILE: campaign, job, Slurm view of the job (node list, partition, exclusive), EAR
prov_job_common() {
    local d="$1" f="$2" v
    prov_kv "$f" campaign_id "$CAMPAIGN_ID" || return 1
    prov_kv "$f" job_name "$JOB_NAME"
    prov_kv "$f" slurm_job_id "${SLURM_JOB_ID:-unavailable(no SLURM_JOB_ID)}"
    prov_kv "$f" cluster "${SLURM_CLUSTER_NAME:-unavailable(no SLURM_CLUSTER_NAME)}"
    prov_kv "$f" profile "$PROFILE_NAME"
    prov_kv "$f" job_start "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    prov_cmd "$d/scontrol.txt" scontrol show job "${SLURM_JOB_ID:-0}"
    prov_kv "$f" node_list "$(prov_field "$d/scontrol.txt" 's/.*[[:space:]]NodeList=\([^[:space:]]*\).*/\1/p')"
    prov_kv "$f" partition "$(prov_field "$d/scontrol.txt" 's/.*[[:space:]]Partition=\([^[:space:]]*\).*/\1/p')"
    v="$(prov_field "$d/scontrol.txt" 's/.*OverSubscribe=\([^[:space:]]*\).*/\1/p')"
    case "$v" in NO) v="yes (OverSubscribe=NO)" ;; unavailable*) ;; *) v="no (OverSubscribe=$v)" ;; esac
    prov_kv "$f" exclusive "$v"
    prov_kv "$f" ear "$EAR_SETTING"
}

# prov_node_hw DIR FILE: CPU, NUMA, memory and clock state of the node this runs on (dumps into DIR)
prov_node_hw() {
    local d="$1" f="$2"
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
}

# prov_job_env DIR FILE: MPI, modules, runtime environment file, pinning / threading environment of the runs
prov_job_env() {
    local d="$1" f="$2" v
    prov_cmd "$d/mpi_version.txt" $MPI_VERSION_CMD
    v="$(grep -v '^[[:space:]]*$' "$d/mpi_version.txt" | head -1)"
    case "$v" in *[Vv]ersion*) ;; *) v="$v $(grep -m1 -i 'version' "$d/mpi_version.txt")" ;; esac
    prov_kv "$f" mpi_version "$v"
    prov_kv "$f" i_mpi_root "${I_MPI_ROOT:-unset}"
    if type module >/dev/null 2>&1; then module list > "$d/modules.txt" 2>&1; v="$(cat "$d/modules.txt")"
    else v="unavailable(no module command)"; echo "$v" > "$d/modules.txt"; fi
    prov_kv "$f" modules_job "$v"
    # the profile's runtime environment file (SuperMUC: RUNTIME_ENV_SH) is sourced by profile_setup, not part of the campaign
    prov_kv "$f" runtime_env_sh "${RUNTIME_ENV_SH:-n/a (profile sources none)}"
    if [ -n "${RUNTIME_ENV_SH:-}" ]; then prov_kv "$f" runtime_env_sh_sha256 "$(prov_sha256 "$RUNTIME_ENV_SH")"
    else prov_kv "$f" runtime_env_sh_sha256 "n/a (profile sources none)"; fi
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

# prov_job DIR: job-level provenance of a one-node job (node, CPU, memory, clock, Slurm, MPI, modules, environment)
# into DIR/job.prov
prov_job() {
    local d="$1" f="$1/job.prov"
    prov_job_common "$d" "$f" || return 1
    prov_kv "$f" node "${SLURMD_NODENAME:-unavailable(no SLURMD_NODENAME)}"
    prov_kv "$f" hostname "$(hostname)"
    prov_node_hw "$d" "$f"
    prov_job_env "$d" "$f"
}

# prov_pin_check FILE RANKS THREADS [THREADS_PER_CORE]: verify the pinning record of a run.
# Sources: Intel MPI I_MPI_DEBUG>=4 table lines "[0] MPI startup(): <rank> <pid> <node> <cpus>" with <cpus> as
# "{0,1,2}", "0,1,2", "0-11" or a mix of lists and ranges (braces optional), and srun --cpu-bind=verbose lines
# "cpu-bind=MASK - <node>, task <rank> <local> [<pid>]: mask 0x<hex> set" (masks of any width). When both are present,
# the Intel table is judged (it is the final pinning, set inside MPI_Init) and the srun masks are informational.
# Ownership is keyed by (node, cpu): the same cpu numbers on different nodes are no overlap. Requires exactly RANKS
# ranks, each pinned to THREADS cpus (or THREADS x THREADS_PER_CORE hardware threads for core binding with SMT), and
# pairwise disjoint sets per node. Prints "OK <source>: ..." / "BAD ..."; returns 0 / 1. Table lines with a rank
# that cannot be parsed unambiguously are quoted in the BAD reason.
prov_pin_check() {
    [ -r "$1" ] || { echo "BAD no pinning record file"; return 1; }
    awk -v R="$2" -v T="$3" -v TPC="${4:-1}" '
        function put(src, r, node, c) { k = ++cnt[src]; er[src, k] = r; en[src, k] = node; ec[src, k] = c }
        function hexset(src, r, node, h,   i, L, v, b) {
            h = tolower(h); sub(/^0x/, "", h); L = length(h)
            for (i = L; i >= 1; i--) {
                v = index("0123456789abcdef", substr(h, i, 1)) - 1
                for (b = 0; b < 4; b++) if (int(v / 2 ^ b) % 2 == 1) put(src, r, node, (L - i) * 4 + b)
            }
        }
        function listset(src, r, node, s,   k, a, parts, lh, lo, hi, c) {
            gsub(/[{} \t]/, "", s); k = split(s, parts, ",")
            for (a = 1; a <= k; a++) {
                if (parts[a] == "") continue
                if (split(parts[a], lh, "-") == 2) { lo = lh[1] + 0; hi = lh[2] + 0 } else { lo = parts[a] + 0; hi = lo }
                for (c = lo; c <= hi; c++) put(src, r, node, c)
            }
        }
        /MPI startup\(\):/ {
            s = $0; sub(/.*MPI startup\(\):[ \t]*/, "", s); nw = split(s, w, /[ \t]+/)
            if (w[1] !~ /^[0-9]+$/) next                       # header, version and environment lines
            set = s; for (i = 1; i <= 3; i++) sub(/^[^ \t]+[ \t]+/, "", set); sub(/[ \t]+$/, "", set)
            if (nw >= 4 && w[2] ~ /^[0-9]+$/ && set ~ /^\{?[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*\}?$/) {
                listset("I", w[1] + 0, w[3], set)
            } else if (raw == "") raw = $0
            next
        }
        /cpu[-_]bind/ && /mask 0x[0-9a-fA-F]+/ && /task[ \t]+[0-9]+/ {
            node = $0; sub(/^.*cpu[-_]bind[^-]*-[ \t]*/, "", node); sub(/,.*/, "", node)
            s = $0; sub(/.*task[ \t]+/, "", s); r = s + 0
            m = $0; sub(/.*mask[ \t]+/, "", m); sub(/[ \t].*/, "", m)
            hexset("S", r, node, m)
        }
        END {
            if (cnt["I"] == 0 && raw != "") {   # an Intel table exists but cannot be judged: never fall back to srun
                print "BAD unparsable I_MPI_DEBUG pin line: " raw; exit 1
            }
            if (cnt["I"] > 0) { src = "I"; name = "I_MPI_DEBUG pin table" (cnt["S"] > 0 ? " (srun masks informational)" : "") }
            else if (cnt["S"] > 0) { src = "S"; name = "srun cpu-bind masks" }
            else {
                print "BAD no parsable I_MPI_DEBUG pin lines or cpu-bind masks" (raw != "" ? "; unparsed: " raw : "")
                exit 1
            }
            ov = ""; nov = 0
            for (k = 1; k <= cnt[src]; k++) {
                r = er[src, k]; node = en[src, k]; c = ec[src, k]
                ranks[r] = 1
                if ((r, node, c) in has) continue
                has[r, node, c] = 1; n[r]++
                if (((node, c) in owner) && owner[node, c] != r) {
                    if (nov < 5) ov = ov " " node ":cpu" c "(r" owner[node, c] ",r" r ")"
                    nov++
                } else owner[node, c] = r
            }
            nr = 0; for (r in ranks) nr++
            if (nr != R) { print "BAD " nr " ranks in the " name ", expected " R (raw != "" ? "; unparsed: " raw : ""); exit 1 }
            for (r in ranks) if (n[r] != T && n[r] != T * TPC) {
                print "BAD rank " r " pinned to " n[r] " cpus, expected " T (TPC > 1 ? " or " T * TPC : "") " (" name ")"; exit 1
            }
            if (nov) { print "BAD overlapping domains (" name "):" ov; exit 1 }
            for (r in ranks) { sz = n[r]; break }
            print "OK " name ": " R " ranks x " sz " cpus, disjoint"
        }' "$1"
}

# prov_mkl_threads_check THREADS LAYER: non-zero if more than one thread is requested but MKL runs sequentially
prov_mkl_threads_check() {
    [ "$1" -le 1 ] && return 0
    [ "$2" != sequential ]
}

# Machine profile: SuperMUC-NG Phase 1 (LRZ). Sourced by make_campaign.py (static part) and by every job (functions).
# bash 3.2 compatible; never use the double-bracket test syntax here: make_campaign.py refuses any profile that still
# contains the placeholder marker (two opening square brackets), so fill every such value before generating a campaign.
# Facts (docs/audit/2026-10-31-gap-audit.md section 4): 2x Xeon Platinum 8174, 48 cores, AVX-512, thin nodes 96 GB
# (90 usable); partitions test (30 min), micro / general / fat (60 h); slurm_setup is mandatory; EAR off -> 1.9 GHz;
# pinning I_MPI_PIN_CELL=core, I_MPI_PIN_DOMAIN=omp:compact.

# ---- static part (read by make_campaign.py) ---------------------------------------------------------------------------
PROFILE_NAME=supermuc
CORES_PER_NODE=48
NODE_MEM_GB=90
# partition:max wall minutes
PARTITIONS="test:30 micro:3600 general:3600 fat:3600"
ACCOUNT="[[LRZ project id (sbatch --account)]]"
# extra #SBATCH directives, one per line (without the '#SBATCH ' prefix)
SBATCH_EXTRA="--get-user-env
[[EAR switch, e.g. --ear=off (verify the option name in doku.lrz.de energy-aware-runtime)]]"
EAR_SETTING="[[EAR state recorded in the logs, e.g. off (fixed 1.9 GHz)]]"
FAKEINTEL_DEFAULT=off

# ---- runtime part (used inside the job on the compute node) ----------------------------------------------------------
RUNTIME_MODULES="[[modules loaded when configure.sh SUPERMUC-NG built the binary, e.g. intel-mpi intel-mkl hdf5/1.10.7-intel21]]"
MPI_VERSION_CMD="mpiexec --version"
TIME_STYLE=gnu                      # /usr/bin/time -v per rank; falls back to sacct MaxRSS (job level) if absent
FAKEINTEL_LIB=""                    # Intel CPU: never preloaded
BLAS_PROBE_CMD="[[blas_bench binary of Task T2 with tiny arguments, e.g. /path/blas_bench.out --probe]]"
MKL_PROBE_REQUIRED=1
PROV_SYSFS_ROOT=""
PROV_PROC_ROOT=""

profile_setup() {
    module load slurm_setup || return 1
    module load $RUNTIME_MODULES || return 1
    export I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_DOMAIN=omp:compact I_MPI_DEBUG=4
    export OMP_PLACES=cores OMP_PROC_BIND=close
    return 0
}

# profile_launch RANKS THREADS CMD ...   (one process per rank; pin map printed by I_MPI_DEBUG=4)
profile_launch() {
    local np="$1"
    shift 2
    mpiexec -n "$np" "$@"
}

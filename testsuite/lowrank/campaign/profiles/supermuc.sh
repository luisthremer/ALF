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
# partition:max wall minutes; the allocation may use only the test partition (one test job runs at a time, at most 3 queued)
PARTITIONS="test:30"
# partition:max nodes per job (packed jobs, description key pack_nodes, are cut at this node count)
PARTITION_MAX_NODES="test:16"
ACCOUNT="pn73xu"
# extra #SBATCH directives, one per line (without the '#SBATCH ' prefix). EAR: the controller verifies that sbatch accepts
# --ear=off (doku.lrz.de energy-aware-runtime) before the first campaign.
SBATCH_EXTRA="--get-user-env
--ear=off"
EAR_SETTING="off (sbatch --ear=off)"
FAKEINTEL_DEFAULT=off

# ---- runtime part (used inside the job on the compute node) ----------------------------------------------------------
# (the hdf5 module that configure.sh loads for the build is not needed at run time: benchmark binaries are built without HDF5)
# environment file of the runtime stack: it does 'module purge; module unload stack; module load stack/22.2.1' (stack/24.5.0 is sticky
# and conflicts); profile_setup sources it and then loads slurm_setup. The modules actually loaded are recorded (module list).
RUNTIME_ENV_SH="/hppfs/work/pn73xu/di38zuw/CLAUDE/env_alf.sh"
MPI_VERSION_CMD="mpiexec --version"
TIME_STYLE=gnu                      # /usr/bin/time -v per rank; falls back to sacct MaxRSS (job level) if absent
FAKEINTEL_LIB=""                    # Intel CPU: never preloaded
BLAS_PROBE_CMD="/hppfs/work/pn73xu/di38zuw/CLAUDE/ALF/builds/blas_bench/blas_bench_mpi.out --probe"
MKL_PROBE_REQUIRED=1
NODE_HOST_CHECK=1                   # packed jobs: a slot aborts if the launcher did not start it on its host
PROV_SYSFS_ROOT=""
PROV_PROC_ROOT=""

profile_setup() {
    [ -r "$RUNTIME_ENV_SH" ] || { echo "profile_setup: $RUNTIME_ENV_SH not readable" >&2; return 1; }
    . "$RUNTIME_ENV_SH" || return 1
    module load slurm_setup || return 1
    # configure.sh (SUPERMUC-NG) unsets FORT_BLOCKSIZE at build time; the runs must not see it either
    unset FORT_BLOCKSIZE
    # packed jobs: under Slurm Intel MPI takes the per-node process counts from the allocation (respects Slurm's placement) and
    # then ignores -hosts / -ppn; 0 makes `mpiexec -hosts H -n N -ppn N` authoritative, so that slot i really runs on host i
    # (every run's rank hosts are checked against its slot host afterwards, NODE_HOST_CHECK)
    export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
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

# profile_launch_on HOST RANKS THREADS CMD ...   (packed jobs: the same on a given host of the allocation; concurrent calls
# for different hosts are independent steps of the job)
profile_launch_on() {
    local host="$1" np="$2"
    shift 3
    mpiexec -hosts "$host" -n "$np" -ppn "$np" "$@"
}

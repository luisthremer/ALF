# Machine profile: Helma CPU partition (NHR@FAU). Sourced by make_campaign.py (static part) and by every job (functions).
# bash 3.2 compatible; never use the double-bracket test syntax here: make_campaign.py refuses any profile that still
# contains the placeholder marker (two opening square brackets), so fill every such value before generating a campaign.
# Facts (docs/audit/2026-10-31-gap-audit.md section 4): 2x AMD EPYC 9965 (Zen5c), 384 cores, 2.25 GHz (boost 3.35),
# SMT off, 768 GB; partitions cpu (24 h), preempt_cpu (1 node, 48 h); srun --cpu-bind=cores, --cpus-per-task, --exclusive.
# MKL on AMD: without fakeintel MKL takes its AVX2 path; with fakeintel (LD_PRELOAD of a library that makes
# mkl_serv_intel_cpu_true return 1, set ONLY in the per-rank wrapper) the AVX-512 path. Modules: configure.sh HELMA.

# ---- static part (read by make_campaign.py) ---------------------------------------------------------------------------
PROFILE_NAME=helma
CORES_PER_NODE=384
NODE_MEM_GB="[[usable memory per node in GB (node has 768 GB)]]"
PARTITIONS="cpu:1440 preempt_cpu:2880"
# partition:max nodes per job (packed jobs, description key pack_nodes)
PARTITION_MAX_NODES="[[max nodes per job and partition, e.g. cpu:N preempt_cpu:1]]"
ACCOUNT="[[NHR project account (sbatch --account), or empty if not required]]"
SBATCH_EXTRA=""
EAR_SETTING="none (no EAR at NHR@FAU)"
FAKEINTEL_DEFAULT=off

# ---- runtime part ------------------------------------------------------------------------------------------------------
RUNTIME_MODULES="intel/2025.3.1 intelmpi/2021.17.0 mkl/2024.2.2"
MPI_VERSION_CMD="mpirun --version"
TIME_STYLE=gnu
FAKEINTEL_LIB="[[absolute path of the fakeintel shared library on Helma]]"
BLAS_PROBE_CMD="[[blas_bench binary of Task T2 with tiny arguments, e.g. /path/blas_bench.out --probe]]"
MKL_PROBE_REQUIRED=1
NODE_HOST_CHECK=1                   # packed jobs: a slot aborts if the launcher did not start it on its host
PROV_SYSFS_ROOT=""
PROV_PROC_ROOT=""

profile_setup() {
    module --force switch gpu-env/2025 cpu-env/2026 || return 1
    module load $RUNTIME_MODULES || return 1
    export I_MPI_PMI_LIBRARY="[[libpmi2.so path for Intel MPI under srun, see doc.nhr.fau.de]]"
    # I_MPI_DEBUG=4: Intel MPI prints its final pin table after MPI_Init; prov_pin_check judges that table and keeps
    # the srun cpu-bind masks (launcher masks before MPI_Init) as informational
    export I_MPI_DEBUG=4
    export OMP_PLACES=cores OMP_PROC_BIND=close
    return 0
}

# profile_launch RANKS THREADS CMD ...   (srun prints the cpu-bind masks with 'verbose')
profile_launch() {
    local np="$1" nt="$2"
    shift 2
    srun -n "$np" --cpus-per-task="$nt" --cpu-bind=cores,verbose "$@"
}

# profile_launch_on HOST RANKS THREADS CMD ...   (packed jobs; UNTESTED on Helma: same srun style as profile_launch, pinned to
# HOST with --nodelist. Concurrent steps on different hosts of one allocation may need --exact (Slurm >= 21.08) or --overlap)
profile_launch_on() {
    local host="$1" np="$2" nt="$3"
    shift 3
    srun --nodelist="$host" -N 1 -n "$np" --cpus-per-task="$nt" --cpu-bind=cores,verbose "$@"
}

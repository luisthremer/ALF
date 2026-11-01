# Machine profile: local test machine (the Mac), used with the PATH shims in tests/shims (srun, sbatch, squeue,
# scontrol, lscpu, numactl, free, taskset, module, blas_probe). Without the shims the Slurm commands do not exist here.
# bash 3.2 compatible; the LOCAL_* variables let the test driver point provenance at fixtures.

# ---- static part (read by make_campaign.py) ---------------------------------------------------------------------------
PROFILE_NAME=local
CORES_PER_NODE="${LOCAL_CORES_PER_NODE:-4}"
NODE_MEM_GB=16
PARTITIONS="local:60"
PARTITION_MAX_NODES="local:4"
ACCOUNT=""
SBATCH_EXTRA=""
EAR_SETTING="none (local machine)"
FAKEINTEL_DEFAULT=off

# ---- runtime part ------------------------------------------------------------------------------------------------------
RUNTIME_MODULES=""
MPI_VERSION_CMD="mpirun --version"
case "$(uname -s)" in Darwin) TIME_STYLE=bsd ;; *) TIME_STYLE=gnu ;; esac
FAKEINTEL_LIB="${LOCAL_FAKEINTEL_LIB:-}"
BLAS_PROBE_CMD="${LOCAL_BLAS_PROBE_CMD:-blas_probe}"
MKL_PROBE_REQUIRED="${LOCAL_MKL_PROBE_REQUIRED:-1}"
NODE_HOST_CHECK="${LOCAL_NODE_HOST_CHECK:-0}"      # profile_launch_on ignores the host here (one machine), so no check
PROV_SYSFS_ROOT="${LOCAL_SYSFS_ROOT:-}"
PROV_PROC_ROOT="${LOCAL_PROC_ROOT:-}"

profile_setup() {
    export VECLIB_MAXIMUM_THREADS="$OMP_NUM_THREADS"
    export OMP_WAIT_POLICY=PASSIVE OMP_PLACES=cores OMP_PROC_BIND=close
    return 0
}

profile_launch() {
    local np="$1" nt="$2"
    shift 2
    srun -n "$np" --cpus-per-task="$nt" --cpu-bind=cores,verbose "$@"
}

# profile_launch_on HOST RANKS THREADS CMD ...   (packed jobs; the local machine has one host: HOST is ignored, but handed to
# the srun shim, which logs it for the tests)
profile_launch_on() {
    local host="$1"
    shift
    SHIM_LAUNCH_HOST="$host" profile_launch "$@"
}

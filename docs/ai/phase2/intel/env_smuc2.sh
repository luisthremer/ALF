# SuperMUC-NG Phase 2 (Intel Max 1550): toolchain and run environment of the gpu-intel checks (P2.4 step 2).
# Module set verified on gpu-port (docs/gpu/platforms/intel_pvc_supermuc_ng2.md, update 2026-11-13).
module load intel/2025.3.0 intel-mpi/2021.15.0 >/dev/null 2>&1
module unload intel-mkl >/dev/null 2>&1
module load intel-mkl/2025.3.0 >/dev/null 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY
export ZE_FLAT_DEVICE_HIERARCHY=FLAT                         # every tile is a device (8 per node)
export LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576   # gpu-port: L = 24 10.6 -> 7.1 s/sweep
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

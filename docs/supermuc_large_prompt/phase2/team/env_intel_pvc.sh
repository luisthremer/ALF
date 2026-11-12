# Source this: toolchain for ALF GPU=intel on SuperMUC-NG Phase 2 (PVC).
# Never pipe 'module load' (a pipe runs it in a subshell and nothing is loaded).
# ifx >= 2025.3 needed: omp_lib of 2025.1/2025.2 lacks omp_get_interop_* (used by the oneMKL backend). HDF5 .mod files come from ifx 2025.1.0 (LRZ build).
module load slurm_setup 2>/dev/null
module load intel/2025.3.0 intel-mpi/2021.15.0
module unload intel-mkl            # login default is 2025.2.0 (would mix with ifx 2025.1)
module load intel-mkl/2025.3.0 hdf5/1.14.5-intel25      # serial HDF5, static Fortran libs
export OMP_TARGET_OFFLOAD=MANDATORY

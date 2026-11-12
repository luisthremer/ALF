#!/bin/bash
#SBATCH --job-name=lp_hello
#SBATCH --output=lp_hello.%j.log
#SBATCH --error=lp_hello.%j.err
#SBATCH --time=00:10:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

echo "== submit/start: $(date -Is)  host=$(hostname)  jobid=$SLURM_JOB_ID"
module load slurm_setup
module load intel/2025.1.0 intel-mpi/2021.15.0 >/dev/null 2>&1; module -t list 2>&1; which ifx sycl-ls mpiexec
echo "== lscpu"; lscpu | egrep "Model name|^CPU\(s\)|Thread|Core|Socket|NUMA node\(s\)"
echo "== mem"; free -g | head -2
echo "== sycl-ls"; sycl-ls 2>&1
echo "== xpu-smi"; xpu-smi discovery 2>&1 | head -40; ls /dev/dri
echo "== env"; env | sort | egrep -i "^(SLURM|ZE_|ONEAPI|I_MPI|SYCL|HOME|SCRATCH|WORK|TMP|LRZ_|DAOS|CPATH|LD_LIB)" | cut -c1-200
echo "== done $(date -Is)"

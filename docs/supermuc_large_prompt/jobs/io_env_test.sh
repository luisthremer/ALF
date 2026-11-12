#!/bin/bash
#SBATCH --job-name=lp_io
#SBATCH --output=lp_io.%j.log
#SBATCH --error=lp_io.%j.err
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# Short I/O and environment probe, 1 node. Writes only into $SCRATCH/large_prompt/io_test (removed at the end).
module load slurm_setup
echo "== host $(hostname) $(date -Is)"
echo "== env (job)"; env | sort | egrep -i "^(HOME|SCRATCH|WORK|work|TMP|TMPDIR|LRZ_|DSS|DAOS|SLURM_JOB_(PART|QOS|ACC)|ZE_|ONEAPI|I_MPI|OMP_|SYCL)" | cut -c1-200
echo "== df"; df -hT $HOME $SCRATCH /tmp /dev/shm 2>&1; mount | egrep " /tmp | /dev/shm | /hppfs | /dss/dsshome1 " | cut -c1-160
echo "== daos/tmp"; ls -ld /tmp /dev/shm; df -h /tmp | tail -1
D=$SCRATCH/large_prompt/io_test; mkdir -p $D
w() { # name dir
  echo "-- write 2 GiB bs=1M on $1"; dd if=/dev/zero of=$2/ddtest bs=1M count=2048 conv=fdatasync 2>&1 | tail -1
  echo "-- read 2 GiB on $1 (after drop not possible; may be cached)"; dd if=$2/ddtest of=/dev/null bs=1M 2>&1 | tail -1
  echo "-- 2000 small files (4k) on $1"; mkdir -p $2/small; ( time (for i in $(seq 1 2000); do head -c 4096 /dev/zero > $2/small/f$i; done; sync) ) 2>&1 | grep real
  rm -rf $2/ddtest $2/small
}
w SCRATCH $D
mkdir -p $HOME/large_prompt_iotest; w HOME $HOME/large_prompt_iotest; rmdir $HOME/large_prompt_iotest
mkdir -p /dev/shm/lp_io; w /dev/shm /dev/shm/lp_io; rmdir /dev/shm/lp_io
rmdir $D
echo "== quota"; (dssusrinfo all 2>&1 | grep -v '^\*' | head -30) ; mmlsquota 2>&1 | head -5
echo "== done $(date -Is)"

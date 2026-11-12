#!/bin/bash
#SBATCH --job-name=p2_fp
#SBATCH --output=p2_fp.%j.log
#SBATCH --error=p2_fp.%j.err
#SBATCH --time=00:10:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
# per-kernel launch and first-private-array counts, small case (L=6, same code path), NEW build aottim and OLD build aottim
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 OMP_NUM_THREADS=4 SEEDS=$L0/src_p2/Scripts_and_Parameters_files/Start/seeds
for v in new old; do
  [ $v = new ] && B=$L0/build_p2_aottim/Prog/ALF.out || B=$L0/build_aottim/Prog/ALF.out
  D=$L0/work_p2fp/$v; rm -rf $D; mkdir -p $D; cd $D; cp $SEEDS seeds
  python3 - <<'PY'
import re,os,sys
sys.argv=['x','T','L=6','mz=F','scheme=submatrix','ndelay=8','nsweep=2','nbin=1','U=5.3d0','symm=T','beta=6.d0','nphi=18']
src=open(os.environ['J']+'/run_cmp_sng.py').read() if 'J' in os.environ else ''
PY
  cp $J/../src_p2/Scripts_and_Parameters_files/Start/parameters . 2>/dev/null
  # reuse the generated L6 parameter file of the correctness runs
  cp $L0/work_p2gg/runs/G_L6_mzF/gpu/parameters parameters
  sed -i -e 's/NSweep *= *[0-9]*/NSweep = 2/' -e 's/NBin *= *[0-9]*/NBin = 1/' parameters
  LIBOMPTARGET_DEBUG=1 $B 2>&1 | awk '/first-private array/{fp++} /Launching target execution/{n=$0; sub(/.*execution /,"",n); sub(/ with.*/,"",n); sub(/__omp_offloading_[0-9a-f]+_[0-9a-f]+_/,"",n); sub(/_l[0-9]+.*/,"",n); L[n]++; F[n]+=fp; fp=0} END{for(k in L) printf "KERN %-45s launches=%d firstpriv=%d\n",k,L[k],F[k]}' | sort -k3 -t= -n -r | head -30 > fp_$v.txt
  echo "== $v"; cat fp_$v.txt
  cd $J
done
echo done

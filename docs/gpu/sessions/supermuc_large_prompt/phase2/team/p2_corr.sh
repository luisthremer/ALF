#!/bin/bash
#SBATCH --job-name=p2_corr
#SBATCH --output=p2_corr.%j.log
#SBATCH --error=p2_corr.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# prototyper p2: (1) new GPU vs CPU same-seed, (2) new GPU vs old GPU bitwise, (3) planted defect vs old GPU, (4) launch counts
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; J=$L0/jobs; cd $J
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0
echo "host $(hostname) $(date -Is)"
echo "######## CHECK 1: new GPU vs CPU (tm_cmp semantics, no lead-defects)"
bash tm_cmp_p2.sh gpu-txt p2 nodefect 2>&1 | tail -60
export SEEDS=$L0/src_p2/Scripts_and_Parameters_files/Start/seeds
OLD=$L0/build_gpu-txt/Prog/ALF.out; NEW=$L0/build_p2_gpu-txt/Prog/ALF.out; DEF=$L0/build_p2_gpu-txt_def/Prog/ALF.out
gg() { # tag gpubin args...   (A = old GPU build, B = given build, both on the GPU, tol 0)
  local tag=$1 bin=$2; shift 2
  W=$L0/work_p2gg CPUBIN=$OLD GPUBIN=$bin CMP_TOL=0 python3 run_cmp_gg.py $tag "$@" 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "
}
echo "######## CHECK 2: new GPU vs old GPU, bitwise (tol 0)"
gg G_L6_mzF  $NEW L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
gg G_L6_mzT  $NEW L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
gg G_L12_mzF $NEW L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72
echo "######## CHECK 3: planted defect (row factor * (1+1e-12) in k_blk_scale) vs old GPU: must FAIL"
gg D_L6_mzF  $DEF L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
gg D_L12_mzF $DEF L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72
echo "######## COUNTS: new build (aottim), L=24 one tile, LIBOMPTARGET_DEBUG=1, 1 vs 2 sweeps"
export OMP_NUM_THREADS=8
B=$L0/build_p2_aottim/Prog/ALF.out; W=$L0/work_p2cnt
for ns in 1 2; do
  D=$W/s$ns; rm -rf $D; mkdir -p $D; cd $D; cp $L0/ref_L24/seeds .
  sed -e "s/NSweep *= *[0-9]*/NSweep = $ns/" -e 's/NBin *= *[0-9]*/NBin = 1/' -e 's/CPU_MAX *= *[0-9.]*/CPU_MAX = 0.0/' $L0/ref_L24/parameters > parameters
  LIBOMPTARGET_DEBUG=1 $B 2>&1 | awk '/for first-private array/{fp++} /Launching target execution/{la++} /Executing a kernel/{ek++} /Submitting firstprivate data/{sf++} /Submitted copy/{cp++} END{print "COUNT fp="fp" launch="la" kernel="ek" fpsubmit="sf" copies="cp}' | tee cnt.txt
  echo "s$ns: $(cat cnt.txt)"
done
echo "done $(date -Is)"

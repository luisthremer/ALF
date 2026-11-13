#!/bin/bash
#SBATCH --job-name=qf_alf
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:50:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
# Task 3: ALF level. NEW = build from git archive HEAD (p8/qrfix/build_*), OLD = p5 build (build_txtnew / build_aottt, QR header before 7e8b6032).
# V7 = LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 (p5lib.sh var V7)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; Q=$L0/p8/qrfix; cd $L0/p5/jobs
source p5lib.sh; export W=$Q/work_q3; mkdir -p $W
OLDT=$C/build_txtnew/Prog/ALF.out; NEWT=$Q/build_txtnew/Prog/ALF.out; OLDA=$C/build_aottt/Prog/ALF.out; NEWA=$Q/build_aottt/Prog/ALF.out
echo "host $(hostname) $(date -Is)  V7 = $(var V7)"
md5sum $OLDT $NEWT $OLDA $NEWA
cmpdirs() { local d=$W/gg/runs/$1 n=0 bad=0 f; for f in $(ls $d/cpu); do case $f in info|run.log|parameters|seeds|*_info|RUNNING|TREE) continue;; esac; n=$((n+1)); cmp -s $d/cpu/$f $d/gpu/$f || { bad=$((bad+1)); echo "  DIFF $f"; }; done; echo "  byte-compare $1: $n files, $bad differ"; }
pgl() { for k in cpu gpu; do echo "  PG $1 $k: $(grep -h 'Precision Green' $W/gg/runs/$1/$k/info)"; done; }
S=$L0/src/Scripts_and_Parameters_files/Start/seeds
gg() { # tag oldbin newbin  [key=value ...]   (env of the caller applies to both runs)
  local tag=$1 ob=$2 nb=$3; shift 3
  W=$W/gg CPUBIN=$ob GPUBIN=$nb CMP_TOL=0 SEEDS=$S ZE_AFFINITY_MASK=0 python3 $J/run_cmp_gg.py $tag "$@" 2>&1 | grep -E "rc=|PASS|FAIL|acceptance" | sed "s/^/  /"; cmpdirs $tag; pgl $tag; }
CFG6="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
CFG12="L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
CFG24="L=24 mz=F scheme=submatrix ndelay=32 nsweep=2 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288"
echo "######## V7, NEW vs OLD same seed, CMP_TOL=0 (expect bitwise identical)"
for c in "L6F $CFG6" "L12F $CFG12" "L24F $CFG24"; do set -- $c; t=$1; shift
  echo "[V7 old-vs-new $t]"; env $(var V7) bash -c "$(declare -f gg cmpdirs pgl); W=$W; C=$C; J=$J; S=$S; gg ab_$t $OLDT $NEWT $*"; done
echo "######## controls (the check must tell right from wrong)"
echo "[control 1: old vs old, V7, L12F: determinism of the run itself, must be identical]"
env $(var V7) bash -c "$(declare -f gg cmpdirs pgl); W=$W; C=$C; J=$J; S=$S; gg ctl_oo $OLDT $OLDT $CFG12"
printf '#!/bin/bash\nexec env ALFSH_DEVQR=0 ALFSH_DEVUNMQR=0 ALFSH_ASYNCBLK=0 "$@"\n' > $W/wrap_mkl.sh; chmod +x $W/wrap_mkl.sh
echo "[control 2: new V7 vs new with oneMKL QR (wrapper switches the device QR off), L12F: different rounding, must DIFFER]"
env $(var V7) WRAP=$W/wrap_mkl.sh bash -c "$(declare -f gg cmpdirs pgl); W=$W; C=$C; J=$J; S=$S; gg ctl_nm $NEWT $NEWT $CFG12"
echo "######## node config (32 ranks, 4/tile, CCS 4, OMPT=3), V7, 3 sweeps, Precision Green of all ranks (previous: mean 1.960957926459005E-014 max 6.536447432744344E-010, job 560724/g2)"
BIN=$OLDA node NODE_old V7 3; pg NODE_old
BIN=$NEWA node NODE_new V7 3; pg NODE_new
echo "gpu-side kernels used (ALFSH_TIMING absent: check the info line only)"; grep -h -E "Precision Green" $W/runs/NODE_old/info $W/runs/NODE_new/info
echo "end $(date -Is)"

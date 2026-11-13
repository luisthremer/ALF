#!/bin/bash
#SBATCH --job-name=mc_gate2
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --exclude=i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14
# kind 2 (SpinSUN, N_FL = 1, N_SUN = 2, no defect) is not reached by the Hubbard-defect cases of mc_gate: gate it here (same binaries as mc_gate).
M=$SCRATCH/mgr_corr; U=$SCRATCH/large_prompt/unified; L0=$SCRATCH/large_prompt/team; J=$U/jobs
REF=$SCRATCH/large_prompt/batched/build_ib_txtnew/Prog/ALF.out; NEW=$M/build_txt/Prog/ALF.out; CPU=$L0/build_cpu-txt/Prog/ALF.out
source $L0/jobs/env_intel_pvc.sh
echo "start $(date -Is) $(hostname)"
export W0=$M/work2; rm -rf $W0; mkdir -p $W0; export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
SB=LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
V7="$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_QRNB=16"
declare -A CASE=( [L6S]="L=6 mz=F nfl=1 nsun=2 h0=0.d0 scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12S]="L=12 mz=F nfl=1 nsun=2 h0=0.d0 scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" )
TASKS=()
add() { TASKS+=("$*"); }
for c in L6S L12S; do
  add a_$c gg $REF $NEW 0 $c
  add b_$c gg $M/wrap/off.sh $M/wrap/on.sh 1e-12 $c
  add con_$c sng $CPU $M/wrap/on.sh 1e-13 $c
  add coff_$c sng $CPU $M/wrap/off.sh 1e-13 $c
  add p_$c gg $M/wrap/off.sh $M/wrap/plant.sh 1e-12 $c
done
stream() { local s=$1 i=0 t
  for t in "${TASKS[@]}"; do
    if [ $((i % 8)) -eq $s ]; then set -- $t; local name=$1 drv=$2 cb=$3 gb=$4 tol=$5 cs=$6
      env $V7 CMP_TOL=$tol W=$W0/$name CPUBIN=$cb GPUBIN=$gb ZE_AFFINITY_MASK=$s python3 $J/run_cmp_$drv.py $name ${CASE[$cs]} > $W0/out_$name.txt 2>&1
    fi
    i=$((i+1)); done; }
for s in 0 1 2 3 4 5 6 7; do stream $s & done
wait
for t in "${TASKS[@]}"; do set -- $t; f=$W0/out_$1.txt; echo "#### $1 ($2, tol $5)"; grep -E "rc=|PASS|FAIL" $f | cut -c1-100; done
echo "---- switch announced in gpu run.log of on/plant runs:"; grep -l "ALFSH_OBSCORR=1" $W0/*/runs/*/gpu/run.log | wc -l
for n in b_L12S p_L6S; do echo "## $n"; grep -vE "rc=" $W0/out_$n.txt | head -20; done
echo "end $(date -Is)"

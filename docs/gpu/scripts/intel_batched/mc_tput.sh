#!/bin/bash -l
#SBATCH -J mc_tput
#SBATCH -A pn93la
#SBATCH -p test
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 00:55:00
#SBATCH -x i20r02c06s02
#SBATCH -o %x.%j.log
#SBATCH -e %x.%j.err
# Observables corr: ALFSH_OBSCORR 0 (OpenMP kernels) vs 1 (SYCL), SAME aottt binary, V7 switches, L=24 3 sweeps, 1 and 4 chains per tile,
# alternating off/on, 2 repetitions; L=32 2 sweeps 1 repetition.  Timer: TEAM_TIMERS (mean over ranks, as ib_tput.sh).
M=$SCRATCH/mgr_corr
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
export W=$M/tput; mkdir -p $W
BIN=$M/build_aottt/Prog/ALF.out
echo "start $(date -Is) $(hostname) $BIN"
ccs() { local c=$1 s=""; for t in 0 1 2 3 4 5 6 7; do s="$s,$t:$c"; done; echo "ZEX_NUMBER_OF_CCS=${s#,}"; }
tim2() { grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); s+=$(NF-6)/n; u+=($(NF-5)+$(NF-4))/n; b+=($(NF-3)+$(NF-2))/n; c+=$(NF-1)/n; o+=$NF/n; r++} END {printf "MEAN %-16s ranks %d  sweep %6.3f  upd %5.3f  stab %5.3f  cgr %5.3f  obs %5.3f s/sweep\n", t, r, s/r, u/r, b/r, c/r, o/r}'; }
run() { local t=$1 cpt=$2 sw=$3; shift 3; local nr=$((8*cpt)) omp=$((12/cpt))
  env $(var V7) ALFSH_QRNB=16 ALFSH_OBSCORR=$sw $(ccs $cpt) OMPT=$omp "$@" python3 $J/run_ref_sng.py $t $BIN $nr $NS 1 > $W/$t.out 2>&1
  tim2 $t; grep -h "Precision Green" $W/runs/$t/info | head -1; }
CMPJ=$SCRATCH/large_prompt/unified/jobs
NS=3
for rep in 1 2; do for cpt in 1 4; do for sw in 0 1; do run L24_c${cpt}_o${sw}_r$rep $cpt $sw; done; done; done
python3 $CMPJ/compare.py $W/runs/L24_c1_o0_r1 $W/runs/L24_c1_o1_r1 --tol 1e-12 | tail -12
NS=2
for cpt in 1 4; do for sw in 0 1; do run L32_c${cpt}_o${sw}_r1 $cpt $sw LSIZE=32; done; done
echo "end $(date -Is)"

#!/bin/bash -l
#SBATCH -J mc5_run
#SBATCH -A pn93la
#SBATCH -p test
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 00:55:00
#SBATCH -x i20r02c06s02
#SBATCH -o %x.%j.log
#SBATCH -e %x.%j.err
# intel-batched step 2, iteration 4: BLKFUSE=5 + prefetch with two block-site buffers (fix) (flush gathers before the scale kernel, scale factors on the fly).  Part 1 gates (8 tiles in parallel),
# part 2 timing (same aottt binary, V7 + ALFSH_OBSCORR=1 in all runs), part 3 counting runs.  Binaries: mc5_build.sh (login node).
M=$SCRATCH/mgr_corr5; U=$SCRATCH/large_prompt/unified; L0=$SCRATCH/large_prompt/team; JJ=$U/jobs
REF=$SCRATCH/large_prompt/batched/build_ib_txtnew/Prog/ALF.out
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
echo "start $(date -Is) $(hostname)"
NEW=$M/build_txt/Prog/ALF.out
export W0=$M/work; rm -rf $W0; mkdir -p $W0; export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
V7g="LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_QRNB=16"
declare -A CASE=( [L6F]="L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12F]="L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
                  [L6T]="L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12T]="L=12 mz=T scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72"
                  [L6S]="L=6 mz=F nfl=1 nsun=2 h0=0.d0 scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18"
                  [L12S]="L=12 mz=F nfl=1 nsun=2 h0=0.d0 scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" )
declare -A TOL=( [L6F]=1e-13 [L12F]=1e-13 [L6T]=1e-10 [L12T]=1e-10 [L6S]=1e-13 [L12S]=1e-13 )
CPU=$L0/build_cpu-txt/Prog/ALF.out
TASKS=(); add() { TASKS+=("$*"); }   # name driver cpubin gpubin tol case
CS="L12F L12T L12S L6F L6T L6S"
for c in $CS; do add a_$c gg $REF $M/wrap/off0.sh 0 $c; done
for c in $CS; do for n in 7 8; do add b${n}_$c gg $M/wrap/off.sh $M/wrap/on$n.sh 1e-12 $c; done; done
for c in $CS; do add c8_$c sng $CPU $M/wrap/on8.sh ${TOL[$c]} $c; done
for n in 8 9; do for c in L6F L6S L12F; do add p${n}_$c gg $M/wrap/off.sh $M/wrap/plant$n.sh 1e-12 $c; done; done
add p7_L6F gg $M/wrap/off.sh $M/wrap/plant7.sh 1e-12 L6F; add p6_L6F gg $M/wrap/off.sh $M/wrap/plant6.sh 1e-12 L6F
add b6_L12S gg $M/wrap/off.sh $M/wrap/on6.sh 1e-12 L12S; add b5_L12T gg $M/wrap/off.sh $M/wrap/on5.sh 1e-12 L12T
echo "gate tasks: ${#TASKS[@]}"
stream() { local s=$1 i=0 t
  for t in "${TASKS[@]}"; do
    if [ $((i % 8)) -eq $s ]; then set -- $t; local name=$1 drv=$2 cb=$3 gb=$4 tol=$5 cs=$6
      env $V7g CMP_TOL=$tol W=$W0/$name CPUBIN=$cb GPUBIN=$gb ZE_AFFINITY_MASK=$s python3 $JJ/run_cmp_$drv.py $name ${CASE[$cs]} > $W0/out_$name.txt 2>&1
    fi
    i=$((i+1)); set -- $s
  done; }
for s in 0 1 2 3 4 5 6 7; do stream $s & done
wait
for t in "${TASKS[@]}"; do set -- $t; f=$W0/out_$1.txt; echo "#### $1 ($2, tol $5)"; grep -E "PASS|FAIL|rc=[1-9]" $f | cut -c1-110; done
echo "---- runs that announced fused / prefetch (gpu run.log):"; for n in 2 4 5; do echo -n "level $n: "; grep -l "ALFSH_BLKFUSE= *$n" $W0/*/runs/*/gpu/run.log | wc -l; done; echo -n "prefetch: "; grep -l "ALFSH_BLKPF=1" $W0/*/runs/*/gpu/run.log | wc -l
echo "---- prefetch hit/miss (b8_L12F, b8_L12T):"; grep -h "block prefetch" $W0/b8_L12F/runs/*/gpu/info $W0/b8_L12T/runs/*/gpu/info
for n in b8_L12T b7_L12F p8_L6F p9_L6F; do echo "## $n"; grep -vE "rc=" $W0/out_$n.txt | head -20; done
echo "gates done $(date -Is)"

# ---------------- timing ----------------
export SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters
export W=$M/tput; mkdir -p $W
BIN=$M/build_aottt/Prog/ALF.out
ccs() { local c=$1 s=""; for t in 0 1 2 3 4 5 6 7; do s="$s,$t:$c"; done; echo "ZEX_NUMBER_OF_CCS=${s#,}"; }
tim2() { grep -h TEAM_TIMERS $W/runs/$1/run.log | awk -v t=$1 '{n=$(NF-7); s+=$(NF-6)/n; u+=($(NF-5)+$(NF-4))/n; b+=($(NF-3)+$(NF-2))/n; c+=$(NF-1)/n; o+=$NF/n; r++} END {printf "MEAN %-16s ranks %d  sweep %6.3f  upd %5.3f  stab %5.3f  cgr %5.3f  obs %5.3f s/sweep\n", t, r, s/r, u/r, b/r, c/r, o/r}'; }
# variants: tag fuse pf
vt() { case $1 in f0) echo "ALFSH_BLKFUSE=0 ALFSH_BLKPF=0";; f2) echo "ALFSH_BLKFUSE=2 ALFSH_BLKPF=0";; f2p) echo "ALFSH_BLKFUSE=2 ALFSH_BLKPF=1";; f4) echo "ALFSH_BLKFUSE=4 ALFSH_BLKPF=0";; f4p) echo "ALFSH_BLKFUSE=4 ALFSH_BLKPF=1";; f5) echo "ALFSH_BLKFUSE=5 ALFSH_BLKPF=0";; f5p) echo "ALFSH_BLKFUSE=5 ALFSH_BLKPF=1";; esac; }
run() { local t=$1 cpt=$2 v=$3; shift 3; local nr=$((8*cpt)) omp=$((12/cpt))
  env $(var V7) ALFSH_QRNB=16 ALFSH_OBSCORR=1 $(vt $v) $(ccs $cpt) OMPT=$omp "$@" python3 $J/run_ref_sng.py $t $BIN $nr $NS 1 > $W/$t.out 2>&1
  tim2 $t; grep -h "Precision Green" $W/runs/$t/info | head -1; grep -h "P3_TIMERS" $W/runs/$t/run.log | head -1 | cut -c1-200; grep -h "block prefetch" $W/runs/$t/run.log | head -1; }
NS=3
for rep in 1 2; do for v in f0 f2p f5p; do run L24_c1_${v}_r$rep 1 $v; done; done
NS=3
for v in f0 f2p f5p; do run L24_c4_${v}_r1 4 $v; done
NS=2
for v in f0 f5p; do run L32_c4_${v}_r1 4 $v LSIZE=32; done
NS=3
for v in f5p; do t=L24_count_$v
  env $(var V7) ALFSH_QRNB=16 ALFSH_OBSCORR=1 $(vt $v) ALFSB_COUNT=1 LIBOMPTARGET_PLUGIN_PROFILE=T OMPT=12 ZE_AFFINITY_MASK=0 python3 $J/run_ref_sng.py $t $BIN 1 $NS 1 > $W/$t.out 2>&1
  tim2 $t; grep -h "alfsb counts" $W/runs/$t/run.log; grep -hE "Kernel Name|^Kernel|DataAlloc|DataRead|DataWrite" $W/runs/$t/run.log | grep -v "__omp_offloading" | head -30
done
echo "end $(date -Is)"

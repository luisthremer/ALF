#!/bin/bash -l
#SBATCH -J cgalf_T1
#SBATCH -A b133ae
#SBATCH -p preempt
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module load nvhpc/25.7 2>/dev/null
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
B=~/cgalf/src_stab/Prog/ALF.out
export ALF_GPU_DEBUG=1
run() { # L nwrap nsweep tag [env...]
  d=~/cgalf/T1/r_$4; rm -rf $d; mkdir -p $d; cp ~/cgalf/E0/L$1/parameters ~/cgalf/E0/seeds $d/; cd $d
  sed -i -E "s/^(NSweep[[:space:]]*=[[:space:]]*)[0-9]+/\1$3/; s/^(Nwrap[[:space:]]*=[[:space:]]*)[0-9]+/\1$2/" parameters
  t0=$(date +%s.%N); env "${@:5}" $B > stdout 2> stderr; rc=$?; t1=$(date +%s.%N)
  pgm=$(grep "Precision Green  Mean" info | awk "{print \$NF}"); pgx=$(grep -A1 "Precision Green  Mean" info | tail -1 | awk "{print \$NF}")
  cpu=$(grep "CPU Time" info | awk "{print \$NF}"); acc=$(grep -m1 "Acceptance  " info | awk "{print \$NF}")
  echo "RESULT $4 L=$1 nwrap=$2 nsweep=$3 rc=$rc wall=$(echo "$t1 - $t0" | bc) cpu=$cpu acc=$acc pg_mean=$pgm pg_max=$pgx"
  grep -h "rd/cgalf-stab\|Precision Green" info stdout 2>/dev/null | sort -u | head -4
  [ $rc -ne 0 ] && tail -5 stderr
}
nvidia-smi -L
# smoke
run 24 10 1 smoke_cqr ALF_GPU_CHOLQR=1; grep -h "rd/cgalf|Cholesky" ~/cgalf/T1/r_smoke_cqr/info
# timing hook per stabilization (WRAPUR block of Nwrap slices + decompose; CGR)
for L in 24 32; do for nw in 5 10 20 40; do
  run $L $nw 1 tim_qr_L${L}_w$nw ALF_GPU_STABTIME=1 ALF_GPU_EXIT=1
  run $L $nw 1 tim_cqr_L${L}_w$nw ALF_GPU_STABTIME=1 ALF_GPU_EXIT=1 ALF_GPU_CHOLQR=1
  grep -h -i "stabiliz\|time per" ~/cgalf/T1/r_tim_qr_L${L}_w$nw/stdout ~/cgalf/T1/r_tim_qr_L${L}_w$nw/info 2>/dev/null | head -8
  grep -h -i "stabiliz\|time per" ~/cgalf/T1/r_tim_cqr_L${L}_w$nw/stdout ~/cgalf/T1/r_tim_cqr_L${L}_w$nw/info 2>/dev/null | head -8
done; done
# Precision Green vs Nwrap, 3 sweeps, same seeds
for L in 24 32; do for nw in 5 10 20 40; do
  run $L $nw 3 pgd_qr_L${L}_w$nw
  run $L $nw 3 pg_cqr_L${L}_w$nw ALF_GPU_CHOLQR=1
done; done
echo ALLDONE

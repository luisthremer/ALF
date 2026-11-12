#!/bin/bash
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname)"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
# Work per sweep: MKL_VERBOSE log of every oneMKL BLAS/LAPACK call (single chain L=24, V6 = sync devqr2+devunmqr so timings are per call),
# 1 vs 2 sweeps; FLOPs per routine and size from the logged dimensions. Device QR calls (not oneMKL) from ALFSH_TIMING counts.
export W=$SCRATCH/large_prompt/team/p6/work_flop
for ns in 1 2; do one F_s$ns V6 $ns MKL_VERBOSE=1 ALFSH_TIMING=1; done
for ns in 1 2; do
  echo "== sweeps $ns: verbose lines $(grep -c MKL_VERBOSE $W/runs/F_s$ns/run.log)"
  grep -h ALFSH_TIMING $W/runs/F_s$ns/run.log
  grep -h -o "MKL_VERBOSE [A-Z_0-9]*([^)]*)[^C]*" $W/runs/F_s$ns/run.log | awk '
   { split($2,a,"("); f=a[1]; sub(/\).*/,"",a[2]); n=split(a[2],p,","); t=$3;
     ms=0; if (t ~ /ms$/) {sub(/ms/,"",t); ms=t} else if (t ~ /us$/) {sub(/us/,"",t); ms=t/1000} else if (t ~ /s$/) {sub(/s/,"",t); ms=t*1000}
     fl=0; if (f ~ /GEMM$/ && n>=5) fl=8*p[3]*p[4]*p[5]; else if (f ~ /TRSM/ && n>=6) fl=(p[1]=="L"? 4*p[5]*p[5]*p[6] : 4*p[5]*p[6]*p[6]);
     key=f" "p[3]"x"p[4]"x"p[5]; c[key]++; F[key]+=fl; T[key]+=ms; tot+=fl; tt+=ms }
   END { for (k in c) printf "%-36s calls %7d  GFLOP %10.1f  ms %10.1f\n", k, c[k], F[k]/1e9, T[k]; printf "TOTAL GFLOP %.1f  logged ms %.1f\n", tot/1e9, tt }' | sort -k5 -n -r | head -25
  echo "  sample: $(grep -m2 MKL_VERBOSE $W/runs/F_s$ns/run.log | cut -c1-200)"
done
tim F_s2
echo "end $(date -Is)"

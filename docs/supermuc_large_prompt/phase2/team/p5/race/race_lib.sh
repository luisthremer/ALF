# sourced by the race jobs: runcase TAG N SEC [K=V ...]  (32 ranks, 4 per tile, CCS 4, OMP 3, as production)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
P=$SCRATCH/large_prompt/team/p5; R=$P/race; mkdir -p $R/out; cd $R
LIB=$R/build/Libraries/Modules/libalfdevqr.so
echo "node $(hostname) job $SLURM_JOB_ID start $(date -Is)"
for d in 0 1 2 3; do timeout 60 xpu-smi health -d $d 2>&1 | grep -E "Status" | sort | uniq -c | sed "s/^/xpu-smi dev $d: /"; done
runcase() {
  local tag=$1 N=$2 sec=$3; shift 3
  local out=$R/out/$tag.$SLURM_JOB_ID.out
  echo "CASE $tag N=$N sec=$sec [$*] start $(date -Is) node $(hostname)"
  ( for kv in "$@"; do export "$kv"; done
    export OMP_NUM_THREADS=3 ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 STRESS_N=$N STRESS_SEC=$sec ALFSH_TIMING=1 ALFSH_QRNB=16 \
           LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 OMP_TARGET_OFFLOAD=MANDATORY
    timeout $((sec + 420)) mpiexec -n 32 $R/t_stress > $out 2>&1; echo "mpiexec rc=$?" )
  echo "  ranks finished: $(grep -c ^STRESS_RESULT $out) of 32; iterations total: $(grep ^STRESS_RESULT $out | awk '{s+=$5} END {print s+0}'); s/iter mean: $(grep ^STRESS_RESULT $out | awk '{s+=$7;n++} END {if(n) print s/n}')"
  echo "  BAD lines: $(grep -c ^BAD $out); ranks with BAD: $(grep ^BAD $out | awk '{print $3}' | sort -u | wc -l); by what: $(grep ^BAD $out | awk '{print $6}' | sort | uniq -c | tr '\n' ' ')"
  echo "  distinct nonzero-BAD tiles (rank mod 8): $(grep ^BAD $out | awk '{print $3%8}' | sort -un | tr '\n' ' ')"
  echo "  STRESS_RESULT nbad bit/1e-8/1e-3 sums: $(grep ^STRESS_RESULT $out | awk '{a+=$10;b+=$11;c+=$12} END {print a+0,b+0,c+0}')   max maxrel LU,G1,G2: $(grep ^STRESS_RESULT $out | awk '{for(i=1;i<=3;i++) if($(NF-3+i)+0>m[i]) m[i]=$(NF-3+i)+0} END {print m[1]+0,m[2]+0,m[3]+0}')"
  echo "  distinct reference hashes over ranks (1 = all ranks bitwise equal): $(grep ^STRESS_HASH $out | awk '{print $5,$6,$7}' | sort -u | wc -l)"
  echo "  device-kernel calls (ALFSH_TIMING sums): geqrf_dev $(grep 'ALFSH_TIMING geqrf_dev' $out | awk '{s+=$4} END {print s+0}') unmqr_dev $(grep 'ALFSH_TIMING unmqr_dev' $out | awk '{s+=$4} END {print s+0}') geqrf(mkl) $(grep 'ALFSH_TIMING geqrf ' $out | awk '{s+=$4} END {print s+0}') unmqr(mkl) $(grep 'ALFSH_TIMING unmqr ' $out | awk '{s+=$4} END {print s+0}')"
  grep ^BAD $out | sort -k5,5n -s | head -12 | sed 's/^/  /'
  grep -i -E "alfsh:|error|fail|abort|killed" $out | grep -v ALFSH_TIMING | sort | uniq -c | head -6 | sed 's/^/  ! /'
  echo "CASE $tag end $(date -Is)"
}

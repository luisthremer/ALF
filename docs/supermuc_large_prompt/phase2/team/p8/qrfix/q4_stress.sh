#!/bin/bash
#SBATCH --job-name=qf_stress
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:50:00
# Delay-sensitised stress test: the existing detector p5/race/t_stress (read-only, binary + sources untouched) with libalfdevqr.so built from the
# header copies (hdr_tau / hdr_v = plants race_tau / race_v, hdr_new = real header) with -DQR2_TESTDELAY; delays through QR2_PDELAY / QR2_DELAY.
echo "start $(date -Is) node $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; Q=$L0/p8/qrfix; R=$L0/p5/race; cd $Q; mkdir -p out
AOT='-fsycl-targets=spir64_gen -Xs "-device pvc"'
for v in tau v new; do
  rm -rf lib_$v; mkdir lib_$v; cp src_head/Libraries/Modules/alf_la_sycl_devqr.cpp src_head/Libraries/Modules/alf_la_sycl_devqr.hpp hdr_$v/alf_la_sycl_qr2.hpp lib_$v/
  ( cd lib_$v && eval icpx -fsycl -O2 -std=c++17 -fPIC -shared $AOT -DQR2_TESTDELAY -I$MKLROOT/include alf_la_sycl_devqr.cpp -o libalfdevqr.so -qmkl -lmkl_sycl_blas > build.log 2>&1; echo "lib_$v rc=$?" ) &
done
wait; ls -l lib_*/libalfdevqr.so
for d in 0 1 2 3; do timeout 60 xpu-smi health -d $d 2>&1 | grep -E "Status" | sort | uniq -c | sed "s/^/xpu-smi dev $d: /"; done
runcase() {   # as p5/race/race_lib.sh, binary from p5/race, output in p8/qrfix/out
  local tag=$1 N=$2 sec=$3; shift 3
  local out=$Q/out/$tag.$SLURM_JOB_ID.out
  echo "CASE $tag N=$N sec=$sec [$*] start $(date -Is)"
  ( for kv in "$@"; do export "$kv"; done
    export OMP_NUM_THREADS=3 ZEX_NUMBER_OF_CCS=0:4,1:4,2:4,3:4,4:4,5:4,6:4,7:4 STRESS_N=$N STRESS_SEC=$sec ALFSH_TIMING=1 ALFSH_QRNB=16 \
           LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 OMP_TARGET_OFFLOAD=MANDATORY ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1
    timeout $((sec + 300)) mpiexec -n 32 $R/t_stress > $out 2>&1; echo "mpiexec rc=$?" )
  echo "  ranks finished: $(grep -c ^STRESS_RESULT $out) of 32; iterations total: $(grep ^STRESS_RESULT $out | awk '{s+=$5} END {print s+0}'); s/iter mean: $(grep ^STRESS_RESULT $out | awk '{s+=$7;n++} END {if(n) print s/n}')"
  echo "  BAD lines: $(grep -c ^BAD $out); ranks with BAD: $(grep ^BAD $out | awk '{print $3}' | sort -u | wc -l); by what: $(grep ^BAD $out | awk '{print $6}' | sort | uniq -c | tr '\n' ' ')"
  echo "  distinct reference hashes over ranks (1 = all ranks bitwise equal): $(grep ^STRESS_HASH $out | awk '{print $5,$6,$7}' | sort -u | wc -l)"
  echo "  device-kernel calls: geqrf_dev $(grep 'ALFSH_TIMING geqrf_dev' $out | awk '{s+=$4} END {print s+0}') unmqr_dev $(grep 'ALFSH_TIMING unmqr_dev' $out | awk '{s+=$5} END {print s+0}')"
  grep ^BAD $out | sort -k5,5n -s | head -3 | sed 's/^/  /'
  grep -i -E "alfsh:|error|fail|abort|killed" $out | grep -v ALFSH_TIMING | sort | uniq -c | head -4 | sed 's/^/  ! /'
  echo "CASE $tag end $(date -Is)"
}
S=${1:-150}
runcase S1_tau_N576_G1_pd2000  576 $S ALFSH_DEVQR_LIB=$Q/lib_tau/libalfdevqr.so ALFDEVQR_G=1 QR2_PDELAY=2000
runcase S2_v_N576_G1_pd2000    576 $S ALFSH_DEVQR_LIB=$Q/lib_v/libalfdevqr.so   ALFDEVQR_G=1 QR2_PDELAY=2000
runcase S3_tau_N576_Gdef_pd2000 576 $S ALFSH_DEVQR_LIB=$Q/lib_tau/libalfdevqr.so QR2_PDELAY=2000
runcase S4_tau_N576_G1_nodelay 576 $S ALFSH_DEVQR_LIB=$Q/lib_tau/libalfdevqr.so ALFDEVQR_G=1
runcase S5_real_N576_G1_pd2000_d3000 576 $S ALFSH_DEVQR_LIB=$Q/lib_new/libalfdevqr.so ALFDEVQR_G=1 QR2_PDELAY=2000 QR2_DELAY=3000
runcase S6_real_N1024_Gdef_pd2000_d3000 1024 $S ALFSH_DEVQR_LIB=$Q/lib_new/libalfdevqr.so QR2_PDELAY=2000 QR2_DELAY=3000
echo "end $(date -Is)"

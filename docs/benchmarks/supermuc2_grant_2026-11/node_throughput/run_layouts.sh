#!/bin/bash
# Node throughput test (owner OK 2026-11-27): one full Phase 2 node, all layouts of one side and one L, one after another.
# usage: [REPEAT=n RUN_TIMEOUT=s] sbatch -J nt_<side>_L<L> -p test|general -t HH:MM:SS run_layouts.sh SIDE L   (REPEAT: every layout n times)
#   SIDE = cpu  (build_cpu-pp, 410c7b04)       layouts: threads 1 2 4 6 8 -> ranks 112 56 28 18 14
#          gpu  (build_pp, 410c7b04)           layouts: walkers/tile 1 2 4 8 (threads 14 7 3 1) + 4/tile with 1 thread
#          gpub (intel-batched d1060305, production run-time switches of RUN_L36_batched)   same GPU layouts
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --nodes=1
#SBATCH --exclude=i20r02c06s02,i20r02c06s06,i20r03c06s06
SIDE=$1; L=$2
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27; J=$A/jobs/nt
source $A/jobs/env_intel_pvc.sh
export I_MPI_PIN=0 OMP_PROC_BIND=close OMP_PLACES=cores
case $SIDE in
  cpu)  EXE=$A/build_cpu-pp/Prog/ALF.out; LAYOUTS="1:1 2:1 4:1 6:1 8:1"; unset OMP_TARGET_OFFLOAD ;;   # T:-
  gpu)  EXE=$A/build_pp/Prog/ALF.out;     LAYOUTS="1:14 2:7 4:3 8:1 4:1" ;;                          # W:T
  gpub) EXE=$A/bin_batched_d1060305/pp/ALF.out; LAYOUTS="1:14 2:7 4:3 8:1 4:1"
        export LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_QRNB=16 \
               ALFSH_OBSCORR=1 ALFSH_BLKFUSE=5 ALFSH_BLKPF=1
        export LD_LIBRARY_PATH=$A/bin_batched_d1060305/pp:$LD_LIBRARY_PATH ;;
  *) echo "unknown side $SIDE"; exit 2 ;;
esac
[ -n "$LAYOUTS_OVERRIDE" ] && LAYOUTS=$LAYOUTS_OVERRIDE   # multi.sh: one layout per node
[ "$SIDE" != cpu ] && export OMP_TARGET_OFFLOAD=MANDATORY ZE_FLAT_DEVICE_HIERARCHY=FLAT ALF_UNSTABLE_THRESHOLD=1e-4
echo "job $SLURM_JOB_ID host $(hostname) side $SIDE L $L exe $EXE md5 $(md5sum < $EXE | cut -c1-32) start $(date -Is)"
SHM=/dev/shm/nt_${SLURM_JOB_ID}_${SLURM_STEP_ID:-0}; mkdir -p $SHM $A/runs
for rep in $(seq ${REPEAT:-1}); do for lay in $LAYOUTS; do
  a=${lay%:*}; b=${lay#*:}
  if [ $SIDE = cpu ]; then export MODE=cpu T=$a; NR=$(( (56 / a) * 2 )); TH=$a; TAG=t${a}; unset W
  else export MODE=gpu W=$a; NR=$(( 8 * a )); TH=$b; TAG=w${a}t${b}; unset T
       c=$(( a < 4 ? a : 4 )); export ZEX_NUMBER_OF_CCS=$(for t in 0 1 2 3 4 5 6 7; do printf '%s:%s,' $t $c; done | sed 's/,$//'); fi
  export NR OMP_NUM_THREADS=$TH MKL_NUM_THREADS=$TH
  R=${SIDE}_L${L}_${TAG}_r${rep}; D=$SHM/$R; mkdir -p $D/Temp_0; cd $D
  python3 $J/make_params.py $J/template_parameters_RUN_L24.txt $L $NR > parameters || exit 3
  cp parameters Temp_0/parameters
  python3 -c "import random; random.seed(1000+$L); print('\n'.join('%12d' % random.randint(1, 2**30) for _ in range($NR)))" > seeds
  echo "== $R ranks $NR threads $TH ${ZEX_NUMBER_OF_CCS:+CCS $ZEX_NUMBER_OF_CCS} start $(date -Is)"
  XP=""
  if [ $SIDE != cpu ]; then xpu-smi dump -d -1 -m 0,1,8,9,11 -i 1 > xpu.txt 2>&1 & XP=$!; fi
  T0=$(date +%s%N)
  timeout ${RUN_TIMEOUT:-4500} mpiexec ${NT_LOCAL:+-bootstrap fork -hosts $(hostname)} -n $NR -ppn $NR $J/pin.sh $EXE > run.log 2>&1 & P=$!
  python3 $J/watch.py $D/Temp_0/data.h5 $P > watch.txt & WP=$!
  wait $P; RC=$?; T1=$(date +%s%N); wait $WP; [ -n "$XP" ] && kill $XP
  echo "side=$SIDE L=$L tag=$TAG rep=$rep ranks=$NR threads=$TH w=${W:-} rc=$RC t0_ns=$T0 t1_ns=$T1 node=$(hostname) job=$SLURM_JOB_ID exe=$EXE" > meta.txt
  echo "   rc=$RC wall $(( (T1 - T0) / 1000000 )) ms, bin writes $(wc -l < watch.txt)"
  grep -A1 -E "Precision Green" Temp_0/info | tail -1; grep -E "Acceptance" Temp_0/info | tail -1
  rm -f Temp_0/confout_* confout_*   # configurations are not needed for timing
  cp -r $D $A/runs/${R}_$SLURM_JOB_ID; rm -rf $D
done; done
rm -rf $SHM; echo "end $(date -Is)"

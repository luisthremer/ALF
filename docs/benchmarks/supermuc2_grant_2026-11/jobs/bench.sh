#!/bin/bash
# DRAFT, NOT SUBMITTED (owner: "don't submit anything yet", 2026-11-27).
# Scaling benchmark of the space-time point defect production binary (ALF 410c7b04, PP build) on SuperMUC-NG Phase 2.
# usage: sbatch --nodes=N --partition=test|general bench.sh L RANKS_PER_TILE TILES_PER_NODE [WALKERS_PER_U]
#   production layout = bench.sh L 4 8 64   (32 ranks per node, 4 per tile, ZEX_NUMBER_OF_CCS 4)
# Parameters: jakob_params (h0 = 1, N_FL = 2, Beta = L, N_Phi = L^2/2, submatrix, N_delay 64, Nwrap 10) with
# NSweep 2, NBin 4, CPU_MAX 0. Timing: bin completion times of Temp_0/data.h5 (mtime poll), bin 1 discarded.
#SBATCH --job-name=gb_scal
#SBATCH --output=gb_scal.%j.log
#SBATCH --error=gb_scal.%j.err
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --exclude=i20r02c06s02,i20r02c06s06,i20r03c06s06
L=$1; RPT=${2:-4}; TPN=${3:-8}; WPU=${4:-64}
G=$SCRATCH/grant_bench_2026-11
EXE=$G/build_pp/Prog/ALF.out
source $G/jobs/env_intel_pvc.sh
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=3
NR=$(( SLURM_NNODES * TPN * RPT ))
[ $WPU -gt $NR ] && WPU=$NR
NU=$(( NR / WPU ))
CCS=$(for t in $(seq 0 $((TPN-1))); do printf '%s:%s,' $t $RPT; done); export ZEX_NUMBER_OF_CCS=${CCS%,}
RUN=$G/runs/L${L}_n${SLURM_NNODES}_t${TPN}_r${RPT}_${SLURM_JOB_ID}; mkdir -p $RUN; cd $RUN
python3 $G/jobs/make_params.py $L $NU $WPU > parameters   # jakob_params text, U values 5.0 + 0.05 k (k < NU)
cp $G/jobs/seeds .
# tile mask per rank: local rank / RPT (only for partial nodes; a full node uses the production mapping)
cat > mask.sh <<'EOF'
#!/bin/bash
[ -n "$PART_NODE" ] && export ZE_AFFINITY_MASK=$(( ${MPI_LOCALRANKID:-0} / RPT ))
exec "$@"
EOF
chmod +x mask.sh; [ $TPN -lt 8 ] && export PART_NODE=1 RPT
( prev=0; while sleep 1; do m=$(stat -c %Y Temp_0/data.h5 2>/dev/null || echo 0)
  [ "$m" != "$prev" ] && { echo "$(date +%s.%N) $m"; prev=$m; }; done > bins.t ) & POLL=$!
date +%s.%N > t_start
mpiexec -n $NR -ppn $(( TPN * RPT )) ./mask.sh $EXE
date +%s.%N > t_end; kill $POLL
grep -A1 'Precision Green' Temp_0/info | tail -2
echo "L=$L nodes=$SLURM_NNODES ranks=$NR U_sets=$NU walkers_per_U=$WPU"; cat bins.t

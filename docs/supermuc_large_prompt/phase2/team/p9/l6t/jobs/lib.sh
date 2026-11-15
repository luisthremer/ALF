# shared helpers (sourced) for the p9 L6T same-seed test (decision 43)
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p9/l6t; J=$P/jobs
CPUBIN=$L0/build_cpu-txt/Prog/ALF.out            # ifx CPU text build (same Prog/ physics; GPU code is ALF_GPU/ALF_TEAM_TIMERS-only)
GPUBIN=$P/build_txt/Prog/ALF.out                  # GPU text build of git archive HEAD 5924614f (job 561175)
S0=$L0/src/Scripts_and_Parameters_files/Start/seeds
SB=LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576
var() { case $1 in
  oneMKL) echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=0 ALFSH_DEVUNMQR=0";;        # = V2
  V7)     echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1";;
  P*)     echo "$SB ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1 ALFSH_DEVQR_LIB=$P/plant/$1/libalfdevqr.so";;   # V7 with a planted device-QR library
  esac; }
sysargs() { case $1 in
  L6T)  echo "L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18";;
  L6F)  echo "L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18";;
  L12F) echo "L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72";;
  esac; }
# seed set k = seeds file rotated so that it starts at line k (k=1: the original file); with 1 MPI rank only the first line is used
mkseed() { local k=$1; mkdir -p $W/seeds; (tail -n +$k $S0; head -n $((k-1)) $S0) > $W/seeds/seeds_$k; }
# concurrent single-rank runs on one node: fork bootstrap (no srun step per run), no pinning, one tile per worker
export I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 OMP_TARGET_OFFLOAD=MANDATORY
runc() { local k=$1 s=$2; W=$W CPUBIN=$CPUBIN GPUBIN=$GPUBIN SEEDS=$W/seeds/seeds_$k python3 $J/run_cmp_sng.py ${s}_${k}_cpu only=cpu $(sysargs $s) | grep -v "rc=0"; }
rung() { local k=$1 s=$2 v=$3 t=$4; env $(var $v) W=$W CPUBIN=$CPUBIN GPUBIN=$GPUBIN SEEDS=$W/seeds/seeds_$k ZE_AFFINITY_MASK=$t python3 $J/run_cmp_sng.py ${s}_${k}_$v only=gpu $(sysargs $s) | grep -v "rc=0"; }
cmpr() { local k=$1 s=$2 v=$3; python3 $J/compare.py $W/runs/${s}_${k}_cpu/cpu $W/runs/${s}_${k}_$v/gpu --tol 0 > $W/runs/${s}_${k}_$v/cmp.txt
  awk -v k=$k -v s=$s -v v=$v '/^acceptance/{a=$NF} /worst rel dev/{w=$5} /MISSING|SHAPE/{bad=1} END{if(bad||w=="")w="inf"; printf "RES %-5s k=%-3s %-7s worst %s acc%s\n", s, k, v, w, a}' $W/runs/${s}_${k}_$v/cmp.txt; }

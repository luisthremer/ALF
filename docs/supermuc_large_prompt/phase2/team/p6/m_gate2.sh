#!/bin/bash
#SBATCH --job-name=p6_gate2
#SBATCH --output=%x.%j.log
#SBATCH --error=%x.%j.err
#SBATCH --time=00:59:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=preempt
#SBATCH --qos=preempt
#SBATCH -x i20r02c06s02
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
echo "start $(date -Is) $(hostname) job $SLURM_JOB_ID"
source $SCRATCH/large_prompt/team/jobs/env_intel_pvc.sh
source $SCRATCH/large_prompt/team/p5/jobs/p5lib.sh
# Gate 2: GPU vs GPU, same binary (p6 txtnew = prototype), CMP_TOL=0 (bitwise): CPUBIN side runs with the env of this shell, GPUBIN side through WRAP (adds ALF_GPU_WB*).
# Each comparison on a single tile; L6F, L12F and L24 (1 sweep). V7 switches on both sides.
P6=$L0/p6; export W=$P6/work_gate2; BIN=$P6/build_txtnew/Prog/ALF.out
export SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds ZE_AFFINITY_MASK=0 CPUBIN=$BIN GPUBIN=$BIN CMP_TOL=0
export ALF_GPU_DEBUG=1   # ALF refuses ALF_GPU_* variables in production runs unless ALF_GPU_DEBUG=1 (allow-list only, no numerical effect)
mkw() { printf '#!/bin/bash\n%s\nexec "$@"\n' "$1" > $P6/w_$2.sh; chmod +x $P6/w_$2.sh; }
cases() { local tag=$1 wrap=$2; shift 2   # remaining args: environment of the reference side
  for a in "L6F L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18" \
           "L12F L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72" \
           "L24F L=24 mz=F scheme=submatrix ndelay=32 nsweep=1 nbin=1 U=5.3d0 symm=T beta=24.d0 nphi=288"; do
    set -- $a; local n=$1; shift
    echo "## $tag $n"
    env $(var V7) "${REFENV[@]}" WRAP=$P6/w_$wrap.sh python3 $J/run_cmp_gg.py ${tag}_$n "$@" 2>&1 | grep -E "worst rel dev|PASS|FAIL|acceptance"
  done; }
mkw 'export ALF_GPU_WB=4 ALF_GPU_WB_MODE=3' m3b4
mkw 'export ALF_GPU_WB=4 ALF_GPU_WB_MODE=1' m1b4
mkw 'export ALF_GPU_WB=2 ALF_GPU_WB_MODE=1' m1b2
mkw 'export ALF_GPU_WB=16 ALF_GPU_WB_MODE=1' m1b16
mkw 'export ALF_GPU_WB=4 ALF_GPU_WB_MODE=2' m2b4
mkw 'export ALF_GPU_WB=8 ALF_GPU_WB_MODE=2' m2b8
mkw 'export ALF_GPU_WB=4 ALF_GPU_WB_MODE=3 ALF_GPU_WB_PLANT=1' plant
mkw ':' none
echo "######## A: control, default vs default (same binary, no switches): must be bitwise"
REFENV=(); cases A_none none
echo "######## B: mode 3 (separate gemms in the batch harness, B=4) vs default: must be bitwise (tests stand-in slots / slot rotation / copies)"
REFENV=(); cases B_m3b4 m3b4
echo "######## C: planted defect (1e-12 on one element of G, B=4 mode 3) vs default: must FAIL"
REFENV=(); cases C_plant plant
echo "######## D: mode 1 (strided batch) B=4 vs mode 1 B=1 (batched mode run alone): bitwise if the strided kernel does not depend on B"
REFENV=(ALF_GPU_WB=1 ALF_GPU_WB_MODE=1); cases D_m1b4_vs_b1 m1b4
REFENV=(ALF_GPU_WB=1 ALF_GPU_WB_MODE=1); cases D_m1b16_vs_b1 m1b16
echo "######## E: mode 2 (wide gemm left) B=4 vs B=2, B=8 vs B=2"
REFENV=(ALF_GPU_WB=2 ALF_GPU_WB_MODE=2); cases E_m2b4_vs_b2 m2b4
REFENV=(ALF_GPU_WB=2 ALF_GPU_WB_MODE=2); cases E_m2b8_vs_b2 m2b8
echo "######## F: mode 1 B=4 vs default (not expected to be bitwise at N=576: oneMKL batched kernel differs); CMP_TOL=0, the worst deviation is the number of interest"
REFENV=(); cases F_m1b4_vs_def m1b4
echo "end $(date -Is)"

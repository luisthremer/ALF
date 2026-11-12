# Gate B: shim geqrf per call inside ALF, switch off (oneMKL) vs on (ALFSH_DEVQR=1), L=24 / L=32, 1 tile, 8 and 1 host threads.
# per call = (total_b - total_a)/(calls_b - calls_a) of two runs with different sweep counts (removes start-up and the JIT of the kernels).
# usage (queue script): NB=32 bash gB.sh
source $SCRATCH/large_prompt/team/p4/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p4; cd $P/jobs
export W=$P/work_gB SEEDS=$L0/ref_L24/seeds PARAMS=$L0/ref_L24/parameters OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
B=${BIN:-$P/build_aotttnew/Prog/ALF.out}
pc() { # tag-a tag-b : per call geqrf and getrf/unmqr
  for r in geqrf geqrf_dev unmqr getrf; do
    python3 - "$W/runs/$1/run.log" "$W/runs/$2/run.log" $r <<'PY'
import re,sys
def g(f,r):
    for l in open(f):
        m=re.match(r"ALFSH_TIMING %-12s calls\s+(\d+)\s+total\s+([\d.]+) s"%r if False else r"ALFSH_TIMING "+r+r"\s+calls\s+(\d+)\s+total\s+([\d.]+) s",l)
        if m: return int(m.group(1)),float(m.group(2))
    return 0,0.0
a=g(sys.argv[1],sys.argv[3]); b=g(sys.argv[2],sys.argv[3])
if b[0]>a[0]: print("   %-6s calls %d  per call %.3f ms (delta)   [run b: %d calls, %.3f s total, %.3f ms mean]"%(sys.argv[3],b[0]-a[0],1e3*(b[1]-a[1])/(b[0]-a[0]),b[0],b[1],1e3*b[1]/b[0]))
PY
  done
}
IFS=";"; CL=(${CASES:-L24 8 2 6;L24 1 2 6;L32 8 1 3;L32 1 1 3}); unset IFS   # CASES: "tag threads sweeps_a sweeps_b;..."
for c in "${CL[@]}"; do set -- $c
  for sw in ${SWS:-off on}; do
    if [ $sw = on ]; then export ALFSH_DEVQR=1 ALFSH_QRNB=${NB:-32}; else unset ALFSH_DEVQR; fi
    if [ $1 = L32 ]; then export LSIZE=32; else unset LSIZE; fi
    for ns in $3 $4; do OMPT=$2 python3 run_ref_sng.py g_$1_t$2_${sw}${NB}_s$ns $B 1 $ns 1 | tail -1; done
    echo "### $1 threads=$2 devqr=$sw nb=${NB:-32}"; pc g_$1_t$2_${sw}${NB}_s$3 g_$1_t$2_${sw}${NB}_s$4
    grep -h "Precision Green" $W/runs/g_$1_t$2_${sw}${NB}_s$4/info | head -1
  done
done
echo "done $(date -Is)"

# usage: pmc.sh OUTDIR SRCPARAMS NSWEEP BIN  -> counter passes with rocprofv3 (per-kernel CSV)
export ALF_GOMP=$HOME/amd/gccsrc/build/.libs; . ~/amd/env.sh
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 GCN_NUM_THREADS=4
export LD_LIBRARY_PATH=/opt/rocm/core-7.14/lib:$LD_LIBRARY_PATH
O=$1; P=$2; NS=$3; B=${4:-$HOME/amd/work/txt_gpu/Prog/ALF.out}
i=0
for grp in "OccupancyPercent MeanOccupancyPerCU" "VALUBusy MemUnitBusy" "FETCH_SIZE WRITE_SIZE" "SQ_WAVES GRBM_GUI_ACTIVE" "L2CacheHit MemUnitStalled"; do
  i=$((i+1)); d=$O/pass$i; rm -rf $d; mkdir -p $d; sed -e "s/^NSweep .*/NSweep = $NS/" $P/parameters > $d/parameters; cp $P/seeds $d/
  (cd $d; rocprofv3 --pmc $grp --kernel-trace -f csv -d prof -o r -- $B > run.log 2>&1; echo "pass $i [$grp] rc=$?")
done

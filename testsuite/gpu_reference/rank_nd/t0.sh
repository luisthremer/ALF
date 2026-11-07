#!/bin/bash
# Ticket 401-B, tier T0, on nividia or intel (from work/t401b/src): same-seed GPU vs CPU path of the SAME binary (ALF_GPU_DISABLE=1),
# standard / delayed / submatrix, Symm T/F, Ltau 0/1, then the scheme comparison on the GPU.
#   usage: t0.sh BIN        (BIN = nv | ifx; gpu lock taken on nividia)
W=$HOME/alf-gpu/work/${T401:-t401b}; R=$W/src/testsuite/gpu_reference/rank_nd; B=$1
if [ -d /opt/nvidia ]; then . ~/alf-gpu/env/gnu.sh; else . ~/alf-gpu/env/intel.sh; fi
. ~/alf-gpu/venv/bin/activate 2>/dev/null; export TMPDIR=$W/tmp; cd $W || exit 1
LK=0; [ -d /opt/nvidia ] && LK=1
for s in 1 0; do for t in 0 1; do for sc in standard delayed submatrix; do
  n=pv_s${s}_t${t}_${sc:0:3}
  [ -d runs/base/$n ] || python3 $R/mkrun.py src runs/base $n Symm=$([ $s = 1 ] && echo True || echo False) Ltau=$t NSweep=5 NBin=20 Sequential_scheme=\"$sc\" > /dev/null 2>&1
  GPULOCK=$LK TAG=${B}_dev $R/run.sh $B $n
  GPULOCK=0 TAG=${B}_cpu $R/run.sh $B $n ALF_GPU_DISABLE=1
  python3 $R/check_401.py same-seed runs/${B}_dev/$n runs/${B}_cpu/$n | tail -3
  grep -h "GPU path\|GPU stab\|GPU update path\|GPU tau" runs/${B}_dev/$n/info | cut -c1-110
done; done; done

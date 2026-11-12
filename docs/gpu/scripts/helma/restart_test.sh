#!/bin/bash -l
# Restart test (production = chain of links): L=8, 2 groups x 2 walkers on 1 GPU. Run 1: NBin=3. Run 2 in the same dir: NBin=6
# (CPU_MAX / NBin edited like prod_job.sh does) -> ALF must resume from confout_* and append bins to data.h5 (expect 6 bins, 3+3).
# Reference: uninterrupted NBin=6 run with the same seeds; the final confout and observables must be identical (bitwise on the same GPU type).
#SBATCH -A b133ae
#SBATCH --export=NONE
#SBATCH -J alfrestart
unset SLURM_EXPORT_ENV
WS=/hnvme/workspace/b133ae29-ALF_GPU; S=$WS/ALF/docs/gpu/scripts/helma; O=$WS/bench/restart_$SLURM_JOB_ID
PY=/hnvme/workspace/b133ae29-CLAUDE/venv_h5/bin/python
python3 $S/setup_pp.py $O/split 8 10 64 0 2 5.3,5.35 --nbin 3 --master rst > /dev/null
python3 $S/setup_pp.py $O/full  8 10 64 0 2 5.3,5.35 --nbin 6 --master rst > /dev/null
$S/run_pp.sh $O/full 2 2 25 > /dev/null 2>&1; echo "full rc=$(cat $O/full/DONE)"
$S/run_pp.sh $O/split 2 2 25 > /dev/null 2>&1; echo "split part 1 rc=$(cat $O/split/DONE)"
cp $O/split/Temp_0/info $O/split/info_part1
for f in $O/split/Temp_*/parameters; do python3 - $f <<'P'
import sys,re
t=open(sys.argv[1]).read(); open(sys.argv[1],"w").write(re.sub(r"(?m)^(NBin\s*=\s*).*", r"\g<1>6", t))
P
done
$S/run_pp.sh $O/split 2 2 25 > /dev/null 2>&1; echo "split part 2 rc=$(cat $O/split/DONE)"
grep -h "initial configuration" $O/split/Temp_0/info | head -3
$PY - $O <<'P'
import sys, h5py, numpy as np
o = sys.argv[1]
for g in (0, 1):
    a, b = h5py.File(f"{o}/full/Temp_{g}/data.h5"), h5py.File(f"{o}/split/Temp_{g}/data.h5")
    for k in ("Ener_scal/obser", "Green_eq/obser", "SpinZ_local/obser"):
        x, y = a[k][()], b[k][()]
        print("group", g, k, "bins full", x.shape[0], "split", y.shape[0], "max|diff|", np.max(np.abs(x - y)) if x.shape == y.shape else "shape differs")
P

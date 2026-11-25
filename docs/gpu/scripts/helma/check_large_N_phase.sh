#!/bin/bash -l
# Regression check for the NaN determinant phase at large N (docs/gpu/KNOWN_ISSUE_DET_PHASE_UNDERFLOW.md).
# usage: sbatch --export=ALF_BIN=/path/to/ALF.out docs/gpu/scripts/helma/check_large_N_phase.sh [L]     (default L = 56)
# Runs L x L, Nt = 64 (Beta = 6.4), 1 walker, 2 bins x 1 sweep, Ltau = 0, defect off, and fails (exit 1) if the acceptance is NaN or
# below 0.30. Expected: about 0.33 with the fix, 0.27 (L = 56 to 60) or 0.01 (L >= 64) without it. About 1 minute at L = 56.
#SBATCH -A b133ae
#SBATCH -p h100
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=32
#SBATCH -t 00:20:00
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
L=${1:-56}
# Slurm runs a copy of this file from its spool directory: take the helper scripts from the original location
if [ -n "$SLURM_JOB_ID" ]; then
  SRC=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/.*Command=\([^ ]*\).*/\1/p' | head -1)
  HERE=$(cd "$(dirname "$SRC")" && pwd)
else
  HERE=$(cd "$(dirname "$0")" && pwd)
fi
: "${ALF_BIN:?set ALF_BIN to the GPU ALF.out}"
OUT=${OUT:-$PWD/check_large_N_L${L}_$SLURM_JOB_ID}
python3 $HERE/setup_pp.py $OUT $L 1 32 0 1 5.0 --nbin 2 --master check_large_N > /dev/null || exit 2
for f in $OUT/parameters $OUT/Temp_0/parameters; do
  sed -i -E 's/^(Beta[[:space:]]*=[[:space:]]*)[^[:space:]]+/\16.4d0/; s/^(Ltau[[:space:]]*=[[:space:]]*).*/\10/; s/^(ham_h0[[:space:]]*=[[:space:]]*).*/\10.d0/; s/^(ham_U2[[:space:]]*=[[:space:]]*).*/\15.0d0/' $f
done
export ALF_BIN RUN_TIMEOUT=1200
$HERE/run_pp.sh $OUT 1 1 0 > /dev/null 2>&1
acc=$(grep -m1 "Acceptance  " $OUT/Temp_0/info | awk '{print $NF}')
echo "L=$L acceptance=$acc"
awk -v a="$acc" 'BEGIN{ if (a == "" || a+0 != a+0 || a+0 < 0.30) exit 1 }' || { echo "FAIL: acceptance missing, NaN or below 0.30"; exit 1; }
echo "OK"

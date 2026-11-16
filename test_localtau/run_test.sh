#!/bin/bash
# Usage: run_test.sh SRC_SIM_DIR TEST_DIR NBINS   (NBINS=0: copy full data.h5)
# Reads SRC only. Reference: Jakob's original script + ana.out. New: ana_hdf5.out only.
source /etc/profile >/dev/null 2>&1
module --force switch gpu-env/2025 cpu-env/2026 >/dev/null 2>&1
module load intel/2025.3.1 mkl/2024.2.2 >/dev/null 2>&1
set -e
SRC=$1; T=$2; NB=$3; NOARG=${4:-0}   # NOARG=1: run ana_hdf5.out without arguments (all observables)
ALF=/hnvme/workspace/b133ae29-stp_data/ALF_Luis
PY=/home/atuin/b133ae/b133ae29/SOFTWARE/micromamba/envs/general_python/bin/python
rm -rf "$T"; mkdir -p "$T/ref" "$T/new"
for d in ref new; do cp "$SRC/parameters" "$T/$d/"; done
if [ "$NB" = 0 ]; then
  cp "$SRC/data.h5" "$T/ref/data.h5"
else
  $PY - "$SRC/data.h5" "$T/ref/data.h5" "$NB" <<'PYEOF'
import h5py, sys
src = h5py.File(sys.argv[1], "r"); nb = int(sys.argv[3])
with h5py.File(sys.argv[2], "w") as f:
    g = f.create_group("SpinZ_localtau"); s = src["SpinZ_localtau"]
    g.attrs.update(s.attrs); s.copy("lattice", g)
    for k in ("obser", "sign"):
        g.create_dataset(k, data=s[k][:nb]); g[k].attrs.update(s[k].attrs)
PYEOF
fi
cp "$T/ref/data.h5" "$T/new/data.h5"
# --- reference: Jakob's workflow
sed "s#^sim_dir = Path(.*#sim_dir = Path('$T/ref')#" $ALF/Jakob_scripts/make_SpinZlocal_tau_orig.py > "$T/orig_patched.py"
( cd "$T/ref" && $PY "$T/orig_patched.py" > py.log && $ALF/Analysis/ana.out SpinZlocal_tau > ana.log )
# --- new: ana_hdf5 only
if [ "$NOARG" = 1 ]; then ( cd "$T/new" && $ALF/Analysis/ana_hdf5.out > ana.log ); else ( cd "$T/new" && $ALF/Analysis/ana_hdf5.out SpinZ_localtau > ana.log ); fi
# --- compare every analysis output
cd "$T/ref"
files=$(find . -type f ! -name data.h5 ! -name parameters ! -name 'SpinZlocal_tau' ! -name 'SpinZlocal_tau_info' ! -name '*.log' | sort)
nf=0; bad=0
for f in $files; do nf=$((nf+1)); cmp -s "$f" "../new/$f" || { echo "DIFF: $f"; bad=$((bad+1)); }; done
extra=$(cd ../new && find . -type f -path '*SpinZlocal*' | sort | comm -13 <(echo "$files") -)
[ -n "$extra" ] && echo "EXTRA in new: $extra"
echo "$(basename $T): compared $nf files, $bad differ  $( [ $bad = 0 ] && [ -z "$extra" ] && [ $nf -gt 0 ] && echo BIT-IDENTICAL || echo FAILED)"
rm -f "$T/ref/SpinZlocal_tau"   # large text intermediate

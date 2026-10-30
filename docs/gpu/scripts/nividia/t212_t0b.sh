#!/bin/bash
# Ticket 212 (development only), nividia: determinism T0 - the same binary twice, bit-identical data.h5 and Precision lines.
# -> dev/logs/t0b.txt ; runs in dev/runs/t0b_<bin>_<cfg>_<a|b>/
W=$HOME/alf-gpu/work/t212/dev; out=$W/logs/t0b.txt; : > $out
. ~/alf-gpu/env/nvhpc.sh
for pair in "nv_def hub6x6_cb_fd" "nv_log hub6x6_cb_fd" "nv_def hub12x12_cb_b10_p" "nv_log hub12x12_cb_b10_p"; do
  set -- $pair
  for t in a b; do echo "$(TAG=t0b_$1_$2_$t bash $W/t212_run.sh $1 $2)" >> $out; done
  A=$W/runs/t0b_$1_$2_a/$2; B=$W/runs/t0b_$1_$2_b/$2
  python3 -I - $A $B >> $out <<'PY'
import sys, h5py, numpy as np, re
a, b = sys.argv[1], sys.argv[2]
def ds(p):
    o = {}
    with h5py.File(p + "/data.h5", "r") as f:
        f.visititems(lambda n, x: o.__setitem__(n, np.asarray(x[()])) if isinstance(x, h5py.Dataset) else None)
    return o
da, db = ds(a), ds(b)
keys = set(da) | set(db)      # symmetric: a dataset missing on either side counts as differing
nd = sum(1 for k in keys if k not in da or k not in db or da[k].tobytes() != db[k].tobytes())
print("  datasets", len(keys), "differing", nd)
pa = [l for l in open(a + "/info") if "Precision" in l]; pb = [l for l in open(b + "/info") if "Precision" in l]
print("  Precision lines", len(pa), "identical" if pa == pb else "DIFFER")
for l in pa[:3]: print("   ", l.strip()[:150])
PY
done
echo done >> $out

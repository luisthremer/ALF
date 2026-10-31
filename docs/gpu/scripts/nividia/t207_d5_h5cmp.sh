#!/bin/bash
# Ticket 207 dev5 (copy of t207_h5cmp.sh), nividia: bit-for-bit comparison of all datasets of two data.h5 files.  usage: t207_d4_h5cmp.sh A.h5 B.h5
. ~/alf-gpu/env/gnu.sh
python3 -I - "$1" "$2" <<'PY'
import sys, h5py, numpy as np
a, b = h5py.File(sys.argv[1], "r"), h5py.File(sys.argv[2], "r")
names = []
a.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) else None)
nb = []
b.visititems(lambda n, o: nb.append(n) if isinstance(o, h5py.Dataset) else None)
bad = [n for n in names if n not in nb or a[n].shape != b[n].shape or a[n][()].tobytes() != b[n][()].tobytes()]
print("datasets A:", len(names), "B:", len(nb), "differing:", len(bad), bad[:5])
sys.exit(1 if bad or sorted(names) != sorted(nb) else 0)
PY

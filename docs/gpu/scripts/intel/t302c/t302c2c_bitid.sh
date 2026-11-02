#!/bin/bash
# Ticket 302c-2 round 3 (development only), on intel: bit-identity of the round-3 reruns (d3e82a36 binaries) of round 2's main runs with
# round 2's (08559f0b binaries, t302c-2b/runs, read only), same seeds and env.  -> stdout (the queue writes logs/bitid.txt)
#   per run: seeds / parameters / env (stamp.txt) equal; qref_manifest.txt (sha256 of every raw dump of the run, written before
#   processing) byte-identical; qref.txt identical after masking the program's wall-time field 'tq' (the gate compares qref.txt lines
#   'timing tq aside' too), and the raw cmp; data.h5 datasets (h5py, every dataset bit for bit); info without the timing lines.
#   usage: t302c2c_bitid.sh TAG...      e.g. t302c2c_bitid.sh m4_def_dev5 m10_log_dev5
W=$HOME/alf-gpu/work/t302c-2c; O=$HOME/alf-gpu/work/t302c-2b
for t in "$@"; do
  case "$t" in m4_*) c=hub12x12_cb_p;; m10_*) c=hub12x12_cb_b10_p;; *) echo "bad tag $t"; exit 1;; esac
  a=$W/runs/$t/$c; b=$O/runs/$t/$c
  echo "== $t: new $a vs round 2 $b"
  for f in seeds parameters; do echo "  $f: $(cmp -s $a/$f $b/$f && echo identical || echo DIFFERENT)"; done
  echo "  env: $(diff <(grep '^env ' $a/stamp.txt) <(grep '^env ' $b/stamp.txt) > /dev/null && echo identical || echo DIFFERENT) ($(grep '^env ' $a/stamp.txt))"
  echo "  binaries: new $(grep '^binary' $a/stamp.txt | cut -c8-23) round 2 $(grep '^binary' $b/stamp.txt | cut -c8-23)"
  echo "  qref_manifest.txt (sha256 of all $(wc -l < $a/qref_manifest.txt) raw dumps): $(cmp -s $a/qref_manifest.txt $b/qref_manifest.txt && echo identical || echo "DIFFERENT ($(diff $a/qref_manifest.txt $b/qref_manifest.txt | grep -c '^<') lines differ)")"
  echo "  qref.txt raw cmp: $(cmp -s $a/qref.txt $b/qref.txt && echo identical || echo differs)"
  m() { sed -E 's/ tq +[0-9.E+-]+ / tq X /' "$1"; }
  n=$(diff <(m $a/qref.txt) <(m $b/qref.txt) | grep -c '^<')
  echo "  qref.txt with tq masked ($(wc -l < $a/qref.txt) lines): $([ "$n" = 0 ] && echo IDENTICAL || echo "DIFFERENT in $n lines")"
  ~/alf-gpu/venv/bin/python3 -I - $a/data.h5 $b/data.h5 <<'PY'
import sys, h5py, numpy as np
A, B = (h5py.File(f, 'r') for f in sys.argv[1:3])
na, nb = [], []
A.visititems(lambda k, v: na.append(k) if isinstance(v, h5py.Dataset) else None)
B.visititems(lambda k, v: nb.append(k) if isinstance(v, h5py.Dataset) else None)
diff = [k for k in na if k in nb and not (A[k].shape == B[k].shape and np.array_equal(np.asarray(A[k][()]).view(np.uint8), np.asarray(B[k][()]).view(np.uint8)))]
print('  data.h5: %d datasets (%d in round 2), names %s, %s' % (len(na), len(nb), 'equal' if sorted(na) == sorted(nb) else 'DIFFERENT',
      'all bit-identical' if not diff and sorted(na) == sorted(nb) else 'DIFFERENT: %s' % ', '.join(diff[:8])))
PY
  echo "  info (lines without 'time'/'Time'/'CPU'/'Wall'): $(diff <(grep -viE 'time|cpu|wall' $a/info) <(grep -viE 'time|cpu|wall' $b/info) > /dev/null && echo identical || echo DIFFERENT)"
done

import os, struct, sys
import numpy as np
a, b, name = sys.argv[1:4]
fa, fb = sorted(os.listdir(a)), sorted(os.listdir(b))
assert fa == fb and fa, 'dump lists differ or are empty'
st = {}
for f in fa:
    A, B = open(os.path.join(a, f), 'rb').read(), open(os.path.join(b, f), 'rb').read()
    assert len(A) == len(B), f
    h = struct.unpack('<13i', A[:52])
    assert h[0] == 20791 and h[9] == 0 and h[10] == 0 and h[11] == 0 and h == struct.unpack('<13i', B[:52]), (f, h)
    nst, nfe, n, hasb = h[4], h[5], h[7], h[12]
    lo = len(A) - 16 * n * n * (4 + 2 * hasb); hi = lo + 64 * n * n          # the four production outputs G00, GT0, GTT, G0T
    same_in = A[:lo] == B[:lo] and A[hi:] == B[hi:]                          # header, input states, B and Bi
    oa = np.frombuffer(A[lo:hi], dtype='<c16'); ob = np.frombuffer(B[lo:hi], dtype='<c16')
    rel = float(np.max(np.abs(ob - oa)) / np.max(np.abs(oa)))
    s = st.setdefault(nfe, dict(n=0, inputs_identical=0, outputs_identical=0, max_rel=0.0, min_rel=float('inf'), nst=set()))
    s['n'] += 1; s['inputs_identical'] += same_in; s['outputs_identical'] += (A[lo:hi] == B[lo:hi]); s['nst'].add(nst)
    s['max_rel'] = max(s['max_rel'], rel); s['min_rel'] = min(s['min_rel'], rel)
for nfe, s in sorted(st.items()):
    print('  QREF %s AT=0:2 vs no fault, flavour %d: dumps %d (nst %s), inputs+B identical %d, outputs identical %d, output max rel diff %.3e min %.3e'
          % (name, nfe, s['n'], ','.join(map(str, sorted(s['nst']))), s['inputs_identical'], s['outputs_identical'], s['max_rel'], s['min_rel']))
ok = st[1]['outputs_identical'] == st[1]['n'] and st[2]['outputs_identical'] == 0 and all(s['inputs_identical'] == s['n'] for s in st.values())
print('  QREF %s flavour-only check: %s' % (name, 'PASSED (flavour 1 bit-identical, every flavour-2 output perturbed)' if ok else 'FAILED'))

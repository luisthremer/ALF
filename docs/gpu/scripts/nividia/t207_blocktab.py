#!/usr/bin/env python3
"""Ticket 207: per stabilization block table of the tau-pass self-test (run with ALF_GPU_TAUTEST_V=1).

usage: t207_blocktab.py RUN.log [NWRAP]
Reads the lines "GPU tau self-test slice N (...): G00, GT0, GTT, G0T a b c d" (device vs host, one line per time slice, taken before
the CGR2_2 of that slice) and prints, for every block of NWRAP slices (a block ends at a stabilization point), the maximum of each of
the four matrices over the block and the slice where the block maximum is reached.
"""
import re
import sys

log = open(sys.argv[1], encoding="utf-8", errors="replace").read()
nwrap = int(sys.argv[2]) if len(sys.argv) > 2 else 10
rows = {}
for m in re.finditer(r"self-test slice\s+(\d+) \(before a CGR2_2 at this slice\): G00, GT0, GTT, G0T\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", log):
    rows[int(m.group(1))] = [float(m.group(k)) for k in range(2, 6)]
print("block (slices)      max G00     max GT0     max GTT     max G0T    (slice of the max over the 4)")
worst = 0.0
for b in range(1, (max(rows) + nwrap - 1) // nwrap + 1):
    sl = [s for s in range(1 + (b - 1) * nwrap, b * nwrap + 1) if s in rows]
    mx = [max(rows[s][k] for s in sl) for k in range(4)]
    best = max(sl, key=lambda s: max(rows[s]))
    worst = max(worst, max(mx))
    print(f"{b:3d} ({sl[0]:3d}-{sl[-1]:3d})   " + "  ".join(f"{x:10.3e}" for x in mx) + f"    {best}")
print(f"overall max {worst:.3e}")

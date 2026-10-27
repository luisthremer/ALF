#!/usr/bin/env python3
"""Ticket 108: print prep specs name:nsweep:nbin for a case list, sized to ~TARGET seconds per run
(sweep time model measured on nividia: 0.13 s * (beta/10) * (L/6)^3.4 for the checkerboard Hubbard model, 1 core).
usage: mkcases108.py TARGET_S NBIN SEEDS "L:beta:U:Nwrap[:dense] ..." """
import sys
target, nbin, seeds, cases = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4].split()
for c in cases:
    p = c.split(":"); L, b, U, w = int(p[0]), float(p[1]), int(p[2]), int(p[3]); dense = len(p) > 4
    t = 0.13 * (b / 10) * (L / 6) ** 3.4 * (1.5 if dense else 1.0) * (1.0 + 2.0 / w)/1.2
    ns = max(2, int(target / t / nbin))
    for s in range(1, seeds + 1):
        print(f"L{L}_b{b:g}_U{U}_w{w}{'_dense' if dense else ''}_s{s}:{ns}:{nbin}")

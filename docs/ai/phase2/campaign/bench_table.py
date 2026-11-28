#!/usr/bin/env python3
"""Fastest scheme per (machine, model, L) from bench/bench_*.tsv (s per sweep = (t6 - t2)/4; negative or missing = unreliable, skipped)."""
import csv, glob, collections
d = collections.OrderedDict()
for f in sorted(glob.glob("bench/bench_*.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        try: v = float(r["s_per_sweep"])
        except ValueError: continue
        if v <= 0: continue
        d.setdefault((r["machine"].replace("-refused", ""), r["model"], int(r["L"])), {})[r["scheme"]] = v
print("| machine | model | L | std | dly16 | sub16 | fastest |\n|---|---|---|---|---|---|---|")
for (m, mo, L), s in sorted(d.items()):
    best = min(s, key=s.get)
    print("| %s | %s | %d | %s | %s | %s | %s |" % (m, mo, L, *["%.3f" % s[k] if k in s else "n/a" for k in ("std", "dly16", "sub16")], best))

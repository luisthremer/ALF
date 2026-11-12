#!/usr/bin/env python3
"""Stream a LIBOMPTARGET_INFO=-1 log (rank0.log) and print, per OpenMP kernel: launches, and the distinct (number of teams, team size) with counts.
usage: parse_info.py LOG   ->  lines  KINFO <launches> <kernel> | <teams> x <teamsize> (<count>) ..."""
import re, sys, collections
rk = re.compile(r"Launching kernel (\S+) with")
rt = re.compile(r"Team sizes = \{(\d+), (\d+), (\d+)\}")
rn = re.compile(r"Number of teams = \{(\d+), (\d+), (\d+)\}")
cnt = collections.defaultdict(collections.Counter); tot = collections.Counter()
name = None; ts = None; nt = None
def fin():
    global ts, nt
    if name and ts and nt:
        cnt[name][(nt, ts)] += 1; tot[name] += 1
    ts = nt = None
for l in open(sys.argv[1], errors="replace"):
    if "Launching kernel" in l:
        fin(); m = rk.search(l); name = m.group(1) if m else None
    elif "Team sizes" in l:
        m = rt.search(l); ts = int(m.group(1)) * int(m.group(2)) * int(m.group(3)) if m else None
    elif "Number of teams" in l:
        m = rn.search(l); nt = int(m.group(1)) * int(m.group(2)) * int(m.group(3)) if m else None
fin()
for k, v in sorted(tot.items(), key=lambda kv: -kv[1]):
    print("KINFO %7d %s | %s" % (v, k, "  ".join("%dx%d(%d)" % (a, b, c) for (a, b), c in cnt[k].most_common(6))))

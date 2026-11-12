#!/usr/bin/env python3
"""P7 analysis helpers: parse LIBOMPTARGET_PLUGIN_PROFILE tables, ALFP7_BLK/ALFP7_MKL lines, TEAM_TIMERS; per-sweep values = (run3 - run1)/2."""
import re, sys, glob, os, collections
def short(n): 
    n = re.sub(r"^__omp_offloading_[0-9a-f]+_[0-9a-f]+_", "", n); return re.sub(r"_mod_mp_", ".", n)
def parse_pp(path):
    names = {}; rows = {}
    for l in open(path, errors="replace"):
        m = re.match(r"Kernel (\d+)\s+: (\S+)", l)
        if m: names["Kernel " + m.group(1)] = short(m.group(2)); continue
        m = re.match(r"(Kernel \d+|DataAlloc|DataRead \(Device to Host\)|DataWrite \(Host to Device\))\s*:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", l)
        if m:
            k = m.group(1); k = names.get(k, k)
            rows[k] = (float(m.group(10)), float(m.group(6)), float(m.group(2)))   # count, device ms total, host ms total
    return rows
def diff(r3, r1, ns=2.0):
    out = {}
    for k in set(r3) | set(r1):
        a = r3.get(k, (0, 0, 0)); b = r1.get(k, (0, 0, 0))
        out[k] = tuple((x - y) / ns for x, y in zip(a, b))
    return out
def parse_tag(path, pat):
    d = collections.OrderedDict()
    for l in open(path, errors="replace"):
        m = re.match(pat, l)
        if m: d[m.group(1).strip()] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return d
def parse_blk(path): return parse_tag(path, r"ALFP7_BLK (.+?)\s+calls\s+(\d+)\s+dev_s\s+([\d.]+)\s+host_s\s+([\d.]+)")
def parse_mkl(path): return parse_tag(path, r"ALFP7_MKL (.+?)\s+calls\s+(\d+)\s+dev_s\s+([\d.]+)\s+host_s\s+([\d.]+)")
def phase(k, cnt):
    if k.startswith("gpu_wrap"): return "U"
    if k.startswith("gpu_obs"): return "O"
    c = round(cnt)
    if c == 144: return "S"
    if c in (98, 96, 52, 46): return "C"
    if c == 242: return "SC"
    return "SC"
def tm(path):
    for l in open(path, errors="replace"):
        if "TEAM_TIMERS" in l:
            v = l.split("]:")[1].split(); n = float(v[0]); t = [float(x) / n for x in v[1:]]
            return dict(n=n, sweep=t[0], upd=t[1] + t[2], stab=t[3] + t[4], cgr=t[5], obs=t[6])
if __name__ == "__main__":
    pass

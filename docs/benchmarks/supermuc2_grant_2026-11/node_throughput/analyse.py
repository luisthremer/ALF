#!/usr/bin/env python3
"""Time per sweep and node throughput from the run directories of run_layouts.sh.

usage: analyse.py RUNS_DIR [RUNS_DIR ...] > results.csv

Per run: watch.txt holds every mtime change of Temp_0/data.h5. The file is created at start-up, then each bin
(NSweep = 1, NBin = 11) is written as a burst of several HDF5 writes (one per observable, ~25 ms apart, measured in
smoke job 566791). So the events form NBIN + 1 = 12 bursts, separated by the 11 largest gaps.
For bins k = 2..11 (bin 1 = first sweep incl. warm-up, dropped):
  t_comp(k) = first write of bin k - last write of bin k-1   (sweep + measurements + bin reduction of all walkers)
  t_io(k)   = last write of bin k  - first write of bin k     (writing one bin)
Production writes one bin per 10 sweeps, so its cost per sweep is t_comp + t_io / 10 (column t_sweep_s).
Walkers synchronise at every bin, so the times are those of the slowest walker.
Throughput per node = ranks * 3600 / t_sweep  [sweeps per node-hour], error from the standard error of the mean.
GPU runs: mean tile utilisation and power from xpu.txt over the timed window.
"""
import csv
import glob
import math
import os
import re
import sys

NBIN = 11
NSWEEP_PROD = 10
COLS = ["side", "L", "tag", "rep", "ranks", "threads", "w", "rc", "n_sweeps", "t_sweep_s", "t_sweep_err_s", "t_comp_s", "t_io_bin_s", "writes_per_bin",
        "throughput_node_h", "throughput_err", "gpu_util_pct", "gpu_power_w", "precision_green_max", "acceptance", "job", "dir"]


def meta(d):
    m = {}
    for tok in open(os.path.join(d, "meta.txt")).read().split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            m[k] = v
    return m


def bins(d):
    """Return [(first, last, n_events)] per burst; burst 0 is the file creation."""
    ev = sorted({int(line.split()[0]) for line in open(os.path.join(d, "watch.txt")) if len(line.split()) == 3})
    if len(ev) < NBIN + 1:
        raise ValueError(f"only {len(ev)} write events, need >= {NBIN + 1}")
    gaps = sorted(range(1, len(ev)), key=lambda i: ev[i] - ev[i - 1])[-NBIN:]   # indices that start a new burst
    cuts = [0] + sorted(gaps) + [len(ev)]
    return [(ev[a], ev[b - 1], b - a) for a, b in zip(cuts, cuts[1:])]


def stats(x):
    m = sum(x) / len(x)
    return m, math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1) / len(x))


def info_value(d, key):
    try:
        txt = open(os.path.join(d, "Temp_0", "info")).read()
    except OSError:
        return ""
    if key == "green":   # "Precision Green  Mean, Max :" -> the Max value wraps to the next line
        m = re.findall(r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s*\n?\s*(\S+)", txt)
        return m[-1][1] if m else ""
    m = re.findall(r"Acceptance\s*:\s*(\S+)", txt)
    return m[-1] if m else ""


def xpu(d, t0, t1):
    """Mean GPU utilisation and power over all devices; xpu-smi dump columns: Timestamp, DeviceId, util, power, ..."""
    f = os.path.join(d, "xpu.txt")
    if not os.path.exists(f):
        return "", ""
    u, p = [], []
    for line in open(f):
        c = [x.strip() for x in line.split(",")]
        if len(c) < 4 or not c[1].isdigit():
            continue
        try:
            u.append(float(c[2])); p.append(float(c[3]))
        except ValueError:
            pass
    if not u:
        return "", ""
    n = len(set(int(line.split(",")[1]) for line in open(f) if line.split(",")[1:2] and line.split(",")[1].strip().isdigit()))
    return f"{sum(u)/len(u):.1f}", f"{sum(p)/len(p)*max(n,1):.0f}"   # power: summed over the 4 cards


def row(d):
    m = meta(d)
    b = bins(d)[1:]                      # drop the file-creation burst -> NBIN bins
    comp = [(b[k][0] - b[k - 1][1]) / 1e9 for k in range(1, len(b))]
    io = [(b[k][1] - b[k][0]) / 1e9 for k in range(1, len(b))]
    sw = [c + i / NSWEEP_PROD for c, i in zip(comp, io)]
    r = {k: "" for k in COLS}
    r.update(side=m["side"], L=m["L"], tag=m["tag"], rep=m.get("rep", "1"), ranks=m["ranks"], threads=m["threads"], w=m.get("w", ""), rc=m["rc"],
             job=m["job"], dir=os.path.basename(d), n_sweeps=len(sw))
    wpb = sorted({n for _, _, n in b[1:]})
    r["writes_per_bin"] = "/".join(map(str, wpb))     # one value expected; several = burst split or merged, check the run
    mean, sem = stats(sw)
    nr = int(m["ranks"])
    thr = nr * 3600 / mean
    r.update(t_sweep_s=f"{mean:.4f}", t_sweep_err_s=f"{sem:.4f}", t_comp_s=f"{stats(comp)[0]:.4f}",
             t_io_bin_s=f"{stats(io)[0]:.4f}", throughput_node_h=f"{thr:.1f}", throughput_err=f"{thr * sem / mean:.1f}")
    r["gpu_util_pct"], r["gpu_power_w"] = xpu(d, int(m["t0_ns"]), int(m["t1_ns"]))
    r["precision_green_max"] = info_value(d, "green")
    r["acceptance"] = info_value(d, "acc")
    return r


def main():
    w = csv.DictWriter(sys.stdout, fieldnames=COLS)
    w.writeheader()
    for top in sys.argv[1:]:
        for d in sorted(glob.glob(os.path.join(top, "*_L*_*"))):
            if os.path.exists(os.path.join(d, "meta.txt")):
                try:
                    w.writerow(row(d))
                except Exception as exc:   # report and go on; one broken run must not hide the others
                    print(f"# {d}: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()

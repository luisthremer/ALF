#!/usr/bin/env python3
"""FAKE results in the format of analyse.py, only to develop the plots before the real runs are in.

usage: make_fake_data.py > results_fake.csv
The model numbers below are invented (rough shapes only). Never use this file for the proposal.
"""
import csv
import random
import sys

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from analyse import COLS  # noqa: E402

random.seed(7)
w = csv.DictWriter(sys.stdout, fieldnames=COLS)
w.writeheader()


def emit(side, L, tag, ranks, threads, wpt, t, util="", power=""):
    gaps = [t * random.gauss(1, 0.02) for _ in range(10)]
    mean = sum(gaps) / 10
    sem = (sum((g - mean) ** 2 for g in gaps) / 9 / 10) ** 0.5
    thr = ranks * 3600 / mean
    w.writerow(dict(side=side, L=L, tag=tag, ranks=ranks, threads=threads, w=wpt, rc=0, n_sweeps=10,
                    t_sweep_s=f"{mean:.4f}", t_sweep_err_s=f"{sem:.4f}", t_comp_s=f"{mean:.4f}", t_io_bin_s="0.1", writes_per_bin="9", throughput_node_h=f"{thr:.1f}",
                    throughput_err=f"{thr * sem / mean:.1f}", gpu_util_pct=util, gpu_power_w=power,
                    precision_green_max="1e-12", acceptance="0.33", job="FAKE", dir="FAKE"))


for L in (16, 24, 32, 36):
    s = (L / 16) ** 3
    for t in (1, 2, 4, 6, 8):
        ranks = (56 // t) * 2
        eff = 0.85 ** (t.bit_length() - 1)          # invented thread efficiency
        emit("cpu", L, f"t{t}", ranks, t, "", 1.2 * s / (t * eff))
    for side, f in (("gpu", 1.0), ("gpub", 0.72)):
        for wpt, th in ((1, 14), (2, 7), (4, 3), (8, 1), (4, 1)):
            base = 0.30 * (L / 16) ** 2.4 * f
            slow = wpt ** 0.55 * (1.08 if (wpt == 4 and th == 1) else 1.0)
            util = min(98, 35 + 15 * wpt ** 0.8)
            emit(side, L, f"w{wpt}t{th}", 8 * wpt, th, wpt, base * slow, f"{util:.0f}", f"{4 * (250 + 3.2 * util):.0f}")

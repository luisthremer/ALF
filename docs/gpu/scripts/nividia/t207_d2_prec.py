#!/usr/bin/env python3
"""Ticket 207 (development only): Precision Green / tau mean (criterion: GPU mean <= 1.5 x the larger CPU mean) and max
(information) of T0 runs.  usage: t207_d2_prec.py RUNS_DIR   (expects the t0c run directories of t207_d2_t0c.sh RUN 1 / cpu2)"""
import sys, os
R = sys.argv[1]
def pr(d):
    t = open(os.path.join(R, d, "info")).read().replace("\n", " ")
    o = {}
    for k, lab in (("green", "Precision Green  Mean, Max :"), ("tau", "Precision tau    Mean, Max :")):
        i = t.index(lab); v = t[i + len(lab):].split()[:2]; o[k] = (float(v[0]), float(v[1]))
    return o
b4, b10 = "hub12x12_cb_p", "hub12x12_cb_b10_p"
cases = [("t0c1_pr_nv_def_b4/" + b4, ["t0c2_pr_nv_none_b4/" + b4, "t0c1_pr_nv_none_pp_b4/" + b4]),
         ("t0c1_pr_nv_def_b10/" + b10, ["t0c2_pr_nv_none_b10/" + b10, "t0c1_pr_nv_none_pp_b10/" + b10]),
         ("t0c1_pr_nv_log_b4/" + b4, ["t0c1_pr_nv_none_log_b4/" + b4]),
         ("t0c1_pr_nv_log_b10/" + b10, ["t0c1_pr_nv_none_log_b10/" + b10]),
         ("t0c1_ctl_pr_ns0_def_b4/" + b4, ["t0c2_pr_nv_none_b4/" + b4, "t0c1_pr_nv_none_pp_b4/" + b4])]
for g, cs in cases:
    G = pr(g); C = [pr(c) for c in cs]
    for k in ("green", "tau"):
        cm = max(c[k][0] for c in C); cx = max(c[k][1] for c in C); r = G[k][0] / cm
        print("%-24s %-5s GPU mean %.4e max %.3e | CPU means %s maxes %s | mean ratio %.3f %s (max ratio %.2f, information)" % (
            g.split("/")[0], k, G[k][0], G[k][1], ["%.4e" % c[k][0] for c in C], ["%.3e" % c[k][1] for c in C], r,
            "PASS" if r <= 1.5 else "FAIL", G[k][1] / cx))

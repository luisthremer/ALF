#!/usr/bin/env python3
"""Rough node-hour / core-hour estimates for the proposed Phase 8 campaigns (planning only, not a benchmark result).

Usage: python3 -I testsuite/lowrank/campaign/estimate_cost.py   (prints a Markdown table)

Model (stated assumptions; uncertainty about a factor 2 either way; a pilot campaign C0 calibrates it):
* Mac reference (project 1, 1 thread, Mz Hubbard Beta = 2, Ltrot = 20): standard 2.56 s/sweep at N = 400.
* Single-thread time per sweep  t1 = c_scheme * 2.56 s * (N/400)^3 * (Ltrot/20) * F_CPU * (N_FL_eff/2), with
  c_standard = 1, c_delayed = c_submatrix = 0.45 (Mac N = 400 ratio; kept although the gain should grow with N),
  F_CPU = 1.4 (one Skylake 8174 core at the fixed 1.9 GHz vs one M1 Pro core; for Zen5c at ~2.3 GHz 1.2).
* T threads per rank (threaded MKL only): speedup T**0.5 (conservative: BLAS-2 parts barely scale).
* A job runs all ranks of one node concurrently (independent chains); wall = setup + n_sweeps * t1 / T**0.5, with
  setup = 1 sweep-equivalent (initial stabilisation) and n_sweeps = 2 (warm-up bin) + 4 (timed); one job per
  (point, layout, repetition) running the schemes one after the other.
* 3D cubic: same model with N = L^3; memory per rank ~ 10 N^2 * 16 B (GR, UDV stacks, work) limits ranks per node.
"""

F_CPU = {"supermuc": 1.4, "helma": 1.2}
CORES = {"supermuc": 48, "helma": 384}
C = {"standard": 1.0, "delayed": 0.45, "submatrix": 0.45}


def t1(N, ltrot, scheme, machine, nfl=2):
    return C[scheme] * 2.56 * (N / 400.0) ** 3 * (ltrot / 20.0) * F_CPU[machine] * (nfl / 2.0)


def job_hours(N, ltrot, schemes, threads, machine, nsweeps=6):
    wall = 0.0
    for s in schemes:
        t = t1(N, ltrot, s, machine) / threads ** 0.5
        wall += t * (1 + nsweeps)
    return wall / 3600.0


def campaign(name, machine, points, layouts, schemes, reps, nd_values=(0,)):
    nh = 0.0
    njobs = 0
    for (N, ltrot) in points:
        for (ranks, threads) in layouts:
            for nd in nd_values:
                nh += reps * job_hours(N, ltrot, schemes, threads, machine)
                njobs += reps
    return dict(name=name, machine=machine, jobs=njobs, node_h=nh, core_h=nh * CORES[machine])


def main():
    sq = lambda Ls, ltrot: [(L * L, ltrot) for L in Ls]
    cu = lambda Ls, ltrot: [(L ** 3, ltrot) for L in Ls]
    S = ("standard", "delayed", "submatrix")
    camps = [
        campaign("C1 SuperMUC 2D Hubbard U=1 and U=4, Beta=2 (Ltrot 20), L=20..60, 3 schemes, auto N_delay, layouts 1x48 2x24 4x12",
                 "supermuc", sq((20, 30, 40, 50, 60), 20) * 2, [(1, 48), (2, 24), (4, 12)], S, 3),
        campaign("C1b SuperMUC 2D Hubbard U=4, Beta=4 (Ltrot 40), L=40, 60, layout 2x24",
                 "supermuc", sq((40, 60), 40), [(2, 24)], S, 3),
        campaign("C2 SuperMUC N_delay scan: L=40, 60, delayed+submatrix, fixed 16..256 (5) + auto + tune, 2x24",
                 "supermuc", sq((40, 60), 20), [(2, 24)], ("delayed", "submatrix"), 3, nd_values=range(7)),
        campaign("C3 SuperMUC 3D cubic Hubbard U=4, Beta=2, L=8, 12, 16, 20 (N up to 8000), 3 schemes, 1x48 (2x24 for L<=16)",
                 "supermuc", cu((8, 12, 16, 20), 20), [(1, 48)], S, 3),
        campaign("C4 SuperMUC tV (rank-2 bonds), L=20, 30, 40, standard vs delayed (+submatrix), 2x24",
                 "supermuc", sq((20, 30, 40), 20), [(2, 24)], S, 3),
        campaign("H1 Helma 2D Hubbard U=4, Beta=2, L=20..60, 3 schemes, auto, layouts 8x48 16x24 2x192, fakeintel on",
                 "helma", sq((20, 30, 40, 50, 60), 20), [(8, 48), (16, 24), (2, 192)], S, 3),
        campaign("H1b Helma fakeintel off control: L=40, 60, 8x48",
                 "helma", sq((40, 60), 20), [(8, 48)], S, 3),
        campaign("H2 Helma 3D cubic L=12, 16, 20, 3 schemes, 8x48",
                 "helma", cu((12, 16, 20), 20), [(8, 48)], S, 3),
        campaign("H3 Helma N_delay scan L=40, 60, delayed+submatrix, 7 values, 8x48",
                 "helma", sq((40, 60), 20), [(8, 48)], ("delayed", "submatrix"), 3, nd_values=range(7)),
    ]
    print("| campaign | jobs | node-h | core-h |")
    print("|---|---:|---:|---:|")
    tot = {"supermuc": [0, 0], "helma": [0, 0]}
    for c in camps:
        print("| %s | %d | %.0f | %.0f |" % (c["name"], c["jobs"], c["node_h"], c["core_h"]))
        tot[c["machine"]][0] += c["node_h"]
        tot[c["machine"]][1] += c["core_h"]
    for m, (nh, ch) in tot.items():
        print("| **total %s** | | **%.0f** | **%.0f** |" % (m, nh, ch))


if __name__ == "__main__":
    main()

#!/usr/bin/env -S python3 -I
"""Collect bench.sh output lines into docs/benchmarks/results.csv and a markdown size-scan table.

Usage: bench_report.py RAW_FILE [--out-dir DIR] [--repeat RAW_REPEAT]
       bench_report.py --ndelay RAW_FILE [--threads RAW_THREADS] [--conclusion MD_FILE] [--out-dir DIR]

RAW_FILE holds bench.sh lines:
"scheme L N N_delay threads sec_per_sweep acceptance nsweep commit date [pgreen_mean pgreen_max pphase_max]"
(commit and date come from the run's `info`; the three precision fields are optional, old files have none).
--repeat: a raw file with repeated runs of points of RAW_FILE (same fields); the Notes quantify the noise.
--ndelay: writes DIR/ndelay_scan.md (sec/sweep and precision against N_delay, per L and scheme) instead of the
size scan; RAW_THREADS (optional) adds the multithreaded comparison; MD_FILE (optional, hand-written) is appended
verbatim as the conclusion (default choice and its justification). Writes DIR/results.csv (scheme,L,N,N_delay,threads,sec_per_sweep,acceptance,nsweep)
and DIR/size_scan.md (sec/sweep per scheme, speedup vs standard per L, provenance, notes), and prints
the table to stdout. DIR defaults to docs/benchmarks. The output is a pure function of RAW_FILE plus the
compiler/flags of the current environment. Exits with an error if the rows come from different commits
or if a (scheme, L, N_delay, threads) combination occurs twice. Warns if the acceptance of the schemes
differs at one L (trajectories are identical with identical seeds, so it must not).
The committed input is docs/benchmarks/raw_size_scan.txt; regenerating from it must leave git clean.
"""
import csv
import os
import subprocess
import sys

ALFDIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SCHEMES = ["standard", "delayed", "submatrix"]
FIELDS = ["scheme", "L", "N", "N_delay", "threads", "sec_per_sweep", "acceptance", "nsweep"]
RAW_FIELDS = FIELDS + ["commit", "date"]
PREC_FIELDS = ["pg_mean", "pg_max", "pp_max"]


def run(cmd):
    return subprocess.run(["bash", "-c", cmd], cwd=ALFDIR, capture_output=True, text=True).stdout.strip()


def read_rows(path, check_commit=True):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if len(t) not in (len(RAW_FIELDS), len(RAW_FIELDS) + len(PREC_FIELDS)):
                sys.exit("bench_report.py: bad line (%d or %d fields expected): %s"
                         % (len(RAW_FIELDS), len(RAW_FIELDS) + len(PREC_FIELDS), line.strip()))
            rows.append(dict(zip(RAW_FIELDS + PREC_FIELDS, t)))
    commits = sorted({r["commit"] for r in rows})
    if check_commit and len(commits) != 1:
        sys.exit("bench_report.py: rows come from different commits: " + ", ".join(commits))
    seen = set()
    for r in rows:
        key = (r["scheme"], r["L"], r["N_delay"], r["threads"])
        if check_commit and key in seen:
            sys.exit("bench_report.py: duplicate row for scheme=%s L=%s N_delay=%s threads=%s" % key)
        seen.add(key)
    rows.sort(key=lambda r: (int(r["L"]), int(r["N_delay"]), int(r["threads"]), SCHEMES.index(r["scheme"])))
    return rows


def table(rows):
    by = {}
    for r in rows:
        by.setdefault((int(r["L"]), int(r["N"]), int(r["N_delay"]), int(r["threads"])), {})[r["scheme"]] = r
    out = ["| L | N | N_delay | threads | NSweep | standard [s/sweep] | delayed [s/sweep] | submatrix [s/sweep] "
           "| delayed speedup | submatrix speedup | acceptance |",
           "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    warn = []
    for (L, N, nd, th), d in sorted(by.items()):
        def t(s):
            return float(d[s]["sec_per_sweep"]) if s in d else None
        def sp(s):
            return "%.2f" % (t("standard") / t(s)) if t("standard") and t(s) else "-"
        def f(s):
            return "%.4g" % t(s) if t(s) else "-"
        accs = sorted({d[s]["acceptance"] for s in d})
        if len(accs) > 1:
            warn.append("WARNING: acceptance differs between schemes at L=%d: %s" % (L, accs))
        nsw = sorted({d[s]["nsweep"] for s in d})
        out.append("| %d | %d | %d | %d | %s | %s | %s | %s | %s | %s | %s |" % (
            L, N, nd, th, "/".join(nsw), f("standard"), f("delayed"), f("submatrix"),
            sp("delayed"), sp("submatrix"), accs[0] if len(accs) == 1 else "/".join(accs)))
    return "\n".join(out), warn


def sec(r):
    return float(r["sec_per_sweep"])


def notes(rows, repeat_rows):
    """Notes of size_scan.md, every number computed from the data."""
    by = {(r["scheme"], int(r["L"])): r for r in rows if r["threads"] == "1"}
    Ls = sorted({L for (_, L) in by})
    out = []
    out.append("Acceptance is identical between the three schemes at every L (same seed, identical Metropolis "
               "trajectories).")
    L0 = Ls[0]
    nd = by[("standard", L0)]["N_delay"]
    acc0 = float(by[("standard", L0)]["acceptance"])
    out.append("N_delay=%s is the maximum, not always the effective delay: the capacity is capped at "
               "min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed at the end of every "
               "slice, so the effective delay is min(N_delay, accepts per slice); at L=%d (%d vertices per slice, "
               "about %.0f accepts) N_delay=%s is %s." % (nd, L0, L0 * L0, acc0 * L0 * L0, nd,
                                                       "never reached" if acc0 * L0 * L0 < int(nd) else "reached"))
    if all(("standard", L) in by and ("delayed", L) in by and ("submatrix", L) in by for L in (14, 16)):
        g = {s: sec(by[(s, 16)]) / sec(by[(s, 14)]) for s in SCHEMES}
        sp = {(s, L): sec(by[("standard", L)]) / sec(by[(s, L)]) for s in ("delayed", "submatrix") for L in Ls}
        out.append("From L=14 to L=16 (N=196 to 256) the time per sweep grows by x%.1f (standard), x%.1f (delayed) "
                   "and x%.1f (submatrix); a pure N^2.4 trend would give x%.1f. "
                   "In relative terms the two low-rank schemes jump more than standard, which is why their "
                   "speedup dips (delayed %.2f -> %.2f, submatrix %.2f -> %.2f from L=14 to L=16)."
                   % (g["standard"], g["delayed"], g["submatrix"], (256 / 196.0) ** 2.4,
                      sp[("delayed", 14)], sp[("delayed", 16)], sp[("submatrix", 14)], sp[("submatrix", 16)]))
        if 18 in Ls:
            import math
            exps = {}
            for s in SCHEMES:
                for a, b in zip(Ls[:-1], Ls[1:]):
                    exps[(s, a, b)] = math.log(sec(by[(s, b)]) / sec(by[(s, a)])) / math.log((b * b) / float(a * a))
            def e(s, a, b):
                return "%.2f" % exps[(s, a, b)]
            recov = all(sp[(s, 18)] > sp[(s, 16)] for s in ("delayed", "submatrix"))
            out.append("L=18 (N=324) tests the guess that the jump at N=256 is a power-of-two / cache effect "
                       "(leading dimension 256). Local exponents d ln t / d ln N, standard / delayed / submatrix: "
                       "L=14-16: %s / %s / %s; L=16-18: %s / %s / %s; L=18-20: %s / %s / %s. Speedups at L=14, 16, "
                       "18, 20: delayed %s, submatrix %s. %s" % (
                           e("standard", 14, 16), e("delayed", 14, 16), e("submatrix", 14, 16),
                           e("standard", 16, 18), e("delayed", 16, 18), e("submatrix", 16, 18),
                           e("standard", 18, 20), e("delayed", 18, 20), e("submatrix", 18, 20),
                           ", ".join("%.2f" % sp[("delayed", L)] for L in (14, 16, 18, 20)),
                           ", ".join("%.2f" % sp[("submatrix", L)] for L in (14, 16, 18, 20)),
                           ("Observation: the speedup of both low-rank schemes recovers at L=18, consistent with a "
                            "localized slowdown at N=256 (an observation, not proof; one run per point)." if recov else
                            "Observation: the speedup does not recover at L=18, so the N=256 hypothesis is not "
                            "supported by this data (one run per point; not a refutation either).")))
    if repeat_rows:
        parts = []
        groups = {}
        for r in repeat_rows:
            groups.setdefault((r["scheme"], int(r["L"])), []).append(sec(r))
        for (sc, L), v in sorted(groups.items()):
            o = by.get((sc, L))
            if o is not None:
                allv = [sec(o)] + v
                parts.append("%s at L=%d: %d runs, %s s/sweep (spread max/min - 1 = %.1f %%)" % (
                    sc, L, len(allv), ", ".join("%.4g" % x for x in allv), 100.0 * (max(allv) / min(allv) - 1)))
        out.append("Repeats (same binary and parameters, run later in the session; raw_repeat.txt): " + "; ".join(parts) +
                   ". This is the only noise estimate (one point); differences of a few percent between delayed and "
                   "submatrix are not significant.")
    else:
        out.append("Single runs: timing noise is not quantified; differences of a few percent are not significant.")
    last = Ls[-1]
    nsw = int(by[("standard", last)]["nsweep"])
    out.append("L=%d has only %d sweeps in total (NSweep=%d, NBin=2)." % (last, 2 * nsw, nsw))
    return out


def provenance(rows, thr_rows=()):
    commits = sorted({r["commit"] for r in list(rows) + list(thr_rows)})
    dates = sorted({r["date"] for r in list(rows) + list(thr_rows)})
    compiler = run("gfortran --version | head -1")
    flags = run("source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1; echo $ALF_FLAGS_PROG")
    flags = " ".join(x for x in flags.split() if not x.startswith("-I"))
    threads = sorted({r["threads"] for r in list(rows) + list(thr_rows)}, key=int)
    return ["## Provenance", "",
            "- Git commit (of the binary, from `info`): %s" % ", ".join(commits),
            "- Compiler: %s" % compiler,
            "- Flags: `%s`" % flags,
            "- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`",
            "- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = %s" % ", ".join(threads),
            "- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)",
            "- Date (START TIME in `info`): %s" % ", ".join(dates), ""]


def ndelay_tables(rows):
    out = []
    warn = []
    Ls = sorted({int(r["L"]) for r in rows})
    for L in Ls:
        sub = [r for r in rows if int(r["L"]) == L]
        std = [r for r in sub if r["scheme"] == "standard"][0]
        accs = sorted({r["acceptance"] for r in sub})
        if len(accs) > 1:
            warn.append("acceptance differs at L=%d: %s" % (L, accs))
        N = int(std["N"])
        out += ["## L=%d (N=%d), NSweep=%s, acceptance %s, about %.0f accepts per slice" % (
            L, N, std["nsweep"], accs[0], float(accs[0]) * N), "",
            "Standard (reference): %.4g s/sweep; Precision Green mean/max %s / %s; Precision Phase max %s." % (
                sec(std), std["pg_mean"], std["pg_max"], std["pp_max"]), "",
            "| N_delay | delayed [s/sweep] | speedup | Green max | max / standard | Green mean | Phase max "
            "| submatrix [s/sweep] | speedup | Green max | max / standard | Green mean | Phase max |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        nds = sorted({int(r["N_delay"]) for r in sub if r["scheme"] != "standard"})
        for nd in nds:
            cells = ["%d" % nd]
            for sc in ("delayed", "submatrix"):
                r = [x for x in sub if x["scheme"] == sc and int(x["N_delay"]) == nd]
                if not r:
                    cells += ["-"] * 6
                    continue
                r = r[0]
                cells += ["%.4g" % sec(r), "%.2f" % (sec(std) / sec(r)), "%.3g" % float(r["pg_max"]),
                          "%.2f" % (float(r["pg_max"]) / float(std["pg_max"])), "%.3g" % float(r["pg_mean"]),
                          "%.3g" % float(r["pp_max"])]
            out.append("| " + " | ".join(cells) + " |")
        out.append("")
    return out, warn


def main_ndelay(argv, outdir):
    rows = read_rows(argv[2])
    for r in rows:
        if "pg_max" not in r:
            sys.exit("bench_report.py: --ndelay needs the precision fields in every line")
    thr = read_rows(argv[argv.index("--threads") + 1], check_commit=False) if "--threads" in argv else []
    out = ["# N_delay scan: delayed and submatrix", "",
           "Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, "
           "NBin=2, default Nwrap=10, N_global=1, same seed for all runs, 1 thread. sec/sweep = `CPU Time` / "
           "(NSweep*NBin), wall time of the sweep phase only (setup excluded). Speedup = standard / scheme. "
           "Precision = `Precision Green Max` / `Mean` and `Precision Phase Max` from `info` (the Green precision "
           "compares G propagated through the sweep with G recomputed at each stabilisation; trajectories are "
           "identical for all runs, so differences are pure round-off). The capacity is "
           "min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed every slice: the effective delay "
           "is min(N_delay, accepts per slice). Produced by `testsuite/lowrank/bench.sh` and "
           "`testsuite/lowrank/bench_report.py --ndelay`; raw data `raw_ndelay_scan.txt`.", ""]
    tabs, warn = ndelay_tables(rows)
    out += tabs
    if warn:
        out += ["## Warnings", ""] + ["- " + w for w in warn] + [""]
    if thr:
        out += ["## Multithreaded comparison (L=20)", "",
                "| scheme | N_delay | threads | NSweep | s/sweep | speedup vs standard, same threads |", "|---|---:|---:|---:|---:|---:|"]
        for r in sorted(thr, key=lambda r: (int(r["threads"]), SCHEMES.index(r["scheme"]))):
            base = [x for x in thr + rows if x["scheme"] == "standard" and x["L"] == r["L"] and x["threads"] == r["threads"]]
            out.append("| %s | %s | %s | %s | %.4g | %s |" % (
                r["scheme"], r["N_delay"], r["threads"], r["nsweep"], sec(r),
                "%.2f" % (sec(base[0]) / sec(r)) if base else "-"))
        out.append("")
    out += provenance(rows, thr)
    if "--conclusion" in argv:
        with open(argv[argv.index("--conclusion") + 1]) as f:
            out += [f.read().rstrip(), ""]
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "ndelay_scan.md"), "w") as f:
        f.write("\n".join(out))
    print("\n".join(tabs))
    for w in warn:
        print("WARNING: " + w, file=sys.stderr)


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    outdir = os.path.join(ALFDIR, "docs", "benchmarks")
    if "--out-dir" in argv:
        outdir = argv[argv.index("--out-dir") + 1]
    if argv[1] == "--ndelay":
        return main_ndelay(argv, outdir)
    rows = read_rows(argv[1])
    repeat_rows = read_rows(argv[argv.index("--repeat") + 1], check_commit=False) if "--repeat" in argv else []
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FIELDS)
        for r in rows:
            w.writerow([r[k] for k in FIELDS])
    tab, warn = table(rows)
    commit = rows[0]["commit"]
    dates = sorted({r["date"] for r in rows})
    compiler = run("gfortran --version | head -1")
    flags = run("source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1; echo $ALF_FLAGS_PROG")
    flags = " ".join(x for x in flags.split() if not x.startswith("-I"))
    threads = sorted({r["threads"] for r in rows})
    md = ["# Size scan: standard vs delayed vs submatrix", "",
          "Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, "
          "LOBS_ST=LOBS_EN=1, NBin=2, default Nwrap=10, N_Global=1, same seed for all schemes. "
          "sec/sweep = `CPU Time` / (NSweep*NBin); NSweep differs per L "
          "(same for the three schemes). Speedup = standard / scheme. Timing = wall time of the sweep phase only "
          "(`CPU Time` in `info`; the timer starts after the initial Green function computation, setup is "
          "excluded). Produced by "
          "`testsuite/lowrank/bench.sh` and `testsuite/lowrank/bench_report.py`.", "",
          tab, "",
          "## Provenance", "",
          "- Git commit (of the binary, from `info`): %s" % commit,
          "- Compiler: %s" % compiler,
          "- Flags: `%s`" % flags,
          "- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`",
          "- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = %s" % ", ".join(threads),
          "- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)",
          "- Date (START TIME in `info`): %s" % ", ".join(dates),
          "- The commit/date of standard L=14 and L=16 come from the `info` of the re-run of those two points "
          "(same binary, same day), not from the original runs.", ""]
    if warn:
        md += ["## Warnings", ""] + ["- " + w for w in warn] + [""]
    else:
        md += ["## Notes", ""] + [("- " + n) for n in notes(rows, repeat_rows)] + [""]
    with open(os.path.join(outdir, "size_scan.md"), "w") as f:
        f.write("\n".join(md))
    print(tab)
    for w in warn:
        print(w, file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)

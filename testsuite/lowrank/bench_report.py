#!/usr/bin/env -S python3 -I
"""Collect bench.sh output lines into docs/benchmarks/results.csv and a markdown size-scan table.

Usage: bench_report.py RAW_FILE [--out-dir DIR]

RAW_FILE holds bench.sh lines:
"scheme L N N_delay threads sec_per_sweep acceptance nsweep commit date" (commit and date come from
the run's `info`). Writes DIR/results.csv (scheme,L,N,N_delay,threads,sec_per_sweep,acceptance,nsweep)
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


def run(cmd):
    return subprocess.run(["bash", "-c", cmd], cwd=ALFDIR, capture_output=True, text=True).stdout.strip()


def read_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if len(t) != len(RAW_FIELDS):
                sys.exit("bench_report.py: bad line (%d fields expected): %s" % (len(RAW_FIELDS), line.strip()))
            rows.append(dict(zip(RAW_FIELDS, t)))
    commits = sorted({r["commit"] for r in rows})
    if len(commits) != 1:
        sys.exit("bench_report.py: rows come from different commits: " + ", ".join(commits))
    seen = set()
    for r in rows:
        key = (r["scheme"], r["L"], r["N_delay"], r["threads"])
        if key in seen:
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


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    outdir = os.path.join(ALFDIR, "docs", "benchmarks")
    if "--out-dir" in argv:
        outdir = argv[argv.index("--out-dir") + 1]
    rows = read_rows(argv[1])
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
          "- Date (START TIME in `info`): %s" % ", ".join(dates), ""]
    if warn:
        md += ["## Warnings", ""] + ["- " + w for w in warn] + [""]
    else:
        md += ["## Notes", "",
               "- Acceptance is identical between the three schemes at every L (same seed, identical Metropolis "
               "trajectories).",
               "- N_delay=32 is the maximum, not always the effective delay: the capacity is capped at "
               "min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed at the end of every "
               "slice, so the effective delay is min(N_delay, accepts per slice); at small L (L=6: 36 vertices "
               "per slice, about 15 accepts) N_delay=32 is never reached.",
               "- All three schemes jump at N=256 (L=16): standard 0.368 -> 1.098, delayed 0.219 -> 0.728, "
               "submatrix 0.210 -> 0.700 s/sweep from L=14 to L=16, steeper than the N^2.4 trend of the other "
               "points, and standard jumps most, so the L=16 speedups (about 1.5-1.6) dip below L=14 (about "
               "1.7). A power-of-two / cache effect (leading dimension 256) is an unverified guess; "
               "L=18 will be added in Task 15 to test it.",
               "- Single runs; the only repeat is standard at L=16 (1.098 vs 1.093 s/sweep, 0.5 %), so timing "
               "noise is otherwise unquantified and differences of a few percent between delayed and "
               "submatrix are not significant. L=20 has only 24 sweeps in total.", ""]
    with open(os.path.join(outdir, "size_scan.md"), "w") as f:
        f.write("\n".join(md))
    print(tab)
    for w in warn:
        print(w, file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)

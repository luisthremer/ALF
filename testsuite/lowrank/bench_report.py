#!/usr/bin/env python3
"""Collect bench.sh output lines into docs/benchmarks/results.csv and a markdown size-scan table.

Usage: bench_report.py RAW_FILE [--out-dir DIR]

RAW_FILE holds bench.sh lines: "scheme L N N_delay threads sec_per_sweep acceptance [nsweep]".
Writes DIR/results.csv (scheme,L,N,N_delay,threads,sec_per_sweep,acceptance,nsweep) and
DIR/size_scan.md (sec/sweep per scheme, speedup vs standard per L, provenance), and prints the
table to stdout. DIR defaults to docs/benchmarks. Warns if the acceptance of the schemes differs
at one L (trajectories are identical with identical seeds, so it must not).
"""
import csv
import datetime
import os
import subprocess
import sys

ALFDIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SCHEMES = ["standard", "delayed", "submatrix"]
FIELDS = ["scheme", "L", "N", "N_delay", "threads", "sec_per_sweep", "acceptance", "nsweep"]


def run(cmd):
    return subprocess.run(cmd, shell=True, cwd=ALFDIR, capture_output=True, text=True).stdout.strip()


def read_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if len(t) not in (7, 8):
                sys.exit("bench_report.py: bad line (7 or 8 fields expected): " + line.strip())
            rows.append(dict(zip(FIELDS, t + [""] * (8 - len(t)))))
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
    commit = run("git rev-parse --short HEAD")
    compiler = run("gfortran --version | head -1")
    flags = run("source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1; echo $ALF_FLAGS_PROG")
    flags = " ".join(x for x in flags.split() if not x.startswith("-I"))
    threads = sorted({r["threads"] for r in rows})
    md = ["# Size scan: standard vs delayed vs submatrix", "",
          "Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, "
          "LOBS_ST=LOBS_EN=1, NBin=2, default Nwrap=10, N_Global=1, same seed for all schemes. "
          "sec/sweep = wall time of the run (`CPU Time` in `info`) / (NSweep*NBin); NSweep differs per L "
          "(same for the three schemes). Speedup = standard / scheme. Produced by "
          "`testsuite/lowrank/bench.sh` and `testsuite/lowrank/bench_report.py`.", "",
          tab, "",
          "## Provenance", "",
          "- Git commit: %s" % commit,
          "- Compiler: %s" % compiler,
          "- Flags: `%s`" % flags,
          "- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`",
          "- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = %s" % ", ".join(threads),
          "- Machine: Apple M1 Pro (8 performance + 2 efficiency cores)",
          "- Date: %s" % datetime.date.today().isoformat(), ""]
    if warn:
        md += ["## Warnings", ""] + ["- " + w for w in warn] + [""]
    with open(os.path.join(outdir, "size_scan.md"), "w") as f:
        f.write("\n".join(md))
    print(tab)
    for w in warn:
        print(w, file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)

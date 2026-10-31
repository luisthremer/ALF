#!/usr/bin/env -S python3 -I
"""Collect bench.sh output lines into docs/benchmarks/results.csv and a markdown size-scan table.

Usage: bench_report.py RAW_FILE [--out-dir DIR] [--repeat RAW_REPEAT]
       bench_report.py --ndelay RAW_FILE --size-scan RAW_SIZE --repeat RAW_REPEAT [--threads RAW_THREADS] [--out-dir DIR]

Regeneration of the committed documents (from the repository root; a missing --repeat/--size-scan/--threads silently
drops the corresponding note or table):
    python3 -I testsuite/lowrank/bench_report.py docs/benchmarks/raw_size_scan.txt --repeat docs/benchmarks/raw_repeat.txt
    python3 -I testsuite/lowrank/bench_report.py --ndelay docs/benchmarks/raw_ndelay_scan.txt \\
        --size-scan docs/benchmarks/raw_size_scan.txt --repeat docs/benchmarks/raw_repeat.txt \\
        --threads docs/benchmarks/raw_threads.txt
The output is a pure function of the raw files plus the compiler/flags of the current environment; regenerating from the
committed raw files must leave git clean.

RAW_FILE holds bench.sh lines:
"scheme L N N_delay threads sec_per_sweep acceptance nsweep commit date [pgreen_mean pgreen_max pphase_max]"
(commit and date come from the run's `info`; the three precision fields are optional, old files have none).
--repeat: a raw file with repeated runs of points of the size scan (same fields; a row may be copied from another raw
file, then it is a run of the same point); the Notes quantify the noise from it.
Size-scan mode writes DIR/results.csv (scheme,L,N,N_delay,threads,sec_per_sweep,acceptance,nsweep) and DIR/size_scan.md
(sec/sweep per scheme, speedup vs standard per L, provenance, notes) and prints the table to stdout. DIR defaults to
docs/benchmarks. Exits with an error if the rows come from different commits or if a (scheme, L, N_delay, threads)
combination occurs twice. Warns if the acceptance of the schemes differs at one L (trajectories are identical with
identical seeds, so it must not).
--ndelay: writes DIR/ndelay_scan.md instead: sec/sweep and precision against N_delay per L and scheme, the
multithreaded table (RAW_THREADS) and the generated "Default N_delay" section (all numbers computed from the rows; the
run-to-run noise comes from --size-scan plus --repeat).
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
            n0, n1 = len(RAW_FIELDS), len(RAW_FIELDS) + len(PREC_FIELDS)
            if len(t) not in (n0, n1, n1 + 2):
                sys.exit("bench_report.py: bad line (%d, %d or %d fields expected): %s" % (n0, n1, n1 + 2, line.strip()))
            rows.append(dict(zip(RAW_FIELDS + PREC_FIELDS + ["model", "keys"], t)))
    commits = sorted({r["commit"] for r in rows})
    if check_commit and len(commits) != 1:
        sys.exit("bench_report.py: rows come from different commits: " + ", ".join(commits))
    models = sorted({r.get("model", "Hubbard") for r in rows})
    if len(models) > 1:   # tables group by (L, N_delay, threads) only: one model per raw file
        sys.exit("bench_report.py: rows of different models in one file: " + ", ".join(models))
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


def fit_exponent(rows, lmax=14):
    """Least-squares slope of ln t vs ln N of the standard rows with L <= lmax (1 thread)."""
    import math
    pts = sorted((int(r["N"]), sec(r)) for r in rows if r["scheme"] == "standard" and r["threads"] == "1"
                 and int(r["L"]) <= lmax and int(r["L"]) >= 8)
    if len(pts) < 2:
        return None
    x = [math.log(a) for a, _ in pts]
    y = [math.log(b) for _, b in pts]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)


def noise_spread(size_rows, repeat_rows):
    """(spread in percent, number of runs, scheme, L) of the most often repeated point, or None."""
    by = {(r["scheme"], int(r["L"])): sec(r) for r in size_rows if r["threads"] == "1"}
    groups = {}
    for r in repeat_rows:
        groups.setdefault((r["scheme"], int(r["L"])), []).append(sec(r))
    best = None
    for (sc, L), v in groups.items():
        if (sc, L) in by:
            allv = [by[(sc, L)]] + v
            if best is None or len(allv) > best[1]:
                best = (100.0 * (max(allv) / min(allv) - 1), len(allv), sc, L)
    return best


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
                   "and x%.1f (submatrix); the N^%.2f trend fitted to standard at L=8-14 would give x%.1f. "
                   "In relative terms the two low-rank schemes jump more than standard, which is why their "
                   "speedup dips (delayed %.2f -> %.2f, submatrix %.2f -> %.2f from L=14 to L=16)."
                   % (g["standard"], g["delayed"], g["submatrix"], fit_exponent(rows), (256 / 196.0) ** fit_exponent(rows),
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
                       "L=14-16: %s / %s / %s; L=16-18: %s / %s / %s; L=18-20: %s / %s / %s. Speedups at L=%s: "
                       "delayed %s, submatrix %s. %s" % (
                           e("standard", 14, 16), e("delayed", 14, 16), e("submatrix", 14, 16),
                           e("standard", 16, 18), e("delayed", 16, 18), e("submatrix", 16, 18),
                           e("standard", 18, 20), e("delayed", 18, 20), e("submatrix", 18, 20),
                           ", ".join(str(L) for L in (14, 16, 18, 20) if L in Ls),
                           ", ".join("%.2f" % sp[("delayed", L)] for L in (14, 16, 18, 20) if L in Ls),
                           ", ".join("%.2f" % sp[("submatrix", L)] for L in (14, 16, 18, 20) if L in Ls),
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
        out.append("Repeats (same binary and parameters, run later in the session; raw_repeat.txt, whose second row is the "
                   "standard L=12 row of the N_delay scan, raw_ndelay_scan.txt, copied there): " + "; ".join(parts) +
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


def pct(x):
    v = 100.0 * x
    return "0.0 %" if abs(v) < 0.05 else "%+.1f %%" % v


def default_section(rows, noise):
    """The generated 'Default N_delay' section of ndelay_scan.md. Every number is computed from the rows."""
    Ls = sorted({int(r["L"]) for r in rows})
    by = {(int(r["L"]), r["scheme"], int(r["N_delay"])): r for r in rows}
    nds = sorted({int(r["N_delay"]) for r in rows if r["scheme"] != "standard"})
    LS = ("delayed", "submatrix")
    D = 32
    tol = noise[0] if noise else 0.0
    best, near = {}, {}
    for L in Ls:
        for sc in LS:
            ts = {nd: sec(by[(L, sc, nd)]) for nd in nds}
            b = min(ts, key=ts.get)
            best[(L, sc)] = b
            near[(L, sc)] = [nd for nd in nds if ts[nd] <= ts[b] * (1 + tol / 100.0)]
    Ns = sorted({int(by[(L, "standard", D)]["N"]) for L in Ls})

    def lst(v):
        return " and ".join(str(x) for x in v) if len(v) <= 2 else ", ".join(str(x) for x in v)

    out = ["## Default N_delay", ""]
    noise_txt = ("the run-to-run noise (%.1f %% spread of %d runs of %s at L=%d)" % (noise[0], noise[1], noise[2], noise[3])
                 if noise else "the run-to-run noise (not quantified in this call: no --size-scan/--repeat given)")
    ok = all(D in near[k] for k in near)
    out.append(
        "**Default N_delay = %d** (scanned at N = %s). The lowest time is at N_delay = (%s); the points within %s of the "
        "best are: %s. %s It was chosen as the single default for N = %d-%d and is untested outside that range." % (
            D, " and ".join(str(n) for n in Ns),
            "; ".join("L=%d %s %d" % (L, sc, best[(L, sc)]) for L in Ls for sc in LS), noise_txt,
            "; ".join("L=%d %s %s" % (L, sc, lst(near[(L, sc)])) for L in Ls for sc in LS),
            ("So %d is best or within noise of the best for both schemes at both sizes." % D) if ok else
            ("%d is NOT within noise of the best for: %s." % (D, ", ".join("L=%d %s" % k for k in near if D not in near[k]))),
            Ns[0], Ns[-1]))
    out.append("")
    # precision
    dm = {(L, sc): float(by[(L, sc, D)]["pg_mean"]) / float(by[(L, "standard", D)]["pg_mean"]) - 1 for L in Ls for sc in LS}
    rise = {}
    mono = {}
    for L in Ls:
        for sc in LS:
            m = [float(by[(L, sc, nd)]["pg_mean"]) for nd in nds]
            rise[(L, sc)] = m[-1] / m[0] - 1
            mono[(L, sc)] = all(b >= a for a, b in zip(m, m[1:]))
    gmax = [float(r["pg_max"]) / float(by[(int(r["L"]), "standard", D)]["pg_max"]) for r in rows if r["scheme"] != "standard"]
    out.append(
        "Precision check: the submatrix Green mean %s with N_delay (%s from N_delay=%d to %d), so precision mildly favours the "
        "lower end; at N_delay=%d the Green mean differs from standard by %s (delayed: %s), so there is no visible degradation "
        "at the chosen value. `Precision Green Max` of every run is within %.1f %% of standard, "
        "presumably set by one event of the (identical) trajectory, so it does not discriminate between N_delay values." % (
            "rises mildly" if all(rise[(L, "submatrix")] > 0 for L in Ls) else "does not rise everywhere",
            ", ".join("%s at L=%d" % (pct(rise[(L, "submatrix")]), L) for L in Ls), nds[0], nds[-1], D,
            ", ".join("%s at L=%d" % (pct(dm[(L, "submatrix")]), L) for L in Ls),
            ", ".join("%s at L=%d" % (pct(dm[(L, "delayed")]), L) for L in Ls), 100.0 * max(abs(g - 1) for g in gmax)))
    out += ["", "Time relative to N_delay=%d (s/sweep at N_delay / s/sweep at %d - 1; the numbers behind the choice):" % (D, D), "",
            "| L | scheme | " + " | ".join(str(nd) for nd in nds) + " |",
            "|---:|---|" + "---:|" * len(nds)]
    for L in Ls:
        for sc in LS:
            out.append("| %d | %s | " % (L, sc) + " | ".join(
                pct(sec(by[(L, sc, nd)]) / sec(by[(L, sc, D)]) - 1) for nd in nds) + " |")
    out += ["", "Interpretation (not measured; no profile against N_delay was taken): beyond the optimum the time rises, which we "
            "attribute to costs that grow with the number k of pending vertices outgrowing the gain in flush efficiency. "
            "For delayed, forming the new column G(:,p) = G0(:,p) - X Y(p,:)^T is O(N k) per accept (BLAS-2); summed over a "
            "slice this grows with k, while the BLAS-3 flush total does not. For submatrix, the per-proposal solve with the "
            "Gamma-tilde^-1 matrix is O(k^2) and the border update per accept is also O(k^2)."]
    # plateau
    parts = []
    for L in Ls:
        r0 = by[(L, "standard", D)]
        A = float(r0["acceptance"]) * int(r0["N"])
        above = [nd for nd in nds if nd >= A]
        if len(above) >= 2:
            a, b = above[-2], above[-1]
            diffs = [abs(sec(by[(L, sc, a)]) / sec(by[(L, sc, b)]) - 1) for sc in LS]
            parts.append(
                "At L=%d (about %.0f accepts per slice) N_delay = %s are all above the cap; N_delay=%d and %d give %s s/sweep "
                "(delayed) and %s (submatrix), %s." % (
                    L, A, lst(above), a, b,
                    "%.4g and %.4g" % (sec(by[(L, "delayed", a)]), sec(by[(L, "delayed", b)])),
                    "%.4g and %.4g" % (sec(by[(L, "submatrix", a)]), sec(by[(L, "submatrix", b)])),
                    ("flat to within %.1f %% (noise %.1f %%): this is the plateau, the effective delay is capped by the "
                     "accepts per slice" % (100 * max(diffs), tol) if max(diffs) * 100 <= 1.5 * tol else
                     "differing by up to %.1f %%, not flat" % (100 * max(diffs))))
            )
        elif len(above) == 1:
            parts.append("At L=%d (about %.0f accepts per slice) only N_delay=%d is above the cap, so the grid has no plateau "
                         "pair there and the rise from %d to the next values happens below the cap." % (L, A, above[0], D))
        else:
            parts.append("At L=%d (about %.0f accepts per slice) no grid point is above the cap." % (L, A))
    out += ["", " ".join(parts), ""]
    return out


def main_ndelay(argv, outdir):
    rows = read_rows(argv[2])
    for r in rows:
        if "pg_max" not in r:
            sys.exit("bench_report.py: --ndelay needs the precision fields in every line")
    thr = read_rows(argv[argv.index("--threads") + 1], check_commit=False) if "--threads" in argv else []
    size_rows = read_rows(argv[argv.index("--size-scan") + 1]) if "--size-scan" in argv else []
    rep_rows = read_rows(argv[argv.index("--repeat") + 1], check_commit=False) if "--repeat" in argv else []
    noise = noise_spread(size_rows, rep_rows)
    out = ["# N_delay scan: delayed and submatrix", "",
           "Mz Hubbard model on an L x L square lattice, U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, "
           "NBin=2, default Nwrap=10, N_global=1, same seed for all runs, 1 thread (except the multithreaded table). sec/sweep = `CPU Time` / "
           "(NSweep*NBin), wall time of the sweep phase only (setup excluded). Speedup = standard / scheme. "
           "Precision = `Precision Green Max` / `Mean` and `Precision Phase Max` from `info` (the Green precision "
           "compares G propagated through the sweep with G recomputed at each stabilisation; trajectories are "
           "identical for all runs, so differences are pure round-off). The capacity is "
           "min(N_delay, vertices per slice[, Ndim for submatrix]) and G is flushed every slice: the effective delay "
           "is min(N_delay, accepts per slice). Produced by `testsuite/lowrank/bench.sh` and "
           "`testsuite/lowrank/bench_report.py --ndelay`; raw data `raw_ndelay_scan.txt`, `raw_threads.txt`. Regeneration: see the "
           "docstring of `bench_report.py`.", ""]
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
    out += default_section(rows, noise)
    out += provenance(rows, thr)
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
          ] + (["- The commit/date of standard L=14 and L=16 come from the `info` of the re-run of those two points "
                "(same binary, same day), not from the original runs (rows without precision fields are from the "
                "original scan)."] if any("pg_max" not in r for r in rows) and
               all(any(r["scheme"] == "standard" and int(r["L"]) == L for r in rows) for L in (14, 16)) else []) + [""]
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

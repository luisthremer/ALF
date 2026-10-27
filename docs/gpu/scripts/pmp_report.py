#!/usr/bin/env python3
"""Summarise folded stack samples (nividia/pmp.py output) into the R1 hot-kernel rows.

usage: pmp_report.py NAME=stacks.txt [NAME=stacks.txt ...]   -> markdown tables on stdout
Each sample is assigned to one kernel row (first matching rule below, i.e. the outermost R1 routine
that owns the work) and its leaf to a BLAS level. Samples whose stack did not unwind to MAIN__
(OpenBLAS assembly kernels without CFI, e.g. dgemm_kernel_ZEN, ztrmm_kernel_*) are spread over the
rows in proportion to the resolved samples of the same BLAS routine (e.g. all stacks containing dgemm_).
"""
import collections
import re
import sys

# (row, predicate on the set of frame names); order matters
RULES = [
    ("CGR2_2 (tau stabilization, 2N×2N)", lambda f: "cgr2_2_mod::cgr2_2" in f),
    ("CGR (equal-time G from UDV)", lambda f: "cgr1_mod::cgr" in f),
    ("UDV decompose (QRP, ZUNGQR, ZTRMM)", lambda f: "udv_state_mod::decompose_udv_state" in f),
    ("Upgrade2 rank-1 update (ZGERU)", lambda f: "upgrade_mod::upgrade2" in f and ("zgeru_" in f or "zgerc_" in f)),
    ("Upgrade2 ratio + vectors (ZGEMV, alloc)", lambda f: "upgrade_mod::upgrade2" in f),
    ("Op_Wrapup/Wrapdo (e^V G e^-V)", lambda f: any(x.startswith("operator_mod::op_wrap") for x in f)),
    ("WRAPUR/WRAPUL B products (eq + tau)", lambda f: "wrapur_mod::wrapur" in f or "wrapul_mod::wrapul" in f),
    ("PROPR/PROPRM1 (tau B products)", lambda f: "tau_m_mod::propr" in f or "tau_m_mod::proprm1" in f),
    ("Hop_mod B_T G B_T^-1 in WRAPGR", lambda f: "wrapgr_mod::wrapgrup" in f or "wrapgr_mod::wrapgrdo" in f),
    ("ObserT (time-displaced measurements)", lambda f: "hamiltonian_main.ham_hubbard_smod::obsert" in f),
    ("Obser (equal-time measurements)", lambda f: "hamiltonian_main.ham_hubbard_smod::obser" in f),
    ("Control_Precision*", lambda f: any(x.startswith("control::control_precision") for x in f)),
    ("tau_m other", lambda f: "tau_m_mod::tau_m" in f),
]
OTHER = "other (setup, I/O, bins, unresolved)"
ROWS = [r for r, _ in RULES] + [OTHER]

LEVELS = [
    ("BLAS-3 (gemm/trsm/trmm kernels)", re.compile(r"^[dz](gemm|trsm|trmm|symm|hemm)_|_kernel_[lrn]_ZEN|^[dz]gemm_|(gemm|trsm|trmm)_(kernel|beta|[io][nt]copy|i[nt]copy|o[nt]copy)")),
    ("BLAS-1/2 (axpy/ger/gemv/scal/nrm2/dot)", re.compile(r"axpy|ger[uc]?_|gemv|scal|nrm2|dot|copy_k|swap|amax")),
    ("LAPACK auxiliary (zlacpy, zlaqp2, zlarf*...)", re.compile(r"^z(la|lar|ung|unm|geq|lap|get)")),
    ("libc / allocation", re.compile(r"^(free|malloc|_int_|memcpy|memset|__memmove|__memset|__memcpy|cfree|calloc)")),
    ("ALF Fortran loops", re.compile(r"::")),
]


def load(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s, n = line.rstrip("\n").rsplit(" ", 1)
            out.append((s.split(";"), int(n)))
    return out


def family(leaf):
    m = re.match(r"([dz][a-z0-9]+?)_(kernel|beta|i[nt]copy|o[nt]copy)", leaf)
    return m.group(1) + "_" if m else None


def level(leaf):
    for name, pat in LEVELS:
        if pat.search(leaf):
            return name
    return "other"


def classify(samples):
    rows, levels = collections.Counter(), collections.Counter()
    by_family = collections.defaultdict(collections.Counter)
    unresolved, total = [], 0
    for frames, n in samples:
        total += n
        levels[level(frames[-1])] += n
        if frames[0] not in ("MAIN__", "main"):
            unresolved.append((frames, n))
            continue
        fs = set(frames)
        row = next((r for r, p in RULES if p(fs)), OTHER)
        rows[row] += n
        for x in fs:
            if re.match(r"^[dz][a-z0-9]+_$", x):
                by_family[x][row] += n
    for frames, n in unresolved:
        fam = family(frames[-1])
        dist = by_family.get(fam)
        if dist:
            tot = sum(dist.values())
            for r, k in dist.items():
                rows[r] += n * k / tot
        else:
            rows[OTHER] += n
    return rows, levels, total, sum(n for _, n in unresolved)


if __name__ == "__main__":
    results = {}
    for arg in sys.argv[1:]:
        name, path = arg.split("=", 1)
        results[name] = classify(load(path))
    names = list(results)
    print("| Kernel row (R1 §2) | " + " | ".join(names) + " |")
    print("|---|" + "---:|" * len(names))
    for r in ROWS:
        print(f"| {r} | " + " | ".join(f"{100 * results[n][0][r] / results[n][2]:.1f}" for n in names) + " |")
    print("| samples (unresolved, redistributed) | " +
          " | ".join(f"{results[n][2]} ({results[n][3]})" for n in names) + " |")
    print()
    print("| Leaf type | " + " | ".join(names) + " |")
    print("|---|" + "---:|" * len(names))
    for lv in [l for l, _ in LEVELS] + ["other"]:
        print(f"| {lv} | " + " | ".join(f"{100 * results[n][1][lv] / results[n][2]:.1f}" for n in names) + " |")

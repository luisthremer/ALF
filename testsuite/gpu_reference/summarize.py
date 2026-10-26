#!/usr/bin/env python3
"""Analyse finished run directories and write a compact result file per run.

usage: summarize.py RUNDIR [RUNDIR ...]
Writes RUNDIR/reference.json: pyALF results (all observables with errors),
precision numbers from the 'info' file, and timing.txt contents.
"""
import json
import re
import sys

from py_alf.analysis import analysis
from py_alf.ana import load_res

INFO_PATTERNS = {
    "precision_green_mean_max": r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)",
    "precision_phase_max": r"Precision Phase, Max\s*:\s*(\S+)",
    "precision_tau_mean_max": r"Precision tau\s+Mean, Max\s*:\s*(\S+)\s+(\S+)",
    "acceptance": r"^\s*Acceptance\s*:\s*(\S+)",
    "cpu_time": r"CPU Time\s*:\s*(\S+)",
}


def parse_info(path):
    out = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    for key, pat in INFO_PATTERNS.items():
        m = re.findall(pat, text, re.M)
        if m:  # last occurrence = last MPI-reduced block
            last = m[-1]
            out[key] = [float(x.replace("D", "E")) for x in (last if isinstance(last, tuple) else (last,))]
    return out


def to_plain(x):
    if hasattr(x, "tolist"):
        x = x.tolist()
    if isinstance(x, complex):
        return [x.real, x.imag]
    if isinstance(x, list):
        return [to_plain(v) for v in x]
    return x


if __name__ == "__main__":
    for d in sys.argv[1:]:
        analysis(d)
        row = load_res([d]).iloc[0]
        res = {k: to_plain(v) for k, v in row.items()}
        res["info"] = parse_info(f"{d}/info")
        try:
            with open(f"{d}/timing.txt", encoding="utf-8") as f:
                res["timing"] = dict(line.strip().split("=", 1) for line in f if "=" in line)
        except FileNotFoundError:
            pass
        with open(f"{d}/reference.json", "w", encoding="utf-8") as f:
            json.dump(res, f, default=str)
        print(d, "Ener", res.get("Ener_scal0"), "+-", res.get("Ener_scal0_err"), res["info"])

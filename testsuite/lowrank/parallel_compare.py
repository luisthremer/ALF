#!/usr/bin/env python3
"""One-line verdict of stat_compare.py for the parallel tier (parallel_tier.sh).

Usage: python3 -I parallel_compare.py DIR_A DIR_B LOGFILE [LABEL]    (A = reference, B = candidate)

Runs stat_compare.compare (the pre-registered rules, unchanged) on the two analysed run directories, writes the full
stat_compare output to LOGFILE and prints one summary line:
    PASS|FAIL[-ARRAY]|ERROR LABEL: scalars max|z|=.. (name) sign z=..; <array>: n=.. frac>2=..%(lim ..%) max|z|=..; ...
Exit codes: 0 PASS; 1 FAIL with a scalar failure (decisive: a bug, no replicate); 3 FAIL of an array rule only
(the pre-registered protocol allows ONE replicate with fresh seeds); 2 unusable input.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stat_compare  # noqa: E402


def main(argv):
    if len(argv) < 4:
        sys.stderr.write(__doc__)
        return 2
    da, db, log = argv[1:4]
    label = argv[4] if len(argv) > 4 else "%s vs %s" % (da, db)
    lines = []
    rc = stat_compare.compare(da, db, out=lines.append)
    with open(log, "w") as f:
        f.write("\n".join(lines) + "\n")
    if rc == 2:
        print("ERROR %s: stat_compare found unusable input, see %s (%s)" % (label, log, lines[-1] if lines else ""))
        return 2

    scal_bad = False
    zmax, zname = 0.0, "-"
    sign_z = None
    arrays = []
    section = "scalar"
    for ln in lines:
        if ln.startswith("array "):
            section = "array"
            continue
        if ln.startswith("precision"):
            section = "prec"
            continue
        t = ln.split()
        if section == "scalar" and t and t[-1] in ("ok", "FAIL") and len(t) >= 7:
            try:
                z = float(t[-2])
            except ValueError:
                continue
            if t[-1] == "FAIL":
                scal_bad = True
            if t[0] == "avg":
                sign_z = z
            elif not math.isfinite(z):
                zmax, zname = math.nan, t[0]
            elif math.isfinite(zmax) and abs(z) > zmax:
                zmax, zname = abs(z), t[0]
        elif section == "array" and t and t[0].endswith(("JK", "JR")) and len(t) >= 10:
            # array unique skip0 nonfin n>2 frac>2 limit max|z| rms z  ok|FAIL|(info)
            if t[0].endswith("JK"):
                arrays.append("%s n=%s frac>2=%s(lim %s) max|z|=%s rms=%s %s" % (
                    t[0], t[1], t[5], t[6], t[7], t[8], "FAIL" if t[-1] == "FAIL" else ""))
    arr_bad = any(a.endswith("FAIL") for a in arrays)
    if rc == 0:
        verdict = "PASS"
    elif scal_bad or not arrays:
        verdict = "FAIL"
    elif arr_bad:
        verdict = "FAIL-ARRAY"
    else:      # non-finite / other failure not attributable to a JK array rule
        verdict = "FAIL"
    wall = ""
    p = os.path.join(db, "wall_seconds")
    if os.path.exists(p):
        wall = " wall=%ss" % open(p).read().strip()
    print("%s %s: scalars max|z|=%.2f (%s) sign z=%s; %s%s" % (
        verdict, label, zmax, zname, "%.2f" % sign_z if sign_z is not None else "-", "; ".join(arrays), wall))
    return {"PASS": 0, "FAIL": 1, "FAIL-ARRAY": 3}[verdict]


if __name__ == "__main__":
    sys.exit(main(sys.argv))

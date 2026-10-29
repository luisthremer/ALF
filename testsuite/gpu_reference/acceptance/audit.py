#!/usr/bin/env python3
"""Audit an implementer's reported numbers against the verdict computed from raw outputs.

usage: audit.py CLAIMS.json VERDICT.json
CLAIMS.json: {"<gate_id>/<index>/<key>": "<number as printed in the report>", ...}, e.g.
  {"same_seed_12/0/overall_max": "3.103e-10", "paired_log/0/sd_ratio": "1.06"}
A claim matches if the recomputed value, rounded to the claim's printed significant digits, equals it.
Prints every claim with MATCH / MISMATCH / MISSING. Exit 0 only if all claims match.
"""
import json
import sys
from decimal import Decimal


def digits(s):
    return len(Decimal(s.lower().replace("e", "E")).as_tuple().digits)


def main(cf, vf):
    claims, v = json.load(open(cf)), json.load(open(vf))
    gates = {g["id"]: g for g in v["gates"]}
    bad = 0
    for key, printed in claims.items():
        gid, idx, field = key.split("/", 2)
        g = gates.get(gid)
        try:
            val = g["result"][int(idx)]
            for f in field.split("."):
                val = val[f]
        except (TypeError, KeyError, IndexError):
            print(f"MISSING  {key} (claimed {printed})")
            bad += 1
            continue
        nd = digits(printed)
        got = f"{float(val):.{max(nd - 1, 0)}e}"
        ok = float(got) == float(f"{float(printed):.{max(nd - 1, 0)}e}")
        print(f"{'MATCH   ' if ok else 'MISMATCH'} {key}: claimed {printed}, recomputed {float(val):.6e}")
        bad += not ok
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))

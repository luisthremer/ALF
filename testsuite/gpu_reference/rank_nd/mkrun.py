#!/usr/bin/env python3
"""Ticket 401: create run directories (parameters + seeds) with pyALF, only_prep.  Run inside the pyALF venv.
usage: mkrun.py ALF_SRC OUTDIR NAME [--ham Hubbard_Plain_Vanilla] [--seed-offset N] key=value ...
Defaults are the 4-site chain of ADR-003 / 401.json: Square L1=4 L2=1, U=4, chem 0, beta 2, dtau 0.1, Nwrap 5.
Key values are parsed as python literals (True/False/numbers) else kept as strings."""
import argparse, ast, os, sys
from py_alf import ALF_source, Simulation

DEFAULT = {"Lattice_type": "Square", "L1": 4, "L2": 1, "ham_T": 1.0, "Ham_chem": 0.0, "ham_U": 4.0, "beta": 2.0, "dtau": 0.1,
           "Nwrap": 5, "Symm": True, "Ltau": 0, "NSweep": 10, "NBin": 10}

def lit(v):
    try:
        return ast.literal_eval(v)
    except Exception:
        return v

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("alf_dir"); p.add_argument("outdir"); p.add_argument("name")
    p.add_argument("--ham", default="Hubbard_Plain_Vanilla")
    p.add_argument("--seed-offset", type=int, default=0)
    p.add_argument("kv", nargs="*")
    a = p.parse_args()
    d = dict(DEFAULT)
    for kv in a.kv:
        k, v = kv.split("=", 1); d[k] = lit(v)
    # keys of the Hubbard model that Plain_Vanilla does not know are dropped by pyALF's parameter check; keep them minimal
    src = ALF_source(alf_dir=os.path.abspath(a.alf_dir))
    sim = Simulation(src, a.ham, d, sim_root=os.path.abspath(a.outdir), sim_dir=a.name, machine="GNU", mpi=False)
    sim.run(only_prep=True)
    sf = os.path.join(os.path.abspath(a.outdir), a.name, "seeds")
    if a.seed_offset:
        lines = open(sf).readlines(); k = a.seed_offset % len(lines)
        open(sf, "w").writelines(lines[k:] + lines[:k])
    print(os.path.join(a.outdir, a.name))

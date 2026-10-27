#!/usr/bin/env python3
"""Ticket 104: create run directories for the pre-pivot QR study.

usage: prep.py OUTDIR ALF_DIR NSWEEP_NBIN_FILE_OR_DEFAULT  (see CASES below)
Writes OUTDIR/<case>/{parameters,seeds}; run_all.sh copies them per build.
Case name: L<L>_b<beta>_U<U>_w<Nwrap>[_dense]
"""
import os
import sys

from py_alf import ALF_source, Simulation

COMMON = {"Lattice_type": "Square", "Model": "Hubbard", "Symm": False, "ham_T": 1.0,
          "ham_chem": 0.0, "Mz": True, "dtau": 0.1, "Ltau": 1, "n_skip": 2,
          "Checkerboard": True}


def case_dict(name, nsweep, nbin):
    parts = name.split("_")
    L = int(parts[0][1:]); beta = float(parts[1][1:]); U = float(parts[2][1:]); nw = int(parts[3][1:])
    d = dict(COMMON, L1=L, L2=L, beta=beta, ham_U=U, Nwrap=nw, NSweep=nsweep, NBin=nbin)
    if "dense" in parts:
        d["Checkerboard"] = False
    return d


if __name__ == "__main__":
    outdir, alf_dir = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    src = ALF_source(alf_dir=alf_dir)
    # arguments after: name:nsweep:nbin
    for spec in sys.argv[3:]:
        name, ns, nb = spec.split(":")
        sim = Simulation(src, "Hubbard", case_dict(name, int(ns), int(nb)), sim_root=outdir,
                         sim_dir=name, machine="GNU", mpi=True)
        sim.run(only_prep=True)
        print(name)

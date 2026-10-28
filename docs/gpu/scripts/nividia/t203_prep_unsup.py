#!/usr/bin/env python3
"""Ticket 203: run directories for models outside the GPU subset: tV (non-diagonal Op_V) and Hubbard Continuous=True
(type-3 fields).  usage: t203_prep_unsup.py OUTDIR ALF_DIR   (pyALF only_prep; seeds taken from the base runs by the caller)"""
import os, sys
from py_alf import ALF_source, Simulation
out, alf = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
src = ALF_source(alf_dir=alf)
common = dict(Lattice_type="Square", L1=4, L2=4, Checkerboard=True, Symm=False, beta=2.0, dtau=0.1, Nwrap=10,
              Ltau=0, NSweep=20, NBin=2)
cases = {"tV_4x4": ("tV", dict(common, Model="tV", ham_V=1.0)),
         "hub4x4_cont": ("Hubbard", dict(common, Model="Hubbard", ham_U=4.0, Mz=True, Continuous=True))}
for name, (ham, d) in cases.items():
    Simulation(src, ham, d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
    print(os.path.join(out, name))

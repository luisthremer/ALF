#!/usr/bin/env python3
"""Create run directories (parameters + seeds) for reference configs.

usage: prep.py OUTDIR CONFIG [CONFIG ...]   (CONFIG names from configs.py, or 'all')
Writes OUTDIR/<config>/{parameters,seeds}. The binary is supplied by run.sh.
"""
import argparse
import os
import sys

from py_alf import ALF_source, Simulation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from configs import CONFIGS  # noqa: E402

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("outdir")
    p.add_argument("configs", nargs="+")
    p.add_argument("--alf-dir", default=os.environ.get("ALF_DIR", "../.."))
    p.add_argument("--nsweep", type=int, help="override NSweep (calibration runs)")
    p.add_argument("--nbin", type=int, help="override NBin (calibration runs)")
    args = p.parse_args()

    names = list(CONFIGS) if args.configs == ["all"] else args.configs
    alf_src = ALF_source(alf_dir=os.path.abspath(args.alf_dir))
    for name in names:
        sim_dict = dict(CONFIGS[name])
        if args.nsweep:
            sim_dict["NSweep"] = args.nsweep
        if args.nbin:
            sim_dict["NBin"] = args.nbin
        sim = Simulation(alf_src, "Hubbard", sim_dict,
                         sim_root=os.path.abspath(args.outdir), sim_dir=name,
                         machine="GNU", mpi=True)
        sim.run(only_prep=True)
        print(os.path.join(args.outdir, name))

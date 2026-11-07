import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from py_alf import ALF_source, Simulation
from ed_cases import CASES
c = CASES[sys.argv[1]]
src = ALF_source(alf_dir=os.path.expanduser("~/alf-gpu/work/t401ed/src"))
d = dict(c["sim_dict"]); d.update(dtau=0.1, NSweep=40, NBin=20)
sim = Simulation(src, c["ham_name"], d, sim_root=os.path.expanduser("~/alf-gpu/work/t401ed/runs"), sim_dir="probe_" + sys.argv[1], machine="GNU", mpi=False)
sim.run(only_prep=True)

# P2.6.1 H100 (Helma, GPU=nvidia) cases: the three models at a larger size, standard / delayed 16, 32 / submatrix 16, plus the
# small P2.6.1 cases and the Hubbard regression.  Runner ../p25_runs/p25_run.py (GPU path vs ALF_GPU_DISABLE=1, same seeds).
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases_p26 import CASES as SMALL, MODELS, c
BIG = {"kondo": {"L1": "8", "L2": "8"}, "tV": {"L1": "12", "L2": "12"}, "sp": {"L1": "12", "L2": "12"}}
SCH = {"std": {}, "dly16": {"Sequential_scheme": '"delayed"', "N_delay": "16"},
       "dly32": {"Sequential_scheme": '"delayed"', "N_delay": "32"}, "sub16": {"Sequential_scheme": '"submatrix"', "N_delay": "16"}}
CASES = list(SMALL)
for m, me in MODELS.items():
    for k, ke in SCH.items():
        e = dict(me); e.update(BIG[m]); e.update(ke); e.update({"NSweep": "10"})
        CASES.append(c("big-%s-%s" % (m, k), "works", e))

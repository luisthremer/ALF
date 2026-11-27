# P2.6 Calc_Fl cases (Hubbard_Plain_Vanilla, flavour 1 reconstructed from flavour 2), the pv-* cases of gpu_e2e.py, for the vendor
# checks with ../p25_runs/p25_run.py (GPU path vs ALF_GPU_DISABLE=1, same seeds): repulsive U at half filling and attractive U off
# half filling x standard / delayed 4 / submatrix 4, plus an 8x8 size and the Hubbard regression.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases_p26 import c, SCHEMES
PV = {"ham_name": '"Hubbard_Plain_Vanilla"'}
NEGU = dict(PV, ham_U="-4.d0", ham_chem="0.5d0")
CASES = []
for k in ("std", "dly4", "sub4"):
    CASES.append(c("pv-half-%s" % k, "works", dict(PV, **SCHEMES[k])))
    CASES.append(c("pv-negU-%s" % k, "works", dict(NEGU, **SCHEMES[k])))
CASES += [c("pv-negU-8x8-sub16", "works", dict(NEGU, L1="8", L2="8", Sequential_scheme='"submatrix"', N_delay="16")),
          c("pv-half-8x8-dly16", "works", dict(PV, L1="8", L2="8", Sequential_scheme='"delayed"', N_delay="16")),
          c("hub-mz-4x4", "works"), c("hub-sub4", "works", SCHEMES["sub4"])]

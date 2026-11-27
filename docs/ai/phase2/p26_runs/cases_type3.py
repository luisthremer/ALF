# P2.6 p26-langevin: field type 3 (continuous HS fields) on the device path, Hubbard `Continuous = T` with sequential updates,
# Mz and SU(2) x standard / delayed / submatrix, 4x4, 8x8 and 12x12 (precision gate: Precision Green mean GPU / CPU <= 2, G10).
# Runner ../p25_runs/p25_run.py (GPU path vs ALF_GPU_DISABLE=1, same seeds).
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases_p26 import c, SCHEMES
CONT = {"Continuous": ".T."}
CASES = []
for mz in (".T.", ".F."):
    tag = "mz" if mz == ".T." else "su2"
    for k in ("std", "dly4", "sub4"):
        CASES.append(c("t3-%s-%s" % (tag, k), "works", dict(CONT, Mz=mz, **SCHEMES[k])))
    CASES.append(c("t3-%s-8x8-dly16" % tag, "works", dict(CONT, Mz=mz, L1="8", L2="8", Sequential_scheme='"delayed"', N_delay="16")))
    CASES.append(c("t3-%s-12x12" % tag, "works", dict(CONT, Mz=mz, L1="12", L2="12", NSweep="2")))
CASES += [c("hub-mz-4x4", "works"), c("hub-sub4", "works", SCHEMES["sub4"])]

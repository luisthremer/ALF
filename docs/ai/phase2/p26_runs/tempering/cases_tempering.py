# P2.6 tempering on the GPU path (TEMPERING without PARALLEL_PARAMS).  Runner ../../p25_runs/p25_run.py --groups N --mpi N
# (Temp_k: ham_U = 4 - 0.5 k, one rank per replica): GPU path vs ALF_GPU_DISABLE=1, same seeds, every Temp_k compared
# (confout byte-identical, other files incl. the tempering acceptance within the e2e rule).  Exchange every sweep.
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from cases_p26 import MODELS
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2", "N_Tempering_frequency": "1", "N_exchange_steps": "1",
     "mpi_per_parameter_set": "1"}
DLY = {"Sequential_scheme": '"delayed"', "N_delay": "4"}
SUB = {"Sequential_scheme": '"submatrix"', "N_delay": "4"}
FAULT = {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT": "temper"}
def c(label, edits=None, expect="works", env=None, note=""):
    e = dict(S); e.update(edits or {})
    return (label, e, expect, env or {}, note)
CASES = [
    c("tp-mz-std"), c("tp-mz-dly4", DLY), c("tp-mz-sub4", SUB),
    c("tp-su2-sub4", dict(SUB, Mz=".F.")),
    c("tp-mz-ltau0", {"Ltau": "0"}),
    c("tp-mz-nex3", {"N_exchange_steps": "3"}),
    c("tp-mz-nodet", {"Tempering_calc_det": ".F."}, "stops",
      note="S0-only weight is wrong for Hubbard U replicas: the CPU path stops too (unstable); both paths must stop alike"),
    c("tp-sp-sub4", dict(MODELS["sp"], **SUB), note="P2.6 model (Spin_Peierls 6x6)"),
    c("tp-mz-8x8-sub16", {"L1": "8", "L2": "8", "Sequential_scheme": '"submatrix"', "N_delay": "16"}, note="bigger (device runs)"),
    c("tp-neg-temper", {}, "FAIL", FAULT, note="negative control: stale nsigma_d after the exchange"),
]

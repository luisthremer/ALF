# P2.6 p26-langevin: Langevin updates on the device path (force pass and storage reset in gpu_langevin_mod), GPU path vs
# ALF_GPU_DISABLE=1, same seeds (runner ../p25_runs/p25_run.py).  The testsuite sets "Langevin" / "Langevin_1"
# (test_branch_parameters.json: Hubbard 6-leg ladder, Beta 4, dt 0.05, Ltau 0 / 1), square lattices Mz / SU(2), the
# Max_Force clipping branch, LRC (own S0 hooks), and 12x12 for the precision gate (Precision Green mean GPU / CPU <= 2, G10).
# The fields are continuous: round-off of the forces enters the next fields and grows with the chain (Langevin dynamics,
# measured on SU(2) 4x4: 1e-17 after one update, 2e-12 after 10, 1e-7 after 20), so confout is compared under the tolerance
# rule (P25_CONFOUT_NUMERIC=1) and the chains are short (NSweep 4 x NBin 2 = 8 updates); long chains: CPU null in the report.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases_p26 import c
LANG = {"Continuous": ".T.", "Langevin": ".T.", "Delta_t_Langevin_HMC": "0.05d0", "NSweep": "4", "NBin": "2"}
NUM = {"P25_CONFOUT_NUMERIC": "1"}
LADDER = dict(LANG, Lattice_type='"N_leg_ladder"', L1="1", L2="6", Beta="4.d0", Dtau="0.1d0")
LRC = {"ham_name": '"LRC"', "Model": '"LRC"', "Langevin": ".T.", "Delta_t_Langevin_HMC": "0.05d0", "NSweep": "4", "NBin": "2"}
CASES = [
    c("lang-ladder", "works", dict(LADDER, Ltau="0"), NUM),
    c("lang-ladder-tau", "works", dict(LADDER, Ltau="1"), NUM),
    c("lang-mz-4x4", "works", dict(LANG, Mz=".T.", Ltau="0"), NUM),
    c("lang-su2-4x4", "works", dict(LANG, Mz=".F.", Ltau="0", NSweep="2"), NUM),   # G9: 4 sweeps amplify round-off to 1.6e-10 (H100), see KNOWN_ISSUES
    c("lang-mz-4x4-tau", "works", dict(LANG, Mz=".T.", Ltau="1"), NUM),
    c("lang-maxforce", "works", dict(LANG, Mz=".T.", Ltau="0", Max_Force="0.05d0"), NUM),
    c("lang-mz-12x12", "works", dict(LANG, Mz=".T.", Ltau="0", L1="12", L2="12", Beta="2.d0", NSweep="2"), NUM),
    c("lang-su2-12x12", "works", dict(LANG, Mz=".F.", Ltau="0", L1="12", L2="12", Beta="2.d0", NSweep="2"), NUM),
    c("lrc-lang-4x4", "works", dict(LRC, Ltau="0"), NUM),
    c("lrc-lang-12x12", "works", dict(LRC, Ltau="0", L1="12", L2="12", Beta="2.d0", NSweep="2"), NUM),
]

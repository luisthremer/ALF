#!/usr/bin/env python3
"""Write production parameter files of the space-time point defect (Jakob's file, L varied, GPU settings).
usage: make_params.py OUTDIR   -> OUTDIR/parameters_L<L>_beta<B>[_MPS]  for the sizes of TUNING_DEFECT.md"""
import os, sys
T = """&VAR_ham_name
ham_name = "Hubbard"
/
&VAR_lattice
L1 = {L}
L2 = {L}
Lattice_type = "Square"
Model = "Hubbard"
/
&VAR_Model_Generic
Checkerboard = .T.
Symm = .T.
N_SUN = 2
N_FL = 2
Phi_X = 0.d0
Phi_Y = 0.d0
Bulk = .T.
N_Phi = {nphi}
Dtau = 0.1d0
Beta = {beta}.d0
Projector = .F.
Theta = 0.1d0
/
&VAR_QMC
Nwrap = {nwrap}
NSweep = 10
NBin = 1
Ltau = 0
LOBS_ST = 0
LOBS_EN = 0
CPU_MAX = 23.8
Propose_S0 = .F.
Global_moves = .F.
N_Global = 1
Global_tau_moves = .F.
N_Global_tau = 1
Sequential = .T.
Nt_sequential_start = 0
Nt_sequential_end = -1
Sequential_scheme = "submatrix"
N_delay = {ndelay}
Langevin = .F.
Delta_t_Langevin_HMC = 0.1
Max_Force = 5.0
HMC = .F.
N_HMC_sweeps = 1
Leapfrog_steps = 0
Amplitude = 1.d0
/
&VAR_errors
n_skip = 100
N_rebin = 1
N_Cov = 0
Extended_Zone = .F.
N_BZ_Zones = 1
/
&VAR_Hubbard
Mz = .F.
Continuous = .F.
ham_T = 1.d0
ham_chem = 0.d0
ham_U = 5.0d0
ham_h0 = 1.d0
ham_T2 = 1.d0
ham_U2 = 4.d0
ham_Tperp = 1.d0
/
"""
# (L, beta): Nwrap 10, N_delay 32 (L <= 28) / 64 (L >= 32); see TUNING_DEFECT.md for memory (24 GB A30 holds L <= 48 at beta = L)
cfg = [(L, L) for L in (12, 16, 20, 24, 28, 32, 36, 40, 48)] + [(48, 24), (48, 12), (56, 28), (64, 32), (64, 16)]
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
for L, b in cfg:
    nd = 32 if L <= 28 else 64
    open(os.path.join(out, f"parameters_L{L}_beta{b}"), "w").write(T.format(L=L, nphi=L*L//2, beta=b, nwrap=10, ndelay=nd))
print(len(cfg), "files in", out)

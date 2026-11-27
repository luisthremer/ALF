# Cases of the P2.5 evidence run, GPU=host, serial build (configure.sh GNU noMPI GPU=host).  See p25_run.py for the columns.
# (label, namelist edits, expectation from GPU_FEATURES.md, extra env, note)
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2"}          # small default: Hubbard 4x4, Beta = 5, Dtau = 0.1 (Ltrot = 50)
def c(label, expect, edits=None, env=None, note=""):
    e = dict(S); e.update(edits or {})
    return (label, e, expect, env or {}, note)

CASES = [
 # --- Hubbard, the supported class
 c("hub-mz-4x4", "works"),
 c("hub-mz-6x6-20sw", "works", {"L1": "6", "L2": "6", "NSweep": "20", "NBin": "5"}, note="Start/parameters as shipped (20 sweeps, 5 bins)"),
 c("hub-su2-6x6-20sw", "works", {"Mz": ".F.", "L1": "6", "L2": "6", "NSweep": "20", "NBin": "5"}),
 c("hub-mz-6x6", "works", {"L1": "6", "L2": "6"}),
 c("hub-su2-4x4", "works", {"Mz": ".F."}),
 c("hub-su2-6x6", "works", {"Mz": ".F.", "L1": "6", "L2": "6"}),
 c("hub-nsun4-su2", "works", {"Mz": ".F.", "N_SUN": "4"}, note="N_SUN row"),
 c("hub-nsun4-mz", "works", {"N_SUN": "4"}, note="N_SUN row"),
 c("hub-complex-phiX", "works", {"Mz": ".F.", "Phi_X": "0.25d0"}, note="complex hopping"),
 c("hub-complex-flux", "works", {"Mz": ".F.", "N_Phi": "1", "Checkerboard": ".F."}, note="complex hopping, flux"),
 c("hub-ladder", "works", {"Lattice_type": '"N_leg_ladder"', "L1": "6", "L2": "2", "ham_T2": "0.d0", "ham_U2": "0.d0", "ham_Tperp": "0.d0"}),
 c("hub-bilayer-square", "works", {"Lattice_type": '"Bilayer_square"', "Checkerboard": ".F.", "Mz": ".F.", "ham_U": "0.d0", "ham_U2": "4.d0"}),
 c("hub-h0-pinning", "works", {"ham_h0": "0.5d0"}),
 c("hub-symm-F", "works", {"Mz": ".F.", "Symm": ".F."}),
 c("hub-checkerboard-F", "works", {"Checkerboard": ".F."}),
 c("hub-ltau0", "works", {"Ltau": "0"}),
 c("hub-lobs-window", "works", {"LOBS_ST": "10", "LOBS_EN": "40"}, note="LOBS_ST/LOBS_EN window"),
 # --- update schemes
 c("sch-delayed-N4", "works", {"Sequential_scheme": '"delayed"', "N_delay": "4"}),
 c("sch-submatrix-N4", "works", {"Sequential_scheme": '"submatrix"', "N_delay": "4"}),
 c("sch-delayed-N16-su2", "works", {"Mz": ".F.", "Sequential_scheme": '"delayed"', "N_delay": "16", "L1": "6", "L2": "6"}),
 c("sch-auto", "refused:auto", {"Sequential_scheme": '"auto"'}),
 c("sch-timed-N-1", "refused:N_delay = -1", {"Sequential_scheme": '"delayed"', "N_delay": "-1"}),
 c("sch-lowrank-cross-slice", "stops", {"Lowrank_cross_slice": ".T."}, note="option does not exist on this line"),
 # --- moves that are refused
 c("mv-global", "refused:global moves", {"Continuous": ".T.", "Global_moves": ".T.", "Global_tau_moves": ".T.", "N_Global": "2", "N_Global_tau": "1", "Nt_sequential_start": "1", "Nt_sequential_end": "6", "Lattice_type": '"N_leg_ladder"', "L1": "1", "L2": "6", "Beta": "4.d0", "ham_Tperp": "1.d0"}, note="continuous fields (discrete Hubbard stops on every path)"),
 c("mv-global-discrete", "stops", {"Global_moves": ".T.", "N_Global": "1"}, note="Hubbard has no discrete global move: model stops, CPU and GPU path alike"),
 c("mv-global-tau", "refused:global moves", {"Global_tau_moves": ".T.", "N_Global_tau": "1"}),
 c("mv-propose-S0", "refused:Propose_S0", {"Propose_S0": ".T."}),
 c("mv-langevin", "refused:Langevin", {"Continuous": ".T.", "Langevin": ".T.", "Sequential": ".F.", "Mz": ".F."}),
 c("mv-hmc", "refused:Langevin", {"Continuous": ".T.", "HMC": ".T.", "Sequential": ".F.", "Leapfrog_steps": "5", "Mz": ".F."}),
 c("mv-mala-seq", "refused:MALA", {"Continuous": ".T.", "Sequential_MALA": ".T.", "Mz": ".F."}),
 c("mv-mala-global", "refused:MALA", {"Continuous": ".T.", "Global_MALA_moves": ".T.", "Sequential": ".F.", "Mz": ".F."}),
 c("mv-mala-global-tau", "refused:MALA", {"Continuous": ".T.", "Global_tau_MALA_moves": ".T.", "Mz": ".F."}),
 c("mv-continuous-seq", "refused:field type", {"Continuous": ".T.", "Mz": ".F."}, note="HS field type 3 with sequential updates"),
 c("proj-hubbard", "works", {"Projector": ".T.", "Theta": "2.d0", "Beta": "2.d0"}),
 # --- stock models (all other than Hubbard)
 c("mod-plain-vanilla", "refused:Calc_Fl", {"ham_name": '"Hubbard_Plain_Vanilla"', "Model": '"Hubbard_Plain_Vanilla"', "Beta": "1.d0", "Dtau": "0.05d0", "Projector": ".F."}),
 c("mod-plain-vanilla-proj", "works", {"ham_name": '"Hubbard_Plain_Vanilla"', "Model": '"Hubbard_Plain_Vanilla"', "Beta": "1.d0", "Dtau": "0.05d0", "Projector": ".T.", "Ltau": "0"}),
 c("mod-kondo", "refused", {"ham_name": '"Kondo"', "Model": '"Kondo"', "Lattice_type": '"Bilayer_square"', "ham_Uf": "1.d0"}),
 c("mod-tV", "refused", {"ham_name": '"tV"', "Model": '"tV"', "N_SUN": "1", "ham_V": "-1.d0", "ham_T2": "0.d0", "ham_Tperp": "0.d0", "Beta": "2.d0", "Ltau": "0"}),
 c("mod-Z2-matter", "refused", {"ham_name": '"Z2_Matter"', "Model": '"Z2_Matter"', "Beta": "2.d0", "Dtau": "0.05d0", "ham_T": "0.d0", "ham_h": "0.d0", "ham_TZ2": "1.d0", "ham_K": "0.d0", "ham_J": "0.d0", "ham_g": "1.d0", "ham_U": "0.d0", "Ltau": "0", "Nwrap": "10"}),
 c("mod-LRC", "refused", {"ham_name": '"LRC"', "Model": '"LRC"', "ham_alpha": "0.d0", "Ltau": "0"}),
 c("mod-spin-peierls", "refused", {"ham_name": '"Spin_Peierls"', "Model": '"Spin_Peierls"', "L1": "6", "L2": "6", "Beta": "2.d0", "ham_U": "1.5d0", "ham_Lambda": "0.1d0", "ham_omega0": "0.5d0", "N_SUN": "2", "N_FL": "1"}),
]

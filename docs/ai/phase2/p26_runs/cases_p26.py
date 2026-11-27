# Cases of P2.6.1 (delayed / submatrix for every model), GPU=host, serial build.  Runner: ../p25_runs/p25_run.py (same columns:
# GPU path ACTIVE vs ALF_GPU_DISABLE=1, same seeds, confout byte-identical, other files within the e2e rule).
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2"}
def c(label, expect, edits=None, env=None, note=""):
    e = dict(S); e.update(edits or {})
    return (label, e, expect, env or {}, note)

MODELS = {
 "kondo": {"ham_name": '"Kondo"', "Model": '"Kondo"', "Lattice_type": '"Bilayer_square"', "ham_Uf": "1.d0"},
 "tV":    {"ham_name": '"tV"', "Model": '"tV"', "N_SUN": "1", "ham_V": "-1.d0", "ham_T2": "0.d0", "ham_Tperp": "0.d0", "Beta": "2.d0", "Ltau": "0"},
 "sp":    {"ham_name": '"Spin_Peierls"', "Model": '"Spin_Peierls"', "L1": "6", "L2": "6", "Beta": "2.d0", "ham_U": "1.5d0", "ham_Lambda": "0.1d0",
           "ham_omega0": "0.5d0", "N_SUN": "2", "N_FL": "1"},
}
SCHEMES = {
 "std":    {},
 "dly4":   {"Sequential_scheme": '"delayed"', "N_delay": "4"},
 "dly16":  {"Sequential_scheme": '"delayed"', "N_delay": "16"},
 "sub4":   {"Sequential_scheme": '"submatrix"', "N_delay": "4"},
}
CASES = []
for m, me in MODELS.items():
    for k, ke in SCHEMES.items():
        e = dict(me); e.update(ke)
        CASES.append(c("%s-%s" % (m, k), "works", e))
# tV with the tau pass (Ltau = 1)
CASES.append(c("tV-std-tau", "works", dict(MODELS["tV"], Ltau="1")))
# Hubbard regression (rank-1 diagonal, must stay as before)
CASES += [c("hub-mz-4x4", "works"),
          c("hub-dly4", "works", SCHEMES["dly4"]),
          c("hub-sub4", "works", SCHEMES["sub4"])]

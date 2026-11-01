# SuperMUC-NG N_delay scan (user decision 2026-11-01: find the correct N_delay and what it depends on).
# For every parameter point: standard once, delayed and submatrix at fixed N_delay 8..256, automatic (0) and automatic+timing (-1);
# 3 repetitions (scheme order rotated per repetition by the generator); one run per node (packed jobs, exclusive nodes), so no two
# runs share a node and no co-runner can bias one N_delay value. Same seed for all runs of a point => same Markov chain, identical
# work per sweep; differences are timing only. Sweep time from the timers (bins >= 2).
# Varied parameters (one at a time around the base point Hubbard U=4, Beta=2 (Ltrot 20), L=40, 2x24):
#   system size L = 20, 30, 40, 50 | time slices Beta = 1, 2, 4 (Ltrot 10, 20, 40) | U = 1, 4, 8 (acceptance) |
#   vertex rank / model: tV (rank-2 bonds) L = 20, 30, 40; Kondo (bilayer, rank-2 J vertices) L = 20 |
#   threads per rank: layouts 1x48, 2x24, 4x12 at the base point.
# sec_per_sweep_est = standard scheme's value at 2x24 (conservative for all schemes); from pilot stages S/A and B2.
_BASE = {"Dtau": "0.1d0", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"}
_TV = {"ham_name": "\"tV\"", "Model": "\"tV\"", "ham_V": "-1.d0", "N_SUN": "1"}
_KONDO = {"ham_name": "\"Kondo\"", "Model": "\"Kondo\"", "Lattice_type": "\"Bilayer_square\""}
_ND = [8, 16, 32, 64, 128, 256, 0, -1]
_SCH = ["standard", "delayed", "submatrix"]

def _pt(name, L, est, nsweep, nbin, U="4.d0", beta="2.d0", nwrap="10", layouts=("2x24",), model="Hubbard", extra=None,
        n_fl=2, setup=60):
    k = dict(_BASE); k.update({"Beta": beta, "Nwrap": nwrap})
    if model == "Hubbard":
        k["ham_U"] = U
    if extra:
        k.update(extra)
    return {"name": name, "model": model, "keys": k, "L": [L], "schemes": _SCH, "n_delay": _ND, "nsweep": nsweep, "nbin": nbin,
            "layouts": list(layouts), "n_fl": n_fl, "sec_per_sweep_est": est, "setup_s_est": setup}

CAMPAIGN = {
    "id": "smuc_ndscan",
    "profile": "supermuc",
    "partition": "test",
    "binary": "/hppfs/work/pn73xu/di38zuw/CLAUDE/alf_bin/ALF_timers_mpi_78bc38ee.out",
    "template": "/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_DELAYED_SUBMATRIX/ALF/runs/campaigns/templates/78bc38ee",
    "fakeintel": "off",
    "repetitions": 3,
    "seed": 790789,
    "max_jobs_in_queue": 3,
    "max_attempts": 2,
    "wall_safety": 1.15,
    "pack_nodes": 16,
    "points": [
        # system size (base U=4, Beta=2)
        _pt("hubL20", 20, 2.0, 20, 4),
        _pt("hubL30", 30, 15.0, 3, 4),
        _pt("hubL40", 40, 70.0, 1, 4, layouts=("2x24", "1x48", "4x12")),
        _pt("hubL50", 50, 260.0, 1, 3, setup=90),
        # time slices at L=40
        _pt("hubL40b1", 40, 36.0, 2, 4, beta="1.d0"),
        _pt("hubL40b4", 40, 130.0, 1, 3, beta="4.d0"),
        # interaction / acceptance at L=40
        _pt("hubL40U1", 40, 70.0, 1, 4, U="1.d0"),
        _pt("hubL40U8", 40, 80.0, 1, 4, U="8.d0", nwrap="5"),
        # vertex rank / model
        _pt("tvL20", 20, 2.0, 20, 4, model="tV", extra=_TV, n_fl=1),
        _pt("tvL30", 30, 15.0, 3, 4, model="tV", extra=_TV, n_fl=1),
        _pt("tvL40", 40, 70.0, 1, 4, model="tV", extra=_TV, n_fl=1),
        _pt("kondoL20", 20, 12.0, 4, 4, model="Kondo", extra=_KONDO, n_fl=2),
    ],
}

# SuperMUC-NG pilot C0, stage A: paired large-N precision + first timings (16 nodes, one run per node, test partition).
# Same seed for every run => the same Markov chain for every scheme (trajectory tier), so Precision Green/Phase are paired;
# fixed sweep counts (no CPU_MAX) so every scheme of a point covers the same sweeps.
# Grid adapted to the 30-min limit (stage S, job 5956035: standard 2x24 needs ~123 s/sweep at L=40 Beta=4 and ~690 s/sweep at
# L=60 Beta=2): L=60 Beta=4 does not fit; L=56 (N=3136, auto n0: delayed 128, submatrix 256) and L=48 (N=2304, n0 128/128)
# at Beta=2 instead. U=8 with Nwrap=5 (the Mac pre-check showed Precision Green max 8.8e-2 for standard with Nwrap=10 at U=8).
_KEYS = {"Dtau": "0.1d0", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"}

def _pt(name, L, U, beta, nwrap, nbin, est, setup, schemes, n_delay):
    k = dict(_KEYS); k.update({"ham_U": U, "Beta": beta, "Nwrap": nwrap})
    return {"name": name, "model": "Hubbard", "keys": k, "L": [L], "schemes": schemes, "n_delay": n_delay,
            "nsweep": 1, "nbin": nbin, "layouts": ["2x24"], "n_fl": 2, "sec_per_sweep_est": est, "setup_s_est": setup}

_ALL = ["standard", "delayed", "submatrix"]
CAMPAIGN = {
    "id": "smuc_c0_A",
    "profile": "supermuc",
    "partition": "test",
    "binary": "/hppfs/work/pn73xu/di38zuw/CLAUDE/alf_bin/ALF_timers_mpi_78bc38ee.out",
    "template": "/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_DELAYED_SUBMATRIX/ALF/runs/campaigns/templates/78bc38ee",
    "fakeintel": "off",
    "pilot": True,
    "repetitions": 1,
    "seed": 790789,
    "max_jobs_in_queue": 1,
    "max_attempts": 1,
    "wall_safety": 1.1,
    "pack_nodes": 16,
    "points": [
        # sec_per_sweep_est = the standard scheme's (largest) value, so every slot's estimate is conservative
        _pt("hubL40b4U4", 40, "4.d0", "4.d0", "10", 8, 125, 60, _ALL, [0, 32]),
        _pt("hubL40b4U8", 40, "8.d0", "4.d0", "5", 6, 150, 60, _ALL, [0, 32]),
        _pt("hubL56b2U4", 56, "4.d0", "2.d0", "10", 2, 470, 120, _ALL, [0]),
        _pt("hubL48b2U8", 48, "8.d0", "2.d0", "5", 4, 260, 90, _ALL, [0]),
    ],
}

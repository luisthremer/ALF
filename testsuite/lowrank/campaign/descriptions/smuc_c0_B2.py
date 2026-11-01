# SuperMUC-NG pilot C0, stage B2: seconds-per-sweep calibration of the later campaign points (10 nodes, one run per node, 2x24).
# Runs stop gracefully via CPU_MAX (cpu_max_h); sec/sweep from the sweeps actually done (collect.py, truncated runs).
# (3D cubic points are not calibrated: this ALF version has no cubic lattice.)
_KEYS = {"Dtau": "0.1d0", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1", "Beta": "2.d0", "Nwrap": "10"}
_TV = {"ham_name": "\"tV\"", "Model": "\"tV\"", "ham_V": "-1.d0", "N_SUN": "1"}

def _pt(name, L, schemes, extra, model="Hubbard", setup=60):
    k = dict(_KEYS); k.update(extra)
    return {"name": name, "model": model, "keys": k, "L": [L], "schemes": schemes, "n_delay": [0],
            "nsweep": 1, "cpu_max_h": 0.2, "layouts": ["2x24"], "n_fl": 2 if model == "Hubbard" else 1, "setup_s_est": setup}

CAMPAIGN = {
    "id": "smuc_c0_B2",
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
        _pt("hubL20b2U4", 20, ["standard", "delayed"], {"ham_U": "4.d0"}),
        _pt("hubL30b2U4", 30, ["standard", "delayed", "submatrix"], {"ham_U": "4.d0"}),
        _pt("hubL50b2U4", 50, ["standard", "delayed", "submatrix"], {"ham_U": "4.d0"}, setup=90),
        _pt("tvL30b2", 30, ["standard", "delayed"], _TV, model="tV"),
    ],
}

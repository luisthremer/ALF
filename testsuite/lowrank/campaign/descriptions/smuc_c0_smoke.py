# SuperMUC-NG pilot, stage S (part of the approved pilot C0): first packed test-partition job (4 nodes, one run per node).
# Purpose: (1) verify the packed launch on the cluster (placement, pinning, node-local provenance, overlapping slots);
# (2) measure seconds per sweep of the two most expensive pilot points with standard and automatic delayed, to size the fixed
# sweep counts of stage A (paired precision runs need the same number of sweeps for every scheme).
# Runs stop gracefully via ALF's CPU_MAX (cpu_max_h); one sweep per bin, so the per-bin timers give seconds per sweep.
CAMPAIGN = {
    "id": "smuc_c0_smoke",
    "profile": "supermuc",
    "partition": "test",
    "binary": "/hppfs/work/pn73xu/di38zuw/CLAUDE/alf_bin/ALF_timers_mpi_78bc38ee.out",
    "template": "/Users/luis/Documents/_WORK/Coding/CLAUDE_ALF_DELAYED_SUBMATRIX/ALF/runs/campaigns/templates/78bc38ee",   # = git show 78bc38ee:Scripts_and_Parameters_files/Start/*
    "fakeintel": "off",
    "pilot": True,
    "repetitions": 1,
    "seed": 790789,
    "max_jobs_in_queue": 1,
    "max_attempts": 1,
    "wall_safety": 1.1,
    "pack_nodes": 4,
    "points": [
        {
            "name": "hubL60b2",
            "model": "Hubbard",
            "keys": {"ham_U": "4.d0", "Beta": "2.d0", "Dtau": "0.1d0", "Nwrap": "10", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"},
            "L": [60],
            "schemes": ["standard", "delayed"],
            "n_delay": [0],
            "nsweep": 1,
            "cpu_max_h": 0.2,
            "layouts": ["2x24"],
            "n_fl": 2,
            "setup_s_est": 120,
        },
        {
            "name": "hubL40b4",
            "model": "Hubbard",
            "keys": {"ham_U": "4.d0", "Beta": "4.d0", "Dtau": "0.1d0", "Nwrap": "10", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"},
            "L": [40],
            "schemes": ["standard", "delayed"],
            "n_delay": [0],
            "nsweep": 1,
            "cpu_max_h": 0.2,
            "layouts": ["2x24"],
            "n_fl": 2,
            "setup_s_est": 60,
        },
    ],
}

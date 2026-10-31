# Mini campaign of the local end-to-end test (run_local_test.sh): Hubbard L = 8, standard and delayed (N_delay 8),
# layouts 1x1 and 2x2 (MPI binary), 2 repetitions. The test driver sets the environment variables below.
import os

CAMPAIGN = {
    "id": os.environ["CAMPAIGN_TEST_ID"],
    "profile": "local",
    "partition": "local",
    "binary": os.environ["CAMPAIGN_TEST_BIN"],
    "fakeintel": os.environ.get("CAMPAIGN_TEST_FAKEINTEL", "off"),
    "repetitions": int(os.environ.get("CAMPAIGN_TEST_REPS", "2")),
    "seed": 790789,
    "max_jobs_in_queue": 4,
    "max_attempts": 2,
    "wall_safety": 1.5,
    "points": [{
        "name": "hub",
        "model": "Hubbard",
        "keys": {"ham_U": "4.d0", "Beta": "2.d0", "Dtau": "0.1d0", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"},
        "L": [8],
        "schemes": os.environ.get("CAMPAIGN_TEST_SCHEMES", "standard delayed").split(),
        "n_delay": [8],
        "nsweep": 20,
        "nbin": 3,
        "layouts": os.environ.get("CAMPAIGN_TEST_LAYOUTS", "1x1 2x2").split(),
        "partial_node": True,
        "n_fl": 2,
        "sec_per_sweep_est": 0.03,
        "setup_s_est": 2,
    }],
}

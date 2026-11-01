# Mini campaign of the packed-job test (run_pack_test.sh): Hubbard L = 6 and 8, standard + delayed + submatrix (N_delay 8),
# layout 1x1 (noMPI binary), 2 repetitions, pack_nodes = 3 on the local profile (partition local: 60 min, 4 nodes).
# The wall ESTIMATES (setup_s_est, not the real run time) are chosen so that the slots are fixed by the first-fit-decreasing
# rule: a run of L = 8 is estimated at 901.8 s, of L = 6 at 301.8 s; a slot holds at most 2080 s (60 min limit, safety 1.5,
# 120 s start-up, 5 min margin) => two L = 8 runs per slot (3 slots for the 6 runs), the six L = 6 runs in a fourth slot
# => 4 slots, cut into jobs of 3 + 1. The test driver sets the environment variables below.
import os

CAMPAIGN = {
    "id": os.environ["CAMPAIGN_TEST_ID"],
    "profile": os.environ.get("CAMPAIGN_TEST_PROFILE", "local"),
    "partition": "local",
    "binary": os.environ["CAMPAIGN_TEST_BIN"],
    "repetitions": int(os.environ.get("CAMPAIGN_TEST_REPS", "2")),
    "seed": 790789,
    "max_jobs_in_queue": 4,
    "max_attempts": 2,
    "wall_safety": 1.5,
    "pack_nodes": int(os.environ.get("CAMPAIGN_TEST_PACK", "3")),
    "points": [{
        "name": "hub",
        "model": "Hubbard",
        "keys": {"ham_U": "4.d0", "Beta": "2.d0", "Dtau": "0.1d0", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"},
        "L": [6, 8],
        "schemes": ["standard", "delayed", "submatrix"],
        "n_delay": [8],
        "nsweep": 20,
        "nbin": 3,
        "layouts": ["1x1"],
        "partial_node": True,
        "n_fl": 2,
        "sec_per_sweep_est": 0.03,
        "setup_s_est": {6: float(os.environ.get("CAMPAIGN_TEST_SETUP6", "300")), 8: float(os.environ.get("CAMPAIGN_TEST_SETUP8", "900"))},
    }],
}

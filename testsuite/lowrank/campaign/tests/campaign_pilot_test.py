# Mini campaign of the pilot / cpu_max_h tests (run_pack_test.sh, section "pilot and cpu_max_h"): Hubbard L = 6, standard +
# delayed (N_delay 8), layout 1x1, ONE repetition, pilot = True, on the local profile packed into 2 slots (the estimate of a run is
# setup_s_est 1100 s + the cpu_max_h budget; two of them exceed the 60 min limit), with cpu_max_h = 0.0008 h (2.88 s):
# ALF ignores nbin (set to 10^7 internally), runs bins until the budget is used and stops gracefully after a bin.
# nbin = 100000 would make the nominal estimate (60000 s) exceed the 60 min limit: only the cpu_max_h cap makes it packable.
# The test driver sets the environment variables below.
import os

CAMPAIGN = {
    "id": os.environ["CAMPAIGN_TEST_ID"],
    "profile": os.environ.get("CAMPAIGN_TEST_PROFILE", "local"),
    "partition": os.environ.get("CAMPAIGN_TEST_PARTITION", "local"),
    "binary": os.environ["CAMPAIGN_TEST_BIN"],
    "repetitions": int(os.environ.get("CAMPAIGN_TEST_REPS", "1")),
    "pilot": os.environ.get("CAMPAIGN_TEST_PILOT", "1") == "1",
    "seed": 790789,
    "max_jobs_in_queue": 4,
    "max_attempts": 2,
    "wall_safety": 1.5,
    "pack_nodes": int(os.environ.get("CAMPAIGN_TEST_PACK", "2")),
    "points": [{
        "name": "hub",
        "model": "Hubbard",
        "keys": {"ham_U": "4.d0", "Beta": "2.d0", "Dtau": "0.1d0", "Ltau": "0", "LOBS_ST": "1", "LOBS_EN": "1"},
        "L": [6],
        "schemes": ["standard", "delayed"],
        "n_delay": [8],
        "nsweep": 20,
        "nbin": 100000,
        "layouts": ["1x1"],
        "partial_node": True,
        "n_fl": 2,
        "sec_per_sweep_est": 0.03,
        "setup_s_est": float(os.environ.get("CAMPAIGN_TEST_SETUP", "1100")),
        "cpu_max_h": float(os.environ.get("CAMPAIGN_TEST_CPU_MAX_H", "0.0008")),
    }],
}

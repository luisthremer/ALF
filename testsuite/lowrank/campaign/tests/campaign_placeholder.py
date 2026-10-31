# A cluster campaign: refused while the profile has placeholders; with a filled copy of the profile (run_local_test.sh)
# it is refused for repetitions < 3 and accepted (--dry-run) with 3.
import os

CAMPAIGN = {
    "id": "placeholder_check",
    "profile": os.environ.get("CAMPAIGN_TEST_PROFILE", "supermuc"),
    "partition": "test",
    "binary": "/nonexistent/ALF.out",
    "repetitions": int(os.environ.get("CAMPAIGN_TEST_REPS", "3")),
    "seed": 1,
    "points": [{"name": "hub", "model": "Hubbard", "keys": {"Beta": "2.d0", "Dtau": "0.1d0"}, "L": [12],
                "schemes": ["standard", "delayed"], "n_delay": [32], "nsweep": 10, "nbin": 3, "layouts": ["48x1"],
                "sec_per_sweep_est": 0.1, "setup_s_est": 5}],
}

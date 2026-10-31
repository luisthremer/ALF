# A campaign for a cluster profile whose placeholders are not filled: make_campaign.py must refuse it (exit 1).
import os

CAMPAIGN = {
    "id": "placeholder_check",
    "profile": os.environ.get("CAMPAIGN_TEST_PROFILE", "supermuc"),
    "partition": "test",
    "binary": "/nonexistent/ALF.out",
    "repetitions": 3,
    "seed": 1,
    "points": [{"name": "hub", "model": "Hubbard", "keys": {"Beta": "2.d0", "Dtau": "0.1d0"}, "L": [12],
                "schemes": ["standard", "delayed"], "n_delay": [32], "nsweep": 10, "nbin": 3, "layouts": ["48x1"],
                "sec_per_sweep_est": 0.1, "setup_s_est": 5}],
}

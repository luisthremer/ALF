# TEMPERING build without PARALLEL_PARAMS (real replica exchange): refused at compile time; two groups, one rank each
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2", "N_exchange_steps": "2", "N_Tempering_frequency": "2"}
CASES = [("temper-mz-4x4", dict(S), "refused:tempering", {}, ""), ("temper-su2-4x4-delayed", dict(S, Mz=".F.", Sequential_scheme='"delayed"', N_delay="4"), "refused:tempering", {}, "")]

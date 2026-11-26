# PARALLEL_PARAMS / TEMPERING builds: two groups (Temp_0: U = 4, Temp_1: U = 3.5), one rank per group
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2", "N_exchange_steps": "2", "N_Tempering_frequency": "2"}
CASES = [("grp-mz-4x4", dict(S), "works", {}, ""), ("grp-su2-4x4-delayed", dict(S, Mz=".F.", Sequential_scheme='"delayed"', N_delay="4"), "works", {}, "")]

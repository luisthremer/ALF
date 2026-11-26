# restart is handled by --restart in p25_run.py (second run per directory starts from confout -> confin)
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2"}
CASES = [("restart-mz", dict(S), "works", {}, "restart"), ("restart-su2-delayed", dict(S, Mz=".F.", Sequential_scheme='"delayed"', N_delay="4"), "works", {}, "restart")]

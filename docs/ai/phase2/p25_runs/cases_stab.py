# STAB1/2/3 builds: refused at compile time for the whole build (no restart-free "works" case exists)
S = {"L1": "4", "L2": "4", "NSweep": "4", "NBin": "2"}
CASES = [("stab-mz-4x4", dict(S), "refused:STAB", {}, ""), ("stab-su2-4x4", dict(S, Mz=".F."), "refused:STAB", {}, "")]

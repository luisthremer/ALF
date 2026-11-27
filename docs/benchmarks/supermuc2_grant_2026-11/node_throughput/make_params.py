#!/usr/bin/env python3
"""Write an ALF parameters file for the node throughput test.

usage: make_params.py TEMPLATE L NRANKS > parameters

TEMPLATE is the production file RUN_L24/Temp_0/parameters (jakob_params, h0 = 1, Beta = L).
Changed here (owner, 2026-11-27): L1 = L2 = L, N_Phi = L^2/2, Beta = 4, Ltau = 0, ham_h0 = 0,
ham_U = 5.25, NSweep = 1, NBin = 11, CPU_MAX = 0, mpi_per_parameter_set = NRANKS (one parameter set).
Mz = .F., Symm = .T., Checkerboard = .T., N_FL = 2, Nwrap 10, submatrix, N_delay 64 stay as in the template.
"""
import re
import sys

tmpl, L, nr = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
new = {
    "L1": str(L), "L2": str(L), "N_Phi": str(L * L // 2), "Beta": "4.d0", "Ltau": "0",
    "ham_h0": "0.d0", "ham_U": "5.25d0", "NSweep": "1", "NBin": "11", "CPU_MAX": "0.0",
    "mpi_per_parameter_set": str(nr),
    # kept, but asserted, because the test depends on them
    "Mz": ".F.", "Symm": ".T.", "Checkerboard": ".T.", "N_FL": "2", "Dtau": "0.1d0",
}
seen = set()
out = []
for line in open(tmpl):
    m = re.match(r"^(\s*)(\w+)(\s*=\s*)(.*)$", line.rstrip("\n"))
    if m and m.group(2) in new:
        out.append(f"{m.group(1)}{m.group(2)}{m.group(3)}{new[m.group(2)]}")
        seen.add(m.group(2))
    else:
        out.append(line.rstrip("\n"))
missing = set(new) - seen
if missing:
    sys.exit(f"make_params.py: keys not found in template: {sorted(missing)}")
print("\n".join(out))

"""Hubbard reference configurations for the GPU port (docs/gpu/BASELINE.md).

Square lattice, half filling, finite T, sequential local updates, Ltau=1.
Each config is one pyALF sim_dict; names encode size and hopping variant.
"""

COMMON = {
    "Lattice_type": "Square",
    "Model": "Hubbard",
    "Symm": False,
    "ham_T": 1.0,
    "ham_U": 4.0,
    "ham_chem": 0.0,
    "Mz": True,
    "beta": 4.0,
    "dtau": 0.1,
    "Nwrap": 10,
    "Ltau": 1,
    "n_skip": 2,
}

# (NSweep, NBin) per rank, sized from measured sweep times on nividia (GNU, 16 ranks)
# so that one run stays below ~12 min; see docs/gpu/BASELINE.md.
SIZES = {
    (6, True): (50, 100), (6, False): (50, 100),
    (12, True): (20, 30), (12, False): (20, 30),
    (20, True): (6, 10), (20, False): (4, 10),
}

CONFIGS = {}
for (L, cb), (nsweep, nbin) in SIZES.items():
    name = f"hub{L}x{L}_{'cb' if cb else 'dense'}"
    CONFIGS[name] = {**COMMON, "L1": L, "L2": L, "Checkerboard": cb,
                     "NSweep": nsweep, "NBin": nbin}

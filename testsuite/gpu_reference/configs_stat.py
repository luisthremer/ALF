"""Short-autocorrelation configurations for validating compare2.py (ticket 106).

8x8 square Hubbard, half filling, beta=2, dtau=0.1, Ltau=1 (tau-resolved data needed for compare2.py).
At U<=4, beta=2 the AF mode decorrelates within about a bin, so many bins give honest errors.
NSweep/NBin are per MPI rank (run with 4 ranks); sized from timing on nividia, see README.md.
Names: s<L>_U<U>_b<beta>_<cb|dense>; 'p' stands for the decimal point.
"""
from configs import COMMON

NSWEEP, NBIN = 40, 250  # ~0.05 s/sweep/rank (GNU, 4 ranks): ~8 min per run, ifx ~1.3x

CONFIGS = {}


def _add(L, U, beta, cb, nsweep=NSWEEP, nbin=NBIN):
    ustr = f"{U:g}".replace(".", "p")
    name = f"s{L}_U{ustr}_b{int(beta)}_{'cb' if cb else 'dense'}"
    CONFIGS[name] = {**COMMON, "L1": L, "L2": L, "ham_U": float(U), "beta": float(beta),
                     "Checkerboard": cb, "NSweep": nsweep, "NBin": nbin}


_add(8, 2.0, 2, True)
_add(8, 4.0, 2, True)
_add(8, 4.0, 2, False)   # must-fail partner of the previous one: cb vs dense (Trotter difference)
_add(8, 2.1, 2, True)    # must-fail partner of U=2: small physical shift

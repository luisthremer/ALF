#!/usr/bin/env python3
"""Short deterministic end-to-end checks of fork features with the normal ALF executable (serial, plain-text output).

Scenario (argument 2):
  schemes   Sequential_scheme = standard / delayed / submatrix on the standard Hubbard model (square L = 4, N_FL = 2,
            no pinning field), same seed. The delayed and submatrix schemes are exact rewrites of the rank-1 update (only
            the summation order differs, docs/gpu/DELAYED_SUBMATRIX.md), so the Markov chain is the same: the final
            configurations must be identical and every observable must agree to rounding (RTOL below).
"""

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile

# Observable agreement between update schemes. All three schemes compute the same G up to rounding (Precision Green
# of the runs ~1e-13 absolute); an observable is a bin average of products of G entries, so scheme differences are a few
# orders of magnitude above machine epsilon. Measured (gfortran 16, arm64, OpenBLAS): max relative deviation 4.5e-13
# (delayed) and 3.4e-13 (submatrix). RTOL = 1e-8 leaves 4 orders of margin for other compilers/BLAS; planted update
# errors of relative size 1e-6 (flush or M-update coefficient 1 + 1e-6) give ~2e-4 and are caught.
RTOL = 1.0e-8
ATOL = 1.0e-10     # for entries that vanish (e.g. imaginary parts ~1e-16)


def parameters(*, scheme="standard", n_delay=4, ltau=0, h0=0.0, sweeps=20, bins=2):
    return f"""&VAR_ham_name
ham_name = "Hubbard"
/
&VAR_lattice
L1 = 4, L2 = 4, Lattice_type = "Square", Model = "Hubbard"
/
&VAR_Model_Generic
Checkerboard = .T., Symm = .T., N_SUN = 2, N_FL = 2, Phi_X = 0.d0, Phi_Y = 0.d0, Bulk = .T., N_Phi = 0
Dtau = 0.1d0, Beta = 2.d0, Projector = .F., Theta = 0.1d0
/
&VAR_QMC
Nwrap = 10, NSweep = {sweeps}, NBin = {bins}, Ltau = {ltau}, LOBS_ST = 0, LOBS_EN = 0, CPU_MAX = 0.0
Sequential = .T., Sequential_scheme = "{scheme}", N_delay = {n_delay}
/
&VAR_errors
n_skip = 0, N_rebin = 1, N_Cov = 0
/
&VAR_Hubbard
Mz = .F., Continuous = .F., ham_T = 1.d0, ham_chem = 0.d0, ham_U = 4.d0, ham_h0 = {h0:.1f}d0
/
"""


def run(binary, directory, **options):
    directory.mkdir(parents=True)
    (directory / "parameters").write_text(parameters(**options))
    (directory / "seeds").write_text("8752143\n")
    with (directory / "stdout").open("w") as log:
        result = subprocess.run([str(binary)], cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                timeout=60, check=False)
    if result.returncode:
        raise AssertionError(f"ALF failed in {directory}:\n" + (directory / "stdout").read_text())
    return directory


def observable_files(directory):
    skip = {"parameters", "seeds", "stdout", "info", "info_git", "confout_0"}
    files = sorted(p.name for p in directory.iterdir() if p.is_file() and p.name not in skip
                   and not p.name.endswith("_info"))
    assert files, f"no observable files in {directory}"
    return files


def numbers(path):
    text = path.read_text()
    assert text.strip(), f"empty observable file: {path}"
    values = [float(t.replace("D", "E").replace("d", "e")) for t in re.split(r"[\s(),]+", text.strip()) if t]
    assert all(math.isfinite(v) for v in values), f"non-finite output: {path}"
    return values


def config(directory):
    return (directory / "confout_0").read_text().split()


def schemes(binary, root):
    ref = run(binary, root / "standard", scheme="standard")
    files = observable_files(ref)
    failures = []
    for scheme in ("delayed", "submatrix"):
        d = run(binary, root / scheme, scheme=scheme)
        info = (d / "info").read_text()
        if f"Sequential scheme: {scheme.upper()}, N_delay = 4" not in info:
            failures.append(f"{scheme}: 'Sequential scheme: {scheme.upper()}, N_delay = 4' missing in info (fallback?)")
        assert observable_files(d) == files, f"{scheme}: different set of observable files"
        if config(d) != config(ref):
            failures.append(f"{scheme}: final configuration differs from the standard scheme")
        worst = (0.0, "")
        for name in files:
            a, b = numbers(ref / name), numbers(d / name)
            assert len(a) == len(b), f"{scheme}/{name}: different lengths"
            for x, y in zip(a, b):
                dev = abs(x - y)
                rel = dev / max(abs(x), abs(y), ATOL / RTOL)
                if rel > worst[0]:
                    worst = (rel, name)
        print(f"{scheme}: max relative deviation from standard {worst[0]:.3e} ({worst[1]}), RTOL {RTOL:.0e}")
        if not worst[0] <= RTOL:
            failures.append(f"{scheme}: {worst[1]} deviates by {worst[0]:.3e} > RTOL {RTOL:.0e}")
    if failures:
        raise AssertionError("\n".join(failures))


def main():
    binary, scenario = Path(sys.argv[1]).resolve(), sys.argv[2]
    with tempfile.TemporaryDirectory(prefix=f"alf-run-{scenario}-") as tmp:
        {"schemes": schemes}[scenario](binary, Path(tmp))
    print(f"PASS: {scenario}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Short deterministic end-to-end checks of fork features with the normal ALF executable (serial, plain-text output).

Scenarios (argument 2):
  schemes   Sequential_scheme = standard / delayed / submatrix on the space-time point defect model (Hubbard, square
            L = 4, ham_h0 = 1, N_FL = 2), same seed. The delayed and submatrix schemes are exact rewrites of the rank-1
            update (only the summation order differs, docs/gpu/DELAYED_SUBMATRIX.md), so the Markov chain is the same:
            the final configurations must be identical and every observable must agree to rounding (RTOL below).
  pinning   Hubbard rule: Ltau = 1 switches the pinning field off. ham_h0 = 1 and ham_h0 = 0 with Ltau = 1 must give
            bit-identical observable files and configurations (same code path after Ham_Set). Negative control: with
            Ltau = 0, ham_h0 = 1 and ham_h0 = 0 must give different observables (the comparison can see the field).
  local     SpinZ_local / SpinZ_localtau of the point defect model (Hubbard, ham_h0 = 1): file layout (bins, unit cells,
            tau grid Ntau = Ltrot + 1 with the run's dtau) and, per bin and site, SpinZ_local = trapezoidal tau average of
            SpinZ_localtau, (sum of the interior points + half the end points) / Ltrot. This holds exactly up to rounding
            because the end points tau = 0, beta are measured half as often and carry weight 2 (Hamiltonian_Hubbard_smod,
            Obser), and SpinZ_localtau is normalized by its visits of Ltrot/2. Tolerance LOCAL_TOL (absolute; values O(0.5),
            21-term sums, measured deviation <= 1e-16).
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
LOCAL_TOL = 1.0e-12


def parameters(*, scheme="standard", n_delay=4, ltau=0, h0=1.0, sweeps=20, bins=2):
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


def pinning(binary, root):
    on = run(binary, root / "ltau1_h1", ltau=1, h0=1.0)
    off = run(binary, root / "ltau1_h0", ltau=1, h0=0.0)
    if "pinning field ham_h0 is set to 0" not in (on / "stdout").read_text():
        raise AssertionError("Ltau = 1, ham_h0 = 1: no 'pinning field ham_h0 is set to 0' message")
    files = observable_files(off)
    assert observable_files(on) == files, "different set of observable files"
    assert any(f.endswith("_tau") for f in files), "Ltau = 1 run wrote no time-displaced observables"
    diff = [f for f in files if (on / f).read_bytes() != (off / f).read_bytes()]
    if config(on) != config(off):
        diff.append("confout_0")
    if diff:
        raise AssertionError(f"Ltau = 1: ham_h0 = 1 and ham_h0 = 0 differ in {diff}")
    print(f"Ltau = 1: ham_h0 = 1 and 0 bit-identical in {len(files)} observable files and confout_0")
    # negative control: without Ltau the field is active and visible in the observables
    on0 = run(binary, root / "ltau0_h1", ltau=0, h0=1.0)
    off0 = run(binary, root / "ltau0_h0", ltau=0, h0=0.0)
    if (on0 / "Ener_scal").read_bytes() == (off0 / "Ener_scal").read_bytes():
        raise AssertionError("negative control: Ltau = 0, ham_h0 = 1 and 0 give the same Ener_scal")
    print("negative control: Ltau = 0, ham_h0 = 1 and 0 differ (Ener_scal)")


def local_bins(path, timed):
    """Bins of an Obser_Latt_Local file: [(header, [((x, y), [values over tau])])]."""
    v = numbers(path)
    out, i = [], 0
    while i < len(v):
        norb, nunit = int(v[i+1]), int(v[i+2])
        if timed:
            header = (v[i], norb, nunit, int(v[i+3]), v[i+4]); nt = header[3]; i += 5
        else:
            header = (v[i], norb, nunit); nt = 1; i += 3
        sites = []
        for _ in range(norb*nunit):
            xy = (v[i], v[i+1]); i += 2
            sites.append((xy, [complex(v[i+2*k], v[i+2*k+1]) for k in range(nt)])); i += 2*nt
        out.append((header, sites))
    return out


def local(binary, root):
    d = run(binary, root / "defect", ltau=0, h0=1.0, bins=3)
    loc = local_bins(d / "SpinZ_local", False)
    tau = local_bins(d / "SpinZ_localtau", True)
    ltrot = 20                                   # Beta / Dtau of parameters()
    assert len(loc) == 3 and len(tau) == 3, f"expected 3 bins, got {len(loc)} and {len(tau)}"
    worst, biggest = 0.0, 0.0
    for (h1, s1), (h2, s2) in zip(loc, tau):
        assert h1[1:] == (1, 16) and h2[1:3] == (1, 16), f"unexpected lattice in headers {h1}, {h2}"
        assert h2[3] == ltrot + 1, f"SpinZ_localtau: Ntau = {h2[3]}, expected Ltrot + 1 = {ltrot + 1}"
        assert abs(h2[4] - 0.1) < 1e-12, f"SpinZ_localtau: dtau = {h2[4]}, expected 0.1"
        assert [xy for xy, _ in s1] == [xy for xy, _ in s2], "different site lists"
        for (_, a), (_, b) in zip(s1, s2):
            trap = (sum(b[1:-1]) + 0.5*(b[0] + b[-1])) / ltrot
            worst = max(worst, abs(a[0] - trap))
            biggest = max(biggest, abs(a[0]))
    print(f"SpinZ_local vs trapezoidal average of SpinZ_localtau: max deviation {worst:.3e} (tol {LOCAL_TOL:.0e}),"
          f" max |SpinZ_local| {biggest:.3f}")
    assert biggest > 1e-3, "SpinZ_local vanishes everywhere: the pinning field is not visible"
    assert worst <= LOCAL_TOL, f"SpinZ_local and SpinZ_localtau inconsistent: {worst:.3e} > {LOCAL_TOL:.0e}"


def main():
    binary, scenario = Path(sys.argv[1]).resolve(), sys.argv[2]
    with tempfile.TemporaryDirectory(prefix=f"alf-run-{scenario}-") as tmp:
        {"schemes": schemes, "pinning": pinning, "local": local}[scenario](binary, Path(tmp))
    print(f"PASS: {scenario}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""End-to-end test of the Gaussian bootstrap in Analysis/Max_SAC.F90.

Synthetic fermionic ("P" channel) data for a known spectrum (two Gaussian
peaks) are written to g_dat, Max_SAC.out is run in temporary directories, and
the following properties are checked:

  1. N_boot = 1 reproduces the original single run (no Aom_ps_boot_* files,
     exactly one MaxEnt run).
  2. N_boot > 1 performs N_boot + 1 runs (samples + unperturbed reference run),
     each started without restoring accumulated histograms from a dump.
  3. The bootstrap samples are independent: the Monte Carlo energy traces of
     all runs (from Max_stoch_log) differ pairwise. With the old fixed-seed
     reset in MaxEnt_stoch, samples 2..N_boot were identical.
  4. Aom_ps_boot_* exist, are finite and have positive errors where the
     spectrum has weight.
  5. The reference spectrum agrees with the bootstrap mean within errors.
  6. The added noise has the size given by the covariance: the chi^2 of the
     bootstrap samples exceeds that of the reference run by O(1) per tau
     point, and the bootstrap error exceeds the Monte Carlo error of a
     single run.

Usage:
    python3 2-maxsac-bootstrap.py [path/to/Max_SAC.out] [--keep] [--mpi N]

--mpi N runs Max_SAC.out with "mpirun -np N" (MPI build needed). Rank 0 then does
the reference run, ranks 1..N-1 the bootstrap samples in MaxSAC_rank_<r>/.

Exit code 0 if all checks pass, 1 otherwise. Only the Python standard library
is needed.
"""

import glob
import math
import os
import random
import shutil
import statistics
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------- settings
BETA = 10.0
NTAU = 41               # tau = 0, 0.25, ..., beta
SIGMA = 1.0e-2          # absolute error of G(tau)
NBINS = 40              # MaxEnt bins per run
NWARM = 10              # warmup bins
N_BOOT = 8
N_ALPHA = 14
NDIS = 200
OM_ST, OM_EN = -5.0, 5.0
SEED = 1234

PARAMETERS = """&VAR_errors
n_skip  = 1
N_rebin = 1
N_Cov   = 0
N_Back  = 1
N_auto  = 0
/
&VAR_Max_Stoch
N_boot     = {n_boot}
N_boot_checkpoint = .F.
Ngamma     = 100
Om_st      = {om_st}
Om_en      = {om_en}
NDis       = {ndis}
Nbins      = {nbins}
Nsweeps    = 30
NWarm      = {nwarm}
N_alpha    = {n_alpha}
alpha_st   = 1.d0
R          = 1.3d0
Checkpoint = .F.
Tolerance  = 100.d0
Stochastic = .true.
/
"""


def spectrum(om):
    """A(om) with int A = pi, i.e. G(0) + G(beta) = 1."""
    g = lambda x, m, s: math.exp(-0.5 * ((x - m) / s) ** 2) / (math.sqrt(2 * math.pi) * s)
    return math.pi * (0.5 * g(om, -1.0, 0.3) + 0.5 * g(om, 1.0, 0.3))


def green(tau, n=8000, om_min=-8.0, om_max=8.0):
    """G(tau) = int dom exp(-tau om) / (pi (1 + exp(-beta om))) A(om), midpoint rule."""
    dom = (om_max - om_min) / n
    tot = 0.0
    for i in range(n):
        om = om_min + (i + 0.5) * dom
        tot += math.exp(-tau * om) / (math.pi * (1.0 + math.exp(-BETA * om))) * spectrum(om)
    return tot * dom


def write_inputs(wdir, sigma, n_boot, seed):
    rng = random.Random(seed)
    tau = [BETA * i / (NTAU - 1) for i in range(NTAU)]
    g = [green(t) + sigma * rng.gauss(0.0, 1.0) for t in tau]
    with open(os.path.join(wdir, "g_dat"), "w") as f:
        f.write(f"{NTAU} 1000 {BETA} 1 P\n")
        for t, x in zip(tau, g):
            f.write(f"{t:.10e} {x:.17e} {sigma:.17e}\n")
    with open(os.path.join(wdir, "parameters"), "w") as f:
        f.write(PARAMETERS.format(n_boot=n_boot, om_st=OM_ST, om_en=OM_EN,
                                  ndis=NDIS, n_alpha=N_ALPHA, nbins=NBINS, nwarm=NWARM))


MPI_NP = 0   # set by --mpi N


def run(exe, wdir):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    cmd = ["mpirun", "-np", str(MPI_NP), exe] if MPI_NP > 1 else [exe]
    res = subprocess.run(cmd, cwd=wdir, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Max_SAC.out failed in {wdir}:\n{res.stdout}\n{res.stderr}")


def energy_traces(wdir):
    """'Alpha, En_m, Acc' lines of all Max_stoch_log files (MPI rank directories first,
    then the current directory, so that the reference run is always the last run)."""
    en, n_fresh = [], 0
    logs = sorted(glob.glob(os.path.join(wdir, "MaxSAC_rank_*", "Max_stoch_log")))
    for log in logs + [os.path.join(wdir, "Max_stoch_log")]:
        with open(log) as f:
            for line in f:
                if line.strip().startswith("Alpha, En_m, Acc"):
                    en.append(float(line.split()[4]))
                elif "No dump data" in line:
                    n_fresh += 1
    return en, n_fresh


def loadtxt(path):
    """Rows of floats, columns transposed: returns list of columns."""
    with open(path) as f:
        rows = [[float(x) for x in line.split()] for line in f if line.strip()]
    return rows


def load_boot(wdir, ia):
    return loadtxt(os.path.join(wdir, f"Aom_ps_boot_{ia}"))


# ------------------------------------------------------------------ checks
failures = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        failures.append(msg)


def main():
    global MPI_NP
    argv = sys.argv[1:]
    if "--mpi" in argv:
        k = argv.index("--mpi")
        MPI_NP = int(argv[k + 1])
        del argv[k:k + 2]
    args = [a for a in argv if a != "--keep"]
    keep = "--keep" in argv
    here = os.path.dirname(os.path.abspath(__file__))
    exe = os.path.abspath(args[0]) if args else \
        os.path.normpath(os.path.join(here, "..", "..", "Analysis", "Max_SAC.out"))
    if not os.access(exe, os.X_OK):
        print(f"Max_SAC.out not found or not executable: {exe}")
        return 1

    root = tempfile.mkdtemp(prefix="maxsac_bootstrap_")
    print(f"Executable: {exe}\nWork dir:   {root}" + (f"\nMPI ranks:  {MPI_NP}" if MPI_NP > 1 else ""))
    ia_last = N_ALPHA - 10
    n_per_run = NBINS * N_ALPHA        # Nbins * N_alpha log lines per MaxEnt run

    # 1. Backward compatibility: N_boot = 1.
    d1 = os.path.join(root, "nboot1"); os.mkdir(d1)
    write_inputs(d1, SIGMA, 1, SEED)
    run(exe, d1)
    en, n_fresh = energy_traces(d1)
    check(len(en) == n_per_run and n_fresh == 1, "N_boot=1: exactly one MaxEnt run")
    check(not os.path.exists(os.path.join(d1, "Aom_ps_boot_1")),
          "N_boot=1: no Aom_ps_boot_* written")

    # 2.-6. Bootstrap.
    d = os.path.join(root, "boot"); os.mkdir(d)
    write_inputs(d, SIGMA, N_BOOT, SEED)
    run(exe, d)
    en, n_fresh = energy_traces(d)
    n_runs = N_BOOT + 1
    check(len(en) == n_runs * n_per_run,
          f"N_boot+1 = {n_runs} MaxEnt runs (found {len(en) / n_per_run:g})")
    check(n_fresh == n_runs, f"every run starts without dump ({n_fresh}/{n_runs})")
    if len(en) == n_runs * n_per_run:
        tr = [en[k * n_per_run:(k + 1) * n_per_run] for k in range(n_runs)]
        same = [(i + 1, j + 1) for i in range(n_runs) for j in range(i + 1, n_runs)
                if tr[i] == tr[j]]
        check(not same, f"all runs differ pairwise (identical: {same})")

        # chi^2 at the largest alpha (last of the N_alpha replicas), after warmup.
        def chi2(run_trace):
            vals = [run_trace[b * N_ALPHA + N_ALPHA - 1] for b in range(NWARM, NBINS)]
            return sum(vals) / len(vals)
        chi2_ref = chi2(tr[-1])
        chi2_boot = sum(chi2(t) for t in tr[:-1]) / N_BOOT
        n_eff = NTAU - 1        # P channel: the tau = beta point is removed
        excess = (chi2_boot - chi2_ref) / n_eff
        print(f"      chi2 reference {chi2_ref:.1f}, bootstrap mean {chi2_boot:.1f}, "
              f"excess per tau point {excess:.2f}")
        check(0.4 < excess < 1.6,
              "added noise has the size of the covariance: chi2 excess per point in (0.4, 1.6)")

    ok = True
    for ia in range(1, ia_last + 1):
        b = load_boot(d, ia)
        ok &= len(b) == NDIS and all(len(r) == 3 and all(math.isfinite(x) for x in r)
                                     for r in b)
    check(ok, f"Aom_ps_boot_1..{ia_last} present, {NDIS}x3, finite")

    b = load_boot(d, ia_last)
    om = [r[0] for r in b]; a_mean = [r[1] for r in b]; a_err = [r[2] for r in b]
    peak = [k for k in range(NDIS) if a_mean[k] > 0.1 * max(a_mean)]
    check(all(a_err[k] > 0 for k in peak), "error > 0 wherever A has weight")

    ref = loadtxt(os.path.join(d, f"Aom_ps_{ia_last}"))
    check(all(abs(ref[k][0] - om[k]) < 1e-8 for k in range(NDIS)),
          "reference and bootstrap grids agree")
    dev = statistics.median(abs(ref[k][1] - a_mean[k]) / a_err[k] for k in peak
                            if a_err[k] > 0) if peak else float("inf")
    check(dev < 3.0, f"reference vs bootstrap mean: median |dA|/err = {dev:.2f} < 3")

    boot_err = sum(a_err)
    mc_err = sum(r[2] for r in ref)
    check(boot_err > mc_err,
          f"bootstrap error ({boot_err:.3g}) exceeds Monte Carlo error of one run ({mc_err:.3g})")

    # 7. Bootstrap of the (F,A) observable (DIDV.dat in the P channel, QFI_ph.dat in PH).
    fb_path = os.path.join(d, "DIDV_boot.dat")
    check(os.path.exists(fb_path), "DIDV_boot.dat written")
    if os.path.exists(fb_path):
        fb = loadtxt(fb_path)
        fr = loadtxt(os.path.join(d, "DIDV.dat"))
        check(len(fb) == N_ALPHA and all(len(r) == 3 and all(math.isfinite(x) for x in r)
                                         for r in fb),
              f"DIDV_boot.dat: {N_ALPHA}x3, finite")
        check(all(abs(fb[k][0] - fr[k][0]) <= 1e-8 * abs(fr[k][0]) for k in range(N_ALPHA)),
              "DIDV_boot.dat and DIDV.dat have the same alphas")
        check(all(r[2] > 0 for r in fb), "(F,A) bootstrap error > 0 for every alpha")
        devf = statistics.median(abs(fr[k][1] - fb[k][1]) / fb[k][2] for k in range(N_ALPHA)
                                 if fb[k][2] > 0)
        check(devf < 3.0, f"reference (F,A) vs bootstrap mean: median |dF|/err = {devf:.2f} < 3")
        bf = sum(r[2] for r in fb); mf = sum(r[2] for r in fr)
        check(bf > mf, f"(F,A) bootstrap error ({bf:.3g}) exceeds Monte Carlo error ({mf:.3g})")

    if keep or failures:
        print(f"Outputs kept in {root}")
    else:
        shutil.rmtree(root)
    print("\nALL CHECKS PASSED" if not failures else f"\n{len(failures)} CHECK(S) FAILED")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())

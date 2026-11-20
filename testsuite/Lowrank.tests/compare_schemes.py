#!/usr/bin/env python3
"""Compare the low-rank sequential schemes with the standard update on fixed seeds.

Written with AI assistance on master_ai; not part of the ALF baseline.

One ALF executable built with the LOWRANK switch runs every scenario (the serial scenarios of
testsuite/test_branch_parameters.json and those of scenarios.py) three times, with Sequential_scheme =
"standard", "delayed" and "submatrix" (N_delay = 0, automatic, and N_delay = 3) and the same seeds. Two levels of
agreement are checked for the delayed and the submatrix run against the standard run:

 1. STRICT: every file confout* (the final field configuration and the random number generator state)
    is byte-identical. The proposals of the fields do not depend on the Green function and all schemes
    draw the same random numbers in the same order, so identical files mean that every accept/reject
    decision was the same, i.e. the Markov chain is the same. Sequential MALA and Global_tau MALA are
    the exception: their proposals contain the force, which is computed from the Green function, so a
    roundoff difference of G enters the field and only the relative tolerance below is required there
    (scenarios with the flag force_in_proposal).
 2. TOLERANCE: every number in every other output file (bins, info excluded) agrees within --tol
    relative (max |a - b| / max(|a|, |b|, 1)); the low-rank schemes sum in a different order.

Scenarios may run serial, on 2 MPI ranks, as parallel tempering with 2 replicas, or with one parameter
file per replica (PARALLEL_PARAMS); which of these a given executable supports is given by --build.
A restart scenario runs two stages: the second starts from confout of the first (renamed to confin).

--n-delay -1 runs the low-rank schemes with the timed automatic N_delay (the candidates switch between
slices and the choice follows the machine timing; the chain must still agree up to roundoff). The
"Auto N_delay:" decision line of each run is printed when timing took place (>= 3 sweeps, at least two candidates).

Usage:
    compare_schemes.py --exe <ALF>/Prog/ALF.out
                       --start <ALF>/Scripts_and_Parameters_files/Start
                       [--build serial|mpi|tempering|parallel_params] [--mpiexec "mpirun ..."]
                       [--work DIR] [--only NAME ...] [--max-bins N] [--tol T] [--quick]
                       [--n-delay 0 3 -1 ...] [--test-models] [--report FILE]
"""
import argparse
import filecmp
import glob
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from baseline_guard import set_parameter  # noqa: E402
from scenarios import SCENARIOS  # noqa: E402

SCHEMES = ("standard", "delayed", "submatrix")
N_DELAYS = (0, 3)  # automatic N_delay and a small fixed one (many flushes)
IGNORED = {"info", "info_git", "parameters", "seeds", "stdout"}
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
MPI_BUILDS = ("mpi", "tempering", "parallel_params")
# scenarios of test_branch_parameters.json that belong to the quick set (one per model and some variants)
QUICK_JSON = ("Hubbard_N_Leg_Ladder", "Hubbard_PAM_2", "Hubbard_Plain_Vanilla", "Kondo", "tV", "LRC",
              "Z2_Matter", "Spin_Peierls")


def fortran_value(value):
    if isinstance(value, bool):
        return ".T." if value else ".F."
    if isinstance(value, str):
        return '"%s"' % value
    return repr(value)


def make_text(template, pars, max_bins=None, nbin=None):
    """Template parameter file with the assignments of `pars`; the second result lists unknown keys.
    The entry "nml" of `pars` (namelist name -> assignments) is appended as further namelists: the
    parameters of a model that has no namelist in the template."""
    text = template
    missing = []
    pars = dict(pars)
    nml = pars.pop("nml", {})
    if "Model" in pars:  # the model name is read in two namelists
        pars["ham_name"] = pars["Model"]
    for key, value in pars.items():
        if key == "NBin":
            if max_bins is not None:
                value = min(value, max_bins)
            if nbin is not None:
                value = nbin
        text, n = set_parameter(text, key, value)
        if n == 0:
            missing.append(key)
    for name, assignments in nml.items():
        text += "\n&%s\n%s/\n" % (name, "".join("%s = %s\n" % (k, fortran_value(v))
                                                 for k, v in assignments.items()))
    return text, missing


def execute(exe, directory, launcher, ranks):
    """Run ALF in `directory`, on `ranks` MPI ranks if the build needs a launcher; return the exit code."""
    cmd = (launcher.split() + ["-np", str(ranks)] if launcher else []) + [exe]
    with open(os.path.join(directory, "stdout"), "a") as out:
        res = subprocess.run(cmd, cwd=directory, stdout=out, stderr=subprocess.STDOUT,
                             env=dict(os.environ, OMP_NUM_THREADS="1"))
    return res.returncode


def prepare(directory, texts, seeds, mpi_per, build):
    """Create the run directory: `parameters` (+ Temp_i/parameters in the tempering builds), seeds."""
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    shutil.copy(seeds, os.path.join(directory, "seeds"))
    top = texts[0]
    if build in ("tempering", "parallel_params"):
        top = set_parameter(top, "mpi_per_parameter_set", mpi_per)[0]
        for i, t in enumerate(texts):
            os.makedirs(os.path.join(directory, "Temp_%d" % i))
            with open(os.path.join(directory, "Temp_%d" % i, "parameters"), "w") as f:
                f.write(set_parameter(t, "mpi_per_parameter_set", mpi_per)[0])
    with open(os.path.join(directory, "parameters"), "w") as f:
        f.write(top)


def restart_files(directory):
    """confout_k -> confin_k (also in Temp_i/): the next run continues the chain."""
    for root, _, files in os.walk(directory):
        for name in files:
            if name.startswith("confout"):
                os.replace(os.path.join(root, name),
                           os.path.join(root, "confin" + name[len("confout"):]))


def find_h5dump():
    """h5dump from the PATH or from the HDF5 installation ALF was built with (ALF_HDF5_DIR, or HDF5/ of the
    source tree, one subdirectory per compiler)."""
    found = shutil.which("h5dump")
    if found:
        return found
    roots = [os.environ.get("ALF_HDF5_DIR", ""), os.path.join(HERE, "..", "..", "HDF5")]
    for root in roots:
        for pattern in ("*/bin/h5dump", "bin/h5dump", "*/*/bin/h5dump"):
            hits = sorted(glob.glob(os.path.join(root, pattern))) if root else []
            if hits:
                return hits[0]
    return None


def h5_text(path):
    """Contents of an HDF5 file as text (h5dump), so that numbers() can read it; None if h5dump is missing."""
    exe = find_h5dump()
    if exe is None:
        return None
    try:
        res = subprocess.run([exe, path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return res.stdout.decode(errors="replace")


def numbers(path):
    if path.endswith(".h5"):
        text = h5_text(path)
        if text is None:
            raise RuntimeError("h5dump is needed to compare HDF5 output, but it was not found")
    else:
        with open(path) as f:
            text = f.read()
    return [float(t.replace("D", "E").replace("d", "e")) for t in NUMBER.findall(text)]


def compare(dir_a, dir_b):
    """Return (confout_identical, max_rel_diff, problem). The relative difference of the observables
    is max|x - y| / max(|x|, |y|, 1) over all numbers of all files except info and parameters. Where the
    configurations are not identical the numbers of confout* are included (MALA scenarios)."""
    worst = 0.0
    conf_same = True
    n_conf = 0
    for root, dirs, files in os.walk(dir_a):
        rel = os.path.relpath(root, dir_a)
        other = os.path.join(dir_b, rel)
        if not os.path.isdir(other):
            return conf_same, worst, "missing directory %s" % rel
        entries = sorted(os.listdir(other))
        if sorted(files + dirs) != entries:
            return conf_same, worst, "different files in %s" % rel
        for name in files:
            a, b = os.path.join(root, name), os.path.join(other, name)
            if name.startswith("confout"):
                n_conf += 1
                if name.endswith(".h5"):  # HDF5 files differ in metadata: compare the contents exactly
                    if numbers(a) != numbers(b):
                        conf_same = False
                elif not filecmp.cmp(a, b, shallow=False):
                    conf_same = False
                    x, y = numbers(a), numbers(b)
                    if len(x) != len(y):
                        return conf_same, worst, "%s has a different number of entries" % os.path.join(rel, name)
                    for u, v in zip(x, y):
                        worst = max(worst, abs(u - v) / max(abs(u), abs(v), 1.0))
                continue
            if name in IGNORED or name.endswith("_info"):
                continue
            x, y = numbers(a), numbers(b)
            if len(x) != len(y):
                return conf_same, worst, "%s has a different number of entries" % os.path.join(rel, name)
            for u, v in zip(x, y):
                worst = max(worst, abs(u - v) / max(abs(u), abs(v), 1.0))
    if n_conf == 0:
        return conf_same, worst, "no confout file written"
    return conf_same, worst, None


def lowrank_fraction(directory):
    """Fraction of the vertices handled by the low-rank kernel, from the info file (1.0 if all)."""
    for root, _, files in os.walk(directory):
        if "info" in files:
            with open(os.path.join(root, "info")) as f:
                for line in f:
                    m = re.search(r": (\d+) of (\d+) vertices have a support larger", line)
                    if m:
                        return 1.0 - int(m.group(1)) / int(m.group(2))
            return 1.0
    return None


def load_scenarios(path):
    """The scenarios of test_branch_parameters.json (lists: tempering replicas) and of scenarios.py."""
    scen = {}
    with open(path) as f:
        for name, (_, pars) in json.load(f).items():
            if isinstance(pars, list):
                scen[name] = dict(pars=pars[0], replicas=pars, mode="tempering")
            else:
                scen[name] = dict(pars=pars, quick=name in QUICK_JSON)
    scen.update(SCENARIOS)
    return scen


# With one replica these models stop in a TEMPERING build of ALF master itself (checked with upstream/master:
# Z2_Matter "Error in field setup", Hubbard_Plain_Vanilla a BLAS argument error)
NO_SINGLE_REPLICA = ("Z2_Matter", "Hubbard_Plain_Vanilla")


def runs_here(sc, build, test_models):
    if sc.get("testmodel") and not test_models:
        return False
    mode = sc.get("mode", "serial")
    if mode == "serial" and build in ("tempering", "parallel_params") \
            and sc["pars"].get("Model") in NO_SINGLE_REPLICA:
        return False
    if mode == "tempering":
        return build in ("tempering", "parallel_params")
    if mode == "params":
        return build == "parallel_params"
    if mode == "mpi":  # one parameter set: a TEMPERING build needs partners to exchange with
        return build in ("mpi", "parallel_params")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="ALF.out built with the LOWRANK switch")
    ap.add_argument("--start", required=True, help="directory with template 'parameters' and 'seeds'")
    ap.add_argument("--scenarios", default=os.path.join(HERE, "..", "test_branch_parameters.json"))
    ap.add_argument("--work", default="lowrank_runs")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--max-bins", type=int, default=None)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--build", default="serial", choices=("serial",) + MPI_BUILDS,
                    help="switches the executable was built with (default: serial)")
    ap.add_argument("--mpiexec", default=None,
                    help='MPI launcher without -np, e.g. "mpirun --oversubscribe"; needed for MPI builds')
    ap.add_argument("--n-delay", type=int, nargs="*", default=list(N_DELAYS),
                    help="values of N_delay for delayed and submatrix (0 = static automatic, -1 = timed automatic)")
    ap.add_argument("--test-models", action="store_true",
                    help="the executable contains the test models of testsuite/Lowrank.tests/models "
                         "(see add_test_models.sh); enables the scenarios that use them")
    ap.add_argument("--quick", action="store_true",
                    help="only a representative subset: scenarios flagged quick and those needing MPI ranks")
    ap.add_argument("--report", default=None, help="write a markdown table of the results")
    args = ap.parse_args()

    build = args.build
    if build in MPI_BUILDS and not args.mpiexec:
        ap.error("MPI builds need --mpiexec")
    launcher = args.mpiexec if build in MPI_BUILDS else None

    with open(os.path.join(args.start, "parameters")) as f:
        template = f.read()
    seeds = os.path.join(args.start, "seeds")
    exe = os.path.abspath(args.exe)

    failed = 0
    rows = []
    for name, sc in load_scenarios(args.scenarios).items():
        if args.only and name not in args.only:
            continue
        pars = sc["pars"]
        if pars.get("Langevin") or pars.get("HMC") or pars.get("Sequential") is False:
            continue  # no sequential updates: nothing to compare
        mode = sc.get("mode", "serial")
        if not runs_here(sc, build, args.test_models):
            if not sc.get("testmodel"):
                print("SKIP  %-34s needs another build than %s" % (name, build))
            continue
        if args.quick and not (sc.get("quick") or mode != "serial"):
            continue
        replicas = sc.get("replicas", [pars])
        mpi_per = sc.get("mpi_per", 2 if mode == "mpi" else 1)
        ranks = mpi_per * (len(replicas) if mode in ("tempering", "params") else 1)
        n_stage = [None]
        if sc.get("restart"):
            n_stage = [sc["restart"], sc["restart"]]

        variants = {"standard": ("standard", 0)}
        for scheme in SCHEMES[1:]:
            for nd in args.n_delay:
                variants[scheme if nd == 0 else "%s-N%d" % (scheme, nd)] = (scheme, nd)
        dirs, rcs, bad = {}, {}, None
        for label, (scheme, nd) in variants.items():
            lr = '\n&VAR_LOWRANK\nSequential_scheme = "%s"\nN_delay = %d\n/\n' % (scheme, nd)
            texts = []
            for rep in replicas:
                text, missing = make_text(template, rep, args.max_bins, n_stage[0])
                if missing:
                    bad = "keys not in the parameter template: %s" % ", ".join(missing)
                texts.append(text + lr)
            if bad:
                break
            dirs[label] = os.path.join(args.work, name, label)
            prepare(dirs[label], texts, seeds, mpi_per, build)
            rcs[label] = execute(exe, dirs[label], launcher, ranks)
            for _ in n_stage[1:]:
                if rcs[label] != 0:
                    break
                restart_files(dirs[label])
                rcs[label] = execute(exe, dirs[label], launcher, ranks)
        if bad:
            print("FAIL  %-34s scenario error: %s" % (name, bad))
            failed += 1
            continue
        if rcs["standard"] != 0:
            print("FAIL  %-34s standard run stops with exit code %d (see %s/stdout)"
                  % (name, rcs["standard"], dirs["standard"]))
            failed += 1
            continue
        frac = lowrank_fraction(dirs["delayed"])
        for scheme in variants:
            if scheme == "standard":
                continue
            if rcs[scheme] != 0:
                print("FAIL  %-34s %-14s exit code %d" % (name, scheme, rcs[scheme]))
                failed += 1
                continue
            same, d, problem = compare(dirs["standard"], dirs[scheme])
            strict = not sc.get("force_in_proposal")
            ok = problem is None and d <= args.tol and (same or not strict)
            failed += 0 if ok else 1
            rows.append((name, mode, scheme, same, d, frac, ok, problem))
            print("%s  %-34s %-14s confout %-9s max rel diff %.2e%s"
                  % ("PASS" if ok else "FAIL", name, scheme,
                     "identical" if same else "DIFFERS", d, "  " + problem if problem else ""))
            if variants[scheme][1] == -1:  # timed automatic N_delay: show what it decided
                with open(os.path.join(dirs[scheme], "info")) as f:
                    lines = [l.strip() for l in f if l.startswith("Auto N_delay")]
                print("      %s" % (lines[0] if lines else "no decision line (fewer than 3 sweeps, or a single candidate)"))
    if args.report:
        with open(args.report, "w") as f:
            f.write("| scenario | mode | scheme | confout identical | max rel diff | low-rank vertices | status |\n")
            f.write("|---|---|---|---|---|---|---|\n")
            for name, mode, scheme, same, d, frac, ok, problem in rows:
                f.write("| %s | %s | %s | %s | %.1e | %s | %s |\n"
                        % (name, mode, scheme, "yes" if same else "no", d,
                           "n/a" if frac is None else "%.0f%%" % (100 * frac), "ok" if ok else "FAIL"))
    print("%d scenario/scheme combinations compared, %d failed" % (len(rows), failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

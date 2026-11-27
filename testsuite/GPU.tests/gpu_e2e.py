#!/usr/bin/env python3
"""End-to-end check of the GPU build (P2.4 steps 1b, 1c, 1d, 1e): Hubbard 4x4, same seeds, for each case

  active    ALF.out of the GPU build, GPU path active (device wraps, local updates and stabilisation; ALF_GPU_SELFTEST=1
            additionally runs the mirror self-test at start-up)
  disabled  the same binary with ALF_GPU_DISABLE=1 (CPU path, mirror not built)
  none      optional (--ref-exe): ALF.out of a GPU=none build (must be byte-identical to 'disabled': CPU path unchanged)

Requirements per case: the info file of the active run says "ACTIVE (mirror, wraps and updates on the device path)" (cases that
must stay inactive: an INACTIVE line with the expected reason, and then every output file byte-identical to 'disabled');
confout_0 byte-identical to 'disabled' (the Markov chain is the CPU one for these short runs); every other output file equal to
'disabled' to max(1e-11, 100 x Precision Green mean, Precision Green max) of the 'disabled' run (|x - y| / max(1, |y|)).  The 1e-11 is the round-off of the device
kernels (e.g. FMA).  Since step 1d the device decomposition (ZGEQRF + ZUNMQR, column-norm pre-pivot) differs from the CPU one (ZGEQP3,
ZUNGQR) in the last bits of G at every stabilisation point, and the wrap instability amplifies that by the factor that also makes the
Precision Green line of the info file (mean deviation of the wrapped from the recomputed G) non-zero: 6e-11 for the Mz-Hubbard run
below (observable differences 1e-10 .. 1e-9), 1e-14 for the SU(2) runs (differences 1e-14).  The criterion ties the tolerance to that
noise floor of the run itself.  (The Precision Green means of 'active' and 'disabled' are printed: with Nwrap >= 25, where the
wrap is unstable, the device wraps and updates of step 1c give a mean 4 - 6x that of the CPU, already without the device stabilisation.)
Negative control (ALF_GPU_FAULT=upd, a deliberately wrong device
update): the observables must then differ by more than 1e-8, i.e. the comparison can see an error.
Step 1e (device observables and tau pass): every "ok" case also needs the info lines "GPU observables (step 1e): DEVICE" and
"GPU tau pass (step 1e): ACTIVE" (Ltau = 1; "OFF, reason: Ltau /= 1" otherwise) and the counter "on-demand G copies ... D->H 0" (no host
reader of G is left); the tau files (Green_tau, ...) are compared like all other files.  NOOBS_CASES are run a third time with
ALF_GPU_NOOBS=1 (host observables, host TAU_M; same binary, same device chain): equal-time and scalar files must agree with 'active' to
TOL_OBS = 1e-13 (summation order only), tau files to the tolerance of the case, and the info of that run must say
"HOST, reason: ALF_GPU_NOOBS=1" and show a D->H count > 0.  Faults ALF_GPU_FAULT=obs and ALF_GPU_FAULT_TAU = prop, blocks must be seen.
Step 1f (CGR phase from the tracked det(U), speedup S1): PHASE_CASES are run three more times, with the source of the CGR phase named in the
environment (ALF_GPU_DEBUG=1, ALF_GPU_CGR_LU = 0 tracked, 1 LU, 2 both), so the checks do not depend on the default of the build:
tracked against LU: confout byte-identical and every other file equal to TOL_PH (the phases differ at rounding level); check mode (2):
every file byte-identical to the tracked run (the tracked phase is the one used), the info file reports the number of CGR calls compared
and max |phase(tracked) - phase(LU)|, which must be below 1e-8 (the run aborts above that); tracked against 'disabled' like 'active'.
Fault ALF_GPU_FAULT_STAB=det in check mode must stop the run.
Cases: Mz-Hubbard (N_FL = 2) standard / delayed / submatrix with N_delay = 4, Ltau = 0, SU(2) Hubbard (Mz = .F., N_FL = 1) standard /
delayed, small, odd and large Nwrap (2, 7, 25, 50 of Ltrot = 50), timed N_delay = -1 (refused).  The run of the build is Ltau = 1
(time-displaced Green function: udvst exported to the host TAU_M) unless the case sets Ltau = 0.  .h5 files are compared with h5diff
if it is on the PATH, otherwise skipped with a note."""
import argparse, filecmp, os, re, shutil, subprocess, sys

TOL = 1.0e-11          # round-off of the device kernels: |x - y| / max(1, |y|) over all numbers of an output file
TOL_OBS = 1.0e-13      # device against host observables of the same chain: summation order only
FAULT_SEEN = 1.0e-8    # a negative control must change the observables by more than this
PREC_RATIO = 2.0       # "prec" cases: Precision Green mean of the device path at most this multiple of the CPU path's (KNOWN_ISSUES G10)
TOL_PH = 1.0e-13       # tracked against LU phase of the CGR (same binary): the phases differ at rounding level (1e-15), observables by 1e-13
PH_ABORT = 1.0e-8      # check mode limit
NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
NUM_WS = re.compile(r"\s*" + NUM.pattern)      # a number with its leading blanks (a minus sign takes the place of a blank)
ACTIVE_TEXT = "ACTIVE (mirror; wraps, updates and stabilisation on the device path)"

# label, edits of the parameter file, env of the active run, expectation ("ok", "okh", "hostobs", "inactive:<text>", "fault" or "abort");
# "okh" = "ok" for models whose observables and tau pass stay on the host (no step-1e line checks);
# "hostobs" = "okh" plus the checks of a Calc_Fl run (update path ACTIVE, observables HOST, tau pass OFF for reconstruction);
# "prec" = "okh" lines, confout byte-identical and Precision Green mean active / disabled <= PREC_RATIO; the other files are compared
# and printed but do not decide (at this size the 100 x Precision Green rule is below the round-off of the run, KNOWN_ISSUES G7)
KONDO = {"ham_name": '"Kondo"', "Model": '"Kondo"', "Lattice_type": '"Bilayer_square"', "ham_Uf": "1.d0"}
TV = {"ham_name": '"tV"', "Model": '"tV"', "N_SUN": "1", "ham_V": "-1.d0", "ham_T2": "0.d0", "ham_Tperp": "0.d0", "Beta": "2.d0", "Ltau": "0"}
SP = {"ham_name": '"Spin_Peierls"', "Model": '"Spin_Peierls"', "L1": "6", "L2": "6", "Beta": "2.d0", "ham_U": "1.5d0", "ham_Lambda": "0.1d0",
      "ham_omega0": "0.5d0", "N_SUN": "2", "N_FL": "1"}
PV = '"Hubbard_Plain_Vanilla"'
GT_RANGE = {"Global_tau_moves": ".T.", "N_Global_tau": "3", "Nt_sequential_start": "3", "Nt_sequential_end": "12"}
Z2 = {"ham_name": '"Z2_Matter"', "Model": '"Z2_Matter"', "Beta": "2.d0", "Ltau": "0"}
CASES = [
    ("mz-standard", {}, {"ALF_GPU_SELFTEST": "1"}, "ok"),
    ("mz-delayed-N4", {"Sequential_scheme": '"delayed"', "N_delay": "4"}, {}, "ok"),
    ("mz-submatrix-N4", {"Sequential_scheme": '"submatrix"', "N_delay": "4"}, {}, "ok"),
    ("su2-standard", {"Mz": ".F."}, {}, "ok"),
    ("su2-delayed-N4", {"Mz": ".F.", "Sequential_scheme": '"delayed"', "N_delay": "4"}, {}, "ok"),
    ("mz-ltau0", {"Ltau": "0"}, {}, "ok"),
    ("su2-nwrap2", {"Mz": ".F.", "Nwrap": "2"}, {}, "ok"),
    ("su2-nwrap7", {"Mz": ".F.", "Nwrap": "7"}, {}, "ok"),
    ("su2-nwrap25", {"Mz": ".F.", "Nwrap": "25"}, {}, "ok"),
    ("su2-nwrap50", {"Mz": ".F.", "Nwrap": "50"}, {}, "ok"),
    # step 1e: measurement at every slice (LOBS 1..Ltrot = 50), pinning field (local observable on the device), bilayer lattice with
    # U /= U2 (interaction energy weights, H7), Symm and Checkerboard off
    ("mz-lobs-all", {"LOBS_ST": "1", "LOBS_EN": "50"}, {}, "ok"),
    ("su2-lobs-all", {"Mz": ".F.", "LOBS_ST": "1", "LOBS_EN": "50"}, {}, "ok"),
    ("mz-h0", {"ham_h0": "1.d0"}, {}, "ok"),
    ("mz-bilayer-honeycomb", {"Lattice_type": '"Bilayer_honeycomb"', "L1": "2", "L2": "2", "ham_U2": "2.d0"}, {}, "ok"),
    ("su2-symm-F", {"Mz": ".F.", "Symm": ".F."}, {}, "ok"),
    ("mz-checkerboard-F", {"Checkerboard": ".F."}, {}, "ok"),
    ("timed-N-1", {"Sequential_scheme": '"delayed"', "N_delay": "-1"}, {}, "inactive:N_delay = -1"),
    ("fault-upd", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT": "upd"}, "fault"),
    # a wrong pivot permutation of the device decomposition (ALF_GPU_FAULT_STAB=perm) makes G wrong at the first stabilisation: the run
    # must stop with an error (expectation "abort": the active run exits with a non-zero return code)
    ("fault-obs", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT": "obs"}, "fault"),
    ("fault-tau-prop", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT_TAU": "prop"}, "fault"),
    ("fault-tau-blocks", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT_TAU": "blocks"}, "fault"),
    ("fault-stab-perm", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT_STAB": "perm"}, "abort"),
    # a flipped sign of the tracked det(U) phase must be stopped by the check mode (compares with the LU phase)
    ("fault-stab-det", {}, {"ALF_GPU_DEBUG": "1", "ALF_GPU_CGR_LU": "2", "ALF_GPU_FAULT_STAB": "det"}, "abort"),
    # P2.6.1: non-diagonal (Kondo, tV), rank-2 and type-4 (Spin_Peierls) vertices on the device path, standard / delayed (rank-r
    # delayed path) / submatrix; the self-tests ALF_GPU_SELFTEST=2 (device wrap kernels against Operator_mod) on two of them
    ("kondo-standard", KONDO, {"ALF_GPU_SELFTEST": "2"}, "okh"),
    ("kondo-delayed-N4", dict(KONDO, Sequential_scheme='"delayed"', N_delay="4"), {}, "okh"),
    ("kondo-submatrix-N4", dict(KONDO, Sequential_scheme='"submatrix"', N_delay="4"), {}, "okh"),
    ("tV-standard", TV, {}, "okh"),
    ("tV-delayed-N4", dict(TV, Sequential_scheme='"delayed"', N_delay="4"), {}, "okh"),
    ("tV-submatrix-N4", dict(TV, Sequential_scheme='"submatrix"', N_delay="4"), {}, "okh"),
    ("sp-standard", SP, {"ALF_GPU_SELFTEST": "2"}, "okh"),
    ("sp-delayed-N4", dict(SP, Sequential_scheme='"delayed"', N_delay="4"), {}, "okh"),
    ("sp-submatrix-N4", dict(SP, Sequential_scheme='"submatrix"', N_delay="4"), {}, "okh"),
    # G10: Spin_Peierls 12x12, standard: the rank-2 accept of the type-4 J vertices in Upgrade2's grouping (device / CPU ~1.1; the
    # Woodbury form G(:,P) W R gave ~30, the same chain)
    ("sp12-precision", dict(SP, L1="12", L2="12", NSweep="4"), {}, "prec"),
    # P2.6 Calc_Fl: Plain_Vanilla reconstructs flavour 1 from flavour 2 (weight, G, GT0/G0T); repulsive U at half filling and attractive
    # U off half filling (the two Calc_Fl triggers of the model). Device updates and stabilisation, host Obser and TAU_M ("hostobs").
    ("pv-half-standard", {"ham_name": PV}, {}, "hostobs"),
    ("pv-half-delayed-N4", {"ham_name": PV, "Sequential_scheme": '"delayed"', "N_delay": "4"}, {}, "hostobs"),
    ("pv-half-submatrix-N4", {"ham_name": PV, "Sequential_scheme": '"submatrix"', "N_delay": "4"}, {}, "hostobs"),
    ("pv-negU-standard", {"ham_name": PV, "ham_U": "-4.d0", "ham_chem": "0.5d0"}, {}, "hostobs"),
    ("pv-negU-delayed-N4", {"ham_name": PV, "ham_U": "-4.d0", "ham_chem": "0.5d0", "Sequential_scheme": '"delayed"',
                            "N_delay": "4"}, {}, "hostobs"),
    ("pv-negU-submatrix-N4", {"ham_name": PV, "ham_U": "-4.d0", "ham_chem": "0.5d0", "Sequential_scheme": '"submatrix"',
                              "N_delay": "4"}, {}, "hostobs"),
    ("pv-projector", {"ham_name": PV, "Projector": ".T."}, {}, "inactive:projector"),
    # P2.6 global: global tau moves on the device path. Hubbard: one random vertex per move (FINAL mode), full sequential range and a
    # partial one (H5 placement of G around it), standard / delayed / submatrix; Z2_Matter: moves of 4 or 5 vertices (INTERMEDIATE
    # mode, device copy of G for a rejected move), the fault ALF_GPU_FAULT=gtau (no G restore after a rejected move) must be seen.
    # Propose_S0 acts on type-1 vertices only: Hubbard (S0 = 1: half of the proposals skip the update, also inside delayed /
    # submatrix blocks) and Z2_Matter (S0 of the matter field).
    ("mz-gtau-N3", {"Global_tau_moves": ".T.", "N_Global_tau": "3", "Nt_sequential_start": "1", "Nt_sequential_end": "16"}, {}, "ok"),
    ("mz-gtau-range", dict(GT_RANGE), {}, "ok"),
    ("mz-gtau-range-delayed-N4", dict(GT_RANGE, Sequential_scheme='"delayed"', N_delay="4"), {}, "ok"),
    ("mz-gtau-range-submatrix-N4", dict(GT_RANGE, Sequential_scheme='"submatrix"', N_delay="4"), {}, "ok"),
    ("z2-gtau", dict(Z2, Global_tau_moves=".T.", N_Global_tau="4"), {}, "okh"),
    ("z2-gtau-delayed-N4", dict(Z2, Global_tau_moves=".T.", N_Global_tau="4", Sequential_scheme='"delayed"', N_delay="4"), {}, "okh"),
    ("fault-gtau", dict(Z2, Global_tau_moves=".T.", N_Global_tau="4"), {"ALF_GPU_DEBUG": "1", "ALF_GPU_FAULT": "gtau"}, "fault"),
    ("mz-propose-s0", {"Propose_S0": ".T."}, {}, "ok"),
    ("mz-propose-s0-delayed-N4", {"Propose_S0": ".T.", "Sequential_scheme": '"delayed"', "N_delay": "4"}, {}, "ok"),
    ("mz-propose-s0-submatrix-N4", {"Propose_S0": ".T.", "Sequential_scheme": '"submatrix"', "N_delay": "4"}, {}, "ok"),
    ("z2-propose-s0", dict(Z2, Propose_S0=".T."), {}, "okh"),
    # ALF_GPU_FAULT=nsig (stale device copy of the field) cannot be seen on the host backend: host and "device" array are the
    # same memory there, so it is not a case here (to be added for the vendor backends, step 2).
]


NOOBS_CASES = ("mz-standard", "su2-standard", "mz-h0", "mz-bilayer-honeycomb", "mz-lobs-all", "su2-symm-F")
# step 1f: Mz-Hubbard (N_FL = 2, real), SU(2) Hubbard (N_FL = 1, complex phases of det U), both with Ltau = 1; large Nwrap; delayed
PHASE_CASES = ("mz-standard", "su2-standard", "su2-nwrap25", "su2-delayed-N4")
PHASE_MODES = (("cgr0", "0", "CGR phase: tracked reflector determinants (no LU)"), ("cgr1", "1", "CGR phase: LU (ZGETRF)"),
               ("cgr2", "2", "CGR phase: tracked and LU, compared at every CGR (ALF_GPU_CGR_LU=2)"))


def prepare(d, start, edits):
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for f in ("parameters", "seeds"):
        shutil.copy(os.path.join(start, f), os.path.join(d, f))
    p = os.path.join(d, "parameters")
    t = open(p).read()
    for key, val in (("L1", 4), ("L2", 4), ("NSweep", 6), ("NBin", 2)):
        t, n = re.subn(r"(?m)^(%s\s*=\s*)\d+" % key, r"\g<1>%d" % val, t, count=1)
        if n != 1:
            sys.exit("cannot set %s in parameters" % key)
    for key, val in edits.items():
        # existing key: replace the value in every namelist that has it (only the model's namelists are read); new key of &VAR_QMC:
        # inserted before Ltau
        t, n = re.subn(r"(?mi)^(%s\s*=\s*)[^\s!]+" % re.escape(key), r"\g<1>%s" % val, t)   # namelist names: any case
        if n == 0:
            t, n = re.subn(r"(?m)^(Ltau\s*=)", "%s = %s\n\\g<1>" % (key, val), t, count=1)
        if n < 1:
            sys.exit("cannot set %s in parameters" % key)
    open(p, "w").write(t)


def run(exe, d, env_extra, mpiexec):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for k in [k for k in env if k.startswith("ALF_GPU_") or k.startswith("ALF_LA_TEST_")]:
        del env[k]
    env.update(env_extra)
    cmd = (mpiexec.split() + [exe]) if mpiexec else [exe]
    with open(os.path.join(d, "stdout"), "w") as out:
        rc = subprocess.call(cmd, cwd=d, env=env, stdout=out, stderr=subprocess.STDOUT)
    return rc


def files(d):
    return [f for f in sorted(os.listdir(d)) if f not in ("info", "info_git", "stdout") and not os.path.isdir(os.path.join(d, f))]


def numdiff(fa, fb):
    """(max abs difference, max of |x-y|/max(1,|y|)) over the numbers of two text files; None if the layout differs."""
    ta, tb = open(fa, errors="replace").read(), open(fb, errors="replace").read()
    if NUM_WS.sub("#", ta) != NUM_WS.sub("#", tb):
        return None
    va = [float(x.replace("D", "E").replace("d", "e")) for x in NUM.findall(ta)]
    vb = [float(x.replace("D", "E").replace("d", "e")) for x in NUM.findall(tb)]
    return (max([abs(x - y) for x, y in zip(va, vb)] + [0.0]),
            max([abs(x - y) / max(1.0, abs(y)) for x, y in zip(va, vb)] + [0.0]))


def green_precision(d):
    """(mean, max) of the 'Precision Green' line of the info file of run directory d, or None"""
    m = re.search(r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)", open(os.path.join(d, "info"), errors="replace").read())
    return (float(m.group(1).replace("D", "E")), float(m.group(2).replace("D", "E"))) if m else None


def device_backend(d):
    """True if the info file of run directory d names a device backend (CUDA, ONEMKL, HIP), False for HOST or a CPU run"""
    m = re.search(r"GPU build: backend (\w+)", open(os.path.join(d, "info"), errors="replace").read())
    return bool(m) and m.group(1) != "HOST"


def compare(a, b, strict, tol=TOL):
    """(files that are not within the rules, {file: (max abs diff, max scaled diff)} of the numeric files).  strict: every file
    byte-identical; else confout_* byte-identical and every other file equal to tol (scaled)."""
    bad, diffs = [], {}
    for f in files(a):
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if not os.path.exists(fb):
            bad.append(f + " (missing)")
        elif f.endswith(".h5"):
            if shutil.which("h5diff"):
                if subprocess.call(["h5diff", fa, fb], stdout=subprocess.DEVNULL) != 0:
                    bad.append(f)
            else:
                print("note: no h5diff, %s not compared" % f)
        elif filecmp.cmp(fa, fb, shallow=False):
            diffs[f] = (0.0, 0.0)
        elif strict or f.startswith("confout"):
            bad.append(f)
        else:
            dd = numdiff(fa, fb)
            if dd is None or dd[1] > (tol(f) if callable(tol) else tol):
                bad.append("%s (%s)" % (f, "layout differs" if dd is None else "max abs diff %.3g, scaled %.3g" % dd))
            diffs[f] = dd
    return bad, diffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--ref-exe", default=None)
    ap.add_argument("--mpiexec", default="")
    a = ap.parse_args()
    fail = 0
    for label, edits, env, expect in CASES:
        dirs = {}
        runs = {"active": (a.exe, env), "disabled": (a.exe, {"ALF_GPU_DISABLE": "1"})}
        if label in NOOBS_CASES:
            runs["noobs"] = (a.exe, dict(env, ALF_GPU_DEBUG="1", ALF_GPU_NOOBS="1"))
        if a.ref_exe:
            runs["none"] = (a.ref_exe, {})
        for lab, (exe, e) in runs.items():
            dirs[lab] = os.path.join(a.work, label, lab)
            prepare(dirs[lab], a.start, edits)
            rc = run(exe, dirs[lab], e, a.mpiexec)
            if lab == "active":
                rc_active = rc
            if rc != 0 and not (expect == "abort" and lab == "active"):
                sys.exit("FAIL: %s: run '%s' exited with %d (see %s/stdout)" % (label, lab, rc, dirs[lab]))
        if expect == "abort":
            ok = rc_active != 0
            print("%-16s active run (fault) exit code %d: %s" % (label, rc_active, "OK (run stopped)" if ok else "FAILED (run NOT stopped)"))
            fail += not ok
            continue
        info = open(os.path.join(dirs["active"], "info"), errors="replace").read()
        out = open(os.path.join(dirs["active"], "stdout"), errors="replace").read()
        dinfo = open(os.path.join(dirs["disabled"], "info"), errors="replace").read()
        bad_case = []
        if "INACTIVE, reason: disabled by ALF_GPU_DISABLE=1" not in dinfo:
            bad_case.append("info of the disabled run has no INACTIVE line")
        if expect.startswith("inactive:"):
            m = re.search(r"INACTIVE, reason: (.*)", info)
            if not m or expect[9:] not in m.group(1):
                bad_case.append("expected INACTIVE with reason containing '%s'" % expect[9:])
            reason = m.group(1).strip() if m else "?"
        else:
            reason = ""
            if ACTIVE_TEXT not in info:
                bad_case.append("info of the active run has no ACTIVE line")
            if "ALF_GPU_SELFTEST" in env and "GPU self-test PASSED" not in out:
                bad_case.append("stdout of the active run has no 'GPU self-test PASSED'")
            if expect in ("ok", "okh", "prec") and "GPU update path: ACTIVE" not in info:
                bad_case.append("info of the active run has no 'GPU update path: ACTIVE' line")
            if expect == "ok":
                tau_on = edits.get("Ltau", "1") != "0"
                if "GPU observables (step 1e): DEVICE" not in info:
                    bad_case.append("info of the active run has no 'GPU observables (step 1e): DEVICE' line")
                want = "GPU tau pass (step 1e): ACTIVE" if tau_on else "GPU tau pass (step 1e): OFF, reason: Ltau /= 1"
                if want not in info:
                    bad_case.append("info of the active run has no '%s' line" % want)
                if "local h0 T" in info and "ham_h0" not in edits:
                    bad_case.append("local observable on the device without ham_h0")
                if "ham_h0" in edits and "local h0 T" not in info:
                    bad_case.append("ham_h0 /= 0 but the info line does not say 'local h0 T'")
                m = re.search(r"on-demand G copies H->D (\d+) and D->H (\d+)", info)
                if not m or m.group(2) != "0":
                    bad_case.append("D->H count of on-demand G copies is not 0: %s" % (m.group(0) if m else "line missing"))
            if expect == "hostobs":
                for want in ("GPU update path: ACTIVE", "GPU observables (step 1e): HOST",
                             "GPU tau pass (step 1e): OFF, reason: Projector or reconstruction"):
                    if want not in info:
                        bad_case.append("info of the active run has no '%s' line" % want)
        if not os.path.exists(os.path.join(dirs["active"], "confout_0")) and not os.path.exists(os.path.join(dirs["active"], "data.h5")):
            bad_case.append("no confout_0 / data.h5 in the active run")
        gp_act, gp_dis = green_precision(dirs["active"]), green_precision(dirs["disabled"])
        # interim rule (KNOWN_ISSUES G7): the observables of two runs of the same chain differ by up to the largest wrap error, so
        # the CPU run's Precision Green max is also accepted; the GPU run's max is not (it would hide a worse device path, G10)
        tol = TOL if gp_dis is None else max(TOL, 100.0 * gp_dis[0], gp_dis[1])
        bad, diffs = compare(dirs["active"], dirs["disabled"], strict=expect.startswith("inactive:"), tol=tol)
        scal = ", ".join("%s %.3g / %.3g" % (f[:-5], diffs[f][0], diffs[f][1]) for f in ("Ener_scal", "Kin_scal", "Pot_scal")
                         if diffs.get(f) is not None)
        if expect == "fault":
            # negative control: the comparison must see the injected error
            seen = [f for f in files(dirs["active"]) if not f.startswith("confout") and (diffs.get(f) is None or diffs.get(f)[1] > FAULT_SEEN)]
            seen += [f for f in bad if not f.startswith("confout")]
            status = "OK (fault seen)" if seen else "FAILED (fault NOT seen)"
            if not seen:
                bad_case.append("fault not visible in the observables")
            print("%-16s active vs disabled: max abs / scaled diff (%s) %s" % (label, scal, status))
        elif expect == "prec":
            ratio = gp_act[0] / gp_dis[0] if gp_act and gp_dis and gp_dis[0] > 0 else float("inf")
            if ratio > PREC_RATIO:
                bad_case.append("Precision Green mean active / disabled %.3g > %.1f" % (ratio, PREC_RATIO))
            if not ("confout_0" in diffs and diffs["confout_0"][0] == 0.0):
                bad_case.append("confout_0 not byte-identical")
            mx = max([v[1] for f, v in diffs.items() if v is not None and not f.startswith("confout")] + [0.0])
            print("%-16s active vs disabled: Precision Green mean active %.2g / disabled %.2g = %.3g (limit %.1f), max scaled diff of the "
                  "files %.3g (not decisive) %s" % (label, gp_act[0] if gp_act else -1, gp_dis[0] if gp_dis else -1, ratio, PREC_RATIO, mx,
                                                     "FAILED" if bad_case else "OK"))
        else:
            status = "OK" if not bad and not bad_case else "FAILED"
            conf = "confout byte-identical" if "confout_0" in diffs and diffs["confout_0"][0] == 0.0 else "confout differs"
            gp = (", Precision Green mean active %.2g / disabled %.2g, tol %.2g" % (gp_act[0], gp_dis[0], tol)) if gp_act and gp_dis else ""
            print("%-16s active vs disabled: %s; max abs / scaled diff (%s)%s %s%s" % (label, conf, scal, gp, status,
                                                                                       (" [reason: %s]" % reason) if reason else ""))
            for b in bad:
                print("    not within the rules: " + b)
            fail += bool(bad)
        for b in bad_case:
            print("    " + b)
        fail += bool(bad_case)
        if "noobs" in dirs:
            ninfo = open(os.path.join(dirs["noobs"], "info"), errors="replace").read()
            tol_active = tol
            badn, dn = compare(dirs["active"], dirs["noobs"], strict=False, tol=lambda f: tol_active if "_tau" in f else TOL_OBS)
            mx = lambda sel: max([v[1] for f, v in dn.items() if v is not None and sel(f) and not f.startswith("confout")] + [0.0])
            m = re.search(r"on-demand G copies H->D (\d+) and D->H (\d+)", ninfo)
            okn = (not badn and "HOST, reason: ALF_GPU_NOOBS=1" in ninfo and m is not None and int(m.group(2)) > 0)
            print("%-16s active vs NOOBS=1: max scaled diff equal-time/scalar files %.3g (tolerance %.0e), tau files %.3g, D->H count %s: %s" %
                  (label, mx(lambda f: "_tau" not in f), TOL_OBS, mx(lambda f: "_tau" in f), m.group(2) if m else "?", "OK" if okn else "FAILED"))
            for b in badn:
                print("    not within the rules: " + b)
            fail += not okn
        if label in PHASE_CASES:
            pdirs = {}
            for lab, m, text in PHASE_MODES:
                pdirs[lab] = os.path.join(a.work, label, lab)
                prepare(pdirs[lab], a.start, edits)
                rc = run(a.exe, pdirs[lab], dict(env, ALF_GPU_DEBUG="1", ALF_GPU_CGR_LU=m), a.mpiexec)
                if rc != 0:
                    sys.exit("FAIL: %s: run '%s' exited with %d (see %s/stdout)" % (label, lab, rc, pdirs[lab]))
            bad_p = []
            for lab, m, text in PHASE_MODES:
                if text not in open(os.path.join(pdirs[lab], "info"), errors="replace").read():
                    bad_p.append("info of %s has no '%s'" % (lab, text))
            bt, dt = compare(pdirs["cgr0"], pdirs["cgr1"], strict=False, tol=TOL_PH)           # tracked vs LU
            bad_p += ["tracked vs LU: " + b for b in bt]
            # check mode == tracked: byte for byte on the host backend; on a device backend two runs of the same binary differ at
            # rounding level (device reductions, KNOWN_ISSUES G5 section), so there confout byte-identical and the rest to TOL_OBS
            dev = device_backend(pdirs["cgr0"])
            bc, _ = compare(pdirs["cgr2"], pdirs["cgr0"], strict=not dev, tol=TOL_OBS)
            bad_p += ["check vs tracked: " + b for b in bc]
            bd, _ = compare(pdirs["cgr0"], dirs["disabled"], strict=False, tol=tol)             # tracked vs CPU path
            bad_p += ["tracked vs disabled: " + b for b in bd]
            m = re.search(r"CGR calls compared (\d+), max \|phase\(tracked\) - phase\(LU\)\|\s+(\S+)",
                          open(os.path.join(pdirs["cgr2"], "info"), errors="replace").read())
            ncmp, pdiff = (int(m.group(1)), float(m.group(2).replace("D", "E"))) if m else (0, float("nan"))
            if not m or ncmp < 1 or not pdiff < PH_ABORT:
                bad_p.append("check mode report missing or out of range")
            obs = max([v[1] for f, v in dt.items() if v is not None and not f.startswith("confout")] + [0.0])
            conf = "confout byte-identical" if dt.get("confout_0", (1, 1))[0] == 0.0 else "confout differs"
            print("%-16s CGR phase tracked vs LU: %s, max scaled diff of the other files %.3g (tol %.0e); check mode: %d CGR calls, "
                  "max phase diff %.3g (limit %.0e); check run == tracked run %s: %s; %s" %
                  (label, conf, obs, TOL_PH, ncmp, pdiff, PH_ABORT, "to 1e-13 (device)" if dev else "byte for byte",
                   "yes" if not bc else "NO", "OK" if not bad_p else "FAILED"))
            for b in bad_p:
                print("    " + b)
            fail += bool(bad_p)
        if a.ref_exe:
            bad_n = compare(dirs["none"], dirs["disabled"], strict=True)[0]
            print("%-16s none   vs disabled: %s" % (label, "byte-identical" if not bad_n else "DIFFER: " + " ".join(bad_n)))
            fail += bool(bad_n)
    print("PASSED" if not fail else "FAILED")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())

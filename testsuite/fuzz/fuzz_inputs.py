#!/usr/bin/env python3
"""P2.8 input fuzzing and size probes for ALF (docs/ai/phase2/ROADMAP.md, P2.8.1 and P2.8.3).

Each case starts from Scripts_and_Parameters_files/Start (small 4x4 Hubbard run), changes one input
(parameters, seeds, confin_0, confin_0.h5, data.h5) and runs ALF.out with a time and memory limit.
The outcome is classified:

  clean    stopped with rc != 0 and ALF's own message (Terminate_on_error / error text)
  rtstop   stopped by a Fortran runtime error that names the problem (gfortran message)
  crash    signal or "Program received signal" (segfault, FPE, abort)
  oom      hit the memory limit
  hang     hit the time limit
  nan      rc 0, but NaN/Inf in info or an observable file
  runs     rc 0 with finite output

A case expects "stop" (bad input) or "run" (valid input). FAIL = a stop case that crashes, hangs,
runs out of memory or runs on; or a run case that does not finish cleanly.

Usage: fuzz_inputs.py ALF_DIR ALF_EXE WORK_DIR [--only CASE,CASE,...] [--hdf5 [--h5fc PATH]] [--timeout S] [--mem-gb G]
"""
import argparse
import math
import os
import re
import resource
import shutil
import signal
import subprocess
import sys

BASE = {"L1": 4, "L2": 4, "Beta": 1.0, "NSweep": 3, "NBin": 2, "Nwrap": 5, "Ltau": 0}
OBS_FILES = ("Ener_scal", "Pot_scal", "Kin_scal", "Part_scal")
NAN_RE = re.compile(r"(?i)\bnan\b|\binfinity\b|\binf\b")


def set_parameter(text, key, value):
    """Set every assignment of `key` (all namelists); value is written verbatim if str."""
    pattern = re.compile(r"^(\s*)%s(\s*)=\s*[^!\n]*" % re.escape(key), re.IGNORECASE | re.MULTILINE)
    new, n = pattern.subn(lambda m: "%s%s%s= %s " % (m.group(1), key, m.group(2), value), text)
    assert n > 0, "no key %s" % key
    return new


def drop_namelist(text, name):
    new, n = re.subn(r"(?ims)^\s*&%s\b.*?^\s*/\s*$" % re.escape(name), "", text)
    assert n == 1, name
    return new


def add_to_namelist(text, name, line):
    new, n = re.subn(r"(?im)^(\s*&%s\b.*)$" % re.escape(name), lambda m: m.group(1) + "\n" + line, text)
    assert n == 1, name
    return new


def base_text(alf):
    with open(os.path.join(alf, "Scripts_and_Parameters_files", "Start", "parameters")) as f:
        text = f.read()
    for k, v in BASE.items():
        text = set_parameter(text, k, repr(v) if isinstance(v, float) else str(v))
    return text


# ---------------------------------------------------------------------------------------------
# Cases: (name, expect, setup) where setup(dir, text, ctx) writes the inputs into dir.
# ctx holds paths of a reference confout_0 (and .h5) from a clean run, made once.
# ---------------------------------------------------------------------------------------------

def p(**kv):
    """Parameter change(s): values are strings written verbatim."""
    def setup(d, text, ctx):
        for k, v in kv.items():
            text = set_parameter(text, k, v)
        write(d, "parameters", text)
    return setup


def ptext(fn):
    def setup(d, text, ctx):
        write(d, "parameters", fn(text))
    return setup


def write(d, name, content, mode="w"):
    with open(os.path.join(d, name), mode) as f:
        f.write(content)


def seeds(content):
    def setup(d, text, ctx):
        write(d, "parameters", text)
        write(d, "seeds", content)
    return setup


def confin(fn, **kv):
    """confin_0 = fn(lines of a good confout_0); optional parameter changes."""
    def setup(d, text, ctx):
        for k, v in kv.items():
            text = set_parameter(text, k, v)
        write(d, "parameters", text)
        src = ctx["confout_cont"] if kv.get("Continuous") == ".T." else ctx["confout"]
        with open(src) as f:
            lines = f.read().splitlines(True)
        out = fn(lines)
        write(d, "confin_0", out, "wb" if isinstance(out, bytes) else "w")
    return setup


def confin_other_L(L):
    def setup(d, text, ctx):
        write(d, "parameters", text)
        shutil.copy(ctx["confout_L%d" % L], os.path.join(d, "confin_0"))
    return setup


def replace_fields(value, every=1):
    def fn(lines):
        out = lines[:ctx_seed_lines]
        for i, l in enumerate(lines[ctx_seed_lines:]):
            out.append((value + "\n") if i % every == 0 else l)
        return "".join(out)
    return fn


def h5(fn):
    """Corrupt a good confout_0.h5 written into confin_0.h5 (HDF5 builds)."""
    def setup(d, text, ctx):
        write(d, "parameters", text)
        dst = os.path.join(d, "confin_0.h5")
        shutil.copy(ctx["confout_h5"], dst)
        fn(dst, ctx)
    return setup


def h5_truncate(nbytes):
    def fn(path, ctx):
        with open(path, "r+b") as f:
            f.truncate(nbytes)
    return fn


def h5_garbage(path, ctx):
    with open(path, "wb") as f:
        f.write(os.urandom(4096))


def h5_other_L(L):
    def fn(path, ctx):
        shutil.copy(ctx["confout_h5_L%d" % L], path)
    return fn


def h5_mut(op):
    """Edit the HDF5 file with h5mut (testsuite/fuzz/h5mut.f90, built with ALF's h5fc)."""
    def fn(path, ctx):
        subprocess.run([ctx["h5mut"], path, op], check=True)
    return fn


def datah5(fn):
    """Restart with a corrupted data.h5 next to a good confin_0.h5."""
    def setup(d, text, ctx):
        write(d, "parameters", text)
        shutil.copy(ctx["confout_h5"], os.path.join(d, "confin_0.h5"))
        shutil.copy(ctx["data_h5"], os.path.join(d, "data.h5"))
        fn(os.path.join(d, "data.h5"))
    return setup


ctx_seed_lines = None  # set after the reference run (length of the seed block in confout_0)

TEXT_CASES = [
    # --- valid controls
    ("ok_base", "run", p()),
    ("ok_crlf", "run", ptext(lambda t: t.replace("\n", "\r\n"))),
    ("ok_tabs", "run", ptext(lambda t: t.replace("            ", "\t"))),
    # --- parameters: missing pieces
    ("par_missing_file", "stop", lambda d, t, c: None),
    ("par_empty_file", "stop", ptext(lambda t: "")),
    ("par_no_VAR_QMC", "stop", ptext(lambda t: drop_namelist(t, "VAR_QMC"))),
    ("par_no_VAR_lattice", "stop", ptext(lambda t: drop_namelist(t, "VAR_lattice"))),
    ("par_no_VAR_Model_Generic", "stop", ptext(lambda t: drop_namelist(t, "VAR_Model_Generic"))),
    ("par_no_VAR_Hubbard", "stop", ptext(lambda t: drop_namelist(t, "VAR_Hubbard"))),
    ("par_no_VAR_ham_name", "stop", ptext(lambda t: drop_namelist(t, "VAR_ham_name"))),
    ("par_truncated_in_QMC", "stop", ptext(lambda t: t[:t.index("NSweep")])),
    ("par_no_slash_QMC", "stop", ptext(lambda t: re.sub(r"(?s)(&VAR_QMC.*?)\n/", r"\1\n", t, count=1))),
    ("par_bom", "stop?", ptext(lambda t: "﻿" + t)),
    # --- parameters: unknown keys, wrong types
    ("par_unknown_key_QMC", "stop", ptext(lambda t: add_to_namelist(t, "VAR_QMC", "Nswep = 3"))),
    ("par_unknown_key_lattice", "stop", ptext(lambda t: add_to_namelist(t, "VAR_lattice", "L3 = 4"))),
    ("par_unknown_key_Hubbard", "stop", ptext(lambda t: add_to_namelist(t, "VAR_Hubbard", "ham_J = 1.0"))),
    ("par_L1_real", "stop", p(L1="4.5")),
    ("par_L1_string", "stop", p(L1='"four"')),
    ("par_L1_huge", "stop", p(L1="99999999999")),
    ("par_Beta_string", "stop", p(Beta="abc")),
    ("par_Mz_int", "stop", p(Mz="7")),
    ("par_NSweep_real", "stop", p(NSweep="1e3")),
    ("par_ham_name_unknown", "stop", p(ham_name='"Hubard"')),
    ("par_ham_name_empty", "stop", p(ham_name='""')),
    ("par_lattice_empty", "stop", p(Lattice_type='""')),
    ("par_dup_key", "run", ptext(lambda t: add_to_namelist(t, "VAR_QMC", "NSweep = 2"))),
    # --- parameters: NaN / Inf / extremes (namelist reads "NaN" and "Inf")
    ("par_U_nan", "stop", p(ham_U="NaN")),
    ("par_U_inf", "stop", p(ham_U="Inf")),
    ("par_U_1e308", "stop", p(ham_U="1e308")),
    ("par_U_1e5", "stop?", p(ham_U="1e5")),
    ("par_T_nan", "stop", p(ham_T="NaN")),
    ("par_T_inf", "stop", p(ham_T="Inf")),
    ("par_T_1e200", "stop", p(ham_T="1e200")),
    ("par_chem_nan", "stop", p(ham_chem="NaN")),
    ("par_chem_1e5", "stop?", p(ham_chem="1e5")),
    ("par_h0_nan", "stop", p(ham_h0="NaN")),
    ("par_PhiX_nan", "stop", p(Phi_X="NaN")),
    ("par_Beta_nan", "stop", p(Beta="NaN")),
    ("par_Beta_inf", "stop", p(Beta="Inf")),
    ("par_Dtau_nan", "stop", p(Dtau="NaN")),
    ("par_Dtau_inf", "stop", p(Dtau="Inf")),
    ("par_Dtau_neg", "stop", p(Dtau="-0.1")),
    ("par_Dtau_1e-300", "stop", p(Dtau="1e-300")),
    ("par_Theta_huge_proj", "stop", p(Projector=".T.", Theta="1e9")),
    ("par_Theta_nan_proj", "stop", p(Projector=".T.", Theta="NaN")),
    ("par_CPU_MAX_nan", "stop", p(CPU_MAX="NaN")),
    ("par_CPU_MAX_neg", "stop?", p(CPU_MAX="-1.0")),
    ("par_NBin_0", "stop", p(NBin="0")),
    ("par_NBin_neg", "stop", p(NBin="-1")),
    ("par_Ltau_2", "stop", p(Ltau="2")),
    ("par_Ltau_neg", "stop", p(Ltau="-1")),
    ("par_N_SUN_neg", "stop", p(N_SUN="-2")),
    ("par_N_FL_0", "stop", p(N_FL="0")),
    ("par_L1_neg", "stop", p(L1="-4")),
    ("par_L2_0", "stop", p(L2="0")),
    ("par_NSweep_huge", "hang?", p(NSweep="2147483647")),
    ("par_LOBS_neg", "stop", p(LOBS_ST="-5", LOBS_EN="-1")),
    ("par_LOBS_beyond", "stop", p(LOBS_ST="50", LOBS_EN="60")),
    ("par_Amplitude_nan_cont", "stop", p(Continuous=".T.", Amplitude="NaN")),
    ("par_Amplitude_inf_cont", "stop", p(Continuous=".T.", Amplitude="Inf")),
    ("par_dtLangevin_nan", "stop", p(Continuous=".T.", Langevin=".T.", Delta_t_Langevin_HMC="NaN")),
    ("par_MaxForce_inf", "stop", p(Continuous=".T.", Langevin=".T.", Max_Force="Inf")),
    ("par_U2_nan", "stop", p(ham_U2="NaN")),
    ("par_Nwrap_huge", "stop?", p(Nwrap="2147483647")),
    ("par_N_Phi_huge", "stop?", p(N_Phi="2147483647")),
    # --- seeds
    ("seeds_missing", "stop", lambda d, t, c: write(d, "parameters", t)),
    ("seeds_empty", "stop", seeds("")),
    ("seeds_garbage", "stop", seeds("hello\nworld\n")),
    ("seeds_real", "stop?", seeds("1.5\n")),
    ("seeds_zero", "run", seeds("0\n")),
    ("seeds_neg", "run", seeds("-12345\n")),
    ("seeds_huge", "stop", seeds("99999999999\n")),
    ("seeds_binary", "stop", seeds("\x00\x01\x02\xff\n")),
    # --- confin_0 (text)
    ("confin_empty", "stop", confin(lambda l: "")),
    ("confin_seed_only", "stop", confin(lambda l: "".join(l[:ctx_seed_lines]))),
    ("confin_short_seed", "stop", confin(lambda l: "".join(l[:2]))),
    ("confin_truncated", "stop", confin(lambda l: "".join(l[:len(l) // 2]))),
    ("confin_garbage", "stop", confin(lambda l: "foo bar\n" * len(l))),
    ("confin_binary", "stop", confin(lambda l: os.urandom(2000))),
    ("confin_field_0", "stop", confin(replace_fields("0"))),
    ("confin_field_3", "stop", confin(replace_fields("3"))),
    ("confin_field_neg2", "run", confin(replace_fields("-2"))),  # Hubbard discrete = type 2 (+-1, +-2)
    ("confin_field_real", "stop", confin(replace_fields("0.5"))),
    ("confin_field_nan", "stop", confin(replace_fields("NaN"))),
    ("confin_extra_lines", "stop", confin(lambda l: "".join(l) + "1\n" * 500)),
    ("confin_seed_nan", "stop", confin(lambda l: "NaN\n" + "".join(l[1:]))),
    ("confin_from_L6", "stop", confin_other_L(6)),
    ("confin_from_L2", "stop", confin_other_L(2)),
    # --- confin_0, continuous fields (type 3)
    ("cont_ok_restart", "run", confin(lambda l: "".join(l), Continuous=".T.")),
    ("cont_field_nan", "stop", confin(replace_fields("NaN"), Continuous=".T.")),
    ("cont_field_inf", "stop", confin(replace_fields("Inf"), Continuous=".T.")),
    ("cont_field_1e308", "stop", confin(replace_fields("1e308"), Continuous=".T.")),
    ("cont_field_1e3", "stop?", confin(replace_fields("1e3"), Continuous=".T.")),
    ("cont_field_nan_one", "stop", confin(replace_fields("NaN", every=10**9), Continuous=".T.")),
    # --- both formats present
    ("confin_and_h5", "stop", lambda d, t, c: (write(d, "parameters", t), shutil.copy(c["confout"], os.path.join(d, "confin_0")),
                                               write(d, "confin_0.h5", "x"))),
    # --- P2.8.1 size probes (memory limit set by the runner)
    ("size_L1_1", "stop", p(L1="1", L2="1")),
    ("size_chain_L2_1", "run", p(L2="1", Lattice_type='"N_leg_ladder"', Checkerboard=".F.")),
    ("size_Ltrot_1", "run", p(Beta="0.1", Nwrap="1")),
    ("size_L_200_memlimit", "stop", p(L1="200", L2="200")),
    ("size_N_SUN_200", "stop?", p(N_SUN="200", Mz=".F.")),
    ("size_Ltrot_1e8", "stop", p(Beta="1e7")),
    ("size_Ltrot_2e9", "stop", p(Beta="2e8")),
]

H5_CASES = [
    ("h5_ok_restart", "run", h5(lambda path, ctx: None)),
    ("h5_empty", "stop", h5(h5_truncate(0))),
    ("h5_truncated", "stop", h5(lambda path, ctx: h5_truncate(os.path.getsize(path) // 2)(path, ctx))),
    ("h5_garbage", "stop", h5(h5_garbage)),
    ("h5_from_L6", "stop", h5(h5_other_L(6))),
    ("h5_from_L2", "stop", h5(h5_other_L(2))),
    ("h5_no_configuration", "stop", h5(h5_mut("del_conf"))),
    ("h5_no_seed", "stop", h5(h5_mut("del_seed"))),
    ("h5_short_seed", "stop", h5(h5_mut("short_seed"))),
    ("h5_field_0", "stop", h5(h5_mut("conf_zero"))),
    ("h5_field_nan", "stop", h5(h5_mut("conf_nan"))),
    ("h5_config_int", "stop?", h5(h5_mut("conf_int"))),
    ("datah5_truncated", "stop", datah5(lambda path: h5_truncate(os.path.getsize(path) // 2)(path, None))),
    ("datah5_garbage", "stop", datah5(lambda path: h5_garbage(path, None))),
]


# ---------------------------------------------------------------------------------------------

def limit(mem_gb):
    def fn():
        b = int(mem_gb * 2**30)
        try:
            resource.setrlimit(resource.RLIMIT_AS, (b, b))
        except (ValueError, OSError):
            pass  # macOS: RLIMIT_AS cannot be set; the memory cases are skipped there (MEM_LIMIT)
    return fn


def mem_limit_works():
    probe = "import resource; resource.setrlimit(resource.RLIMIT_AS, (2**31, 2**31))"
    return subprocess.run([sys.executable, "-c", probe], capture_output=True).returncode == 0


MEM_CASES = ("size_L_200_memlimit", "size_Ltrot_1e8", "size_Ltrot_2e9")  # allocate > 20 GB without the limit


def run(exe, d, timeout, mem_gb):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with open(os.path.join(d, "stdout"), "w") as so, open(os.path.join(d, "stderr"), "w") as se:
        proc = subprocess.Popen([exe], cwd=d, stdout=so, stderr=se, env=env, preexec_fn=limit(mem_gb),
                                start_new_session=True)
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            return "hang", None
    return None, rc


def has_nan(d):
    """NaN/Inf in an observable file or in the precision/acceptance lines of info (not the parameter echo)."""
    for name in ("info",) + OBS_FILES:
        path = os.path.join(d, name)
        if os.path.exists(path):
            with open(path, errors="replace") as f:
                for line in f:
                    if name == "info" and not re.search(r"Precision|Acceptance", line):
                        continue
                    if NAN_RE.search(line):
                        return name
    return None


def frozen(d):
    """Acceptance exactly 0: no move was accepted (the chain did not move)."""
    path = os.path.join(d, "info")
    if os.path.exists(path):
        with open(path, errors="replace") as f:
            for line in f:
                m = re.match(r"\s*Acceptance\s*:\s*(\S+)", line)
                if m and float(m.group(1)) == 0.0:
                    return True
    return False


def classify(d, timeout, mem_gb, exe):
    kind, rc = run(exe, d, timeout, mem_gb)
    with open(os.path.join(d, "stderr"), errors="replace") as f:
        err = f.read()
    with open(os.path.join(d, "stdout"), errors="replace") as f:
        out = f.read()
    text = err + out
    if kind == "hang":
        return "hang", None, text
    if re.search(r"Error allocating \d+ bytes", text) and rc not in (None, 0) and rc > 0:
        return "memstop", rc, text  # gfortran stops with file, line and size: a clean stop
    if re.search(r"(?i)memory limit|cannot allocate|allocation would exceed|Out of memory|std::bad_alloc", text):
        return "oom", rc, text
    if rc < 0 or "Program received signal" in text or rc in (134, 136, 139):
        return "crash", rc, text
    if rc != 0:
        if "Fortran runtime error" in text:
            return "rtstop", rc, text
        return "clean", rc, text
    bad = has_nan(d)
    if bad:
        return "nan", rc, "NaN/Inf in " + bad
    if frozen(d):
        return "frozen", rc, "acceptance 0"
    return "runs", rc, text


def verdict(expect, kind):
    if expect == "run":
        return "ok" if kind == "runs" else "FAIL"
    if expect == "stop":
        return "ok" if kind in ("clean", "rtstop", "memstop") else "FAIL"
    if expect == "hang?":
        return "ok" if kind in ("clean", "rtstop", "hang") else "FAIL"
    # "stop?": either outcome is acceptable as long as it is not a crash, NaN, hang or OOM
    return "ok" if kind in ("clean", "rtstop", "runs") else "FAIL"


def first_message(text):
    for line in text.splitlines():
        s = line.strip()
        if s and not s.startswith(("Note:", "ERROR STOP", "Setting N_FL")) and not re.match(r"^#\d", s) and "ALF usually" not in s:
            if re.search(r"(?i)error|warning|stop|abort|signal|invalid|unable|cannot|not|NaN|limit", s):
                return s[:110]
    return ""


def reference(alf, exe, work, args, hdf5):
    """Clean runs whose confout files the confin cases corrupt."""
    global ctx_seed_lines
    ctx = {}
    for tag, changes in (("", {}), ("_cont", {"Continuous": ".T."}), ("_L6", {"L1": "6", "L2": "6"}),
                         ("_L2", {"L1": "2", "L2": "2"})):
        d = os.path.join(work, "ref" + tag)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        text = base_text(alf)
        for k, v in changes.items():
            text = set_parameter(text, k, v)
        write(d, "parameters", text)
        shutil.copy(os.path.join(alf, "Scripts_and_Parameters_files", "Start", "seeds"), d)
        kind, rc, _ = classify(d, args.timeout, args.mem_gb, exe)
        assert kind == "runs", "reference run %s: %s" % (d, kind)
        if hdf5:
            ctx["confout_h5" + tag] = os.path.join(d, "confout_0.h5")
            if tag == "":
                ctx["data_h5"] = os.path.join(d, "data.h5")
        else:
            ctx["confout" + tag] = os.path.join(d, "confout_0")
    if not hdf5:
        # the seed block = the lines before the first field line (fields are +-1 for Ising)
        with open(ctx["confout"]) as f:
            lines = f.read().splitlines()
        n_fields = len(lines)
        k = 0
        while lines[k].strip() not in ("1", "-1"):
            k += 1
        ctx_seed_lines = k
        assert n_fields > k
    return ctx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("alf")
    ap.add_argument("exe")
    ap.add_argument("work")
    ap.add_argument("--only", default="", help="comma-separated case names")
    ap.add_argument("--hdf5", action="store_true", help="ALF_EXE is an HDF5 build: run the HDF5 cases")
    ap.add_argument("--h5fc", default="", help="h5fc of the HDF5 that ALF was built with (default: PATH, then ALF_DIR/HDF5)")
    ap.add_argument("--timeout", type=float, default=60)
    ap.add_argument("--mem-gb", type=float, default=2.0)
    args = ap.parse_args()
    exe = os.path.abspath(args.exe)
    args.work = os.path.abspath(args.work)
    os.makedirs(args.work, exist_ok=True)
    ctx = reference(args.alf, exe, args.work, args, args.hdf5)
    if args.hdf5:
        h5fc = args.h5fc or shutil.which("h5fc") or \
            [os.path.join(r, "h5fc") for r, _, f in os.walk(os.path.abspath(os.path.join(args.alf, "HDF5"))) if "h5fc" in f][0]
        ctx["h5mut"] = os.path.join(args.work, "h5mut")
        subprocess.run([h5fc, "-o", ctx["h5mut"], os.path.join(os.path.dirname(os.path.abspath(__file__)), "h5mut.f90")],
                       check=True, cwd=args.work)
    cases = H5_CASES if args.hdf5 else TEXT_CASES
    if args.only:
        only = args.only.split(",")
        unknown = set(only) - {c[0] for c in cases}
        assert not unknown, "unknown cases: %s" % sorted(unknown)
        cases = [c for c in cases if c[0] in only]
    n_fail = 0
    mem_ok = mem_limit_works()
    for name, expect, setup in cases:
        if name in MEM_CASES and not mem_ok:
            print("%-26s %-6s %-7s skipped: the memory limit cannot be set on this system" % (name, expect, "skip"))
            continue
        d = os.path.join(args.work, name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        if name != "seeds_missing":
            shutil.copy(os.path.join(args.alf, "Scripts_and_Parameters_files", "Start", "seeds"), d)
        setup(d, base_text(args.alf), ctx)
        kind, rc, text = classify(d, args.timeout, args.mem_gb, exe)
        v = verdict(expect, kind)
        n_fail += v == "FAIL"
        msg = text if kind in ("nan", "frozen") else first_message(text)
        print("%-26s %-6s %-7s rc=%-4s %-4s %s" % (name, expect, kind, rc, v, msg), flush=True)
    print("FAILED: %d" % n_fail)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

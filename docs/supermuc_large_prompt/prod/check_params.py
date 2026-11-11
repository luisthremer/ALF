#!/usr/bin/env python3
"""Check the parameter files of a PP run directory against the reference parameter set (e.g. docs/supermuc_large_prompt/ref_L24/parameters).

usage: check_params.py RUNDIR REFFILE        (python3 standard library only; exit 1 on any violation)

Files checked: RUNDIR/parameters and every RUNDIR/Temp_<i>/parameters.
Rules
  * key = value lines of Fortran namelists, keys case-insensitive, whitespace normalised, trailing ! comments dropped.
  * Allowed to differ from the reference: L1, L2, Beta, N_Phi, ham_U, NSweep, N_delay, CPU_MAX, NBin.
  * The namelist &VAR_TEMP may be present in addition (it is not in the reference).
  * Everything else must be textually identical. 1.d0 vs 1.0d0 counts as a difference (reported as "same number, other spelling").
  * Per Temp_<i> (and the top-level file): N_FL=2, ham_h0=1.d0, Mz=.F., L1=L2=L, Beta=L.d0, N_Phi=L*L/2, NSweep=10, Nwrap=10.
  * Across the Temp_<i>: only ham_U may differ between them, and the ham_U values must be distinct.
"""
import glob, os, re, sys

ALLOWED = {"l1", "l2", "beta", "n_phi", "ham_u", "nsweep", "n_delay", "cpu_max", "nbin"}
OPTIONAL_GROUPS = {"var_temp"}


def parse(path):
    """-> dict {(group, key): value}, group/key lower case; a key given twice in one group is recorded as ('!dup', 'group.key')."""
    out, group = {}, ""
    with open(path) as fh:
        for raw in fh:
            res, q = [], None
            for ch in raw.strip():          # drop comments outside of quotes
                if q:
                    res.append(ch)
                    if ch == q:
                        q = None
                elif ch in "\"'":
                    q = ch
                    res.append(ch)
                elif ch == "!":
                    break
                else:
                    res.append(ch)
            line = "".join(res).strip()
            if not line:
                continue
            if line.startswith("&"):
                group = line[1:].split()[0].lower()
            elif line == "/":
                group = ""
            elif "=" in line:
                k, v = line.split("=", 1)
                key = (group, k.strip().lower())
                if key in out:
                    out[("!dup", key[0] + "." + key[1])] = ""
                out[key] = " ".join(v.split())
    return out


def num(v):
    try:
        return float(re.sub(r"[dD]", "e", v.strip()))
    except (ValueError, AttributeError):
        return None


def logical(v):
    s = (v or "").strip().lower().strip(".")
    return {"t": True, "true": True, "f": False, "false": False}.get(s)


def get(p, key):
    for (g, k), v in p.items():
        if k == key and g != "!dup":
            return v
    return None


def table(title, head, data):
    print("\n### %s\n" % title)
    print("| " + " | ".join(head) + " |")
    print("|" + "|".join("---" for _ in head) + "|")
    for r in data:
        print("| " + " | ".join(str(x) for x in r) + " |")
    if not data:
        print("| (none) |" + " |" * (len(head) - 1))


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    rundir, reffile = sys.argv[1], sys.argv[2]
    ref = parse(reffile)

    def tnum(p):
        m = re.search(r"Temp_(\d+)", p)
        return int(m.group(1)) if m else 0
    temps = sorted(glob.glob(os.path.join(rundir, "Temp_*", "parameters")), key=tnum)
    files = [os.path.join(rundir, "parameters")] + temps
    rows, asserts, params = [], [], {}      # rows: differences to the reference; asserts: fixed checks
    nviol = 0
    if not os.path.exists(files[0]):
        print("VIOLATION: %s missing" % files[0])
        return 1
    if not temps:
        rows.append(("-", "-", "-", "-", "VIOLATION: no Temp_*/parameters"))
        nviol += 1

    for f in files:
        p = parse(f)
        params[f] = p
        short = os.path.relpath(f, rundir)
        for (g, k) in sorted(p):
            if g == "!dup":
                rows.append((short, k, "-", "-", "VIOLATION: key given twice"))
                nviol += 1
        for key in sorted(set(ref) | set(p)):
            g, k = key
            if g == "!dup":
                continue
            rv, fv = ref.get(key), p.get(key)
            name = "%s.%s" % (g, k)
            if g in OPTIONAL_GROUPS:
                if fv is not None:
                    rows.append((short, name, "(not in ref)", fv, "ok (&VAR_TEMP block added)"))
            elif rv is None:
                rows.append((short, name, "(absent)", fv, "VIOLATION: extra key"))
                nviol += 1
            elif fv is None:
                rows.append((short, name, rv, "(absent)", "VIOLATION: key missing"))
                nviol += 1
            elif rv != fv:
                if k in ALLOWED:
                    rows.append((short, name, rv, fv, "ok (allowed to differ)"))
                elif num(rv) is not None and num(rv) == num(fv):
                    rows.append((short, name, rv, fv, "VIOLATION: same number, other spelling"))
                    nviol += 1
                else:
                    rows.append((short, name, rv, fv, "VIOLATION: differs"))
                    nviol += 1

        # fixed assertions
        L, L2 = num(get(p, "l1")), num(get(p, "l2"))
        Li = int(L) if L is not None else None
        checks = [
            ("L1 == L2", "equal", L is not None and L == L2, "%s / %s" % (get(p, "l1"), get(p, "l2"))),
            ("N_FL", "2", num(get(p, "n_fl")) == 2, get(p, "n_fl")),
            ("ham_h0", "1.d0", num(get(p, "ham_h0")) == 1.0, get(p, "ham_h0")),
            ("Mz", ".F.", logical(get(p, "mz")) is False, get(p, "mz")),
            ("Beta", "L.d0 (= %s)" % ("%d.d0" % Li if Li else "?"), L is not None and num(get(p, "beta")) == L, get(p, "beta")),
            ("N_Phi", "L*L/2 (= %s)" % (Li * Li // 2 if Li else "?"), Li is not None and num(get(p, "n_phi")) == Li * Li // 2, get(p, "n_phi")),
            ("NSweep", "10", num(get(p, "nsweep")) == 10, get(p, "nsweep")),
            ("Nwrap", "10", num(get(p, "nwrap")) == 10, get(p, "nwrap")),
        ]
        for what, exp, ok, found in checks:
            asserts.append((short, what, exp, found, "ok" if ok else "VIOLATION"))
            nviol += 0 if ok else 1

    # consistency across the sets and with the top-level file
    if len(temps) > 1:
        base, us = params[temps[0]], []
        for f in temps:
            p = params[f]
            us.append(get(p, "ham_u"))
            for key in sorted(set(base) | set(p)):
                if key[1] != "ham_u" and key[0] != "!dup" and base.get(key) != p.get(key):
                    rows.append((os.path.relpath(f, rundir), "%s.%s" % key, base.get(key), p.get(key), "VIOLATION: differs from Temp_0"))
                    nviol += 1
        if len(set(us)) != len(us):
            rows.append(("Temp_*", "ham_u", "distinct", ",".join(str(u) for u in us), "VIOLATION: ham_U values not distinct"))
            nviol += 1
    top = params[files[0]]
    if temps:
        for key in ("l1", "nsweep", "n_delay", "cpu_max", "nbin"):
            if get(top, key) != get(params[temps[0]], key):
                rows.append(("parameters vs Temp_0", key, get(top, key), get(params[temps[0]], key), "VIOLATION: top-level and Temp_0 differ"))
                nviol += 1
        if get(top, "mpi_per_parameter_set") is None:
            rows.append(("parameters", "var_temp.mpi_per_parameter_set", "-", "(absent)", "VIOLATION: &VAR_TEMP missing in top-level file"))
            nviol += 1

    print("## check_params: %s against %s" % (rundir, reffile))
    print("%d file(s): %s" % (len(files), ", ".join(os.path.relpath(f, rundir) for f in files[:3])) + (" ..." if len(files) > 3 else ""))
    okrows = [r for r in rows if r[4].startswith("ok")]
    rows = [r for r in rows if not r[4].startswith("ok")]
    print("\n### Allowed differences to the reference (informational)\n")
    for f in files:
        short = os.path.relpath(f, rundir)
        mine = ["%s: %s -> %s" % (r[1].split(".")[-1], r[2], r[3]) for r in okrows if r[0] == short and "VAR_TEMP" not in r[4]]
        print("- %s: %s%s" % (short, "; ".join(mine) if mine else "none", " (+ &VAR_TEMP block)" if any(r[0] == short and "VAR_TEMP" in r[4] for r in okrows) else ""))
    table("Violations of the reference comparison", ["file", "key", "reference", "found", "verdict"], rows)
    bad = [a for a in asserts if a[4] != "ok"]
    table("Fixed assertions (%d checks, %d failed)" % (len(asserts), len(bad)), ["file", "check", "expected", "found", "verdict"], bad)
    print("\n**%s** (%d violation(s))" % ("FAIL" if nviol else "PASS", nviol))
    return 1 if nviol else 0


if __name__ == "__main__":
    sys.exit(main())

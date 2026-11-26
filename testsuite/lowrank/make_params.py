#!/usr/bin/env python3
"""Create an ALF run directory from the Start template and override parameters.

Usage: make_params.py [--template DIR] OUT_DIR KEY=VALUE ...

Copies Scripts_and_Parameters_files/Start/{parameters,seeds} (or DIR/{parameters,seeds},
e.g. the Start directory of the baseline worktree) to OUT_DIR and sets
each KEY = VALUE in every namelist that contains KEY (case-insensitive match; the
trailing '! comment' of the line is kept). A KEY found in no namelist (for example
Sequential_scheme before it exists in the template) is appended to &VAR_QMC.
SEED=<int> replaces the first seed in 'seeds' instead.

The test-only model LR_Jz_Test has no namelist in the Start template: if a key sets ham_name to it, the namelist
&VAR_LR_Jz_Test (defaults of Hamiltonian_LR_Jz_Test_smod.F90) is appended to the copied parameters file first.
"""
import os
import re
import shutil
import sys

START = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", "Scripts_and_Parameters_files", "Start")


EXTRA_NAMELISTS = {
    "LR_Jz_Test": """
&VAR_LR_Jz_Test            !! Variables for the test-only model LR_Jz_Test
ham_T      = 1.d0           ! Hopping parameter
ham_U      = 4.d0           ! Hubbard interaction (Mz)
ham_Jz     = 1.d0           ! Ising coupling on nearest-neighbour bonds
ham_chem   = 0.d0           ! Chemical potential
Dtau       = 0.1d0          ! Thereby Ltrot=Beta/dtau
Beta       = 2.d0           ! Inverse temperature
/
""",
}


def set_key(lines, key, value):
    """Set key=value in all namelists containing key; append to &VAR_QMC if absent."""
    pat = re.compile(r"^(\s*)(" + re.escape(key) + r")(\s*=\s*)([^!\n]*?)(\s*)(!.*)?$", re.I)
    found = False
    for i, line in enumerate(lines):
        m = pat.match(line.rstrip("\n"))
        if m:
            found = True
            comment = m.group(6)
            tail = ("  " + comment) if comment else ""
            lines[i] = "%s%s%s%s%s\n" % (m.group(1), m.group(2), m.group(3), value, tail)
    if not found:
        start = None
        for i, line in enumerate(lines):
            if re.match(r"^\s*&VAR_QMC\b", line, re.I):
                start = i
            elif start is not None and line.strip() == "/":
                lines.insert(i, "%s = %s\n" % (key, value))
                return
        sys.exit("make_params.py: &VAR_QMC namelist not found, cannot append " + key)


def main(argv):
    template = START
    if len(argv) >= 3 and argv[1] == "--template":
        template = argv[2]
        argv = argv[:1] + argv[3:]
    if len(argv) < 2:
        sys.exit(__doc__)
    out = argv[1]
    os.makedirs(out, exist_ok=True)
    for name in ("parameters", "seeds"):
        shutil.copyfile(os.path.join(template, name), os.path.join(out, name))
    with open(os.path.join(out, "parameters")) as f:
        params = f.readlines()
    with open(os.path.join(out, "seeds")) as f:
        seeds = f.readlines()
    for arg in argv[2:]:
        if "=" in arg and arg.split("=", 1)[0].lower() == "ham_name":
            name = arg.split("=", 1)[1].strip().strip("\"'")
            nml = EXTRA_NAMELISTS.get(name)
            if nml and not any(re.match(r"^\s*&VAR_" + name + r"\b", l, re.I) for l in params):
                params.extend(nml.splitlines(True))
    for arg in argv[2:]:
        if "=" not in arg:
            sys.exit("make_params.py: expected KEY=VALUE, got " + arg)
        key, value = arg.split("=", 1)
        if key.upper() == "SEED":
            seeds[0] = "%12d\n" % int(value)
        else:
            set_key(params, key, value)
    with open(os.path.join(out, "parameters"), "w") as f:
        f.writelines(params)
    with open(os.path.join(out, "seeds"), "w") as f:
        f.writelines(seeds)


if __name__ == "__main__":
    main(sys.argv)

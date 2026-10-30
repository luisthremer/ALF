#!/bin/bash
# Ticket 212 (development only), nividia: run directories (parameters + dev seeds) for configs of acceptance/212.json
# in dev/runs/base.   usage: t212_prep.sh CONFIG [CONFIG ...]
W=$HOME/alf-gpu/work/t212/dev; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
: "${HOME:?}"
mkdir -p "$W/runs"; cd $W/src/testsuite/gpu_reference || exit 1
python3 -I - "$W/runs/base" "$W/src" "$W/src/testsuite/gpu_reference/acceptance/212.json" "$@" <<'PY'
import hashlib, json, os, sys
sys.path.insert(0, os.getcwd())
from py_alf import ALF_source, Simulation
from configs import CONFIGS
from configs_stat import CONFIGS as CS
CONFIGS = {**CONFIGS, **CS}
out, alf, accf, names = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3], sys.argv[4:]
acc = json.load(open(accf))["configs"]
src = ALF_source(alf_dir=alf)
for name in names:
    d = {**CONFIGS[acc[name]["base"]], **acc[name].get("overrides", {})}
    Simulation(src, "Hubbard", d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
    h = hashlib.sha256(open(os.path.join(out, name, "parameters"), "rb").read()).hexdigest()
    print(name, "parameters sha256", h, "MATCH" if h == acc[name]["parameters_sha256"] else "differs from 212.json")
PY

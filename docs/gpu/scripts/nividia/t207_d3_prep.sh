#!/bin/bash
# Ticket 207 (dev3: diagnosis of the default-path paired spread, development only), nividia: run directories (parameters + committed dev seeds) for configs of
# acceptance/207.draft.json plus the extra Ltau = 0 twins below, into dev3/runs/<SET>.  SET = base (dev seeds) or devN (the dev
# seeds rotated by 37*N lines, as acc_driver.py --who agent maps accN).   usage: t207_d2_prep.sh SET CONFIG [CONFIG ...]
W=$HOME/alf-gpu/work/t207/dev3; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
SET=$1; shift
case "$SET" in base|dev[0-9]*) ;; *) echo "bad set $SET"; exit 1;; esac
: "${HOME:?}"
mkdir -p "$W/runs"; cd $W/src/testsuite/gpu_reference || exit 1
~/alf-gpu/venv/bin/python3 -I - "$W/runs/$SET" "$W/src" "$W/src/testsuite/gpu_reference/acceptance/207.draft.json" "$SET" "$@" <<'PY'
import hashlib, json, os, sys
sys.path.insert(0, os.getcwd())
from py_alf import ALF_source, Simulation
from configs import CONFIGS
from configs_stat import CONFIGS as CS
CONFIGS = {**CONFIGS, **CS}
out, alf, accf, sset, names = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5:]
acc = json.load(open(accf))["configs"]
# Ltau = 0 twins (the tau pass leaves the sweep unchanged) and other dev-only configs
acc["s8_U4_b2_cb_t0"] = {"base": "s8_U4_b2_cb", "overrides": {"NSweep": 15, "NBin": 120}}   # compare2 T0 length (<= 15 min)
acc["hub12x12_cb_b10_s"] = {"base": "hub12x12_cb", "overrides": {"beta": 10.0, "NSweep": 10, "NBin": 3}}   # short determinism runs
for k in list(acc):
    acc[k + "_lt0"] = {"base": acc[k]["base"], "overrides": {**acc[k].get("overrides", {}), "Ltau": 0}}
src = ALF_source(alf_dir=alf)
for name in names:
    d = {**CONFIGS[acc[name]["base"]], **acc[name].get("overrides", {})}
    Simulation(src, "Hubbard", d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
    sf = os.path.join(out, name, "seeds")
    if sset != "base":
        lines = open(sf).read().splitlines()
        k = (37 * int(sset[3:] or 1)) % len(lines)
        open(sf, "w").write("\n".join(lines[k:] + lines[:k]) + "\n")
    h = hashlib.sha256(open(os.path.join(out, name, "parameters"), "rb").read()).hexdigest()
    hs = hashlib.sha256(open(sf, "rb").read()).hexdigest()
    print(sset, name, "parameters sha256", h, "seeds sha256", hs[:16])
PY

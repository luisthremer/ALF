#!/bin/bash
# Ticket 302c-2 (development only), on intel: run directories (parameters + committed dev seeds) for configs of
# acceptance/207.draft.json into t302c-2/runs/<SET>, exactly as t207_3_prep.sh / t207_d5_prep.sh on nividia (same seed rotation, so
# devN here has the same seeds file as devN there).  SET = base or devN (dev seeds rotated by 37*N lines).
#   usage: t302c2_prep.sh SET CONFIG [CONFIG ...]      e.g. t302c2_prep.sh dev1 hub12x12_cb_p hub12x12_cb_b10_p
W=${T302C2_W:-$HOME/alf-gpu/work/t302c-2}   # round 2: T302C2_W=$HOME/alf-gpu/work/t302c-2b (only these two)
case "$W" in "$HOME/alf-gpu/work/t302c-2"|"$HOME/alf-gpu/work/t302c-2b") ;; *) echo "refusing work dir $W"; exit 1;; esac
:; export TMPDIR=$W/tmp
SET=$1; shift
case "$SET" in base|dev[0-9]*) ;; *) echo "bad set $SET"; exit 1;; esac
: "${HOME:?}"
mkdir -p "$W/runs" $TMPDIR; cd $W/src/testsuite/gpu_reference || exit 1
~/alf-gpu/venv/bin/python3 -I - "$W/runs/$SET" "$W/src" "$W/src/testsuite/gpu_reference/acceptance/207.draft.json" "$SET" "$@" <<'PY'
import hashlib, json, os, sys
sys.path.insert(0, os.getcwd())
from py_alf import ALF_source, Simulation
from configs import CONFIGS
from configs_stat import CONFIGS as CS
CONFIGS = {**CONFIGS, **CS}
out, alf, accf, sset, names = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5:]
acc = json.load(open(accf))["configs"]
# round 2 (T0 as t207_3_prep.sh on nividia): short 6x6 / beta = 10 configs and the Ltau = 0 twins (dev-only, not in 207.draft.json)
acc["hub12x12_cb_b10_s"] = {"base": "hub12x12_cb", "overrides": {"beta": 10.0, "NSweep": 10, "NBin": 3}}
acc["hub6x6_cb_s"] = {"base": acc["hub6x6_cb"]["base"], "overrides": {**acc["hub6x6_cb"].get("overrides", {}), "NSweep": 10, "NBin": 5}}
for k in list(acc):
    acc.setdefault(k + "_lt0", {"base": acc[k]["base"], "overrides": {**acc[k].get("overrides", {}), "Ltau": 0}})
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
    ok = ('dev-only config' if 'parameters_sha256' not in acc[name] else
          'parameters_sha256 matches 207.draft.json' if h == acc[name]['parameters_sha256'] else 'PARAMETERS DIFFER FROM 207.draft.json')
    print(sset, name, "parameters sha256", h, "seeds sha256", hs[:16], ok)
PY

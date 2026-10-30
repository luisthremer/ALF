#!/bin/bash
# Ticket 207, nividia: additional run directories next to runs/base (does not touch the existing ones):
#   first-divergence runs (1 sweep per bin) for the Symm and complex-hopping (flux) paths, and the dense hopping.
#   usage: t207_prep2.sh [NBIN]     -> runs/base/hub6x6_cb_symm_fd, hub6x6_cb_flux_fd, hub6x6_dense_fd
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
cd $W/src_new/testsuite/gpu_reference || exit 1
: "${HOME:?}"
python3 -I - "$W/runs/base" "$W/src_new" "${1:-60}" <<'PY'
import os, sys
sys.path.insert(0, os.getcwd())
from py_alf import ALF_source, Simulation
from configs import CONFIGS
out, alf, nb = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), int(sys.argv[3])
src = ALF_source(alf_dir=alf)
for name, base, extra in (("hub6x6_cb_symm_fd", "hub6x6_cb_symm", dict()),
                          ("hub6x6_cb_flux_fd", "hub6x6_cb", dict(N_Phi=1, Phi_X=0.3)),
                          ("hub6x6_dense_fd", "hub6x6_dense", dict())):
    d = {**CONFIGS[base], "NSweep": 1, "NBin": nb, **extra}
    Simulation(src, "Hubbard", d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
    print(os.path.join(out, name))
PY

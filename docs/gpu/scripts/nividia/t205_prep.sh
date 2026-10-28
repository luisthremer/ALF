#!/bin/bash
# Ticket 205, nividia: run directories (parameters + seeds) in work/t205/runs/base.
#   standard configs of testsuite/gpu_reference (cb 6/12/20/32, dense 6/12/20/32, mz off, symm) plus two with magnetic flux
#   (complex hopping: complex checkerboard bonds and a complex dense element, which no gpu_reference config has).
#   usage: t205_prep.sh [NSWEEP] [NBIN]
W=$HOME/alf-gpu/work/t205; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
cd $W/src_new/testsuite/gpu_reference || exit 1
: "${HOME:?}"
rm -rf -- "$W/runs/base"; mkdir -p "$W/runs"
python3 prep.py $W/runs/base hub6x6_cb hub6x6_dense hub12x12_cb hub6x6_cb_mzoff hub6x6_cb_symm hub20x20_cb hub12x12_dense hub20x20_dense hub32x32_cb hub32x32_dense \
    --nsweep ${1:-20} --nbin ${2:-2} --alf-dir $W/src_new
python3 -I - "$W/runs/base" "$W/src_new" "${1:-20}" "${2:-2}" <<'PY'
import os, sys
sys.path.insert(0, os.getcwd())
from py_alf import ALF_source, Simulation
from configs import CONFIGS
out, alf, ns, nb = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
src = ALF_source(alf_dir=alf)
for name, base, extra in (("hub6x6_cb_flux", "hub6x6_cb", dict(N_Phi=1, Phi_X=0.3)),
                          ("hub6x6_dense_flux", "hub6x6_dense", dict(N_Phi=1, Phi_X=0.3))):
    d = dict(CONFIGS[base], NSweep=ns, NBin=nb, **extra)
    Simulation(src, "Hubbard", d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
    print(os.path.join(out, name))
PY

#!/bin/bash
# Ticket 207, nividia: run directories (parameters + seeds) in work/t207/runs/base.
#   standard configs of testsuite/gpu_reference (cb 6/12/20/32, dense 6/12/20/32, mz off, symm) plus two with magnetic flux
#   (complex hopping: complex checkerboard bonds and a complex dense element, which no gpu_reference config has).
#   usage: t207_prep.sh [NSWEEP] [NBIN]
W=$HOME/alf-gpu/work/t207; export TMPDIR=$W/tmp; . ~/alf-gpu/env/gnu.sh
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
# ticket 206/207: beta = 10 cases (more stabilization blocks, larger scales): Nwrap = 10 (the 104 regime C/B at U = 4) and
# Nwrap = 5 (regime A, the comfortable one); NSweep/NBin as for the beta = 4 configs of the same size unless overridden below
for name, base, extra in (("hub6x6_cb_flux", "hub6x6_cb", dict(N_Phi=1, Phi_X=0.3)),
                          ("hub6x6_dense_flux", "hub6x6_dense", dict(N_Phi=1, Phi_X=0.3)),
                          ("hub6x6_cb_b10", "hub6x6_cb", dict(beta=10.0)),
                          ("hub6x6_dense_b10", "hub6x6_dense", dict(beta=10.0)),
                          ("hub12x12_cb_b10", "hub12x12_cb", dict(beta=10.0)),
                          ("hub12x12_cb_b10_w5", "hub12x12_cb", dict(beta=10.0, Nwrap=5)),
                          ("hub12x12_dense_b10", "hub12x12_dense", dict(beta=10.0)),
                          ("hub20x20_cb_b10", "hub20x20_cb", dict(beta=10.0)),
                          # first-divergence runs: one sweep per bin
                          ("hub6x6_cb_fd", "hub6x6_cb", dict(NSweep=1, NBin=200)),
                          ("hub12x12_cb_fd", "hub12x12_cb", dict(NSweep=1, NBin=300)),
                          ("hub12x12_cb_b10_fd", "hub12x12_cb", dict(beta=10.0, NSweep=1, NBin=150)),
                          # precision runs (sizes chosen for <= 15 min per run on the A30: 1.2 s/sweep at beta = 4, 3.6 s/sweep at beta = 10, GPU path)
                          ("hub12x12_cb_p", "hub12x12_cb", dict(NSweep=20, NBin=20)),
                          ("hub12x12_cb_b10_p", "hub12x12_cb", dict(beta=10.0, NSweep=20, NBin=8)),
                          ("hub12x12_cb_b10_w5_p", "hub12x12_cb", dict(beta=10.0, Nwrap=5, NSweep=20, NBin=8))):
    d = {**CONFIGS[base], "NSweep": ns, "NBin": nb, **extra}
    Simulation(src, "Hubbard", d, sim_root=out, sim_dir=name, machine="GNU", mpi=True).run(only_prep=True)
    print(os.path.join(out, name))
PY

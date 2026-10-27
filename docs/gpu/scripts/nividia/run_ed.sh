#!/bin/bash
# usage: run_ed.sh BINARY_TAG TEST_NAME   (e.g. gnu-DEFAULT hubbard_finT_chain)
# Runs upstream testsuite/test_vs_ed with a prebuilt binary: upstream prep_test.py + analysis.py,
# simulation step replaced by pyALF only_prep + mpirun -n 2 (as CI's n_mpi=2), all dtau values in parallel.
tag=$1; test=$2; t=${tag%%-*}
. ~/alf-gpu/env/$t.sh
src=~/alf-gpu/src/ALF-$t; out=~/alf-gpu/runs/ed_$tag; rm -rf $out; mkdir -p ~/alf-gpu/runs
cd ~/alf-gpu/runs && python $src/testsuite/test_vs_ed/prep_test.py --specs_file $src/testsuite/test_vs_ed/test_specs.yaml --output_dir $out >/dev/null
d=$out/${test}_Bookworm; cd $d
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# HPC-X Open MPI 5 binds each mpirun -n 2 to cores 0,1 by default: five concurrent runs would share 2 cores
export PRTE_MCA_hwloc_default_binding_policy=none OMPI_MCA_hwloc_base_binding_policy=none
start=$(date +%s)
for i in */; do i=${i%/}
  ( cd $i && python - <<P && mpirun -n 2 ~/alf-gpu/bin/ALF-$tag.out > run.log 2>&1 ) &
import yaml; from py_alf import Simulation, ALF_source
s = yaml.safe_load(open("spec.yaml"))
Simulation(ALF_source(alf_dir="$src"), s["ham_name"], s["sim_dict"], sim_root=".", sim_dir=".", machine="GNU", mpi=True).run(only_prep=True)
P
done; wait
echo "sim wall $(( $(date +%s) - start )) s"
python $src/testsuite/test_vs_ed/analysis.py 2>&1 | grep -vE "^Analyzing|^ *$|Warning|warn" | tail -8

#!/bin/bash
# Adds the test-only Hamiltonians of testsuite/Lowrank.tests/models to a source tree, before "make":
#   Lowrank.tests/add_test_models.sh <ALF root>
# The files are copied to Prog/Hamiltonians and registered through Prog/Hamiltonians.list.d (the
# directory is read by the build, Prog/Hamiltonians.list itself stays unchanged). Then run
# compare_schemes.py with --test-models.
set -e
root=${1:?usage: add_test_models.sh <ALF root>}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$root/Prog/Hamiltonians.list.d"
: > "$root/Prog/Hamiltonians.list.d/lowrank_test_models"
for f in "$here"/models/Hamiltonian_*_smod.F90; do
  name=$(basename "$f" _smod.F90); name=${name#Hamiltonian_}
  cp "$f" "$root/Prog/Hamiltonians/"
  echo "$name" >> "$root/Prog/Hamiltonians.list.d/lowrank_test_models"
done

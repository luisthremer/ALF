#!/bin/bash
# usage: rr.sh 'remote cmd'   (filters the module noise)
timeout 90 ssh supermuc2 "$1" 2>&1 | grep -v -E "^Unloading|WARNING: Module|ERROR: Unload|EAR|SLURM_EAR|spack/release|^ +export SLURM|^$"

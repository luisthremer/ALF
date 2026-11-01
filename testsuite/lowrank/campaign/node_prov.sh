#!/bin/bash
# node_prov.sh OUTDIR REQUESTED_HOST   node-level provenance of the node this script runs on (bash 3.2, no python).
# Started by a packed job (job_template_packed.sh) once per host through the profile's launcher
# (profile_launch_on HOST 1 1 node_prov.sh ...), so it runs ON that host: lscpu, numactl, free, the cpufreq clock state.
# Writes OUTDIR/node.prov (KEY=VALUE) and the raw dumps lscpu.txt, numactl.txt, free.txt next to it.
# node / hostname are the OS hostname of the node the script really ran on (not SLURMD_NODENAME: through mpiexec that
# variable may be inherited from the batch host); the job compares it with the requested host.
# PROV_SYSFS_ROOT / PROV_PROC_ROOT come from the job's environment (the local test points them at fixtures).
out="$1"
want="$2"
. "$(cd "$(dirname "$0")" && pwd)/provenance.sh"
f="$out/node.prov"
prov_kv "$f" slot_host_requested "$want"
prov_kv "$f" node "$(hostname)"
prov_kv "$f" hostname "$(hostname)"
prov_kv "$f" slurmd_nodename_env "${SLURMD_NODENAME:-unset}"
prov_node_hw "$out" "$f"

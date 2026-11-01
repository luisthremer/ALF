#!/bin/bash
# node_prov.sh --freq RUN_PROV TAG   the clock state of this node (prov_freq, keys *_TAG) and freq_host_TAG into RUN_PROV
# node_prov.sh OUTDIR REQUESTED_HOST   node-level provenance of the node this script runs on (bash 3.2, no python).
# Started by a packed job (job_template_packed.sh) once per host through the profile's launcher
# (profile_launch_on HOST 1 1 node_prov.sh ...), so it runs ON that host: lscpu, numactl, free, the cpufreq clock state.
# Writes OUTDIR/node.prov (KEY=VALUE) and the raw dumps lscpu.txt, numactl.txt, free.txt next to it.
# node / hostname are the OS hostname of the node the script really ran on (not SLURMD_NODENAME: through mpiexec that
# variable may be inherited from the batch host); the job compares it with the requested host.
# PROV_SYSFS_ROOT / PROV_PROC_ROOT come from the job's environment (the local test points them at fixtures).
. "$(cd "$(dirname "$0")" && pwd)/provenance.sh"
if [ "$1" = --freq ]; then
    prov_freq "$2" "$3"
    prov_kv "$2" "freq_host_$3" "$(hostname)" || exit 1
    exit 0
fi
out="$1"
want="$2"
f="$out/node.prov"
prov_kv "$f" slot_host_requested "$want"
prov_kv "$f" node "$(hostname)"
prov_kv "$f" hostname "$(hostname)"
prov_kv "$f" slurmd_nodename_env "${SLURMD_NODENAME:-unset}"
prov_node_hw "$out" "$f"

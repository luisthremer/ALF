#!/bin/bash
#SBATCH --job-name=gb_topo
#SBATCH --output=gb_topo.%j.log
#SBATCH --time=00:02:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn93la
#SBATCH --partition=test
#SBATCH --nodes=1
module load slurm_setup >/dev/null 2>&1
hostname; lscpu; echo ====NUMACTL; numactl -H 2>&1 || cat /sys/devices/system/node/node*/cpulist
echo ====XPU; xpu-smi discovery 2>&1 | head -60; xpu-smi topology -m 2>&1 | head -40
echo ====PCI_NUMA; for d in /sys/class/drm/card*/device; do echo "$d numa=$(cat $d/numa_node 2>/dev/null) $(cat $d/device 2>/dev/null)"; done

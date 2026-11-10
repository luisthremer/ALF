#!/bin/bash
# usage: run_case.sh BIN NR NB NS TAG   (env passes through). L=24 production parameters, CPU_MAX=0 (fixed bins)
BIN=$1; NR=$2; NB=$3; NS=$4; TAG=$5
D=$RUNROOT/$TAG; rm -rf $D; mkdir -p $D; cd $D
sed -e 's/CPU_MAX              = 14.31/CPU_MAX              = 0.0/' -e "s/NBin                 = 100000/NBin                 = $NB/" -e "s/NSweep               = 50/NSweep               = $NS/" -e "s/N_delay              = 32/N_delay              = ${ND:-32}/" -e "s/L1           = 24/L1           = ${LL:-24}/" -e "s/L2           = 24/L2           = ${LL:-24}/" -e "s/N_Phi        = 288/N_Phi        = $(( ${LL:-24} * ${LL:-24} / 2 ))/" -e "s/Beta         = 24.d0/Beta         = ${LL:-24}.d0/" ~/alf-gpu/params_L24 > parameters
seq 2001 $((2000+NR)) > seeds
t0=$(date +%s.%N)
mpiexec -n $NR $BIN > run.log 2>&1; rc=$?
t1=$(date +%s.%N)
ct=$(grep -m1 "CPU Time" info | awk '{print $NF}')
echo "RESULT $TAG nr=$NR bins=$NB sweeps/bin=$NS rc=$rc wall=$(echo "$t1 - $t0" | bc) cputime_info=$ct sweeps_per_chain=$((NB*NS)) total_sweeps=$((NR*NB*NS))"

#!/bin/bash
# usage: prof.sh LABEL L NSWEEP NBIN scheme [nsys]
LABEL=$1; L=$2; NS=$3; NB=$4; sch=$5; NSYS=$6
. ~/alf-gpu/env/nvhpc.sh
export TMPDIR=~/alf-gpu/work/dlsub/tmp OMP_TARGET_OFFLOAD=MANDATORY
SRC=~/alf-gpu/work/dlsub/src
d=~/alf-gpu/work/dlsub/p_${LABEL}_$sch; rm -rf $d; mkdir -p $d
cp $SRC/Scripts_and_Parameters_files/Start/{parameters,seeds} $d/
sed -i "s/^Sequential_scheme .*/Sequential_scheme   = \"$sch\"/; s/^NSweep .*/NSweep = $NS/; s/^NBin .*/NBin = $NB/; s/^L1 .*/L1 = $L/; s/^L2 .*/L2 = $L/" $d/parameters
cd $d
s=$(date +%s.%N)
if [ -n "$NSYS" ]; then nsys profile -o prof --force-overwrite true --trace=cuda,openmp --stats=false $SRC/Prog/ALF.out > out.log 2>&1; nsys stats --report cuda_api_sum,cuda_gpu_kern_sum,cuda_gpu_mem_time_sum --format csv prof.nsys-rep > stats.csv 2>/dev/null; else $SRC/Prog/ALF.out > out.log 2>&1; fi
e=$(date +%s.%N); echo "$sch wall=$(echo "$e - $s" | bc)"

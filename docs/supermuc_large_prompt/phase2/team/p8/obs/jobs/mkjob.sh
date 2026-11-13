#!/bin/bash
# usage: mkjob.sh NAME TIME < body  -> writes jobs/NAME.sh (header + body)
n=$1; t=$2; d=$(dirname $0)
{ head -1 $d/hdr.txt; echo "#SBATCH --job-name=p8_$n"; echo "#SBATCH --time=$t"; tail -n +2 $d/hdr.txt | grep -v -- '--time'; cat; } > $d/$n.sh

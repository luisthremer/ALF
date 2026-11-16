#!/bin/bash
#SBATCH -J localtau_test
#SBATCH -p cpu
#SBATCH -A b133ae
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -t 03:00:00
#SBATCH -o /hnvme/workspace/b133ae29-stp_data/ALF_Luis/test_localtau/job_%j.log
cd /hnvme/workspace/b133ae29-stp_data/ALF_Luis/test_localtau
D=/hnvme/workspace/b133ae29-stp_data/DATA/JAKOB/CPU_NEW
lscpu | grep "Model name"
# fast: 120 bins, every L
./run_test.sh $D/QFI/L4/Temp_0  $PWD/fast/L4  120
for L in 8 12 16 20 24; do ./run_test.sh $D/FSS/L$L/Temp_0 $PWD/fast/L$L 120; done
# full size
./run_test.sh $D/QFI/L4/Temp_0      $PWD/full/QFI_L4      0
./run_test.sh $D/FSS/L8/Temp_0      $PWD/full/FSS_L8      0
./run_test.sh $D/FSS/L24/Temp_0     $PWD/full/FSS_L24     0
# full size, ana_hdf5.out without arguments (all observables)
./run_test.sh $D/Spectral/L20/Temp_0 $PWD/full/Spectral_L20_noarg 0 1
echo JOB DONE

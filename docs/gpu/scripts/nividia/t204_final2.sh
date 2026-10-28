#!/bin/bash
# Ticket 204, nividia: rebuild nv_def / nv_log from the final source, self-tests of all configs, negative control vexp, timings.
W=$HOME/alf-gpu/work/t204; cd $W; rm -f logs/final2_done
bash t204_build.sh src_new nv_def nvhpc PGI noMPI HDF5 GPU=nvidia > logs/g_nv.txt 2>&1
bash t204_build.sh src_new nv_log nvhpc PGI noMPI HDF5 LOG GPU=nvidia > logs/g_log.txt 2>&1
C6="hub6x6_cb hub6x6_dense hub6x6_cb_flux hub6x6_dense_flux hub6x6_cb_mzoff hub6x6_cb_symm hub12x12_cb"
bash t204_selftest_all.sh nv_def $C6; bash t204_selftest_all.sh nv_log $C6
( bash t204_selftest_all.sh nv_def2 hub20x20_cb >/dev/null 2>&1 ) &
cp bin/ALF_nv_def.out bin/ALF_nv_def2.out
TAG=nc_vexp bash t204_run.sh nv_def hub6x6_cb ALF_GPU_SELFTEST=2 ALF_GPU_FAULT=vexp ALF_GPU_EXIT=1 > logs/negctl2.txt
bash t204_timing_all.sh nv_def final hub12x12_cb hub12x12_dense hub20x20_cb hub20x20_dense hub32x32_cb hub32x32_dense
wait
echo done > logs/final2_done

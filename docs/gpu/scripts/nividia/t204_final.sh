#!/bin/bash
# Ticket 204, nividia: final binaries from src_new, quick self-tests, identity run GPU=nvidia vs GPU=none, GNU MPI nm check, ctest.
W=$HOME/alf-gpu/work/t204; cd $W
bash t204_build.sh src_new nv_def nvhpc PGI noMPI HDF5 GPU=nvidia > logs/f_nv.txt 2>&1
bash t204_build.sh src_new nv_none nvhpc PGI noMPI HDF5 GPU=none > logs/f_none.txt 2>&1
bash t204_selftest_all.sh nv_def hub6x6_cb hub6x6_dense hub6x6_cb_flux hub6x6_dense_flux hub6x6_cb_mzoff hub6x6_cb_symm hub12x12_cb
for e in hop hopc vexp nan order; do c=hub6x6_cb; [ $e = hopc ] && c=hub6x6_cb_flux; TAG=nc_$e bash t204_run.sh nv_def $c ALF_GPU_SELFTEST=2 ALF_GPU_FAULT=$e ALF_GPU_EXIT=1 >> logs/negctl.txt; done
bash t204_run.sh nv_none hub6x6_cb > logs/ident.txt; bash t204_run.sh nv_def hub6x6_cb >> logs/ident.txt
bash t204_build.sh src_head gnu_head gnu GNU MPI HDF5 > logs/f_gh.txt 2>&1
bash t204_build.sh src_new gnu_new gnu GNU MPI HDF5 > logs/f_gn.txt 2>&1
bash t204_ctest.sh > logs/f_ctest.txt 2>&1
echo done > logs/final_done

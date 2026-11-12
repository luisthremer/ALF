## Per-kernel table, V7, 1 rank on tile 0, per sweep (C: jobs 560831 / 560845 / 560846)
| # | phase | kernel / op | calls/sweep | device s/sweep | mean us | launch (teams x threads) |
|--|--|--|--|--|--|--|
| 1 | SC | devqr geqrf (sync call) | 192 | 0.8220 | 4281.3 |  |
| 2 | C | getrf oneMKL (hybrid, sync call) | 96 | 0.3785 | 3942.7 |  |
| 3 | U | sycl gemm m576 | 34560 | 0.3119 | 9.0 |  |
| 4 | O | omp gpu_obs.corr__l437 | 480 | 0.3004 | 625.8 | 1x1(480) |
| 5 | SC | devqr unmqr (sync call) | 192 | 0.2625 | 1367.2 |  |
| 6 | O | omp gpu_obs.corr__l473 | 480 | 0.1662 | 346.3 | 1x1(480) |
| 7 | O | omp gpu_obs.corr__l450 | 480 | 0.1531 | 319.0 | 512x1024(480) |
| 8 | O | mkl zgemm 576x576x576 tNN [gpu_obs_mod_mp_symm_apply_] | 960 | 0.1049 | 109.2 |  |
| 9 | O | mkl zgemm 576x576x576 tNN [gpu_obs_mod_mp_symm_apply_] | 960 | 0.1013 | 105.5 |  |
| 10 | U | sycl scale | 8640 | 0.0933 | 10.8 |  |
| 11 | O | omp gpu_obs.corr__l418 | 480 | 0.0897 | 186.9 | 512x1024(480) |
| 12 | U | sycl fgather | 17280 | 0.0838 | 4.8 |  |
| 13 | O | omp gpu_obs.corr__l386 | 480 | 0.0680 | 141.6 | 512x1024(480) |
| 14 | C | omp gpu_stab.k_cgr_phase__l1232 | 96 | 0.0585 | 609.7 | 1x64(98) |
| 15 | U | mkl zgemm 576x576x576 tNN [gpu_wrap_mod_mp_gemm_left_] | 480 | 0.0530 | 110.4 |  |
| 16 | U | mkl zgemm 576x576x576 tNN [gpu_wrap_mod_mp_hopd_left_] | 480 | 0.0529 | 110.3 |  |
| 17 | U | mkl zgemm 576x576x576 tNN [gpu_wrap_mod_mp_hopd_right_] | 480 | 0.0529 | 110.1 |  |
| 18 | U | mkl zgemm 576x576x576 tNN [gpu_wrap_mod_mp_hopd_left_] | 480 | 0.0518 | 107.8 |  |
| 19 | U | mkl zgemm 576x576x576 tNN [gpu_wrap_mod_mp_hopd_right_] | 480 | 0.0517 | 107.6 |  |
| 20 | C | omp gpu_stab.k_prec__l1327 | 96 | 0.0508 | 528.7 | 4x1024(96) |
| 21 | SC | omp alf_la_om.k_colnorms__l845 | 192 | 0.0444 | 231.3 | 288x2(242) |
| 22 | U | sycl gather | 8640 | 0.0430 | 5.0 |  |
| 23 | C | mkl ztrsm 576x576 LUNN [gpu_stab_mod_mp_cgr_solve_] | 50 | 0.0386 | 771.0 |  |
| 24 | C | omp gpu_stab.k_det_phase__l1207 | 96 | 0.0337 | 350.6 | 1x64(98) |
| 25 | U | sycl copy | 43200 | 0.0290 | 0.7 |  |
| 26 | C | mkl ztrsm 576x576 RUCN [gpu_stab_mod_mp_cgr_solve_] | 46 | 0.0279 | 606.5 |  |
| 27 | S | omp gpu_stab.k_phase__l932 | 96 | 0.0256 | 267.1 | 1x64(144) |
| 28 | U | omp gpu_wrap.k_vslice_r__l980 | 960 | 0.0174 | 18.1 | 1152x64(1440) |
| 29 | O | omp gpu_obs.k_cpy__l207 | 2880 | 0.0169 | 5.9 | 5184x64(2880) |
| 30 | O | omp gpu_obs.gpu_obs_measure__l284 | 480 | 0.0162 | 33.8 | 24x1024(480) |

device busy per phase from event timestamps (omp kernels + sycl block ops + dispatched BLAS), excluding QR/getrf sync calls: {'U': 0.914, 'O': 1.047, 'C': 0.262, 'SC': 0.056, 'S': 0.068}

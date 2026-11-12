### Top 20 kernels / ops, V7, one rank on tile 0, per sweep (C; jobs 560831 pp, 560845 off, 560847 callsites)

| # | phase | kernel / op | calls/sweep | device s/sweep | mean us | launch size |
|--|--|--|--|--|--|--|
| 1 | S+C | devqr geqrf call (sync; panel 1 work-group 60%, apply 40%) | 192 | 0.822 | 4281.3 | panel: 1 WG x 512 WI; apply: ~70 WG x 256 WI |
| 2 | C | getrf oneMKL call (hybrid host/device, sync) | 96 | 0.379 | 3942.7 | oneMKL, partly host |
| 3 | U | mkl zgemm 576x576x576 tNN | 2880 | 0.315 | 109.4 | oneMKL kernel(s), full tile |
| 4 | U | sycl gemm m576 | 34560 | 0.312 | 9.0 | oneMKL 576x576xk, k=2..23 (mostly 10-14) |
| 5 | O | omp gpu_obs:corr l437 | 480 | 0.300 | 625.8 | 1x1(480) |
| 6 | S+C | devqr unmqr call (sync) | 192 | 0.262 | 1367.2 | ~72 WG x 256 WI |
| 7 | O | mkl zgemm 576x576x576 tNN | 1920 | 0.206 | 107.3 | oneMKL kernel(s), full tile |
| 8 | O | omp gpu_obs:corr l473 | 480 | 0.166 | 346.3 | 1x1(480) |
| 9 | O | omp gpu_obs:corr l450 | 480 | 0.153 | 319.0 | 512x1024(480) |
| 10 | U | sycl scale | 8640 | 0.093 | 10.8 | global 2xbx576, mean 36,864 WI |
| 11 | O | omp gpu_obs:corr l418 | 480 | 0.090 | 186.9 | 512x1024(480) |
| 12 | U | sycl fgather | 17280 | 0.084 | 4.8 | global kx576, mean 6,045 WI |
| 13 | O | omp gpu_obs:corr l386 | 480 | 0.068 | 141.6 | 512x1024(480) |
| 14 | C | omp gpu_stab:k_cgr_phase l1232 | 96 | 0.059 | 609.7 | 1x64(98) |
| 15 | C | omp gpu_stab:k_prec l1327 | 96 | 0.051 | 528.7 | 4x1024(96) |
| 16 | SC | omp alf_la_om:k_colnorms l845 | 192 | 0.044 | 231.3 | 288x2(242) |
| 17 | U | sycl gather | 8640 | 0.043 | 5.0 | global 2x32x32 = 2048 WI (~64 sub-groups) |
| 18 | C | mkl ztrsm 576x576 LUNN | 50 | 0.039 | 771.0 | oneMKL kernel(s), full tile |
| 19 | C | omp gpu_stab:k_det_phase l1207 | 96 | 0.034 | 350.6 | 1x64(98) |
| 20 | U | sycl copy | 43200 | 0.029 | 0.7 | DMA copy (mean 13.7 kB) |

### Per phase, top 6 each

**phase U**

| kernel / op | calls/sweep | dev s/sweep | mean us | launch |
|--|--|--|--|--|
| mkl zgemm 576x576x576 tNN | 2880 | 0.315 | 109.4 | oneMKL kernel(s), full tile |
| sycl gemm m576 | 34560 | 0.312 | 9.0 | oneMKL 576x576xk, k=2..23 (mostly 10-14) |
| sycl scale | 8640 | 0.093 | 10.8 | global 2xbx576, mean 36,864 WI |
| sycl fgather | 17280 | 0.084 | 4.8 | global kx576, mean 6,045 WI |
| sycl gather | 8640 | 0.043 | 5.0 | global 2x32x32 = 2048 WI (~64 sub-groups) |
| sycl copy | 43200 | 0.029 | 0.7 | DMA copy (mean 13.7 kB) |

**phase S**

| kernel / op | calls/sweep | dev s/sweep | mean us | launch |
|--|--|--|--|--|
| devqr geqrf call (sync; panel 1 work-group 60%, apply 40%) | 192 | 0.822 | 4281.3 | panel: 1 WG x 512 WI; apply: ~70 WG x 256 WI |
| devqr unmqr call (sync) | 192 | 0.262 | 1367.2 | ~72 WG x 256 WI |
| omp gpu_stab:k_phase l932 | 96 | 0.026 | 267.1 | 1x64(144) |
| mkl zgemm 576x576x576 tNN | 144 | 0.017 | 119.1 | oneMKL kernel(s), full tile |
| mkl zgemm 576x576x576 tCN | 96 | 0.013 | 137.5 | oneMKL kernel(s), full tile |
| mkl zgemm 576x576x576 tNC | 48 | 0.005 | 111.5 | oneMKL kernel(s), full tile |

**phase C**

| kernel / op | calls/sweep | dev s/sweep | mean us | launch |
|--|--|--|--|--|
| devqr geqrf call (sync; panel 1 work-group 60%, apply 40%) | 192 | 0.822 | 4281.3 | panel: 1 WG x 512 WI; apply: ~70 WG x 256 WI |
| getrf oneMKL call (hybrid host/device, sync) | 96 | 0.379 | 3942.7 | oneMKL, partly host |
| devqr unmqr call (sync) | 192 | 0.262 | 1367.2 | ~72 WG x 256 WI |
| omp gpu_stab:k_cgr_phase l1232 | 96 | 0.059 | 609.7 | 1x64(98) |
| omp gpu_stab:k_prec l1327 | 96 | 0.051 | 528.7 | 4x1024(96) |
| mkl ztrsm 576x576 LUNN | 50 | 0.039 | 771.0 | oneMKL kernel(s), full tile |

**phase O**

| kernel / op | calls/sweep | dev s/sweep | mean us | launch |
|--|--|--|--|--|
| omp gpu_obs:corr l437 | 480 | 0.300 | 625.8 | 1x1(480) |
| mkl zgemm 576x576x576 tNN | 1920 | 0.206 | 107.3 | oneMKL kernel(s), full tile |
| omp gpu_obs:corr l473 | 480 | 0.166 | 346.3 | 1x1(480) |
| omp gpu_obs:corr l450 | 480 | 0.153 | 319.0 | 512x1024(480) |
| omp gpu_obs:corr l418 | 480 | 0.090 | 186.9 | 512x1024(480) |
| omp gpu_obs:corr l386 | 480 | 0.068 | 141.6 | 512x1024(480) |

device time from event timestamps per phase [s/sweep] (omp kernels + sycl block ops + oneMKL BLAS, QR/getrf calls excluded): {'O': 1.047, 'C': 0.262, 'SC': 0.056, 'U': 0.914, 'S': 0.068}
SC (shared stab/cgr kernels) split ~50/50 by call count.

### Node layout (32 ranks, 4/tile, CCS 4), V7: omp kernels, mean over the 32 ranks vs one rank alone (C: job 560841 vs 560831)

| kernel | calls/sweep | alone ms/sweep | 4/tile ms/sweep | slow-down | launch |
|--|--|--|--|--|--|
| gpu_obs.corr__l450 | 480 | 153.1 | 590.9 | 3.86x | 512x1024(480) |
| gpu_obs.corr__l437 | 480 | 300.4 | 316.2 | 1.05x | 1x1(480) |
| gpu_obs.corr__l418 | 480 | 89.7 | 297.5 | 3.32x | 512x1024(480) |
| gpu_obs.corr__l386 | 480 | 68.0 | 238.3 | 3.51x | 512x1024(480) |
| gpu_obs.corr__l473 | 480 | 166.2 | 178.6 | 1.07x | 1x1(480) |
| gpu_stab.k_cgr_phase__l1232 | 96 | 58.5 | 61.6 | 1.05x | 1x64(98) |
| gpu_stab.k_prec__l1327 | 96 | 50.8 | 59.7 | 1.18x | 4x1024(96) |
| alf_la_om.k_colnorms__l845 | 192 | 44.4 | 48.3 | 1.09x | 288x2(242) |
| gpu_obs.k_cpy__l207 | 2880 | 16.9 | 44.3 | 2.62x | 5184x64(2880) |
| gpu_wrap.k_cpy_n__l308 | 2400 | 16.0 | 38.4 | 2.39x | 5184x64(2400) |
| gpu_wrap.k_vslice_r__l980 | 960 | 17.4 | 38.1 | 2.19x | 1152x64(1440) |
| gpu_stab.k_det_phase__l1207 | 96 | 33.7 | 35.2 | 1.05x | 1x64(98) |

omp-kernel device time per phase and sweep, per rank at 4/tile: {'O': 1.718, 'C': 0.186, 'SC': 0.077, 'U': 0.088, 'S': 0.042}

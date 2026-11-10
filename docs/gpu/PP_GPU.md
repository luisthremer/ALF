# Parallel parameters (PP) together with the GPU path

Branch `gpu-delayed-submatrix-defect-production-mpi-pp` (from `...-production-mpi`). Build: `. ./configure.sh PGI PP GPU=nvidia HDF5 NO-INTERACTIVE`
(`PP` = `-DMPI -DTEMPERING -DPARALLEL_PARAMS`).

Change: `gpu_state_mod.F90` refused every `TEMPERING` build; it now refuses only real tempering (`TEMPERING` without `PARALLEL_PARAMS`).
With PP the groups are independent runs (own `Temp_<i>/parameters`, own device mirror per rank, no exchange step: `main.F90` skips
`Exchange_Step` under `PARALLEL_PARAMS`), so nothing else changed. Output goes to `Temp_<i>/` (`data.h5`, `info`, `confout_*`).

Setup of a raster: `scripts/pp/prepare_sim.py` (from `luis_playground`, `Scripts_and_Parameters_files/Start_PP`) expands `{a, b, c}` lists in
`parameters` into `Temp_<i>/`. It needs `ham_name = {"Hubbard"}` as a list, otherwise the `var_hubbard` patch is dropped (one `Temp_0` without the
varied value). `mpi_per_parameter_set` in `&VAR_TEMP` is NOT set by the script; set it by hand (`pp_setup.py` does).
Ranks = `mpi_per_parameter_set` x number of parameter sets.

## Measured 2026-11-10 (nividia, A30, CUDA MPS, pi-flux Hubbard L = 24, Beta = 24, N_FL = 1, submatrix, 3 U values)

Wall time of 6 sweeps including ~5-8 s start-up (`scripts/pp/pp_bench.sh`):

| ranks | MPS thread share | `N_delay` | wall (s) |
|---|---|---|---|
| 3 (1 walker per U) | 100 % | 32 | 29.1 |
| 9 (3 walkers per U) | 100 % | 32 | 64.4 |
| 9 | 50 % | 32 | 59.8 |
| 9 | 25 % | 16 / 32 / 64 / 128 | 65.8 / 57.0 / 53.6 / 52.9 |

One chain at 3 ranks: ~4.1 s per sweep; at 9 ranks (25 %, `N_delay` 64): ~7.6 s per sweep. Three walkers per U give about 1.5x the sweeps per U per
second of one walker per U, at the price of 9 host cores (project rule: 4) and a 3x slower single chain. All three groups report
GPU path ACTIVE and device measurement ACTIVE; Precision Green 4.5e-15.

MPS thread share, repeated (9 ranks, `N_delay` 64, 10 sweeps, wall in s, two repeats each): 12 %: 96.5 / 96.4; 25 %: 86.9 / 86.9;
50 %: 91.2 / 91.1; 100 %: 98.3 / 98.6. The cap is per client (9 x 25 % oversubscribes the device); it reduces contention between the ranks and is
12 % faster than no cap. Production run: `NSweep = 10` per bin, `N_delay` 64, 25 %, `CPU_MAX` 23.8 h (`scripts/pp/pp_prod.sh`).

## Profile of the production run and optimisation round (2026-11-10, branch `...-production-mpi-pp-opt`)

nsys of all 9 ranks (L = 24, `N_delay` 64, MPS 25 %, 150 s window; `scripts/pp/pp_prod_prof.sh`, `union.py`): at least one kernel is running during 99.97 % of
the window and on average 5.6 kernels run at the same time, so the device is the limit, not the host. Time of rank 0 in kernels: flush GEMMs ~47 %,
stabilisation (QR, LU, trsm, ...) ~25 %, measurement kernels ~17 %, block kernels ~6 %, wraps ~5 %; 72 % of the host time is `cuStreamSynchronize`.
One rank spends ~4300 blocks (`N_delay` 64) and ~44 000 syncs per sweep.

Changes: (1) the scalar density sums in `corr` ran as a serial loop in one GPU thread (760 us per call, 7 % of the kernel time); now a parallel
reduction. Results identical to round-off (SpinZ_eq 1e-16, Ener_scal 1e-13 vs the baseline) but NO speed-up (88.9 s vs 88.9 s: the kernel occupies one
thread, the device throughput is what limits). (2) `ALF_GPU_OBS_SKIP=1,3` skips Green_eq and Den_eq on the device (they stay zero); the env guard requires
`ALF_GPU_DEBUG=1` for it, a production switch is still to be decided.

Wall time for 10 sweeps, 9 ranks, 25 %, `N_delay` 64 (incl. ~8 s start-up): baseline 88.9 / 89.7 s; opt 88.9; opt + skip Green/Den 82.5; opt + `Nwrap` 20 78.5;
opt + skip + `Nwrap` 20 71.8 (-19 %, -21 % without start-up). Precision Green mean / max, `Nwrap` 10 -> 20: U = 3: 4e-15 / 1.6e-11 -> 9e-14 / 4.4e-10;
U = 5.25: 2e-14 / 2.6e-11 -> 1.5e-13 / 2.2e-10; U = 8: 1.5e-13 / 5e-10 -> 7.6e-13 / 4.8e-9. Observables after 10 sweeps agree with `Nwrap` 10 to <= 4e-12 relative.

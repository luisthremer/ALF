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

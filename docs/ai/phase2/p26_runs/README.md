# P2.6.1 evidence: delayed / submatrix updates for every model (GPU=host, melmac)

Build: `configure.sh GNU noMPI GPU=host`, branch `cm/dev/p26-models`. Base for "unchanged" checks: `Code-modularity-GPU` 5c0af0a4.

| File | What |
|---|---|
| `cases_p26.py` | Kondo (2-layer 4x4), tV (4x4, Beta 2), Spin_Peierls (6x6, type-4 J vertices) x standard / delayed N_delay 4, 16 / submatrix 4; tV with Ltau = 1; Hubbard regression |
| `results_gpu_host.tsv` | `../p25_runs/p25_run.py` on these cases: GPU path ACTIVE against ALF_GPU_DISABLE=1 (CPU path), same seeds: confout byte-identical, every other file within the e2e rule (100 x Precision Green mean, at least 1e-11) |
| `k4.py`, `k4_schemes.txt` | gate K4: standard vs delayed vs submatrix on the GPU path, bin for bin (threshold 1e-9); confout identical |
| `plant.sh`, `planted_defects.txt` | planted defects (built by hand, never committed); each must fail a check: drop N_Type 2, wrong type-4 factor index, W transposed, no sqrt in the type-4 Phi, X/Y of the delayed path not wrapped (Y) |

Not in these files (reported in the commit messages / the report to the manager): K1 (`ALF_GPU_SELFTEST=2`, device kernels against
Operator_mod: Kondo and Spin_Peierls 0, tV <= 4.6e-17), Hubbard outputs byte-identical between the old and the new GPU=host binary
(hub-mz-4x4, hub-dly4, hub-sub4) and the ctest of the GPU=host build.

## H100 (Helma, GPU=nvidia, nvfortran 25.7, job 1467539, 12 min)

Build of p26-7 (cce147cf; p26-8 changes only an argument list) and of the base 5c0af0a4. `cases_h100.py`: the small cases above plus
`big-*` (Kondo 8x8 bilayer Ndim 128, tV 12x12, Spin_Peierls 12x12; NSweep 10, NBin 2). `results_h100.tsv`, `k4_h100.txt`.

- 24 of 28 cases OK (GPU path vs CPU path of the same binary, confout identical everywhere).  `big-kondo-*` (all four schemes):
  Dimer_tau differs by 1.26e-11 against the floor 1e-11 of the e2e rule (Precision Green mean 7e-14), confout identical: round-off of
  the device tau pass at this size, KNOWN_ISSUES G9 (C).
- K4 on the H100: schemes bin for bin <= 3.3e-11, confout identical.
- Hubbard, old vs new binary: confout identical; Kin_scal (hub-mz-4x4) and Ener_scal (hub-dly4) differ by one ulp (3.6e-15 abs):
  the device team reductions of gpu_obs_mod (not changed here); not confirmed by a repeat run of the old binary.
- Timings (s, GPU path / CPU path single core, whole run): big-kondo std 23.9 / 10.1, delayed 16 31.0, delayed 32 31.8, submatrix
  16 25.4; big-tV std 39.8 / 6.2, delayed 16 61.9, delayed 32 62.3, submatrix 16 40.0; big-sp std 34.6 / 13.5, delayed 16 47.6,
  delayed 32 49.4, submatrix 16 35.6.  At these sizes the GPU path is 2.4 - 10x slower than one CPU core and the delayed device path
  is slower than standard (per-vertex launches and transfers dominate; Hubbard 4x4 is 2x slower too).  Scaling with Ndim: P2.6.2.

## p26-langevin (part A: field type 3 and Langevin force pass on the device), 2026-11-27/28

Base `Code-modularity-GPU` bdc30bf9; code commits rebased from `mac/p26-langevin-wip` c8300128. Cases: `cases_langevin.py` (10), `cases_type3.py` (10),
`cases_p26.py`, `cases_calcfl.py` (together 28), runner `../p25_runs/p25_run.py` (`P25_CONFOUT_NUMERIC=1`: confout under the tolerance rule, the fields are continuous).
Precision Green ratio: `pgratio.py` (GPU mean / CPU mean per case).

| Machine | Result | Precision Green ratio |
|---|---|---|
| GPU=host (melmac) | ctest 159/159; 48/48 cases OK (Langevin max scaled diff 1.6e-12) | 0.82 - 1.40 |
| H100 (Helma, nvhpc 25.7, job 1475700) | 48/48 OK | 0.84 - 1.84, except lang-mz-12x12 2.31 (8.9e-13 vs 3.9e-13) and lrc-lang-12x12 2.68 (1.2e-14 vs 4.6e-15): ratios of numbers near machine precision, not G10 (KNOWN_ISSUES G7) |
| MI210 (amd, gfortran amdgcn + HIP) | 48/48 OK | 0.60 - 1.22 (lang-mz-12x12 0.60, lrc-lang-12x12 0.68) |
| GPU=none | byte-identical to the base on the 16 bitwise cases (`bw_gpu_amd.py`) | n/a |

`lang-su2-4x4` with NSweep 4 failed once on both H100 (Den_eq, Ener_scal 1.5e-10) and MI210 (Pot_scal 5.3e-11 against 1e-11). CPU-only nulls of the same binary
(Nwrap 10 -> 5; H100 job 1473891: Den_eq 3.6e-10 against GPU 1.6e-10; MI210: Pot_scal 4.4e-11) are as large: round-off amplification of Langevin dynamics, KNOWN_ISSUES G9.
The case now runs NSweep 2.

Outside the agreed area (named in the commit bodies): `qmc_setup_mod.F90` (one line, `Call gpu_lang_register()`), `gpu_obs_mod.F90` (device observables refused
for Langevin with Ltau = 1: Tau_m measures on the host), `p25_run.py`. HMC and MALA stay refused (part B). GPU=intel refuses type 3 and Langevin.

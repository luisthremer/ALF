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

## Max 1550 (SuperMUC-NG Phase 2, GPU=intel, ifx 2025.3 + oneMKL, job 566852)

AOT build of Code-modularity-GPU d1e0dce6 (Intel switches at their defaults), `cases_h100.py` spread over the 8 tiles
(`../intel/job_p26.sbatch`), the same runner and rule. `results_max1550.tsv`; CPU floor and Kondo bisection in `floor_max1550.txt`
(`../intel/floor.py`, `../intel/xb.py`, melmac ifx). Hubbard 4x4 sanity cases OK.

- 19 of 28 OK; confout identical in all 28. The 9 others exceed the e2e rule at round-off level:
  - `big-kondo-*` (all four schemes): Dimer_tau 6.61e-10 (Den_tau, Greenf_tau 3.9e-11), tol 1e-11. Not round-off of the device: the CPU
    floor (ifx, MKL_CBWR AUTO vs COMPATIBLE) is 6.95e-12, but a CPU-only build with `-DQR_PREPIVOT` against the default CPU build gives
    5.99e-10 and ifx GPU=host 3.28e-10. The column-norm pre-pivoted QR of the GPU path (instead of ZGEQP3) is the cause, not oneMKL or
    the device; how far it moves Dimer_tau depends on the BLAS (H100 1.26e-11). KNOWN_ISSUES G9 (C), family of G7.
  - `big-tV` delayed 16 / 32, submatrix 16: 2.4e-11 .. 5.5e-11 against tol 1.6e-11 .. 2.1e-11; the CPU floor itself is 1.7e-11 .. 3.0e-11
    (above its own tol for dly16 and sub16).
  - `sp-std` 1.38e-10 (tol 9e-11; floor 4.3e-11), `big-sp-sub16` 5.64e-10 (tol 4.6e-10; floor 6.0e-10).
  - Class C: the tolerance rule (100 x Precision Green mean, at least 1e-11) is at or below the CPU's own rounding noise for these
    sizes; it needs the size scaling of G7 / G9.
- Timings (s, GPU path / CPU path single core of the same binary, whole run): big-kondo std 41.0 / 15.1, delayed 16 59.5, delayed 32
  63.0, submatrix 16 41.2; big-tV std 77.6 / 15.4, delayed 16 138.1, delayed 32 131.6, submatrix 16 85.6; big-sp std 60.6 / 19.2,
  delayed 16 101.6, delayed 32 106.2, submatrix 16 63.5; small cases 5 - 24 s against 2 - 3 s. At these test sizes the GPU path is
  2.7 - 5.6x (standard, submatrix) and up to 8.3x (delayed) slower than one CPU core (Max 1550 tile, 8 runs per node; as on the H100), and delayed is slower than standard and submatrix.
  Input for P2.6.2: the crossover size per model.

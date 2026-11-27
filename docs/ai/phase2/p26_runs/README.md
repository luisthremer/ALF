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
(hub-mz-4x4, hub-dly4, hub-sub4), ctest of the GPU=host build, and the H100 runs.

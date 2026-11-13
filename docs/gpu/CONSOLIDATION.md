# How `gpu-port-unified` was put together (2026-11-13)

`gpu-port-unified` starts at `gpu-delayed-submatrix-defect-production-mpi` (78409721), the last
commit shared by the three vendor lines, and merges every line with `--no-ff`. No branch was
changed or deleted; the history of each line is intact below its merge commit.

## Merged

| Order | Branch (tip) | Brings | Conflicts |
|---|---|---|---|
| 1 | `gpu-delayed-submatrix-defect-production-mpi-pp-helma` (cc473abc) | parallel tempering on the GPU path, `HELMAGPU` machine, Helma scripts and report | none |
| 2 | `gpu-delayed-submatrix-defect-production-mpi-pp-opt` (810a9122) | parallel reduction of the `acc0` density sums, opt-in `ALF_GPU_OBS_SKIP`, PP profiling scripts; plus two uncommitted script files from its worktree (commit 8f83cbdb) | none |
| 3 | `production-mpi-intel-unstable-guard` (97fce720, local tip) | `ALF_UNSTABLE_THRESHOLD`, SuperMUC production scripts and ledger; contains `production-mpi-intel-large-prompt` | none |
| 4 | `production-mpi-intel-large-prompt-team-combined` (7e8b6032) | device Householder QR (SYCL), async block path, explicit-shape kernels, team measurements; contains `-team`, `-team-async`, `-team-sycqr` | none |
| 5 | `intel-supermuc-ng2` (c3a91f72) | first SuperMUC-NG Phase 2 session: build scripts, AOT wrapper, layout benchmarks | `alf_la_mod.F90`, two jobfiles |
| 6 | `gpu-delayed-submatrix-defect-production-mpi-amd` (156a0b01) | `GPU=amd`: HIP backend, HIP kernels, device QR, AMD notes and patches | `gpu_lowrank_mod.F90`, `gpu_obs_mod.F90`, `configure.sh` |
| 7 | `origin/luis_playground` (b86df3ad), on branch `gpu-port-unified-maxent` | Gaussian bootstrap of the stochastic MaxEnt (`N_boot`, `Aom_ps_boot_n`) in `Analysis/Max_SAC.F90`, `maxent_stoch_mod`, test `2-maxsac-bootstrap`, `Start_PP` template; scratch files removed after the merge; documented in `Documentation/maxent.tex` (Sec. "Statistical errors of the spectrum: Gaussian bootstrap") and `files.tex` | none |

Contained in these and therefore not merged separately: `master`, `JAKOB_BACKUP`,
`gpu-delayed-submatrix-defect`, `-fullgpu`, `-production`, `-production-mpi`, `-production-mpi-pp`.
`production-mpi-intel-unstable-guard-asyncblk-test` is only a merge of 3 and 4 at older tips.

## How the conflicts were resolved

- **`Libraries/Modules/alf_la_mod.F90` (5):** both Intel sessions declared `omp_get_interop_int/ptr`
  through `bind(C)` (ifx <= 2025.3 lacks them in `omp_lib`). The large-prompt version is kept; it is
  the one in every production Intel binary.
- **SuperMUC jobfiles (5):** the large-prompt versions keep their names; the ng2 versions are
  `Scripts_and_Parameters_files/JobfileSuperMUC-NG-Phase2-GPU{,-PP}-ng2.sh`.
- **`Prog/gpu_lowrank_mod.F90` (6):** the Intel kernels are the base (explicit-shape dummy
  arguments, fused `k_blk_scale`, async block code). The HIP calls (`hk_flush_gather`,
  `hk_blk_gather`, `hk_blk_end`, `blk_flush(..., gathered)`) are re-applied on them under
  `ALF_LA_HIP`. On HIP the block end still uses the fused `hk_blk_end` kernel instead of `k_blk_scale`.
- **`Prog/gpu_obs_mod.F90` (6):** the AMD line had rewritten common code (chunked `corr`, real and
  imaginary parts reduced separately) for the libgomp-GCN runtime, and later replaced it on HIP by
  the `hk_*` kernels. That rewrite would have changed the summation order on NVIDIA and Intel
  without a measurement there. The unified file keeps the pp-opt code for every non-HIP backend
  and calls `hk_scalars`, `hk_local`, `hk_corr`, `hk_copyz` under `ALF_LA_HIP`. The chunked `corr`
  remains on the amd branch (see SPEEDUPS.md).
- **`configure.sh` (6):** usage text lists `HELMAGPU` and `GPU=amd`.

Changes made after the merges (commit 609a6347):
- `alf_exit(code)` in `runtime_error_mod` replaces `stop`, `error stop` (rejected by the CI lint
  `lint_stopcommands`) and the non-standard `call exit` (does not link with `gfortran -std=f2008`)
  in six places. Exit codes are unchanged.
- `la_backend_name` returns `'hip'` for `GPU=amd`.

## Not merged, on purpose

| Branch | Reason |
|---|---|
| `gpu-delayed-submatrix-defect-production-old-Ntau` (8438cf0f) | Variant with the older tau labelling of `SpinZ_localtau`, kept for comparison with older data. It differs from production in 4 files and must not flow into production. |
| `origin/gpu-port-207-4`, `-207-harness`, `-207-wip`, `-212-hardening`, `-302b-docs`, `origin/feature/delayed-submatrix*` | Early ticket branches with a few commits that never reached the port. |

## Tree changes after the merge

- `docs/gpu/README.md` is the entry page. Platform pages are in `docs/gpu/platforms/`
  (`HELMA_H100.md`, `SUPERMUC_NG2.md` and `docs/gpu/amd/` moved there).
- `docs/supermuc_large_prompt/` moved to `docs/gpu/sessions/supermuc_large_prompt/`; path
  references in scripts and documents were updated. Production runs that started from the older
  branches use their own checkouts and are not affected.
- New: `RUNTIME_SWITCHES.md`, `SPEEDUPS.md`, this file.

## Differences from upstream ALF (`master` 521ddfca) in a `GPU=none` build

Checked on 2026-11-13 with gfortran, L = 4, same parameters and seeds. With `LOBS_ST = 1` set
explicitly every upstream feature tested is bitwise identical to upstream (Hubbard `Mz` T/F,
`Symm`, `Checkerboard`, `Ltau` 0/1, projector, `N_SUN`, continuous fields, `Propose_S0`, global
moves, Langevin, HMC, MALA, Kondo, tV, Plain Vanilla). Two intended differences remain:

| Difference | Where | Effect |
|---|---|---|
| Default `LOBS_ST` is 0 (upstream: 1) | `Prog/QMC_runtime_var_mod.F90`, commit 0b63e5d0 | The down sweep also measures tau = 0: 2 Ltrot instead of 2 Ltrot - 1 equal-time measurements per sweep, for every model at finite temperature with sequential updates. Results differ statistically from upstream unless `LOBS_ST = 1` is set. |
| `ham_h0` is the space-time point defect | `Prog/Hamiltonians/Hamiltonian_Hubbard_smod.F90` | Upstream: static pinning field on the origin, orbital 1, all time slices. Here: a field on one time slice (`g_t(1)`) at site `Latt%N/2 - L/2`, all orbitals, with the observables `SpinZ_local` and `SpinZ_localtau`. The site formula assumes a square lattice. |

New rule (2026-11-13): with `Ltau = 1` the pinning field is switched off (`ham_h0 = 0`); the
`info` file says so.

The delayed and submatrix schemes reproduce the standard update to 1e-12 with identical
trajectories (with and without the defect). `SpinZ_localtau` agrees with an exact U = 0
calculation to 4.7e-14 for all tau.

Known defects of the point-defect observable outside the production setting (finite temperature,
sequential updates, default `LOBS_ST`), found in the same check. `Ham_Set` now **refuses** these
combinations with `ham_h0 /= 0` (the normalisation itself is unchanged):

1. `SpinZ_localtau` is NaN if slice `Ltrot/2` is never measured (for example `LOBS_ST = 15`,
   `Ltrot = 20`): the normalisation counter is only incremented at that slice.
2. With `Projector = .T.` the endpoint weights are wrong (first column never measured, last
   doubled) and the pin sits inside the projection region.
3. With Langevin or HMC updates the endpoint weights are wrong (tau = 0 column zero, tau = beta
   doubled).

A proper fix (count every measurement, weight only end points that are visited half as often)
needs the same change in the host code, `gpu_obs_mod.F90` and the HIP kernel, and is left open.

## Verification of the unified branch (2026-11-13)

Method: the unified build against the pre-merge tip of the same backend and against a CPU build
of the same compiler, same parameters and seeds. Commits tested: 609a6347, 4919ce0d, 05d64293
(these differ only in declarations, comments and one preprocessor line).

| Backend | Machine, compiler | Build | GPU vs same-compiler CPU | Unified vs pre-merge tip | Timing |
|---|---|---|---|---|---|
| `GPU=none` | macOS arm64, gfortran 16 | pass | n/a | bitwise identical to cc473abc (L=4, 6; Mz T/F; delayed, submatrix) | n/a |
| `GPU=none` vs upstream | same | pass | n/a | see "Differences from upstream" above | n/a |
| `GPU=amd` | MI210, gfortran 14 GCN + HIP | pass | 1.15e-14 (L=6), 5.0e-14 (L=12), acceptance identical; same with `ALF_HIP_ASYNC=0`, `ALF_HIP_QR=0` | bitwise identical to 156a0b01 at L=6, 12, 16, 24 | L=24 about 4.5 s per sweep, unchanged |
| `GPU=intel` | Max 1550 (SuperMUC-NG Phase 2), ifx 2025.3 | pass (JIT and AOT) | default switches: 2.2e-15 (L6F), 1.5e-14 (L12F), 5.5e-11 (L6T). V7 switches: 2.3e-15, 1.6e-14, 1.36e-10 (L6T, same value as before the merge; tolerance 1e-10) | Green functions, scalars, acceptance and configurations bitwise identical to the p5 build; `Den_eq`, `SpinT_eq`, `SpinZ_eq` differ by at most 3.3e-16 (the `acc0` reduction from pp-opt) | L=24, one chain, V7: 4.92 to 4.98 s per sweep before, 4.42 to 4.48 after (observables 1.18 to 0.74 s) |
| `GPU=nvidia` | A30, nvfortran 26.9 | pass (GPU and CPU reference) | **not run**: the card was occupied by a production run | **not run** | **not run** |
| `GPU=nvidia` | H100 (Helma), nvhpc 25.7 | pass (GPU and CPU reference of 609a6347, GPU of cc473abc and 4919ce0d) | **not run yet**: jobs 1315946 and 1315963 (preempt) queued behind the production chains (account GPU limit) | **not run yet** | **not run yet** |

Helma results will appear in `/hnvme/workspace/b133ae29-ALF_GPU/unified_check/logs/job_1315946.out` and
`job2_1315963.out` (GPU vs CPU at L=6 and 8, old vs new GPU binary bitwise at L=6 and 12, timing at L=24).

V7 switches: `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2
ALFSH_DEVUNMQR=1 ALFSH_QRNB=16`. Intel jobs 560852 to 560860, about 0.85 node-hours. Run
directories: `amd:~/amd/work/uni1`, SuperMUC `$SCRATCH/large_prompt/unified`,
`nividia:~/alf-gpu/work/uni1` (binaries built, `run_cmp.py` ready to run when the card is free).

The checks were run before the later commits on this branch, which change the Hubbard
Hamiltonian set-up only (the `Ltau = 1` rule and the refusal of unsupported defect settings),
tests and documentation.

Added after these checks and verified separately on 2026-11-13 (commit 819a96ec on `speedup/unified`):
the LU-free determinant phase in `dev_cgr`, enabled for oneMKL only. Max 1550: gates against CPU
unchanged (2.2e-15, 1.5e-14, 5.5e-11), phase equal to the LU phase to 5.8e-15, L=24 single chain
4.45 to 4.04 s per sweep with the V7 switches. MI210 with the same commit: bitwise identical to
05d64293 at L=6, 12, 24 (HIP keeps the LU). Note on the Intel runs: a node has 8 flat devices,
`ZE_AFFINITY_MASK` 0 to 7; larger values crash or run silently on the host.

An independent review of the merge compared, per backend, the code the preprocessor selects with
that backend's own tip; the two defects it found (a preprocessor directive not in column 1, eight
scalar-buffer dummies declared too short) are fixed in 4919ce0d.

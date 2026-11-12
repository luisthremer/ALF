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
| `origin/luis_playground` (b86df3ad) | MaxEnt bootstrap work in `Analysis/` on an older upstream base; not part of the GPU port. |
| `origin/gpu-port-207-4`, `-207-harness`, `-207-wip`, `-212-hardening`, `-302b-docs`, `origin/feature/delayed-submatrix*` | Early ticket branches with a few commits that never reached the port. |

## Tree changes after the merge

- `docs/gpu/README.md` is the entry page. Platform pages are in `docs/gpu/platforms/`
  (`HELMA_H100.md`, `SUPERMUC_NG2.md` and `docs/gpu/amd/` moved there).
- `docs/supermuc_large_prompt/` moved to `docs/gpu/sessions/supermuc_large_prompt/`; path
  references in scripts and documents were updated. Production runs that started from the older
  branches use their own checkouts and are not affected.
- New: `RUNTIME_SWITCHES.md`, `SPEEDUPS.md`, this file.

## Verification

See the table in the section "Verification of the unified branch" below; it is filled from the
runs of 2026-11-13.

VERIFICATION_TABLE_PLACEHOLDER

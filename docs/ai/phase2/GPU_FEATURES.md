# P2.5 — GPU feature table (DRAFT)

Draft of 2026-11-26, written read-only on `Code-modularity-GPU` (tip c9ba049c: P2.4 steps 1a-1c merged). **Nothing was built or run
for this draft.** Spec: ROADMAP.md section P2.5. It drives P2.6 (first: delayed/submatrix for every model).

## How to read it

Cell values:
- **works** — the device (or the stated host-assisted path) runs it.
- **host** — GPU path active, the named part runs on the host (host-assisted).
- **refused** — `gpu_supported` refuses; `info` says `GPU path (gpu_state): INACTIVE, reason: ...` and the whole run uses the CPU path
  (never a silent wrong answer).
- **stops** — the model itself stops the run (all builds).
- **untested** — nothing in the sources settles it.
- **n/a** — no code for it on that line.

Provenance of each column (kept apart on purpose):
- **gpu-port (doc)**: copied from `fork/gpu-port:docs/gpu/{FEATURES,NOT_IMPLEMENTED,VERIFICATION}.md` and `docs/gpu/platforms/*`
  (state 2026-11-13, commit b7bf0bb3 of `gpu-port-unified`). *Not re-checked here.* "works" there means "the code does it"; whether a
  vendor has *run* it is in section 3. A vendor cell reads `works*` where the doc says the code path exists but no run on that vendor is documented.
- **new line (read)**: what the code on `Code-modularity-GPU` does, read from `Prog/gpu_state_mod.F90` `gpu_supported` (lines 575-691) and the
  modules present in `Prog/` (only `gpu_state_mod`, `gpu_wrap_mod`, `gpu_update_mod`, `gpu_lowrank_mod` exist; there is no `gpu_stab_mod`,
  `gpu_obs_mod`, `gpu_tau_mod` yet). No run by this draft. The new line has **no vendor backend yet** (`GPU=host` only; `GPU=nvidia|intel|amd` is
  P2.4 step 2), so there is one new-line column, "GPU=host". When step 2 lands, a vendor column is added.
- A new-line cell ends in `(1d)` / `(1e)` when the feature belongs to a step that is not merged yet; the cell then describes the plan
  ([P2_4_GPU_PLAN.md](P2_4_GPU_PLAN.md)), not code.

Test evidence on the new line that exists (listed, not run here): ctest `50-gpu-state`, `51-gpu-wrap`, `gpu-state-e2e` (`testsuite/GPU.tests/`,
GPU=host vs CPU path of the same binary, Hubbard 4x4), CI job `gpu_host` (`.github/workflows/tests.yml:209`).

## 1. Feature table

Backends of gpu-port: NVIDIA = nvfortran OpenMP offload + cuBLAS/cuSOLVER; Intel = ifx + oneMKL + SYCL shim; AMD = gfortran amdgcn + hipBLAS + HIP kernels.

| Feature | NVIDIA (gpu-port doc) | Intel (doc) | AMD (doc) | New line, GPU=host (read) | Test / run / doc |
|---|---|---|---|---|---|
| **Models** | | | | | |
| Hubbard, `N_FL` 1 or 2, field types 1, 2, rank-1 diagonal vertices | works | works | works | works (1b, 1c); observables (1e) | gpu-port `FEATURES.md` row 1; new: `gpu_supported`, `gpu-state-e2e` |
| Other models (Kondo, tV, Spin_Peierls, Z2_Matter, LR_Jz_Test, LRC, Plain_Vanilla) | host-assisted only if all `Op_V` pass the structure tests, else refused; stock models **not run** | same | same | refused by decision D3 only through the structure tests (same as gpu-port); stock models not run | `NOT_IMPLEMENTED.md` §2 notes ("not verified per model"); P2.6 |
| Non-diagonal `Op_V` | refused | refused | refused | refused (`Op_V(n,nf) is not diagonal`) | `NOT_IMPLEMENTED.md` row 12; new: `gpu_state_mod.F90:631` |
| HS field types 3, 4 (continuous, LRC) | refused | refused | refused | refused (`field type k (only 1 and 2 supported)`) | row 13; `:636` |
| Rank > 1 vertices (`N_non_zero` > 1) | refused | refused | refused | refused | row 15 (P2.6 first item, ticket 407 / "401-A") |
| Field type differs between flavours | refused | refused | refused | refused | row 14 |
| Time-dependent coupling of `Op_V` (`g_t`) | refused | refused | refused | refused | row 16 |
| `Calc_Fl` / weight reconstruction (`ham%weight_reconstruction`, Plain_Vanilla) | refused | refused | refused | refused | row 10; `:615` |
| **Hopping, Trotter** | | | | | |
| Checkerboard hopping | works (dense matrix for `Ndim` <= 800) | works | works | works (1c, group kernels) | `FEATURES.md`; `51-gpu-wrap` self-test |
| Dense hopping | works (GEMM) | works | works | works (1c) | same |
| Symmetric Trotter (`Symm`) | works | works | works | works (wrap, 1c) / observables (1e) | `gpu_wrap_mod` `gpu_symm` |
| Time-dependent hopping `Op_T` with `g_t` (space-time point defect) | works | works | works | works in the mirror (`td_*` tables); defect model itself not on this line (P2.6, "defect model") | `TIMEDEP_HOPPING.md` (gpu-port); new: `gpu_state_mod.F90:150-158` |
| Time-dependent hopping, static element after a time-dependent one | refused | refused | refused | refused | row 20 |
| N_SUN | untested (not refused) | untested | untested | untested (not refused) | `NOT_IMPLEMENTED.md` §2 notes |
| Complex phases / complex hopping | works (complex `Op_T`, ZHEMM path) | works | works | works (complex elements, `hop_Cval`) | `gpu_state_mod.F90` header; real vs complex by model not run for non-Hubbard |
| **Update schemes** | | | | | |
| Sequential local updates, rank 1 | works | works | works | works (1c) | `FEATURES.md`; `gpu-state-e2e` |
| Delayed updates (`Sequential_scheme` delayed, `N_delay`) | works | works (async block opt-in) | works (`hk_blk_end`) | works (1c); `N_delay` fixed at first use | `gpu_lowrank_mod.F90`; `CONSOLIDATION.md` merge 6 |
| Submatrix updates | works (matrix `M` on host) | works | works | works (1c), `M` on host | same |
| Delayed/submatrix for non-Hubbard models, rank > 1, non-diagonal vertices | refused (no device code) | refused | refused | refused | **P2.6 first**; ticket 407 |
| `Sequential_scheme = auto` | n/a (not a gpu-port option) | n/a | n/a | refused | `:607` |
| Timed `N_delay = -1` with delayed/submatrix | untested | untested | untested | refused | `:612` |
| `Lowrank_cross_slice` | refused | refused | refused | n/a (CM has no such option) | gpu-port row 9; CG-ALF T6 |
| Global moves, `Global_tau` (`N_Global_tau` > 0) | refused | refused | refused | refused | row 6; `:596` |
| Langevin, HMC | refused | refused | refused | refused | row 5; `:593` |
| MALA (sequential, global, global-tau) | refused | refused | refused | refused | row 8; `:602` |
| `Propose_S0` | refused | refused | refused | refused | row 7; `:599` |
| **Stabilisation** | | | | | |
| Default (pre-pivoted QR) | works | works (oneMKL QR) | works (own HIP QR m <= 2304) | **1d** (device stab on host fallback), not merged | `FEATURES.md`; P2_4 plan step 1d |
| `LOG` | host-assisted (`ZGEQP3` on host) | same | same | **1d** | `gpu_stab_mod` (gpu-port) |
| STAB1, STAB2, STAB3 | refused (compile-time, whole build) | refused | refused | refused (compile-time) | row 3; `:587` |
| `QRREF` | untested | untested | untested | untested | `NOT_IMPLEMENTED.md` §2 notes |
| Determinant-phase underflow, N >~ 2950 | fixed only on `fix/det-phase-underflow`, not on `gpu-port` | same | AMD HIP line on its own branch | **1d must carry the fix** (class A, KNOWN_ISSUES DP-GPU) | KNOWN_ISSUES.md; plan step 1d |
| **Observables** | | | | | |
| Equal-time, Hubbard | works (device) | works | works (`hk_*`) | host `Obser` until **1e** | `gpu_obs_mod` (gpu-port) |
| Equal-time, other models | host-assisted | host-assisted | host-assisted | host (CPU `Obser`); no device obs by D3 | `main.F90:1033,1109` (gpu-port) |
| Custom user observables | host-assisted | host-assisted | host-assisted | host | same |
| Time-displaced, `Ltau = 1`, Hubbard | works | code only; **not run** on Max 1550 | code only; **not run** | host `tau_m` until **1e** | `gpu_tau_mod`; `FEATURES.md` run evidence |
| Defect local observables `SpinZ_local(tau)` | works | works | works | not on this line (defect = P2.6) | `TIMEDEP_HOPPING.md` §3 |
| Entanglement (Renyi) | CPU code only | CPU only | CPU only | CPU only | `FEATURES.md` row |
| **Projector, MPI, I/O** | | | | | |
| Projector / trial wave function (`Projector = .T.`) | refused | refused | refused | refused | row 4; `:590` |
| MPI independent walkers | works (rank `r` -> device `r mod n`, MPS) | works (one rank per tile) | **untested** (no MPI run documented) | untested on GPU=host (CPU MPI is covered by the normal MPI CI jobs, not by `gpu_host`, which is `noMPI`) | `MPI_GPU.md`; CI `tests.yml:231` |
| `PARALLEL_PARAMS` (PP) | works | works | works (same code) | untested (mirror per rank expected to work; `gpu_state_mod.F90` comment) | `PP_GPU.md` (gpu-port) |
| Tempering with state exchange (`TEMPERING` without PP) | refused (compile-time) | refused | refused | refused (compile-time); GPU ctest skipped for TEMPERING builds | row 2; `testsuite/GPU.tests/CMakeLists.txt` |
| HDF5 output | works, used in production (H100, A30) | works (production) | **untested** | untested (host-side code, not a GPU path; no `gpu_host` HDF5 job known) | `platforms/helma_h100.md` §3, `intel_pvc_supermuc_ng2.md` |
| Restart (`confin`/`confout`) | works (H100, 9 bins) | untested in docs | untested | untested | `helma_h100.md` §3 |
| Rank-to-GPU binding | wrapper `CUDA_VISIBLE_DEVICES` | `ZE_AFFINITY_MASK` | untested | n/a (no device) | `NOT_IMPLEMENTED.md` §5.5 (known limitation, no check in code) |
| `ALF_GPU_DISABLE=1` | refused-by-request | same | same | same (row 1) | `:577` |
| No offload device | refused | refused | refused | allowed **only** in GPU=host (D2) | row 11; P2_4 plan §5 |

## 2. What decides P2.6 (reading of the table)

1. Everything that is not a Hubbard-like rank-1 diagonal vertex falls back to the CPU: first, delayed/submatrix with rank > 1 and non-diagonal vertices
   for Kondo, tV, Spin_Peierls (ROADMAP P2.6, ticket 407).
2. Not tested by `gpu_supported`, so neither refused nor shown to work: `QRREF`, `N_SUN`, `LOBS_ST/LOBS_EN` windows, HDF5, MPI count, `Symm = F` and
   `Mz = T` on AMD. These need a test each before they can be called "works".
3. Vendor columns cannot be filled for the new line until P2.4 step 2. Until then the gpu-port columns are the only vendor evidence, and they are
   Hubbard / defect-model runs only.

## 3. Run evidence in the gpu-port docs (not repeated on the new line)

| Backend | Where | What | Scope |
|---|---|---|---|
| AMD MI210 | `CONSOLIDATION.md` (gpu-port) | GPU vs same-compiler CPU 1.15e-14 (L=6), 5.0e-14 (L=12); L=6..24 | Hubbard defect, `Symm = T`, `N_FL = 2` |
| Intel Max 1550 | same, `platforms/intel_pvc_supermuc_ng2.md` | 2.2e-15 (L6F) .. 5.5e-11 (L6T); L6T with V7 switches 1.36e-10 against tolerance 1e-10 (noise-limited, open) | Hubbard defect, L=6..24; 8..32 rank runs |
| NVIDIA A30 | `BENCHMARKS.md`, `MPI_GPU.md`, `PP_GPU.md` | production and acceptance; unified-commit comparison **not run** | Hubbard, defect, MPI, PP |
| NVIDIA H100 (Helma) | `platforms/helma_h100.md` | L=8 5.2e-13, restart, benchmarks to L=32; nvhpc 24.11 aborts at L=8 (`la_qrpp` NaN), use 25.7 | Hubbard defect, PP, MPS |
| New line, GPU=host | `testsuite/GPU.tests`, `gpu-state-e2e` | GPU path vs CPU path of the same binary, Hubbard 4x4; plan probe: observables to ~1e-11 relative | Hubbard 4x4, no vendor |

Not run on any backend in the sources: any model other than Hubbard on the device; `Ltau = 1` on Max 1550 and MI210; MPI on AMD.

## 4. Open items for the owner / manager

- Whether the vendor columns for the new line should be filled by re-running gpu-port-style checks at step 2, or only by the new acceptance suite.
- Stock non-Hubbard models: which refusal each gets was never run; a cheap `GPU=host` run per stock model (read-only on code) would fill the first
  "untested" cells today. Proposed as the first P2.5 follow-up.
- Status of 1d/1e cells must be refreshed at each merge (ROADMAP P2.5/P2.6).

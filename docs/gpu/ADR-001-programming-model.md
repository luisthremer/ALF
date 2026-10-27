# ADR-001 – Programming model for the ALF GPU port

- **Status:** accepted (2026-10-27, owner at the Phase 1 gate).
- **Inputs:** HARD.md, PLAN.md §1, tickets 101 (A30 kernel benchmark), 103 (ifx/oneMKL), 104 (pre-pivoted QR),
  105 (nvfortran CPU build), 106 (statistical gate), 102 (OpenACC launch overhead), 108 (STABLOG pivoting).

## Context
- One Fortran source tree must run on NVIDIA (nvfortran, A30 here; H100-class later) and on Intel PVC (ifx + oneMKL,
  SuperMUC-NG Phase 2), and `GPU=none` must stay upstream ALF (H10).
- The work splits into two kinds (101):
  1. **Bulk work** – stabilization QR, dense B_T products, fused element-wise wraps, measurements: 5–40× one core at
     N ≥ 400 on the device. Library calls plus a few simple kernels.
  2. **Per-site work** – the sequential single-site update, 2·L·N proposals per sweep: bound by launch latency
     (OpenMP target region 12–15 µs, 7–9 µs with `nowait`; cuBLAS 5 µs) and by the D→H→D round trip of the host
     decision (rank-1 + decision ≈ 28–34 µs for N ≤ 400).
- ifx has no OpenACC. nvfortran supports OpenMP target, OpenACC, CUDA Fortran and `do concurrent` offload.

## Decision
1. **Directives: OpenMP 5.x target offload, one conservative subset for both compilers** — `target enter/exit data`,
   `target update`, `target teams distribute parallel do` on explicit-shape contiguous arrays, `use_device_addr`,
   `nowait`/`depend` for asynchrony, `interop` (ifx) only inside the backend layer. No `declare mapper`, USM, `loop`
   construct, polymorphism, device allocation or type-bound procedures in device code. Verified with nvfortran 26.9
   (101) and ifx 2025.3 JIT + host fallback (103).
2. **Vendor libraries through one backend layer `alf_la`** (ADR-002): cuBLAS/cuSOLVER via `bind(C)` (CUDA), oneMKL
   `!$omp dispatch` for BLAS plus a small C++ SYCL shim for the LAPACK factorizations (ONEMKL), the unchanged
   BLAS/LAPACK calls (HOST, the `GPU=none` path and test oracle). No other module knows the backend.
3. **Configure-time selection** `GPU=none|nvidia|intel` (PLAN §1); `none` adds no flag and no object code.
4. **Stabilization QR on the device: column-norm pre-pivot + unpivoted `geqrf`/`ungqr`/`ormqr`** for the default
   stabilization in decompose, CGR and CGR2_2 (104: precision same order as ZGEQP3). **STABLOG (108):** no
   GPU-friendly pivoting meets the precision criteria unconditionally (static sort and blocked re-ranking fail; a
   two-pass sort comes close but has outliers up to 4.6× in max precision and ≈ 1.4× higher means). So under STABLOG the
   CGR (N×N) and CGR2_2 (2N×2N) factorizations run as **host ZGEQP3** with the matrix copied D→H and the factors back,
   i.e. upstream numerics. Cost per call on one core: 14 ms (N=400), 186 ms (N=1024), plus 2.6 / 16.8 MB each way;
   the transfers go into the PLAN.md transfer table. STABLOG's own decompose already pre-pivots and stays on the device.
   `GPU=none` keeps ZGEQP3 everywhere; the GPU default path's pre-pivot is approved by the owner (Phase 1 gate), and
   no further STABLOG pivot study is planned (revisit only if profiling shows the host QR dominates).
5. **Per-site update, MVP (Phase 2):** ratio and rank-1 update on the device, decision on the host, so that the GPU
   build reproduces the CPU trajectory flip by flip for the same compiler and seeds (validation level 2).
   The update code is organised as `propose → ratio → decide → apply/defer → flush` so that a device-side decision
   and delayed/submatrix updates replace pieces, not the structure. The device-side decision (below) is
   not in the MVP: with OpenMP it gains nothing (102); it becomes useful together with delayed updates (403/404).
6. **OpenACC: not adopted (102).** On the A30 the only model difference is the asynchronous launch (OpenMP `nowait`
   ≈ 10 µs, OpenACC `async` ≈ 2.6 µs). It matters only for a device-side per-site decision, where OpenACC gives
   6.5–8.8 µs per site vs 15–18 µs for the host decision (N ≤ 400); OpenMP gains nothing there. Even then the single-site
   update alone stays ≈ 2× slower than the 16-core node at N = 144, while rank-k updates cost 0.2–0.4 µs per update
   (101). The lever is therefore algorithmic (delayed/submatrix updates, 401–403), which amortizes launches in either
   model; a second directive model for one vendor is not worth it. Revisit in 404 only as a small NVIDIA-only module
   for the per-site kernels if profiling after 403 still shows launch-bound per-site work.

### Device-side decision without changing the Markov chain (design note for 5)
Per site, ALF draws random numbers in a fixed order: the proposal (`nsigma%flip`: none for type 1, `nranf(3)` for
type 2, `ranf_wrap` for type 3), the `T0_proposal` test (`Prog/Wrapgr_mod.F90:133`, always one `ranf_wrap`), and the
accept test in `Upgrade2` (`Prog/upgrade_mod.F90:222`, one `ranf_wrap`, only if the proposal test passed).
`ranf_wrap` is the intrinsic `random_number` (`Libraries/Modules/random_wrap_mod.F90`), whose state can be saved and
restored (`random_seed(get/put)`). So for one time slice the host can save the state, pre-generate the maximal
number of draws, let the device consume them in order with a device-side counter, read back the count once, restore
the state and skip that many. The trajectory is then identical to the CPU code, with one round trip per slice
instead of per site.
Preconditions, checked at setup (otherwise the host decision is used):
- `ham%S0` ≡ 1 and `Propose_S0 = .false.` – `S0` is user code that may depend on neighbouring fields, which change
  during the slice; it cannot run on the device.
- no `ham%weight_reconstruction` (user code in `Upgrade2`);
- field types 1 and 2 with a diagonal Op_V (H4 MVP subset).
The Hubbard models used for validation satisfy these.

## Alternatives considered
| Option | Why not (now) |
|---|---|
| OpenACC (nvfortran) | NVIDIA only; would need a second directive set. 102: 4× cheaper async launch, but no gain on synchronous launches, cuBLAS or the host decision. |
| CUDA Fortran | NVIDIA only; same objection, more invasive. |
| `do concurrent` offload (`-stdpar`) | Relies on unified memory, no explicit residency or transfer control, library interop unclear; ifx support for offloaded `do concurrent` is newer than our pinned versions. |
| SYCL/C++ or Kokkos rewrite of the hot path | Large rewrite of Fortran with class-based Hamiltonians; violates "one source tree" in spirit; the C++ is limited to the oneMKL shim. |
| MAGMA for pivoted QR on NVIDIA | Only NVIDIA; option for 404 if the STABLOG host QRCP turns out to dominate on NVIDIA. |

## Consequences
- Device code stays small and conservative; vendor-specific code is confined to `alf_la` (+ the SYCL shim).
- Per-site latency remains the main performance risk (R2) until the device-side decision and delayed/submatrix
  updates (401–404) land; BENCHMARKS.md will report the crossover honestly.
- Every change is compiled with nvfortran and with ifx (`spir64` JIT, host fallback) on nividia; PVC behaviour is
  verified by the owner's runs (Phase 3).

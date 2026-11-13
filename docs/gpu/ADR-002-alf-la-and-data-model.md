# ADR-002 – `alf_la` backend interface and device data model

- **Status:** accepted (2026-10-27, owner at the Phase 1 gate).
- **Depends on:** ADR-001. Prototype: `gpu_bench/src/bench_la.F90` (101, 103).
- **Note 2026-11-13:** the decision text below is kept as accepted. Where it names things that were never built or have since
  changed, the line carries a bracketed note and the section "State 2026-11-13" at the end lists the corrections.

## Context
- ALF's hot path uses class-based objects (`class(ham_base)`, polymorphic hopping containers, `Type(Operator)` with
  pointer components that may be MPI shared-memory windows, `class(UDV_State)`), and allocates inside inner loops (H4).
  None of this can live on the device in both compilers.
- 101/103 showed that a thin layer of plain-array linear-algebra routines over cuBLAS/cuSOLVER, oneMKL dispatch and a
  small SYCL shim works with device-resident arrays from `use_device_addr`, and that kernels and benchmark sections
  never need to know the backend.

## Decision

### 1. `alf_la`: one module, plain arrays, backend chosen by a macro
- `Libraries/Modules/alf_la_mod.F90` (new). Backends: `ALF_LA_HOST` (default; the exact BLAS/LAPACK calls ALF makes
  today), `ALF_LA_CUDA`, `ALF_LA_ONEMKL` (+ `alf_la_sycl_shim.cpp`). Selected by configure (`GPU=none|nvidia|intel`).
  [2026-11-13: a fourth backend `ALF_LA_HIP` (`GPU=amd`) exists; `ALF_LA_ONEMKL` and `ALF_LA_HIP` both define `ALF_LA_OM`.]
- Arguments are explicit-shape arrays plus leading dimensions, in LAPACK order. With a device backend the caller
  guarantees the arrays are device-resident (mapped by `gpu_state`); `alf_la` obtains device addresses itself
  (`target data use_device_addr`). No array descriptors, no derived types cross the interface.
- Routines (initial set, from 101 and the profile): `la_init/finalize/sync`; `la_gemm`, `la_gemm_cr`/`la_gemm_rc`
  [2026-11-13: `la_gemm_rc` was never built, no ALF call site; only `la_gemm_cr` exists]
  (complex × real via DGEMM on the interleaved view; dense B_T), `la_gemm_batched` (strided, flavours), `la_geru`,
  `la_geqrf`, `la_ungqr`, `la_ormqr` (with `_lwork` queries), `la_getrf`, `la_getri` (CGR determinant/inverse paths;
  **as built (202): `la_getri` is not provided, its only call sites are cold; the module has `la_getrs` and `la_trsm`
  instead** — see `alf_la_mod.F90` public list; documentation note of 2026-10-29, no code change),
  `la_qrpp` = column-norm pre-pivot + `geqrf` returning a ZGEQP3-compatible `IPVT` (104 semantics); `la_geqp3_host` for
  STABLOG's CGR/CGR2_2 (108): copies the matrix to the host, calls ZGEQP3, returns factors and `IPVT` to the device
  (a documented transfer).
- Error handling: every routine returns or checks a device `info`; failures stop with `Terminate_on_error`.
- cuBLAS scalars: alpha/beta passed via 16-byte-aligned `c_double_complex` module variables (101 crash finding).
- Streams/queues: library calls must run on the queue the OpenMP runtime uses for the caller's kernels, so that
  kernels and library calls are ordered without host waits. `la_sync` is the only blocking call.
    [2026-11-13: not true for every backend, see "State 2026-11-13" item 4.]
  - **CUDA (amended 2026-10-27, ticket 202, measured with nvfortran 26.9 on the A30):** a stream created in
    `la_init` (as first decided here) is *not* ordered against nvfortran's OpenMP target regions, and the legacy NULL
    stream only appeared to work. Synchronous target regions run on the per-thread default stream, which nvfortran
    returns as `acc_get_cuda_stream(acc_async_sync)`; `nowait` regions run on a separate non-blocking stream.
    `!$omp interop` is not implemented in nvfortran 26.9. So `la_init` sets the cuBLAS/cuSOLVER handles to the
    `acc_async_sync` stream (runtime query only, no OpenACC directives).
  - Consequences: (1) `nowait`/`depend` regions are **not** ordered against `alf_la` calls. Any such use needs a
    `la_sync`, or a later ticket that moves library calls onto the `nowait` stream. (2) `alf_la` is single-threaded.
    (3) The test `t_order` (RAW/WAR/WAW without host waits) is the canary; re-run it whenever the compiler or driver
    changes. (4) `-DALF_LA_SYNC_ALWAYS` brackets every call with `la_sync`, as a safe fallback.
  - ONEMKL: the OpenMP interop `sycl::queue` (ifx implements `interop`), to be verified in 302 with the same test.

### 2. `gpu_state`: the flat device mirror
`Prog/gpu_state_mod.F90` (new, compiled only with `-DALF_GPU`). Allocated and mapped once at setup
(`target enter data`), freed at the end; no device allocation in loops.

| Array | Shape | Built from | Notes |
|---|---|---|---|
| `G` [as built: `G_d`] | (Ndim, Ndim, N_FL) | `GR` | device-resident for the whole run |
| `cb_idx`, `cb_exp` [no such arrays: built as `hop_*`, see below] | per family: bond site pairs; 2×2 exponentials for ±1, ±½ | `Hop_mod` checkerboard containers | product order identical to the CPU (H5) |
| `T_exp` [no such array: dense hopping is `hop_Rval`/`hop_Cval`/`hop_At`/`hop_Ac`] | (Ndim, Ndim, N_FL, 4): exp(±T), exp(±T/2), real | dense `Op_T` | dense hopping |
| `V_site`, `V_exp` [as built: `V_meta`, `V_type`, `V_P`, `V_exp`] | per Op_V: `P(1:N_non_zero)`; `E_exp` table per field value | `Op_V(n,nf)` (`diag`, type 1/2) | MVP subset; setup refuses others (falls back to CPU) |
| `nsigma_d` | (N_op, Ltrot) | `nsigma%f` | updated on the device (decision on device) or per site (host decision) |
| `U/D/V` work, `udvst_d` [as built: `udv_U/V/D/L` (stack), `udw_*` (left/right states)] | per NSTM, N_FL_eff | `udvst` | stabilization stack, scales/logs for STABLOG |
| `GT0, G0T, GTT` [as built: `GT0_d, G0T_d, GTT_d`, plus `G00_d` in `gpu_tau_mod`] | (Ndim, Ndim, N_FL) | τ pass | only if Ltau = 1 |
| scratch | preallocated, size from setup | – | QR work, tau, ipiv, rank-1 vectors |
| `rnd_d`, counter [never built: the accept decision and the random numbers stay on the host] | (max draws per slice) | host `ranf_wrap` | device-side decision only (ADR-001 §5 note); not in the MVP (102) |

- Host objects stay the source of truth for setup and for everything outside the GPU path; `gpu_state` copies
  from them once. [2026-11-13: for the Green function this changed with ticket 208-1: `G_d` is device-authoritative, see below.] A setup check (`gpu_supported(Op_V, Op_T, ham)`) decides whether the run uses the GPU path; if not,
  the run proceeds on the CPU and says so in `info`.
- Upstream routines get at most a guarded early branch (`#ifdef ALF_GPU ... if (gpu_active) then; call gpu_X; return`)
  into GPU modules (H10).

### 3. Transfers
Only the transfers in PLAN.md §1 "Planned host↔device transfers" are allowed; each has a one-line comment in the code
naming its reason. Changes to the table need an ADR update. With a later device-side decision (403/404) the per-site
round trip becomes one per time slice (counter + accepted-flip phase).

## Alternatives considered
| Option | Why not |
|---|---|
| Map ALF's derived types with `declare mapper` / deep copy | ifx support only since 2025.2 and buggy (R3 §2); pointer components into MPI windows; polymorphism. |
| Unified/shared memory (managed) | Hides transfers, unpredictable page migration, contradicts the "every transfer documented" rule; PVC performance unclear. |
| Call cuBLAS / oneMKL directly from ALF routines | Spreads vendor code across `Prog/`; the HOST oracle and `GPU=none` identity would get harder to guarantee. |
| hipfort/MAGMA-style portability library | Extra dependency on SuperMUC; MAGMA remains an option only for an NVIDIA device QRCP (404). |

## Consequences
- The HOST backend makes every `alf_la` call testable on the CPU against the current code, bit for bit.
- Adding a backend means implementing `alf_la` only (103 did this for ONEMKL in the prototype).
- Custom Hamiltonians outside the supported subset run on the CPU; the GPU path covers the predefined Hubbard-type
  models first (FEATURES.md lists coverage).

## State 2026-11-13 (corrections; the text above is the decision as accepted on 2026-10-27)

Verified against the code of branch `gpu-port`. Nothing in the decision is reversed; the items below are what the code does and
what the text above gets wrong or leaves out.

1. **Array names (section 2).** `cb_idx`, `cb_exp`, `T_exp`, `V_site` and `rnd_d` do not exist. What `Prog/gpu_state_mod.F90` builds:
   hopping `hop_start`, `hop_meta`, `hop_P`, `hop_Rval`, `hop_Cval` (flat, per flavour, per element; dense elements also `hop_At`, `hop_Ac`),
   checkerboard groups `hop_gstart`, `hop_gptr`, `hop_big`, `hop_aux`, per-site operator lists `vr_ptr`, `vr_ent`, time-dependent hopping
   `td_grp`, `td_n`; Op_V `V_meta`, `V_type`, `V_P`, `V_exp`; fields `nsigma_d`, `nsigma_t`; Green function `G_d`; stabilization states
   `udv_U/V/D/L` and `udw_U/V/D/L`; τ pass `GT0_d`, `G0T_d`, `GTT_d` (and `G00_d` in `gpu_tau_mod`); scratch `scr_c`, `scr_vec`, `scr_r`,
   `scr_work`, `scr_ipiv`. There is no random-number array on the device: the accept/reject decision and the random numbers are on the
   host (rank-1 update: the ratio block goes to the host per site; delayed/submatrix: `gpu_lowrank_mod`). The exact layout is the header
   of `gpu_state_mod.F90`.
2. **Source of truth for G.** Since ticket 208-1 the device copy `G_d` is authoritative during the sweep. `gpu_update_mod` (WRAPGRUP/DO)
   copies the host `GR` to `G_d` only when the device copy is not current and never copies back; `gpu_stab_cgr` with `lazy = .true.`
   writes the new G to `G_d` only. The flags `g_dev_ok(nf)`, `g_host_ok(nf)` (`gpu_state_mod`) say which copy is current (at least one is);
   host readers call `gpu_g_host(GR)` (D->H of the stale flavours), host writers `gpu_g_host_wrote()` (`gpu_update_mod.F90`). The host
   objects remain the source of truth for the operators, the fields at setup and the UDV states of the host path.
3. **Backends (section 1).** There are four: `ALF_LA_HOST`, `ALF_LA_CUDA` (`GPU=nvidia`), `ALF_LA_ONEMKL` (`GPU=intel`, module
   `alf_la_om_mod` in `alf_la_mod.F90`, `alf_la_sycl_shim.cpp`) and `ALF_LA_HIP` (`GPU=amd`, `alf_la_hip_mod.F90`, `alf_la_hip_shim.cpp`:
   hipBLAS, hipSOLVER). `ALF_LA_OM` (derived inside `alf_la_mod.F90`) is defined for both ONEMKL and HIP: they share the `om_*` call layer,
   and under `ALF_LA_OM` the macro `LA_GUARD` stops the run if a public routine without device implementation is reached. No routine lacks
   one today (every public `la_*` has an `om_*` twin). Not built: `la_gemm_rc`, `la_getri` (see the bracketed notes above).
4. **Blocking calls (section 1, streams).** `la_sync` is not the only blocking call. The oneMKL shim waits for every call (the async block
   path of `libalfsb.so` is the exception, one wait per block). The hipSOLVER calls block (`hipDeviceSynchronize` and the read-back of the
   solver info inside the call), as does `om_sync`; the hipBLAS calls, the own HIP kernels and the own QR only mark their work as pending
   when the deferred synchronisation is on (`ALF_HIP_ASYNC`, default on only with the patched libgomp GCN plugin), and the plugin
   synchronises before the next OpenMP launch or copy. On CUDA the statement of the ADR holds (same stream as the synchronous regions).
5. **Own kernels not mentioned in the ADR.** Device QR and Q-apply written for the backends where the vendor library was slow:
   `Libraries/Modules/alf_hip_qr.hip.cpp` (HIP, `ALF_HIP_QR=0` switches back to hipSOLVER), `alf_la_sycl_devqr.hpp` and `alf_la_sycl_qr2.hpp`
   (oneAPI, behind `libalfdevqr.so`, opt-in `ALFSH_DEVQR`, `ALFSH_DEVUNMQR`). Small hot kernels in HIP: `alf_hip_kernels.hip.cpp` behind
   `alf_hipk_mod.F90` (copies, block gather/scale, correlations, scalar observables, stabilization helpers).
   `libalfsb.so` (`alf_la_sycl_blk.cpp`, built by `Libraries/Modules/Makefile`, linked by `configure.sh` for `GPU=intel`) holds the SYCL
   kernels of the asynchronous submatrix block path (opt-in `ALFSH_ASYNCBLK`). Details and switches: `RUNTIME_SWITCHES.md`.

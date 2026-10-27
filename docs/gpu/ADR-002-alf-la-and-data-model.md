# ADR-002 – `alf_la` backend interface and device data model

- **Status:** draft (2026-10-27). Items marked *[102]* or *[108]* wait for those tickets.
- **Depends on:** ADR-001. Prototype: `gpu_bench/src/bench_la.F90` (101, 103).

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
- Arguments are explicit-shape arrays plus leading dimensions, in LAPACK order. With a device backend the caller
  guarantees the arrays are device-resident (mapped by `gpu_state`); `alf_la` obtains device addresses itself
  (`target data use_device_addr`). No array descriptors, no derived types cross the interface.
- Routines (initial set, from 101 and the profile): `la_init/finalize/sync`; `la_gemm`, `la_gemm_cr`/`la_gemm_rc`
  (complex × real via DGEMM on the interleaved view; dense B_T), `la_gemm_batched` (strided, flavours), `la_geru`,
  `la_geqrf`, `la_ungqr`, `la_ormqr` (with `_lwork` queries), `la_getrf`, `la_getri` (CGR determinant/inverse paths),
  `la_qrpp` = column-norm pre-pivot + `geqrf` returning a ZGEQP3-compatible `IPVT` (104 semantics); *[108]* the STABLOG
  variant or a documented host-QRCP route.
- Error handling: every routine returns or checks a device `info`; failures stop with `Terminate_on_error`.
- cuBLAS scalars: alpha/beta passed via 16-byte-aligned `c_double_complex` module variables (101 crash finding).
- Streams/queues: one in-order stream (CUDA) / the OpenMP interop `sycl::queue` (ONEMKL) per rank, created in
  `la_init`; kernels and library calls are ordered on it. `la_sync` is the only blocking call.

### 2. `gpu_state`: the flat device mirror
`Prog/gpu_state_mod.F90` (new, compiled only with `-DALF_GPU`). Allocated and mapped once at setup
(`target enter data`), freed at the end; no device allocation in loops.

| Array | Shape | Built from | Notes |
|---|---|---|---|
| `G` | (Ndim, Ndim, N_FL) | `GR` | device-resident for the whole run |
| `cb_idx`, `cb_exp` | per family: bond site pairs; 2×2 exponentials for ±1, ±½ | `Hop_mod` checkerboard containers | product order identical to the CPU (H5) |
| `T_exp` | (Ndim, Ndim, N_FL, 4): exp(±T), exp(±T/2), real | dense `Op_T` | dense hopping |
| `V_site`, `V_exp` | per Op_V: `P(1:N_non_zero)`; `E_exp` table per field value | `Op_V(n,nf)` (`diag`, type 1/2) | MVP subset; setup refuses others (falls back to CPU) |
| `nsigma_d` | (N_op, Ltrot) | `nsigma%f` | updated on the device (decision on device) or per site (host decision) |
| `U/D/V` work, `udvst_d` | per NSTM, N_FL_eff | `udvst` | stabilization stack, scales/logs for STABLOG |
| `GT0, G0T, GTT` | (Ndim, Ndim, N_FL) | τ pass | only if Ltau = 1 |
| scratch | preallocated, size from setup | – | QR work, tau, ipiv, rank-1 vectors |
| `rnd_d`, counter | (max draws per slice) | host `ranf_wrap` | device-side decision only (ADR-001 §5 note); not in the MVP (102) |

- Host objects stay the source of truth for setup and for everything outside the GPU path; `gpu_state` copies
  from them once. A setup check (`gpu_supported(Op_V, Op_T, ham)`) decides whether the run uses the GPU path; if not,
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
| hipfort/MAGMA-style portability library | Extra dependency on SuperMUC; MAGMA remains an option only for an NVIDIA device QRCP *[108]*. |

## Consequences
- The HOST backend makes every `alf_la` call testable on the CPU against the current code, bit for bit.
- Adding a backend means implementing `alf_la` only (103 did this for ONEMKL in the prototype).
- Custom Hamiltonians outside the supported subset run on the CPU; the GPU path covers the predefined Hubbard-type
  models first (FEATURES.md lists coverage).

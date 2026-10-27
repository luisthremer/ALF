# CPU profile (Phase 0, step 4) – nividia, 2026-10-27

Where one DQMC sweep spends its time, measured on the upstream code, mapped onto the hot-kernel
table of R1 §2. Raw stacks and run directories: owner's laptop `ALF_GPU_data/phase0_profile/`
(27 MB) and nividia `~/alf-gpu/runs/prof/`.

## Method

`perf` (`perf_event_paranoid=4`) and VTune (`yama ptrace_scope=1`) are blocked for unprivileged
users, and we don't change system settings. Instead:

- **Stack sampling under gdb** (`scripts/nividia/pmp.py`). `gdb-oneapi` from oneAPI launches ALF as
  its child, which `ptrace_scope=1` allows, and stops it every 50 ms (≈16–19 samples/s). The full
  stack of every sample is recorded, inline frames included, through OpenBLAS down to the leaf.
  This gives *inclusive* time per ALF routine, which gprof cannot give for BLAS time.
- **Binary** (`scripts/nividia/profile.sh build`): the reference GNU flags (`-O3 -ffast-math`,
  MPI, HDF5) plus `-g`. Objects are not instrumented. The binary is linked with `-pg` and the
  static OpenBLAS 0.3.26 (pthread, BLAS+LAPACK, Zen kernels). The default stabilization is used.
- **Runs:** 1 MPI rank, 1 OpenMP/BLAS thread, configs from `testsuite/gpu_reference/configs.py`
  (Hubbard square, U=4, β=4, Δτ=0.1 → L=40, Nwrap=10, Mz, Ltau=1, N=L²; F=2 flavours).
  Sweeps per run: 12×12: 80, 20×20 cb: 10, 20×20 dense: 6, 32×32: 2.
- **Report** (`scripts/pmp_report.py`): each sample goes to the first matching row (the outermost
  R1 routine that owns the work, e.g. a ZGEMM inside CGR2_2 counts as CGR2_2). Some OpenBLAS
  assembly kernels have no unwind info (`dgemm_kernel_ZEN`, `ztrmm_kernel_*`), so their stacks end
  at the leaf. Those samples are spread over rows in proportion to the resolved samples of the same
  BLAS routine. That affects 3 % of cb samples and 22–38 % of dense samples, almost all of them
  the DGEMM inside ZLARCM/ZLACRM, i.e. the dense hopping.
- **Cross-checks:**
  - A plain gprof flat profile (12×12 cb) attributes 18.3 % self time to `zaxpy_kernel_4`, the
    inner loop of OpenBLAS ZGERU. The sampler gives 17.2 % for the ZGERU row.
  - An instrumented build (`profile.sh build cg`) gives exact call counts (below).
- **Caveats:**
  - The 12×12 and 20×20 runs shared the node with the Intel test_vs_ed run (10 of 16 cores busy)
    and with each other. Shares are robust; absolute times are indicative only.
  - Sampling adds ≈5–10 % wall time.
  - The 32×32 runs include setup (initial WRAPUL/CGR, dense `Hop_mod_init`) spread over only 2
    sweeps. That is the 5 % "other" in 32×32 dense.

## Result: share of wall time per R1 kernel row (%)

| Kernel row (R1 §2) | 12² cb | 20² cb | 32² cb | 12² dense | 20² dense | 32² dense |
|---|---:|---:|---:|---:|---:|---:|
| CGR2_2 (τ stabilization, 2N×2N) | 27.7 | 29.4 | 23.6 | 19.1 | 16.5 | 14.4 |
| CGR (equal-time G from UDV) | 12.7 | 15.0 | 11.2 | 7.4 | 8.8 | 7.1 |
| UDV decompose (QRP, ZUNGQR, ZTRMM) | 6.6 | 11.0 | 7.5 | 5.5 | 6.7 | 4.9 |
| **Upgrade2 rank-1 update (ZGERU)** | 17.2 | 23.5 | **33.3** | 9.7 | 15.2 | 20.1 |
| Upgrade2 ratio + vectors (ZGEMV, alloc) | 3.2 | 1.0 | 0.4 | 1.7 | 0.6 | 0.3 |
| Op_Wrapup/Wrapdo (e^V G e^-V) | 2.1 | 0.8 | 1.2 | 1.4 | 0.4 | 0.6 |
| WRAPUR/WRAPUL B products (eq + τ) | 4.5 | 2.0 | 5.8 | 14.7 | 12.0 | 13.4 |
| PROPR/PROPRM1 (τ B products) | 4.6 | 2.0 | 4.7 | 12.4 | 13.6 | 14.3 |
| Hop_mod B_T G B_T⁻¹ in WRAPGR | 3.0 | 1.6 | 2.8 | 14.5 | 16.8 | 14.3 |
| ObserT (time-displaced measurements) | 5.0 | 3.9 | 1.9 | 4.3 | 2.2 | 1.2 |
| Obser (equal-time measurements) | 11.9 | 7.6 | 6.2 | 7.4 | 4.1 | 3.9 |
| Control_Precision* | 0.2 | 1.1 | 0.3 | 0.1 | 0.2 | 0.1 |
| tau_m other | 0.6 | 0.3 | 0.1 | 0.9 | 0.4 | 0.1 |
| other (setup, I/O, bins) | 0.9 | 1.0 | 1.0 | 1.0 | 2.5 | 5.3 |
| samples (unresolved, redistributed) | 909 (27) | 1530 (34) | 7393 (51) | 1223 (351) | 1570 (602) | 11281 (3434) |
| wall s / sweep, 1 rank (incl. setup, sampling) | 0.60 | 8.0 | 193 | 0.80 | 13.7 | 294 |

Aggregates:

| | 12² cb | 20² cb | 32² cb | 12² dense | 20² dense | 32² dense |
|---|---:|---:|---:|---:|---:|---:|
| Whole τ pass (`TAU_M`, Ltau=1) | 41.3 | 40.2 | 34.3 | 40.8 | 36.5 | 35.3 |
| Stabilization (CGR + CGR2_2 + decompose) | 47.0 | 55.4 | 42.3 | 32.0 | 32.0 | 26.4 |
| Local updates (Upgrade2 + Op_Wrap) | 22.5 | 25.3 | 34.9 | 12.8 | 16.2 | 21.0 |
| B_T/B_V products (WRAPUR/UL + PROPR + Hop in WRAPGR) | 12.1 | 5.6 | 13.3 | 41.6 | 42.4 | 42.0 |
| Measurements (Obser + ObserT) | 16.9 | 11.5 | 8.1 | 11.7 | 6.3 | 5.1 |

By leaf type (where the CPU actually executes):

| Leaf | 12² cb | 20² cb | 32² cb | 12² dense | 20² dense | 32² dense |
|---|---:|---:|---:|---:|---:|---:|
| BLAS-3 kernels (gemm/trsm/trmm) | 23.2 | 38.3 | 27.7 | 49.9 | 61.3 | 50.4 |
| BLAS-1/2 (axpy/ger/gemv/scal/nrm2) | 34.9 | 36.8 | 48.6 | 20.7 | 23.5 | 29.5 |
| LAPACK auxiliary (zlacpy, zlaqp2, zlarf…) | 4.1 | 2.7 | 1.5 | 9.3 | 4.5 | 8.2 |
| ALF Fortran loops (checkerboard ZDSLSYMM, observables) | 31.4 | 20.3 | 21.1 | 16.3 | 9.1 | 10.0 |
| libc / allocation | 1.8 | 0.2 | 0.1 | 0.7 | 0.3 | 0.0 |

## Call counts per sweep (instrumented build, 12×12 cb, 80 sweeps, L=40, N=144, F=2, Nwrap=10)

| Routine | measured / sweep | R1 §2 formula | |
|---|---:|---|---|
| `upgrade2` | 11 520 | 2·L·N_op·F = 23 040 | **flavour loop is inside Upgrade2: 2·L·N_op calls**; per-flavour work as R1 says |
| `op_wrapup` + `op_wrapdo` | 46 080 | 4·L·N_op·F | ✓ |
| `cgr` | 16 | 2·(L/Nw)·F | ✓ |
| `cgr2_2` | 8 | (L/Nw)·F | ✓ |
| `decompose_udv_state` | 24 | 2·(L/Nw)·F = 16 | **+8 from the WRAPUR inside TAU_M**, which R1's table omits |
| `wrapgrup` / `wrapgrdo` | 40 / 40 | L | ✓ |
| `obser` / `obsert` | 79 / 41 | 2L / L | ✓ (2L−1, L+1) |
| `propr` + `proprm1` | 160 | 4·L·F | 4·L calls (flavour loop inside) |
| `fields_make` (allocating) | 138 395 | – | small arrays allocated per Op_mmult/Op_Wrap call (R1 §4); libc time ≤ 1.8 % |

## Findings

1. **The single rank-1 update is the largest kernel at large N and its share grows with N**:
   17 → 24 → 33 % (cb) and 10 → 15 → 20 % (dense) for N = 144 → 400 → 1024. In flops it is only
   ≈5 of ≈30 F·L·N³ (R1 §2), but it is memory-bound BLAS-2: 2·L·N ZGERUs per sweep, each streaming
   the full 16·N² bytes of G (16 MB at N=1024, beyond the per-CCX L3). On a GPU it stays bandwidth-
   bound and serialised (one launch per accepted flip). This is the case for delayed or submatrix
   updates: confirmed as the dominant post-port bottleneck.
2. **Stabilization is the largest group** (42–55 % cb, 26–32 % dense), and CGR2_2 (τ) is the
   biggest single stabilization kernel. Much of it is QRP (`ZGEQP3` → `zlaqps`/`zlaqp2`, BLAS-2
   column-norm updates) and ZUNMQR/ZTRSM. QRP has no good batched GPU equivalent; this is a design
   question for the ADR.
3. **Ltau=1 makes a sweep 1.5–1.7× as expensive** (the τ pass is 34–41 % of the total), not
   "roughly doubles" (R1 §6). The share falls with N because the rank-1 updates grow faster.
4. **Checkerboard vs dense:** checkerboard B products are cheap O(N²) Fortran loops (`ZDSLSYMM`,
   2–6 %). Dense B products are DGEMM via ZLARCM/ZLACRM, about 42 % at every N. Dense costs
   1.5–1.7× checkerboard per sweep at N ≥ 400. On a GPU, dense B products become efficient GEMMs,
   while the checkerboard's ~2N bond updates of 2×2 need a custom kernel (or a dense or sparse
   formulation).
5. **Measurements are not negligible at small N** (17 % at N=144, 8 % at N=1024 cb). They are O(N²)
   scalar loops with scatter-adds into lattice bins, so they must run on the device too (or G must
   be copied back 3L times per sweep).
6. Allocation overhead (R1 §3/§4) is ≤ 1.8 % on the CPU. On a GPU, per-call allocation would matter
   much more; avoid it in device code.
7. Every ALF routine loops over flavours outside BLAS (F=2 here). These are independent matrices
   → batch or stream them on the GPU.

## How to reproduce (on nividia)
```sh
cd ~/alf-gpu   # profile.sh, pmp.py copied from docs/gpu/scripts/nividia/
./profile.sh build                                  # bin/ALF-gnu-PROF.out  (./profile.sh build cg: call counts)
./profile.sh run PROF hub20x20_cb 10 1 pmp          # -> runs/prof/PROF_hub20x20_cb/stacks.txt
python docs/gpu/scripts/pmp_report.py 20cb=.../stacks.txt   # tables above
```

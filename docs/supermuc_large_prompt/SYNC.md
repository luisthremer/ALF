# Phase 2: sync time and GPU use (running log; every number marked C = confirmed on system)

## Part 4 measurements (2026-11-11, test partition, ~0.5 node-h so far)

### Per-primitive latency, one process, one tile (C; `phase2/latbench.F90`, `jobs/latbench.sh`)
| operation | us/op |
|---|---|
| empty `target` region (synchronous) | 6.1 |
| empty `target nowait` + one `taskwait` at the end | 7.2-8.3 |
| `target update to`, 16 B / 1 KB | 1.0-1.3 |
| `target update from`, 16 B | 2.8 |
| `target update from`, 512 KB | 350 (226 with `LIBOMPTARGET_LEVEL_ZERO_COMMAND_MODE=async_ordered`) |
| small kernel (1024 elements) | 9-12 |
| dispatch ZGEMM 1024x32x32 (includes the implicit wait) | 17.7 |
| dispatch ZGEMM 1024x1024x32 / 1024^3 | 55 / 466 |
With `ZEX_NUMBER_OF_CCS=0:4` a single process is not faster (large kernels are slower: 70 us vs 27 us for a 16 MB streaming kernel).
The assumed "100-170 us per synchronous operation" is NOT confirmed: one synchronous operation costs 6-18 us. 76 000 operations x ~10 us = ~0.8 s per sweep, not the 55-60 % of a 8-15 s sweep.

### Where a sweep spends its time (C; VTune hotspots with call stacks, 1 rank, L=24 reference parameters, 3 sweeps, baseline settings; `phase2/vtune_topdown_L24_1rank.txt`)
Share of all samples (the run includes 19 s of mirror set-up, so the sweep part is ~72 % of the total):
| call tree | % of total |
|---|---|
| gpu_stab_cgr (Green function recomputation) | 25.0 |
| gpu_stab_wrapul + wrapur (UDV decompositions) | 26.5 |
| of both: `geqrf` (oneMKL hybrid: panel on the host) | 37.8 |
| of both: `ormqr` | ~6.3 |
| of both: `getrf` | 3.7 |
| measurement (gpu_obs; 8 % of it waits for ONE long correlation kernel) | 11.1 |
| wrapgrdo (delayed submatrix update, incl. blk_end/blk_flush) | 9.1 |
So QR + ormqr + getrf ~ 48 % of the total = ~two thirds of the sweep. The block updates, which contain almost all of the small synchronous operations, are ~9 % (13 % of the sweep).
Even removing all sync waits in the block path could therefore gain at most ~10 %, below the +-15 % run-to-run noise.

### CholeskyQR2 prototype on the device (C; `phase2/cholqr_bench.F90`, N = 1024, one tile)
herk 0.51 ms, potrf 5.7 ms (host-assisted), trsm 2.06 ms, gemm 0.46 ms; full CholeskyQR2 (2 passes + R product) 16.7 ms, independent of host threads.
Orthogonality |Q^H Q - 1|max: 2.7e-15 for column-scaled test matrices up to condition 1e8, 1.6e3 at 1e10 (as expected: fails above ~1/sqrt(eps)), so a run-time orthogonality check with fall-back to geqrf is mandatory.
For comparison geqrf alone costs 96 ms (1 host thread) / 31 ms (8 threads), qrpp+ungqr 185 / 67 ms (t_timing, FINDINGS section 7).
Open: the `info` argument of the dispatched potrf produced oneMKL exceptions ("Exception caught") in the benchmark; it has to be handled correctly (device info slot) before use.

## Owner decisions, 2026-11-11 (relayed by the "Supermuc-NG Phase 2 optimization comparison" session after the user's strategy talk)
1. The Part 4 diagnosis is accepted: the QR (oneMKL hybrid geqrf, host panel) is the main lever; the earlier "100-170 us per sync" estimate was wrong.
2. No CholeskyQR2: no change of the algorithm. The same-seed criterion (Mz=F <= 1e-13) stays unchanged. No cholqr sub-branch.
3. No STABLOG (its CGR does a host ZGEQP3 with full N^2 copies, gpu_stab_mod.F90:63/171).
4. oneMKL is closed source, so the host panel cannot be removed there: own device-resident blocked Householder QR, ported from reference LAPACK zgeqrf (BSD): zgeqr2 panel + zlarft + zlarfb, level-3 parts as dispatched oneMKL gemm on the device.
5. Prototype first (standalone in phase2/, no ALF change) at N = 576, 784, 1024 against oneMKL geqrf (1 and 8 host threads); variant A one kernel per column step, variant B one fused kernel per panel (nb = 16/32/64). Cheap tests in parallel: geqrf_batch, host threads per rank.
6. Gate: integrate into ALF (behind a run-time switch, Intel path only) only if the device QR is >= 2x faster than geqrf with 8 threads. getrf (13 ms, host-assisted) is the next candidate afterwards.

## Part 5: device-resident blocked Householder QR, prototype benchmark (C; `phase2/devqr_bench.F90`, `jobs/devqr.sh`, job 560130, ~0.1 node-h)
Port of LAPACK zgeqrf (zgeqr2 panel + zlarft + zlarfb; level-3 parts = dispatched oneMKL zgemm), standalone, OpenMP offload (ifx 2025.3), one tile, random complex matrix.
Variant A = several kernels per column step (host drives the columns), variant B = one fused kernel per panel (one work-group, 1024 threads, OpenMP worksharing loops with barriers).
| N | host MKL zgeqrf 1 thr / 8 thr | variant A nb=16 / 32 / 64 | variant B nb=16 / 32 / 64 |
|---|---|---|---|
| 576 | 15.3 / 4.3 ms | 60.7 / 75.5 / 120.9 ms | 55.3 / 59.7 / 98.1 ms |
| 784 | 36.4 / 8.1 ms | 83.1 / 104.3 / 175.1 ms | 76.2 / 83.4 / 146.4 ms |
| 1024 | 82.6 / 16.7 ms | 109.3 / 136.2 / 222.4 ms | 100.7 / 109.5 / 187.4 ms |
- Variant A is correct: max rel. deviation of A (R and V) from host zgeqrf 1.6e-15 .. 4e-15 (tau 5e-14 .. 9e-14), i.e. same numerics as LAPACK. Variant B has a data-race/logic bug (deviation O(1)); its timing is still valid for the cost estimate.
- Time per column step: A ~107 us, B ~98 us (N=1024, nb=16). The cost is NOT the launch count: B has one launch per panel, but every OpenMP worksharing construct / reduction / barrier inside the fused kernel costs ~13 us in ifx offload (6-7 per column step).
- **Gate (device QR >= 2x faster than geqrf with 8 host threads, 31 ms from t_timing) FAILS by a factor 3-6** for both OpenMP variants. Even host-only MKL zgeqrf (no GPU at all) with 8 threads is 16.7 ms at N=1024, faster than oneMKL's hybrid (31 ms).
- The host-MKL reference with 8 threads in this program disagreed with the 1-thread result (the comparison columns for the 8-thread job are meaningless; only the 1-thread job is used for correctness). The dispatched device `zgeqrf` call in this program returned in 0.05 ms and did nothing useful (not used).

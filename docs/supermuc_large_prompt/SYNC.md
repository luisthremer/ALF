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

# PARKING

Parked work for the GPU delayed/submatrix branch (`gpu-port-delayed-submatrix-speedup`), written 2026-11-06.
Nothing here has been run yet. The simulation list at the end is to be run later exactly as written.

Column `x`: the feature or optimisation behind the item is NOT implemented yet (code has to be written first), or (item 14) there is no
access yet. A blank means the code exists and the item only needs testing, measuring or documentation.

## Task numbering (global priority, 1 = highest)

Group letters: A correctness and verification, B model and feature coverage, C performance at small and medium N,
D performance at large N, E Intel and portability, F project housekeeping. The number in brackets is the item's number in the
previous (per-group) list.

| # | x | Group (prev) | Item |
|---:|:-:|---|---|
| 1 |  | A (1) | Acceptance testing of the new schemes: `test_vs_ed` on the GPU for delayed and submatrix, statistical comparison against CPU references with error bars, precision numbers |
| 2 |  | A (2) | `GPU=none` regression gate (211) re-run after the delayed/submatrix merge: identical behaviour and speed to upstream |
| 3 |  | B (15) | Run every predefined Hamiltonian on the GPU to confirm which are refused and why (today only a code-reading guess) |
| 4 |  | A (3) | Cross-compiler statistical agreement (B580 vs A30 trajectories diverge; never compared statistically) |
| 5 |  | A (4) | Wider coverage of the new code: dense hopping, other lattices, bilayer, field type 1, `Ltau = 0`, `Symm = .F.`, fallback paths |
| 6 |  | A (5) | Understand the nvfortran `-O2`/`-O3` miscompile of `blk_end` (pinned to `-O1`); minimal bug report; check ifx |
| 7 |  | A (6) | Rebase onto current master (MALA changes) |
| 8 |  | B (16) | Update the stale `FEATURES.md` (still the end of Phase 0) |
| 9 | x | C (8) | Several chains per GPU (404): the main throughput lever, GPU mostly idle per chain |
| 10 | x | D (35) | Sparse `Symm` instead of dense A·M·B GEMMs: about -37 s per sweep at N = 4096, plus about 240·N² bytes of memory |
| 11 |  | D (36) | Larger `N_delay` at large N (about 128, scanned): update 17 s -> 6–7 s at N = 4096; then host O(k²) algebra needs BLAS or a device version |
| 12 | x | C (9) | Stabilization latency (about 1 ms per small QR): batch flavours, separate streams, device-side pivot sort |
| 13 |  | E (30) | Validation on Ponte Vecchio (SuperMUC-NG Phase 2); ticket 306 (host vs device per routine on PVC) |
| 14 | x | E (31) | Helma and SuperMUC builds (NVHPC missing on julia2; only built on `nividia` and `intel` so far) |
| 15 | x | B (17) | Rank > 1 and non-diagonal vertices (407): opens tV, Kondo, Spin_Peierls |
| 16 | x | B (18) | Flavour-symmetry reconstruction (`Calc_Fl`): opens Plain_Vanilla, Spin_Peierls |
| 17 | x | C (10) | Block mode for the delayed scheme (about 2× slower than submatrix now) |
| 18 |  | D (37) | Real 56×56 run to test the extrapolated scaling and memory limit (fits stop at N = 1024 plus one 48×48 point) |
| 19 | x | B (19) | Device-side observables for models other than Hubbard; entanglement; custom observables on the device |
| 20 | x | B (20) | MPI together with GPU (405); HDF5 output with GPU untested |
| 21 | x | E (33) | Unified accuracy gates: pre-registered acceptance (207-6/7, 302c) predates the new schemes |
| 22 | x | C (11) | Per-bin output (`Pr_obs`): 0.3 s per bin at 12×12, 1.3 s at 20×20 |
| 23 |  | A (7) | Remove the temporary instrumentation in `main.F90` (`gpu_tacc`, extra `info` lines) |
| 24 | x | F (43) | Validation suite (209), final benchmarks (210) and the owner gate redone for the new code |
| 25 | x | D (38) | Memory footprint (about 1650·N² bytes at β = 5): trim dense `Symm`, saved copies, workspaces |
| 26 |  | D (39) | Parameter-level choices (your decision): `Nwrap` 10 -> 20, τ pass every 4th sweep, or `Ltau = 0`. Together with 10 and 11: about 3× to 7.6× at N = 4096 |
| 27 |  | E (32) | Understand the B580's FP64 weakness (stabilization, τ pass, dense products 3–8× slower) and the 0.45 ms fixed oneMKL QR cost |
| 28 | x | C (12) | Remaining launch overhead: merge flavours into one launch, CUDA graphs, async kernels for the hottest loops |
| 29 | x | B (21) | HS field types 3 and 4 (continuous, complex) |
| 30 | x | B (22) | Global-in-τ and global space-time moves |
| 31 | x | B (23) | `Propose_S0` |
| 32 | x | C (13) | Device-side accept decision (406) with pre-drawn random numbers |
| 33 | x | B (24) | Time-dependent couplings and hoppings |
| 34 | x | B (25) | STAB1/2/3 and QRREF stabilization (default and LOG only today) |
| 35 |  | D (40) | Measure the dense-hop threshold (Ndim <= 800 was a guess) |
| 36 | x | C (14) | Start-up (about 0.4–0.6 s; GPU mirror setup 0.25–0.3 s) |
| 37 |  | F (44) | Documentation for the new parameters and modules beyond `DELAYED_SUBMATRIX.md` |
| 38 |  | E (34) | Check whether ifx needs an equivalent of the nvfortran-only `!pgi$r opt=1` |
| 39 | x | B (26) | Cross-slice delay (`Lowrank_cross_slice`) on the GPU (CPU only today) |
| 40 | x | D (41) | τ-pass `CGR2_2` (2N×2N QR): exploit the block structure, 1.5–2× possible (research) |
| 41 | x | B (27) | Parallel tempering and parallel parameters |
| 42 | x | B (28) | Projector (T = 0), including its τ pass |
| 43 | x | D (42) | Cholesky-QR variants (about 2× cheaper; stability doubtful, experiment only) |
| 44 | x | B (29) | Langevin and HMC (a separate port: forces, solvers) |

## Parked (owner decision 2026-11-06)

Not worked on now. Marked as parked at the owner's request; nothing was removed from the list above.

- **#9 several chains per GPU** and **#20 MPI together with GPU**: parked. They raise throughput, not the speed of the algorithm. Owner: "ich will nicht das MPI machen, sondern wirklich einfach bar algorithm speedup".
- **Algorithm speed-ups discussed 2026-11-06**, all parked for now, in the order recommended:
  1. #10 sparse `Symm` (about -37 s per sweep at N = 4096, about 240·N² bytes less memory)
  2. #11 larger `N_delay` at large N (update 17 s -> 6-7 s at N = 4096; needs BLAS or a device version of the host O(k²) algebra)
  3. #26 parameter choices (`Nwrap` 20, τ pass every 4th sweep, `Ltau = 0`): the owner's decision; together with 1 and 2 about 3x to 7.6x at N = 4096
  4. #40 `CGR2_2` block structure (1.5-2x, research)
  5. #12 stabilization latency (batching, streams, device pivot sort)
  6. #28 launch overhead and #32 device-side accept decision
  7. #17 block mode for the delayed scheme

## Simulation time estimates (full version, to be run later exactly as written)

These are estimates of machine time for the simulations only. They do not include writing or debugging the code, or queue waits on shared systems. They are scaled from the speeds I measured, listed below.

**Measured speeds** (per sweep, β = 5, Mz-Hubbard):

| Lattice | A30 (submatrix) | B580 (submatrix) | CPU, 1 core |
|---|---:|---:|---:|
| 6×6 | about 0.2 s (rough, from earlier runs) | about 0.4 s (rough) | about 0.05 s (rough) |
| 12×12 | 0.36 s | about 0.7 s | about 0.6 s (baseline) |
| 20×20 | 1.0 s | about 5 s | about 10 s |
| 32×32 | 4.6 s | not measured | about 190 s (baseline) |
| 48×48 | 31 s | does not fit | far too slow |

The 6×6 numbers are rough. The B580 numbers are from 12×12 and 20×20 runs. The CPU 32×32 value is from the project baseline, taken on a loaded node.

**Assumption for statistical tests:** I assumed about 50,000 sweeps per case to compare with exact diagonalisation (ED) at the 1e-3 level. That is my guess, not taken from the test suite, so these items could be off by a factor of 2–3.

| # | x | What needs simulations | Machine time |
|---:|:-:|---|---|
| 1 |  | `test_vs_ed` for delayed and submatrix, on A30 and B580 (about 4 cases at 6×6, 50k sweeps each) | A30 about 3 h per case, B580 about 6 h per case. About 12 h on A30 plus 24 h on B580 if run one at a time; about 12–24 h wall if the two machines run in parallel |
| 2 |  | `GPU=none` regression gate: unit tests, `compare_dirs`, CPU speed old vs new | about 1 h |
| 3 |  | Short GPU run of each predefined Hamiltonian, on both GPUs | about 15 min of runs; the work is writing the parameter files |
| 4 |  | Cross-compiler statistics: CPU vs A30 vs B580 at 6×6 or 12×12, about 20k sweeps each | about 4 GPU-h, about 2 h wall |
| 5 |  | About 20 extra configurations (dense hopping, other lattices, bilayer, field type 1, `Ltau = 0`, `Symm = .F.`) | 1–2 h |
| 6 |  | Bisecting the nvfortran miscompile: builds (about 30 s each) plus short runs | 2–3 h |
| 7 |  | Re-validation after rebasing: item 2 plus a few smoke runs | 1–2 h |
| 9 | x | Multi-chain throughput: 1, 2 and 4 chains per GPU at 20×20 and 32×32 (MPS needs setup) | 1–2 h on the A30 |
| 10 | x | Sparse `Symm`: dense-vs-sparse check, then 20×20, 32×32 and 48×48 timing (1–3 sweeps each) | about 30 min |
| 11 |  | `N_delay` scan: 5 values × 20×20, 32×32, 48×48 (2 sweeps each) | 30–40 min |
| 12 | x | Stabilization-latency variants (batching, streams, device pivot): timing runs | about 1 h per variant |
| 13 |  | PVC validation on SuperMUC Phase 2: correctness plus a scaling run | 2–4 node-hours, plus queue time (unknown) |
| 14 | x | First runs on Helma and SuperMUC after the builds work | about 30 min per system |
| 15 | x | Rank > 1 vertices: ED validation on tV and Kondo at 4×4, tens of thousands of sweeps | 2–4 GPU-h per model |
| 16 | x | Flavour reconstruction: validation on Plain_Vanilla (and Spin_Peierls) | 2–3 GPU-h per model |
| 17 | x | Delayed block mode: benchmark plus a small validation | about 2 h |
| 18 |  | The 56×56 run, plus 58×58 and 60×60 to find the real memory limit (about 3 sweeps each) | about 1 h in total |
| 19 | x | Device observables for other models: validation per model | 1–2 GPU-h per model |
| 20 | x | MPI together with GPU: scaling on SuperMUC (8 chains per node) | 1–2 node-h plus queue |
| 21 | x | Unified accuracy gates: calibration plus the pre-registered run | 8–10 GPU-h (the previous tickets budgeted about 5 GPU-h each) |
| 22 | x | `Pr_obs` profiling | about 30 min |
| 23 |  | Smoke re-run after removing instrumentation | about 5 min |
| 24 | x | Validation suite (209) and benchmarks (210): CPU 1-core runs up to 32×32 (about 190 s per sweep, few sweeps), GPU all sizes | 3–5 h, plus 209-1 (counted in item 1) |
| 25 | x | Peak-memory measurements at several N | about 30 min |
| 26 |  | `Nwrap` 20 and τ-frequency studies: precision plus observables bias at 12×12, about 20k sweeps × 3 settings | 3–6 GPU-h |
| 27 |  | B580 FP64 and oneMKL QR micro-benchmarks | about 1 h |
| 28 | x | Launch-overhead experiments (CUDA graphs, merged flavours, async) | 1–2 h per technique |
| 29 | x | Field types 3 and 4: ED validation | 2–4 GPU-h |
| 30 | x | Global moves: validation and performance | 3–5 GPU-h |
| 31 | x | `Propose_S0` validation | 1–2 GPU-h |
| 32 | x | Device accept decision: RNG-order test and benchmark | about 2 h |
| 33 | x | Time-dependent couplings and hoppings validation | 1–2 GPU-h |
| 34 | x | STAB1/2/3 and QRREF validation | 1–2 h |
| 35 |  | Dense-hop threshold: N = 400 to 1600, dense vs sparse | about 1 h |
| 36 | x | Start-up profiling | about 15 min |
| 38 |  | Rebuild and rerun on B580 to see whether ifx has the same bug | about 30 min |
| 39 | x | Cross-slice delay on the GPU: validation | 2–3 GPU-h |
| 40 | x | `CGR2_2` block-structure experiments: micro-benchmarks at N = 1024–4096, plus a precision study | about 0.5 h per variant, a few GPU-h of precision runs |
| 41 | x | Parallel tempering on multi-GPU: correctness at small size, scaling on SuperMUC | about 4–8 h plus 1–2 node-h |
| 42 | x | Projector: T = 0 ED validation | 2–4 GPU-h |
| 43 | x | Cholesky-QR: stability and precision experiments | 3–4 h |
| 44 | x | Langevin/HMC: validation | 4–8 GPU-h |

Items 8 and 37 (documentation) need no simulation.

The biggest simulation costs are the ED validations (items 1, 15, 21, 29, 30, 44), where the A30 and B580 time per case is mostly set by the B580's slower speed. The queue-bound items on SuperMUC, items 13, 20 and 41, may take days of waiting regardless of the machine time.

All "GPU-h" estimates assume one chain on one GPU. MPS with 2–4 chains (item 9) could divide the validation items by about 2–3 once it works.


## Results of the first minimal test round (2026-11-06)

Branches: `gpu-port-delayed-submatrix-speedup` (fixes) and `gpu-port-delayed-submatrix-speedup-master` (the same plus the master merge).
`x` as above (not implemented, or no access). Tests that failed on the first attempt and were fixed: see the Result column of items 1, 2, 3, 5, 7.

| # | Group | x | Short description | Result |
|---:|---|:-:|---|---|
| 1 | A Correctness |  | Acceptance test of new schemes vs ED | PASS. 4-site chain vs ED: A30 10k sweeps -1.479±0.007 (0.96σ), B580 -1.470±0.007 (0.41σ), CPU 20k sweeps -1.477±0.005 (0.87σ); delayed = submatrix = standard bin for bin |
| 2 | A Correctness |  | GPU=none regression gate | PASS. Old vs new CPU: 14/14 observable files bit-identical, 12x12 not slower (3.6 s vs 2.8 s); unit tests 53/53 (before and after the master merge) |
| 3 | B Coverage |  | Every Hamiltonian on the GPU | PASS. Only Hubbard runs on the device; Plain_Vanilla (flavour reconstruction), Kondo/Z2_Matter/Spin_Peierls/tV/LR_Jz_Test (non-diagonal vertices), LRC (field type 3) are refused with a reason and run on the CPU path. tV aborts at Nwrap 10 on every build (numerics), fine at Nwrap 4 |
| 4 | A Correctness |  | CPU vs A30 vs B580 statistics | PASS. 6x6, 200 sweeps: energy -29.64±0.11 / -29.52±0.08 / -29.67±0.14; kinetic and potential parts differ up to 2σ (20 bins only) |
| 5 | A Correctness |  | 24 extra configurations | PASS after 2 fixes. Found: dense hopping + Symm crashed on the device (hop auxiliaries held 2 of 4 matrices), cross-slice key silently ignored (now refused). Nwrap 20: schemes differ 1e-7 on CPU too (precision-limited). Global moves refused (Hubbard has none) |
| 6 | A Correctness |  | nvfortran miscompile bisect | DONE. Reproducible at -O2 and above even without -fast; -O1 correct; vect/unroll/inline/LRE/PRE off do not help; root pass not isolated; ifx unaffected |
| 7 | A Correctness |  | Merge with master | DONE. New branch gpu-port-delayed-submatrix-speedup-master: 5 conflicts resolved (MALA vertices use the standard path; GPU refuses MALA); builds with GNU, nvfortran, ifx; CPU and GPU smoke OK; unit tests 53/53 |
| 8 | B Coverage |  | Update stale FEATURES.md | DONE (2026-11-06). Rewritten from the actual state (what `gpu_supported` refuses, what was validated) |
| 9 | C Perf small/medium N | x | Several chains per GPU | not run (not implemented) |
| 10 | D Perf large N | x | Sparse Symm | not run (not implemented) |
| 11 | D Perf large N |  | Larger N_delay scan | DONE. Best N_delay 32-64; default 32 is within 6% (20x20) and 16% (32x32) of the best update time; 128 is slower (host O(k²)) |
| 12 | C Perf small/medium N | x | Stabilization latency | not run (not implemented) |
| 13 | E Intel/portability |  | PVC validation | PARTIAL. No PVC access; on the B580: all 3 schemes correct at 12x12 and 20x20, bin-identical |
| 14 | E Intel/portability | x | Helma / SuperMUC builds | not run: no access yet |
| 15 | B Coverage | x | Rank > 1 / non-diagonal vertices | not run (not implemented) |
| 16 | B Coverage | x | Flavour reconstruction | not run (not implemented) |
| 17 | C Perf small/medium N | x | Delayed block mode | not run (not implemented) |
| 18 | D Perf large N |  | 56x56 run and memory limit | DONE. 56x56: sweep about 73 s (predicted 70 s); 60x60 peaks at 21,265 MiB, 62x62 at 23,967 MiB (fits, 192 s); 63x63 cannot fit |
| 19 | B Coverage | x | Device observables, other models | not run (not implemented) |
| 20 | B Coverage | x | MPI together with GPU | not run (not implemented) |
| 21 | E Intel/portability | x | Unified accuracy gates | not run (not implemented) |
| 22 | C Perf small/medium N | x | Per-bin output (Pr_obs) | not run (not implemented); measured 38 s per bin at 56x56 |
| 23 | A Correctness |  | Instrumentation behind a flag | DONE. Everything is under #ifdef TIMING; builds with and without it pass on A30 and B580 |
| 24 | F Housekeeping | x | Validation suite and benchmarks | not run (not implemented) |
| 25 | D Perf large N | x | Memory footprint | not run (not implemented); memory model confirmed within 1% |
| 26 | D Perf large N |  | Nwrap / parameter choices | DONE. Nwrap 10: Green error 1e-10; Nwrap 20: 2e-4; Nwrap 30: aborts. Not a free speed-up |
| 27 | E Intel/portability |  | B580 FP64 and QR micro-benchmarks | DONE. Raw FP64 5.1 TFlop/s (A30 5.0), FP32 22 vs 9.9; zgemm 0.18/0.70/0.89 TFlop/s at N 144/400/1024, zgeqrf 1.0/7.4/47 ms. Corrects my earlier "emulated FP64": the matrix work is slow, not raw FP64 |
| 28 | C Perf small/medium N | x | Launch overhead | not run (not implemented) |
| 29 | B Coverage | x | Field types 3 and 4 | not run (not implemented) |
| 30 | B Coverage | x | Global moves | not run (not implemented) |
| 31 | B Coverage | x | Propose_S0 | not run (not implemented) |
| 32 | C Perf small/medium N | x | Device accept decision | not run (not implemented) |
| 33 | B Coverage | x | Time-dependent couplings | not run (not implemented) |
| 34 | B Coverage | x | STAB1/2/3, QRREF | not run (not implemented) |
| 35 | D Perf large N |  | Dense-hop threshold | not run |
| 36 | C Perf small/medium N | x | Start-up | not run (not implemented); measured 7.7 s at 56x56 |
| 37 | F Housekeeping |  | Docs for new parameters | not run (no simulation) |
| 38 | E Intel/portability |  | ifx equivalent of the -O1 pin | ANSWERED. ifx -O3 on the B580 needs no pin in any run so far |
| 39 | B Coverage | x | Cross-slice delay on the GPU | not run (not implemented); GPU path now refuses the key |
| 40 | D Perf large N | x | CGR2_2 block structure | not run (not implemented) |
| 41 | B Coverage | x | Parallel tempering | not run (not implemented) |
| 42 | B Coverage | x | Projector | not run (not implemented) |
| 43 | D Perf large N | x | Cholesky-QR | not run (not implemented) |
| 44 | B Coverage | x | Langevin / HMC | not run (not implemented) |

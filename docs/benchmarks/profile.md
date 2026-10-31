# Profile at L=20: standard vs submatrix

Mz Hubbard, L=20 (N=400), U=4, Beta=2, Dtau=0.1 (Ltrot=20), Ltau=0, LOBS_ST=LOBS_EN=1, N_delay=32, 1 thread.

## Setup

- Tool: macOS `sample <pid> 20 -file ...` (1 ms interval requested, about 12.5k main-thread samples per 20 s window),
  attached after the ~24 s setup, i.e. 29 s after launch, so the whole window lies in the sweep phase
  (NSweep=25, 50 sweeps per run; the run was killed after sampling). The binary is the benchmark binary (commit 26b354c7).
- Compiler: GNU Fortran (Homebrew GCC 16.2.0) 16.2.0
- Flags: `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`
- BLAS/LAPACK: Apple Accelerate (vecLib), linked as `-llapack -lblas`
- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1 (the sample shows a single thread)
- Machine: Apple M1 Pro (8 performance + 2 efficiency cores), macOS 27.2.
- Categorisation: `testsuite/lowrank/profile_samples.py` (self time of each frame is attributed to the innermost frame
  of its stack that belongs to a category; frames inside `libBLAS` count for the calling ALF routine). Categories:
  - flush ZGEMM: `Lowrank_flush`, `Submatrix_factors` (and their ZGEMM);
  - rank-1 updates in standard: `ZGERU` below `Upgrade2`;
  - ratio / border / propose: other `Upgrade2` time, `Upgrade_lowrank`/`Lowrank_*` outside the flush, and the self time of `Wrapgr*`
    (random numbers, S0, flips);
  - hopping/interaction wraps: `Hop_mod_*`, `Op_mmult*`, `Op_wrapup/do` (the `ZSCAL`), `RealExpOpt_*`/`ZDSLSYMM`;
  - stabilisation: `CGR`, `QDRP_decompose` (ZGEQP3), `Decompose_UDV_state`, `WRAPUR`/`WRAPUL` self time, `det_c`;
  - measurements: `Obser`/Hamiltonian observables; other: remainder.
- The raw `sample` files are in `runs/bench/` (not committed); the output of `profile_samples.py` on them (the sample counts behind the
  table below) is committed as `raw_profile_tables.txt`.

## Time fractions in the sweep phase

Seconds per sweep use the N_delay=32 runs of the scan (standard 2.592 s, submatrix 1.112 s, 1 thread).

| category | standard fraction | standard [s/sweep] | submatrix fraction | submatrix [s/sweep] |
|---|---:|---:|---:|---:|
| flush ZGEMM (`Lowrank_flush` / `Submatrix_factors`) | - | - | 19.0 % | 0.211 |
| rank-1 updates in standard (`Upgrade2` -> `ZGERU`) | 63.1 % | 1.636 | - | - |
| ratio / border / propose (incl. `Wrapgr` self time) | 1.4 % | 0.036 | 2.2 % | 0.025 |
| hopping / interaction wraps | 8.5 % | 0.220 | 18.5 % | 0.206 |
| stabilisation (QR/UDV, `CGR`, `WRAPUR`/`WRAPUL`) | 26.7 % | 0.692 | 59.5 % | 0.662 |
| measurements | 0.3 % | 0.008 | 0.7 % | 0.008 |
| other | 0.0 % | 0.000 | 0.1 % | 0.001 |
| samples | 12513 | | 12582 | |

Observations.

- In standard, the rank-1 `ZGERU` update of the full N x N Green function is 63 % of the sweep. Submatrix replaces it by the
  flush ZGEMM (0.21 s) plus negligible ratio/border work (0.025 s): the update cost drops from 1.64 s to about 0.24 s per sweep, a factor 7.
- Stabilisation (pivoted QR `ZGEQP3`, dominated by the level-2 `ZLAQPS`/`ZGEMV`) and the wraps cost the same absolute time in both
  schemes (about 0.69/0.66 s and 0.22/0.21 s per sweep); with the update now cheap they are 78 % of the submatrix sweep. The remaining
  time is therefore in stabilisation and wraps, not in the local update: even a free update would give at most
  2.592 / (0.692 + 0.220 + 0.036 + 0.008) = 2.7x over standard at this size, and submatrix reaches 2.3x.
- Measurements are below 1 % with LOBS_ST=LOBS_EN=1 (one measurement slice); a run with measurements on every slice would shift weight to them.

## Setup (not timed, about 24 s at L=20)

Short look, two 12 s samples of the standard run (3-15 s and 15-27 s after launch): the setup is almost entirely
`Hamiltonian_main::validate_ham_variables` -> `Op_V_is_symmetric` -> `_gfortran_matmul_c8` (gfortran's built-in complex
`matmul`, not BLAS): 100 % of the first window and 78 % of the second (the rest of that window is already the sweep phase,
the setup ends about 24.4 s after launch). The initial Green function computation is negligible in comparison.
This is outside the timed region and unrelated to the local update, but it dominates short runs at large N.

## Multithreaded comparison

THREADS=8 at L=20, N_delay=32 (see `ndelay_scan.md`, section "Multithreaded comparison"; raw `raw_threads.txt`):
standard 2.596, delayed 1.185, submatrix 1.121 s/sweep, versus 2.592 / 1.172 / 1.112 with one thread (N_delay scan, `raw_ndelay_scan.txt`). There is no gain from
8 threads for any scheme (the reason was not investigated; the dominant calls are 400 x 400 or smaller BLAS/LAPACK calls, which vecLib may not
thread), so the speedups at 8 threads are the same as at 1 thread: 2.19 (delayed) and 2.32 (submatrix).

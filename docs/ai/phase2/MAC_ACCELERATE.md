# macOS BLAS/LAPACK backend: Accelerate (AMX) vs OpenBLAS vs reference

Date 2026-11-27, Apple M1 Pro (8 P + 2 E cores), gfortran 16.2 (Homebrew), serial build (`GNU noMPI`).

## Step 0: what a Mac build linked before

`configure.sh` used `-llapack -lblas` for GNU. On this Mac the SDK resolves these to Apple's Accelerate
(`otool -L Prog/ALF.out`: `/System/Library/Frameworks/Accelerate.framework/.../libLAPACK.dylib`, `libBLAS.dylib`),
not to Homebrew. So Mac results were already Accelerate results (cf. the ZGEMM rounding remark in G2).

## Change

`configure.sh`, GNU branch, `uname -s = Darwin` only (Linux and every other machine setting untouched):

| `ALF_MAC_BLAS` | link line |
|---|---|
| `accelerate` (default) | `-framework Accelerate -fopenmp` |
| `openblas` | Homebrew OpenBLAS (`brew --prefix openblas`) |
| `reference` | Homebrew reference LAPACK/BLAS 3.12 |

Usage: `ALF_MAC_BLAS=openblas source configure.sh GNU noMPI NO-INTERACTIVE`. The default is bitwise the same library as before.
`ACCELERATE_NEW_LAPACK` is not used: it only renames symbols for the C headers (`zgemm$NEWLAPACK`);
a Fortran call needs a `bind(C, name=...)` interface for every routine in `linalg_mod`. A hand-written test call
to `zgemm$NEWLAPACK` links and works, so this is possible later, but it changes the interface layer and was not in scope.

## Benchmark (s/sweep, whole run / NSweep incl. start-up, 2 repetitions)

Parameters: Scripts_and_Parameters_files/Start defaults (Dtau 0.1, Beta 5, Nwrap 10, U 4), Hubbard 16x16 (N = 256),
Kondo 12x12 bilayer (N = 288), 6 sweeps, 1 bin, Ltau 0, one run at a time. The Mac was shared with other
sessions (load 4-6), so timings carry about 5 % noise. Script: `mac_blas_bench.py`.

| case | build | threads | s/sweep (rep1, rep2) | Precision Green mean, max | rc |
|---|---|---|---|---|---|
| Hubbard16 | accelerate | 1 | 3.24, 3.21 | 6.0963572621226763E-011, 1.0308763863193723E-006 | 0 |
| Hubbard16 | openblas | 1 | 2.93, 2.88 | 5.1567862433111035E-011, 5.1814205900141408E-007 | 0 |
| Hubbard16 | reference | 1 | 3.67, 3.62 | 7.3957307032449402E-011, 3.1992095873151216E-007 | 0 |
| Hubbard16 | accelerate | 4 | 3.42, 3.51 | 6.0963572621226763E-011, 1.0308763863193723E-006 | 0 |
| Hubbard16 | openblas | 4 | 2.99, 2.91 | 4.8685058074097817E-011, 5.4854870255860444E-007 | 0 |
| Hubbard16 | reference | 4 | 3.60, 3.60 | 7.3957307032449402E-011, 3.1992095873151216E-007 | 0 |
| Kondo12 | accelerate | 1 | 1.75, 1.75 | 4.1913057662566035E-014, 2.5732702135624391E-011 | 0 |
| Kondo12 | openblas | 1 | 1.50, 1.50 | 4.2286010149497845E-014, 1.1442732371777858E-011 | 0 |
| Kondo12 | reference | 1 | 0.08, 0.08 | abort | 3 |
| Kondo12 | accelerate | 4 | 1.83, 1.83 | 4.1913057662566035E-014, 2.5732702135624391E-011 | 0 |
| Kondo12 | openblas | 4 | 1.36, 1.36 | 4.2462360608009003E-014, 1.1494402739233216E-011 | 0 |
| Kondo12 | reference | 4 | 0.08, 0.08 | abort | 3 |

## Reading

- At these sizes (N = 256-288) **Accelerate is not faster**: OpenBLAS is 10-20 % faster per sweep, reference BLAS 10 % slower
  (Hubbard). The matrices are too small for AMX to pay off; the cost is dominated by QR/LAPACK and small GEMMs.
- Threads do not help for either library at this N (4 threads: Accelerate 3 % slower, OpenBLAS 0-9 % faster).
  Accelerate's result is bitwise independent of the thread count here; OpenBLAS changes at the 1e-12 level.
- Large GEMM does profit from AMX: ZGEMM N = 1024 gives 233 GFLOP/s (Accelerate, 1 thread) and 368 (4 threads)
  against 175 / 178 (OpenBLAS). It would matter only for much larger N than DQMC lattices of this size use.
- `Precision Green` is the same order for all three (mean 5e-11 to 7e-11 at the Hubbard 16 start-up default,
  4e-14 for Kondo 12). The Hubbard 16 maximum of 1e-7 to 1e-6 is the Nwrap = 10 default, not the library.
- **Reference LAPACK 3.12 (Homebrew) aborts the Kondo 12x12 run** (G difference 107 > 10, return code 3, 1 and 4 threads),
  where Accelerate and OpenBLAS run clean. Cause not investigated; recorded as KNOWN_ISSUES R1 (class C: nobody runs reference BLAS on the Mac).

Verdict: no leverage row (no measured gain from Accelerate). The default stays Accelerate so that Mac results do not change;
`ALF_MAC_BLAS=openblas` is the faster choice for these sizes if the owner wants a Mac production default.

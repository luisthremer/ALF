# Minimal test round scripts (2026-11-06)

The scripts that produced the results in `PARKING.md` and `docs/gpu/DELAYED_SUBMATRIX.md`. They were written for the two test machines
(`nividia`: A30, GNU/NVHPC; `intel`: Arc B580, ifx) and expect the layout `~/alf-gpu/work/dlsub/{src,srcT,srcC,srcP,srcPG,srcU}`
(build trees) with the scripts and `.sed` files copied into `~/alf-gpu/work/dlsub/`. Paths and host names are hard-coded; adapt them.

- `bld.sh HOST VARIANT`: syncs the working tree and builds (`src` = plain, `srcT` = `-DTIMING`). `build.sh`: the older single-variant version.
- `rc.sh MACH LABEL SEDFILE [HAM]`: runs one case (`MACH` = `nv`, `intel` or `cpu`). Environment: `SCHEME`, `VAR` (build tree), `PSRC`
  (parameter template of an older tree), `TMO` (timeout). A `.sed` file holds the parameter substitutions (quotes must be kept in the file).
- `ana.py FILE [REF] [SKIP]`: mean and standard error of a scalar observable (optionally deviation from a reference value).
  `cmp.py A B`: maximum relative difference between two observable files.
- `prof.sh`, `iprof.sh`, `runall.sh`: older timing runners (A30, B580, three schemes). `fma.F90`: FP64/FP32 throughput (same OpenMP code on
  both GPUs). `gq.cpp`: oneMKL `zgemm`/`zgeqrf` on the B580.

| PARKING item | Scripts |
|---|---|
| 1 | `t1.sed` (2,000 sweeps), `t1mid.sed` (10,000), `t1big.sed` (20,000, CPU), run with `rc.sh` and `ana.py` against ED energy -1.47261997 |
| 2 | `t2.sed`, `t2b.sed` (old `srcP` vs new `srcC` on the CPU) |
| 3 | `t3.sed`, `t3run.sh MACH` |
| 4 | `t4.sed` |
| 5 | `t5run.sh MACH` (22 configurations); global-move cases are `t5g.sed`, `t5g2.sed` (Hubbard has no global moves: expected failure) |
| 6 | `test6.sh`, `t6.sed` |
| 7 | merge check: `mm_a.sed`, `mm_b.sed` (without and with `Sequential_MALA`) |
| 11 | `t11.sh` |
| 13 | `t13a.sed`, `t13b.sed` |
| 18 | `t18.sed`, `t18_56.sed`, `t18_60.sed`, `t18_62.sed` |
| 26 | `t26_10.sed`, `t26_20.sed`, `t26_30.sed` |
| 27 | `fma.F90`, `gq.cpp` |

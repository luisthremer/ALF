# Ticket 104 study scripts (pre-pivot QR vs ZGEQP3)

Used on nividia in `~/alf-gpu/work/t104/` (GNU 13.3, `. ~/alf-gpu/env/gnu.sh`).

- Builds (MPI+HDF5, `make all`, copy `Prog/ALF.out` to `bin/ALF-<name>.out`):
  `def-up`, `def-pp` (default stabilization, without/with `ALF_FLAGS_EXT=-DQR_PREPIVOT`),
  `log-up`, `log-pp` (`configure.sh GNU MPI LOG HDF5`), `def-orig` (unmodified Prog sources).
- `prep.py OUTDIR ALF_DIR name:nsweep:nbin ...` writes parameters+seeds (pyALF only_prep);
  case name `L<L>_b<beta>_U<U>_w<Nwrap>[_dense]`.
- `run_all.sh PREPDIR RUNDIR JOBS case:build ...` copies the same parameters+seeds into
  `RUNDIR/<case>/<build>/`, runs one MPI rank, deletes `confout_*`.
- `analyze.py RUNDIR` prints the markdown table (precision from `info`, pyALF E and S(pi,pi),
  first bin where E/Kin/Pot/S_zz(R) differ by more than 1e-8 relative).
  S(pi,pi) is taken as the maximum of the SpinZ_eq K-space result (half filling: peak at (pi,pi)).

## Ticket 108 (STABLOG pivot strategies; work dir `~/alf-gpu/work/t108/`)
`gen_variants.py` inserted the `QR_PREPIVOT == 2/3` blocks (`QR_PPNB` = block size for 3) into the Prog sources;
`build108.sh up|orig|v1|v2|v3nN`, `ctest108.sh NAME FLAGS`, `stage108.sh` (prep + run), `reseed108.sh` (pyALF writes
the same seeds file for all cases; this makes real extra seeds), `run108.sh`, `post108.py` (E/S to res.json, delete bulky
data), `analyze108.py RUNDIR up v1 v2 ...`, `caselist108.sh`, `mkcases108.py`, `test_variants.F90` (standalone check of the blocks
against ZGEQP3), `bench_variants.F90` (host timings). Report: `docs/gpu/reports/108-stablog-pivot.md`.

## Ticket 108 (STABLOG pivot strategies; work dir `~/alf-gpu/work/t108/`)
`gen_variants.py` inserted the `QR_PREPIVOT == 2/3` blocks (`QR_PPNB` = block size for 3) into the Prog sources;
`build108.sh up|orig|v1|v2|v3nN`, `ctest108.sh NAME FLAGS`, `stage108.sh` (prep + run), `reseed108.sh` (pyALF writes
the same seeds file for all cases; this makes real extra seeds), `run108.sh`, `post108.py` (E/S to res.json, delete bulky
data), `analyze108.py RUNDIR up v1 v2 ...`, `caselist108.sh`, `mkcases108.py`, `test_variants.F90` (standalone check of the blocks
against ZGEQP3), `bench_variants.F90` (host timings). Report: `docs/gpu/reports/108-stablog-pivot.md`.

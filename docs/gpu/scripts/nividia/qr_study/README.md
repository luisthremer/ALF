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

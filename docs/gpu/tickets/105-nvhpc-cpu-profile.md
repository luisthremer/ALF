# 105 – NVHPC (nvfortran) CPU build profile and baseline

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev
- **Depends on:** –
- **Status:** done (2026-10-27, PGI profile works unchanged with NVHPC 26.9; ctest 39/39; test_vs_ed PASS −6.5e-4, 1.52σ, run by manager)

## Goal
Upstream ALF builds and passes its tests with nvfortran 26.9 on the CPU (`GPU=none` equivalent). It
produces single-rank Hubbard references, which later become the **same-compiler, same-seed CPU
reference** for the nvfortran GPU build (HARD.md item 7).

## Context
- `configure.sh` has a `PGI` profile using `pgfortran` (`configure.sh:217-221` flags, ~`:400-415`
  selection, HDF5 via `set_hdf5_flags pgcc pgfortran pgc++`). NVHPC 26.9 still ships `pgfortran`.
  First try the existing `PGI` profile unchanged with `. ~/alf-gpu/env/nvhpc.sh`
  (nvfortran/pgfortran, NVHPC BLAS/LAPACK, HPC-X OpenMPI `mpifort` on PATH).
- **Only if needed**, add a new MACHINE `NVHPC` (nvfortran / nvc / nvc++ names, the same flag set as
  PGI unless a flag is broken) as a new `case` branch. Do not change any other profile's behaviour.
  The GPU switch is ticket 201's job; don't add it here.
- configure.sh self-builds HDF5 1.14.6 into `HDF5/<compiler_version>` when asked for `HDF5`. This
  takes a few minutes the first time.
- Build/test recipe (CLAUDE.md, `docs/gpu/BASELINE.md`): `. ./configure.sh <MACHINE> devel noMPI
  NO-INTERACTIVE`, `make -j4 all`, then `testsuite/`: `mkdir -p tests && cd tests && cmake .. &&
  cmake --build . -j4 && ctest -j4`. Production: `. ./configure.sh <MACHINE> MPI HDF5 NO-INTERACTIVE`.
- References: `testsuite/gpu_reference/` on nividia at `~/alf-gpu/gpu_reference`: `prep.py OUTDIR
  CONFIG --alf-dir SRC [--nsweep --nbin]`, `run.sh BIN NRANKS DIR`, `summarize.py DIR`. pyALF needs
  the venv (sourced by `env/nvhpc.sh`). pyALF's `machine=` argument in `prep.py` is only used for
  prep and is fine as "GNU".

## Files allowed
- `configure.sh` (only a new `NVHPC` case/branch if the PGI profile doesn't work, or a minimal
  fix to the PGI branch that keeps its behaviour for pgfortran users; explain which)
- `docs/gpu/scripts/nividia/build_nvhpc.sh` (new), `docs/gpu/scripts/nividia/README.md` if present

## nividia rules
- Work dir `~/alf-gpu/work/t105/`: clone `fork` there (`git clone -b gpu-port
  https://github.com/luisthremer/ALF.git`, public, https works) or rsync from the laptop, and apply
  your configure.sh change there.
- **Budget: 4 CPU cores** (`make -j4`, ≤ 4 MPI ranks). Other tickets run concurrently. No GPU.
- **Do not run test_vs_ed.** It needs 10 ranks; the manager runs it after you report the binary
  path.
- Runs ≤ 15 min each, no leftovers, install nothing.

## Acceptance
- Unit tests: `ctest` result table for nvfortran devel noMPI. Every failure is investigated and
  reported, not patched by editing tests (a compiler bug or ALF portability bug is a finding, put it
  in the report).
- Production build `MPI HDF5` succeeds. Binary copied to `~/alf-gpu/bin/ALF-nvhpc-DEFAULT.out`, and a
  `noMPI HDF5` build to `~/alf-gpu/bin/ALF-nvhpc-noMPI.out` (the GPU MVP is single-rank and no-MPI;
  check whether `noMPI` + `HDF5` is a valid combination and report).
- Single-rank references with the noMPI binary (if valid, else MPI with 1 rank): `hub6x6_cb` and
  `hub12x12_cb` with `--nsweep 20 --nbin 20` (adjust down to stay ≤ 12 min each), summarized
  (`reference.json`), stored in `~/alf-gpu/runs/ref_nvhpc-1rank/`.
- Report:
  - s/sweep vs GNU (BASELINE.md: GNU 12×12 cb 0.61 s/sweep with 16 ranks; measure a GNU 1-rank
    12×12 cb run with the same sweeps for a fair comparison: `~/alf-gpu/bin/ALF-gnu-DEFAULT.out`);
  - Precision Green/tau/phase;
  - energy per site, compared with GNU within errors. Different RNG streams, so a statistical
    comparison, not bin by bin.
- `GPU=none` safety: `git diff configure.sh` touches only the new or PGI branch.

## Out of scope
GPU flags, test_vs_ed (manager), changes to Prog/.

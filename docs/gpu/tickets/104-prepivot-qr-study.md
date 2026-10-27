# 104 – CPU accuracy study: pre-pivoted unpivoted QR in CGR and CGR2_2

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** –
- **Status:** open

## Goal
Find out on the CPU whether replacing true QR with column pivoting (`ZGEQP3`) by
"sort columns by norm once, then unpivoted QR" in **CGR** and **CGR2_2** keeps ALF's numerical
precision, for both the default and the STABLOG stabilization. The GPU port depends on this answer,
because no GPU library has pivoted QR (HARD.md item 2).

## Context
- `Prog/QDRP_decompose_mod.F90:60-100`: `QDRP_decompose` calls `ZGEQP3` (pivoting where
  `IPVT(j) = 0`, column fixed where `IPVT(j) ≠ 0`), then normalises the rows of R by D = |R_ii|.
- `Prog/udv_state_mod.F90:503-542` (STABLOG branch of `decompose`) already pre-pivots: it sorts
  columns by (log) norm with `ZSWAP`, then calls `QDRP_decompose` with `IPVT = 1` (LAPACK pivoting
  off; ZGEQP3 with all columns fixed is unpivoted QR). Use it as the model.
- **CGR** `Prog/cgr1_mod.F90:303-312`: `IPVT = 0` → full QRCP on `TPUP` for every STAB variant. The
  pivots `IPVT` are used afterwards (inverse permutation of the result). **CGR2_2**
  `Prog/cgr2_2_mod.F90:346, 405, 420`: `IPVT = 0` → full QRCP on the 2N×2N `HLPB1`, IPVT used in
  `solve_extended_System`.
- Read how each routine uses `IPVT` after the decomposition. The pre-pivot must produce an `IPVT` with
  exactly ZGEQP3's output semantics (column j of A·P is column IPVT(j) of A), so the downstream code
  stays untouched. Watch for the phase/determinant bookkeeping (`Pivot_Phase`, `DET_C` of the
  permutation) in each routine. A permutation changes the sign of the determinant.
- Implementation: behind a **test-only macro `-DQR_PREPIVOT`** in those two routines, default off.
  - Without the macro the compiled code must be unchanged (same source lines inside `#else`).
  - Pre-pivot = compute column 2-norms of the matrix as it is passed to `QDRP_decompose`, stable-sort
    descending, permute columns (copy or `ZSWAP` sequence), call `QDRP_decompose` with `IPVT = 1`, then
    set `IPVT` to the permutation.
  - Pass the flag via `ALF_FLAGS_EXT="-DQR_PREPIVOT"` before `source configure.sh`.
- ALF's own measure: the `info` file prints `Precision Green Mean, Max` (wrapped vs recomputed G at
  every stabilization), `Precision tau Mean, Max` and `Precision Phase, Max`. Reference values: GNU
  12×12 β=4: ~3e-11 mean, ~3e-6 max (`docs/gpu/BASELINE.md`).
- Test harness: `testsuite/gpu_reference/` (`prep.py`, `run.sh`, `summarize.py`, `configs.py`) on
  nividia under `~/alf-gpu/gpu_reference`. Make your own small parameter sets with pyALF (venv in
  `~/alf-gpu/env/gnu.sh`), e.g. via a script in your work dir; do not change `configs.py`.

## Study design
- Builds (GNU 13.3 + OpenBLAS, MPI + HDF5, from a copy of `~/alf-gpu/src/ALF-gnu`, which already has
  HDF5 built): {default, STABLOG} × {upstream, `-DQR_PREPIVOT`} = 4 binaries.
- Unit tests: run `testsuite/` (`configure.sh GNU devel noMPI` + `ALF_FLAGS_EXT=-DQR_PREPIVOT`; `make all`;
  ctest). In particular `Prog.tests` cgr/cgr2_2 tests must pass with the macro on. Report all.
- Physics runs: Hubbard square lattice, 1 MPI rank each, `Ltau=1`, `Mz=True`, `Symm=False`,
  `Checkerboard=True` (plus one dense check), `dtau=0.1`:
  - grid: L = 6 and 8; β ∈ {4, 10, 20, 40}; U ∈ {4, 8}; Nwrap ∈ {10, 20} (larger Nwrap is a
    harder stabilization test). Trim it if the runtime budget requires; keep at least β up to 40 at
    U=8 for L=8 with Nwrap 10 and 20.
  - each run ≤ ~10 min, enough sweeps for stable precision statistics (precision is per
    stabilization step, so a few hundred stabilizations suffice).
- For each (parameters, stabilization): upstream vs pre-pivot **with identical seeds**. Record:
  1. Precision Green mean/max, Precision tau mean/max, Precision Phase max, both builds;
  2. the first bin where observables (E, Kin, Pot, S_zz(π,π)) differ by more than 1e-8 relative
     (trajectory divergence), or "identical";
  3. E and S_zz(π,π) with errors for both builds.

## Files allowed
- `Prog/cgr1_mod.F90`, `Prog/cgr2_2_mod.F90` (only code inside `#if defined(QR_PREPIVOT)` blocks; no
  other change)
- `docs/gpu/scripts/nividia/qr_study/**` (your run scripts, a short README)
- `docs/gpu/reports/104-prepivot-qr.md` (results table + conclusion)

## nividia rules
Work dir `~/alf-gpu/work/t104/`. **Budget: 6 CPU cores** (≤ 6 concurrent single-rank runs). Other
tickets run concurrently. No GPU. Runs ≤ 15 min each; delete large outputs (`confout*`, `*.h5` not
needed for the report) afterwards. Install nothing.

## Acceptance
- Upstream build (macro off): unit tests pass, and a 6×6 β=4 run is trajectory-identical to the
  same build from unmodified source (proves the macro-off path is unchanged).
- Macro on: unit tests pass (report any failure with numbers; don't edit tests).
- `docs/gpu/reports/104-prepivot-qr.md`:
  - the full table described above;
  - a verdict per stabilization scheme: pre-pivot precision "same order" (max within 10× of upstream
    and mean within 3×), "worse", or "better";
  - where it degrades, at which β/U/Nwrap.
- The numerics-reviewer reviews the diff (permutation semantics, determinant phase).

## Out of scope
GPU code, changing the default behaviour, other STAB variants (STAB1/2/3, QRREF), projector code.

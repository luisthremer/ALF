# 106 – Statistical comparison redesign

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev (Python work)
- **Depends on:** –
- **Status:** done (2026-10-27; all acceptance comparisons re-run by the manager; verdict rule tightened, see Resolution)

## Goal
Write a CPU-vs-GPU (or compiler-vs-compiler) statistical comparison that passes for independent,
correct runs and fails for a real physics difference. It replaces the provisional
`testsuite/gpu_reference/compare.py` as the observable-level gate.

## Context
- Read `docs/gpu/BASELINE.md`, "Cross-checks and what they taught us". The present `compare.py`
  (full-array χ² over all R/k/τ entries) fails GNU vs ifx on 12×12 correlations, for two reasons:
  1. The slow AF mode gives integrated autocorrelation times of several bins, so ALF's
     naive-binning errors are underestimated (28 bins of 20 sweeps).
  2. All correlation entries move coherently with that mode, so a χ² that assumes independent
     entries is inflated.
  Scalars pass.
- ALF averages each bin over all MPI ranks (one value per bin in the output). Per-bin data is in
  `data.h5` (HDF5 builds); pyALF (`py_alf.ana`, `py_alf.analysis`) in the venv
  (`. ~/alf-gpu/env/gnu.sh`) can read bins. Look at how `summarize.py` and pyALF's `ana` read
  `data.h5`, and use pyALF's per-bin readers if they exist (`py_alf.ana.read_scal`, `read_latt` or
  similar); otherwise read `data.h5` with h5py.
- Existing data on nividia:
  - `~/alf-gpu/runs/ref_gnu-DEFAULT/<config>/` (6×6, 12×12, 20×20, cb and dense, 16 ranks);
  - `~/alf-gpu/runs/ref_gnu-LOG/hub12x12_cb`, `~/alf-gpu/runs/ref_intel-DEFAULT/hub12x12_cb`;
  - binaries `~/alf-gpu/bin/ALF-{gnu,intel}-DEFAULT.out`.
  Configs: `testsuite/gpu_reference/configs.py` (square Hubbard U=4, β=4, Ltau=1, Mz).

## Design (to implement)
1. **Curated observables**, a small and nearly independent set, each a scalar per bin:
   E, Kin, Pot, Part (density); S_zz(k) at (π,π) and (π/2,π/2); Den(k) at (π,π); Green(k) equal-time
   at (0,0), (π/2,π/2), (π,0); τ-resolved: Green(k=(π,0), τ=β/2) and S_zz((π,π), τ=β/2); and
   χ_zz(π,π) = ∫dτ S_zz((π,π), τ) if available from the τ data. Justify any change to this list in
   the README. Map k-points robustly via the lattice information in `data.h5`, not by assuming an
   index order (the earlier attempt got the site order wrong).
2. **Errors:** mean ± error from rebinning. Compute the error at rebin factors 1, 2, 4, … up to
   ≥ 8 effective bins; take the plateau (or the largest factor with ≥ 8 bins) and report τ_int
   estimates. Flag an observable as "under-resolved" if the error has not plateaued.
3. **Test per observable:** z = (a − b)/√(σa² + σb²) with a Student-t correction for the small
   effective number of bins (or an equivalent). Fail if any |z| exceeds the threshold for the chosen
   family-wise false-alarm rate (Bonferroni/Šidák over the curated set at overall α = 0.01, i.e. per
   test ≈ 0.01/m). Also report Σz² vs χ²_m as information.
4. **Full-array diagnostics** (informational, never gating): the old per-array χ²/n, for the record.
5. CLI: `compare2.py RUNDIR_A RUNDIR_B [--alpha 0.01] [--json out.json]`. It prints a table and
   exits 0/1. It reads run directories (needs `data.h5`), not `reference.json`.

## Validation of the comparison itself
- **Short-autocorrelation reference set:** add configs to a new file
  `testsuite/gpu_reference/configs_stat.py` (don't change `configs.py`). Choose parameters where the
  AF mode decorrelates fast, e.g. 8×8, U=2, β=2 and 8×8, U=4, β=2, many bins (≥ 100). Run them with
  GNU and with ifx (different RNG streams), 4 ranks each, ≤ 12 min per run.
- **Must PASS:**
  1. GNU vs ifx on each short-autocorrelation config;
  2. GNU vs GNU with different seeds on one config (re-run with other seeds);
  3. GNU default vs GNU LOG 12×12 cb (identical trajectories → z ≈ 0).
- **Must FAIL (power):**
  1. cb vs dense on a short-autocorrelation config (Trotter difference; BASELINE shows ΔE/N ≈ 0.006
     at β=4, U=4, Δτ=0.1);
  2. U=2 vs U=2.1 (or a similarly small, physical shift).
  If a must-fail case passes, report the resolving power and suggest more statistics.
- **Report** for GNU vs ifx 12×12 cb (the old failing case): what the new tool says, and whether the
  plateau check flags S_zz as under-resolved. That is the expected outcome.

## Files allowed
- `testsuite/gpu_reference/compare2.py` (new), `testsuite/gpu_reference/configs_stat.py` (new),
  `testsuite/gpu_reference/README.md` (new or extend), `testsuite/gpu_reference/compare.py`
  (docstring only: point to compare2.py)
- `prep.py` only if it needs an option to import `configs_stat` (keep the CLI backward compatible)

## nividia rules
Work dir `~/alf-gpu/work/t106/` for new runs; read the existing `~/alf-gpu/runs/ref_*` but don't
modify them. **Budget: 4 CPU cores** (one 4-rank run at a time). No GPU. Runs ≤ 15 min each. Clean
up `confout*`. Install nothing (numpy, scipy, h5py, pandas and pyALF are already in the venv).

## Acceptance
- All must-pass cases pass and all must-fail cases fail, shown as a table in the report with z
  values and effective bins.
- `compare2.py --help` documents the method. README explains the curated set, the error method and
  the thresholds.
- Runtime and data volume of the new reference runs reported.

## Out of scope
GPU runs, changes to ALF Fortran, changing pyALF.

## Resolution (manager, 2026-10-27)
- Re-ran all comparisons with the laptop `compare2.py`: must-pass P1a/P1b/P2/P3 PASS, must-fail F1 (E t=−27.8) and
  F2 FAIL.
- **Verdict rule tightened** (stricter than delivered): PASS now requires that no observable is UNDER-RESOLVED. Under
  the delivered rule the old 12×12 GNU vs ifx case returned PASS with all 11 observables under-resolved. It now
  returns INCONCLUSIVE (exit 2). Same-seed runs with identical trajectories (means equal to 1e-9 relative;
  measured ≤ 1.2e-12) PASS without an error estimate.
- Consequence: the Phase 0 12×12/20×20 references (≤ 30 bins) cannot certify a PASS. The observable-level gate for
  the GPU uses the `configs_stat.py` set (250 bins); the old references remain for timing and same-seed checks.

# Known issues

Open and recently fixed problems of this branch line (`Code-modularity` and the branches cut from
it). Every entry: where it applies, symptom, cause, status, and the branch or commit of the fix.
New findings of the bug hunt (docs/ai/phase2/ROADMAP.md, P2.8) are added here first.

## Determinant phase underflow for N >~ 3000 (L >~ 54) in CGR / DET_C

- **Applies to:** this branch and upstream master (`4c690b04`): the CPU code. Not fixed here.
- **Status:** found on the GPU port (entry below), expected on the CPU path, **not yet confirmed
  in a CPU run** (runs at L = 56 / 64 are in progress, 2026-11-25).
- **Cause:** `DET_C` (`Libraries/Modules/mymats_mod.F90`) multiplies the raw LU diagonal,
  `det_C = det_C*mat(i,i)`. For a unitary matrix |det| = 1, but the running product underflows
  double precision for large N. `CGR` (`Prog/cgr1_mod.F90`) forms the phase as
  `PHASE = CONJG(DET_C(RHS, N_size)); PHASE = PHASE/ABS(PHASE)` with RHS = U_R^dagger U_L (default,
  STAB3 and LOG builds; STAB1/STAB2 call `DET_C` in the same way), so the phase can become
  0/0 = NaN. The Metropolis weight is then NaN and every proposal is rejected, with return code 0,
  a finite G and NaN observables. `Control_PrecisionP` does not notice (comparisons with NaN are
  false). The same function is used by `Op_set` (`Prog/Operator_mod.F90`, normalisation det U = 1),
  the likely cause of the documented `Op_set` failure at Ndim = 4096.
- **Evidence:** on the GPU line, measured (entry below). Numerically, the LU diagonal of a random
  unitary matrix (numpy/scipy, `lu_factor`) has a running product with minimum ln|P_k| = -676
  (N = 2000), -1024 (N = 3000), -1367 (N = 4000); the limit is -708, so the raw product is exactly
  0 for N = 3000 and 4000.
- **Who is affected:** runs with N >~ 2950 (square lattice L >~ 54). L <= 48 is not affected.
  Do not use this branch for N >~ 2950 before the fix.
- **Fix (not done, owner's decision):** a phase-only product that normalises every factor
  (as in `fix/det-phase-underflow` on the GPU line) plus a NaN/unit-modulus check of the phase.
  `DET_C` itself must stay as it is where the magnitude of the determinant matters. Whether and
  how this goes to the upstream developers is the owner's decision.
- **See:** branch `fix/det-phase-underflow` (GPU line), `docs/gpu/KNOWN_ISSUE_DET_PHASE_UNDERFLOW.md`.

## GPU: NaN determinant phase for N >~ 3000 (L >~ 54)

- **Applies to:** every GPU line derived from `gpu-port` (CUDA, Intel/oneMKL, AMD), including
  `gpu-port-unified`, `gpu-delayed-submatrix-defect-production-*` and the `feat/*` branches cut from
  them. Not to `Code-modularity` itself (no GPU code yet); it matters as soon as GPU code is merged
  (ROADMAP P2.4).
- **Symptom:** return code 0 and no warning, but almost all proposals are rejected (acceptance
  0.27 at L = 56, 0.01 at L = 64, physical value about 0.33) and the observables in `data.h5` are
  NaN. Precision Green stays near 1e-15.
- **Cause:** `k_det_phase` (`Prog/gpu_stab_mod.F90` and its HIP twin) formed the phase of the LU
  factor from the raw product of the N diagonal elements; the running product underflows double
  precision for N >~ 2950, so the phase is 0/0 = NaN and the Metropolis weight is NaN.
- **Status:** fixed on branch `fix/det-phase-underflow`, commits `0b3ffba4` (fix) and `c4baba2f` (docs) (2026-11-25): every
  factor normalised, plus `Control_Check_Phase` after every device CGR and a NaN test in
  `Control_PrecisionP`. Verified on one H100 (L = 48 to 96); HIP, Intel and Ltau = 1 not rerun.
  Details: `docs/gpu/KNOWN_ISSUE_DET_PHASE_UNDERFLOW.md` on that branch; check script
  `docs/gpu/scripts/helma/check_large_N_phase.sh`.
- **To do:** merge the fix into every GPU branch before it is ported (P2.4); a GPU branch that has
  not merged it must not be used for N >~ 2950.

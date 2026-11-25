# Known issue (fixed on `fix/det-phase-underflow`): NaN determinant phase for N >~ 3000

Found 2026-11-25. Applies to every GPU line derived from `gpu-port` (CUDA, Intel/oneMKL, AMD) with a square lattice of
N = L^2 >~ 2950 sites (L >~ 54), including the production branches `gpu-delayed-submatrix-defect-production-*`, `gpu-port-unified`
and the `feat/*` branches cut from them.

## Symptom

The run finishes with return code 0 and prints no warning, but most or all Monte Carlo proposals are rejected and the observables in
`data.h5` are NaN. Measured acceptance after 2 sweeps from a random start (physical value about 0.33):

| L | N | acceptance |
|---|---|---|
| 48 | 2304 | 0.329 |
| 54 | 2916 | 0.329 |
| 56 to 60 | 3136 to 3600 | 0.27 |
| 64, 96 | 4096, 9216 | 0.01 |

## Cause

`k_det_phase` (`Prog/gpu_stab_mod.F90`, and its HIP twin in `Libraries/Modules/alf_hip_kernels.hip.cpp`) formed the phase of the
LU factor from the raw product of its N diagonal elements, `det = prod_i A(i,i)`, followed by `det/abs(det)`. The product is only
reached when `dev_cgr` has no tracked determinants (first CGR after a host upload). The final |det| is 1, but the running product
dips to about exp(-0.24 N), which underflows double precision (exp(-708)) at N >~ 2950. |det| = 0 gives phase = 0/0 = NaN, and
the Metropolis weight W = S0 T0 |Re(Phase R)/Re(Phase)| is then NaN, so `Weight > ranf` is false for every proposal until the next
stabilisation. The Green function stays finite and the Precision Green check stays near 1e-15, so nothing else notices.

## Who is affected

Only runs with N >~ 2950. Runs with L <= 48 (N <= 2304) are not affected (running minimum of ln|det| about -551, limit -708). An
affected run is easy to recognise: acceptance in `info` far below 0.3, or NaN in `Ener_scal/obser` of `data.h5`.

## Fix

One line (every factor normalised), plus a guard so the failure can never be silent again:

- `Prog/gpu_stab_mod.F90`, `k_det_phase`: `det = det*(A(i,i)/abs(A(i,i)))`
- `Libraries/Modules/alf_hip_kernels.hip.cpp`, `k_det_phase`: same (not built or run, no AMD machine was available)
- `Prog/control_mod.F90`: `Control_Check_Phase` (called after every device CGR) and a NaN test in `Control_PrecisionP` abort the run
  if the phase is NaN or not of unit modulus.

Verified on one H100 (nvhpc 25.7, 2 sweeps, 1 walker, U = 5, Nt = 64): acceptance 0.330 (L = 56), 0.329 (L = 60), 0.329 (L = 64),
0.328 (L = 96). Ltau = 1, the HIP and the Intel paths were not rerun. Script: `docs/gpu/scripts/helma/check_large_N_phase.sh`.

## What a branch that has not merged the fix should do

1. Merge or cherry-pick `fix/det-phase-underflow` as soon as convenient.
2. Until then add this to the branch's `KNOWN_ISSUES.md` (create the file if it does not exist):

```
## GPU: NaN determinant phase for N >~ 3000 (L >~ 54)
Not fixed on this branch. Runs with N >~ 2950 reject almost all proposals (acceptance 0.27 at L = 56, 0.01 at L = 64) and write
NaN observables, with return code 0. L <= 48 is not affected. Cause: k_det_phase multiplies the raw LU diagonal (underflow).
Fix: branch fix/det-phase-underflow (gpu-port), see docs/gpu/KNOWN_ISSUE_DET_PHASE_UNDERFLOW.md.
Do not use this branch for N >~ 2950 before merging the fix.
```

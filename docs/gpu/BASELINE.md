# CPU baseline (Phase 0, step 3) – nividia, 2026-10-26

Upstream source at commit `515c8de0` (branch `gpu-port`, only docs on top of upstream `master` f81b5cbc).
Machine and toolchains: [MACHINES.md](MACHINES.md). Full results (reference.json, info, parameters,
res.pkl) are on nividia under `~/alf-gpu/runs/ref_*` and on the owner's laptop under
`ALF_GPU_data/phase0_cpu_reference/` (≈36 MB, not committed).

## Builds

| Build | Command | Notes |
|---|---|---|
| GNU unit-test | `configure.sh GNU devel noMPI` + `make all` | gfortran 13.3, OpenBLAS (pthread) |
| ifx unit-test | `configure.sh IntelLLVM devel noMPI` + `make all` | ifx 2025.3.3, MKL 2026.2 |
| GNU prod | `configure.sh GNU MPI [LOG] HDF5` | Intel MPI 2021.19 `mpif90` wrapping gfortran (`~/alf-gpu/bin/mpifort`), HDF5 1.14.6 self-built |
| ifx prod | `configure.sh IntelLLVM MPI [LOG] HDF5` | `mpiifort -fc=ifx`, `-xHost` on AMD Zen 3 works |

Build time (`make -j8 all`, Prog+Libraries+Analysis): GNU 27 s, ifx 37 s.

## Unit tests (`testsuite/`, ctest)

| Compiler | Result |
|---|---|
| GNU 13.3 devel noMPI | **39/39 pass** (4.5 s) |
| ifx 2025.3 devel noMPI | **39/39 pass** (4.5 s) |

Needs `make all` (not `make program`): `Analysis.tests/1-maxent-kernel-overflow` uses `Analysis/` objects.

## test_vs_ed (upstream `prep_test.py` + `analysis.py`, 2 MPI ranks per Δτ)

| Test | Build | Extrapolated E | ED | Δ | Δ/σ | Limits | Result |
|---|---|---|---|---|---|---|---|
| hubbard_finT_chain | GNU MPI HDF5 | −1.47289 ± 0.00042 | −1.47261997 | −2.7e-4 | 0.64 | 1e-3, 3σ | **PASS** (679 s) |
| hubbard_finT_chain | ifx MPI HDF5 | −1.47326 ± 0.00044 | −1.47261997 | −6.4e-4 | 1.48 | 1e-3, 3σ | **PASS** (1167 s, rerun 2026-10-27 with owner's OK) |
| hubbard_finT_chain | nvfortran 26.9 (PGI profile) MPI HDF5 | −1.47327 ± 0.00043 | −1.47261997 | −6.5e-4 | 1.52 | 1e-3, 3σ | **PASS** (1315 s; HPC-X Open MPI needs `PRTE_MCA_hwloc_default_binding_policy=none`) |
| hubbard_T0_chain | – | | | | | | not run (projector out of scope) |

The ifx run needs ~20 min on this machine (ifx is ~1.3× slower than GNU on this AMD host); the first
attempt was killed by a 1000 s timeout, the rerun used 1500 s.

## Hubbard reference runs (GNU, default stabilization, 16 MPI ranks × 1 thread)

Square lattice, U=4, t=1, μ=0 (half filling), β=4, Δτ=0.1 (L_trot=40), Nwrap=10, Mz, Symm=F, Ltau=1,
n_skip=2 bins. Configs: `testsuite/gpu_reference/configs.py`. Per site: E, Kin, Pot; S_zz(π,π) = max_k SpinZ_eqK.
Precision columns: mean / max as printed in `info`. Precision Phase = 0 everywhere (real weights).
s/sweep = wall / (NBin·NSweep) per rank, i.e. one sweep of one rank while 16 ranks share the node.

| Config | Build | E | E/N | Kin/N | Pot/N | n | S_zz(π,π) | Prec. G | Prec. τ | Acc. | bins×sweeps | wall s | s/sweep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hub6x6_cb | gnu-DEFAULT | −29.688 ± 0.011 | −0.8247 | −1.3355 ± 0.0002 | 0.5109 ± 0.0003 | 1.0000 | 3.29 ± 0.03 | 2.4e-11 / 3.1e-06 | 5.1e-12 / 1.3e-06 | 0.430 | 100×50 | 145 | 0.029 |
| hub6x6_dense | gnu-DEFAULT | −29.463 ± 0.008 | −0.8184 | −1.3297 ± 0.0002 | 0.5113 ± 0.0002 | 1.0000 | 3.25 ± 0.04 | 2.6e-11 / 3.9e-06 | 5.8e-12 / 8.7e-07 | 0.430 | 100×50 | 168 | 0.034 |
| hub12x12_cb | gnu-DEFAULT | −118.843 ± 0.037 | −0.8253 | −1.3545 ± 0.0003 | 0.5292 ± 0.0002 | 1.0000 | 3.86 ± 0.15 | 3.4e-11 / 3.0e-06 | 4.5e-12 / 2.3e-06 | 0.432 | 30×20 | 365 | 0.608 |
| hub12x12_dense | gnu-DEFAULT | −117.943 ± 0.036 | −0.8190 | −1.3481 ± 0.0004 | 0.5290 ± 0.0003 | 1.0000 | 3.71 ± 0.11 | 3.5e-11 / 5.3e-06 | 5.4e-12 / 7.2e-06 | 0.432 | 30×20 | 516 | 0.860 |
| hub20x20_cb | gnu-DEFAULT | −329.810 ± 0.105 | −0.8245 | −1.3544 ± 0.0005 | 0.5299 ± 0.0004 | 1.0000 | 3.03 ± 0.15 | 2.2e-11 / 1.2e-06 | 3.2e-12 / 8.3e-07 | 0.432 | 10×6 | 546 | 9.10 |
| hub20x20_dense | gnu-DEFAULT | −327.183 ± 0.193 | −0.8180 | −1.3480 ± 0.0008 | 0.5300 ± 0.0005 | 1.0000 | 3.94 ± 0.24 | 2.2e-11 / 2.0e-06 | 3.7e-12 / 1.8e-06 | 0.432 | 10×4 | 593 | 14.83 |
| hub12x12_cb | gnu-LOG | −118.843 ± 0.037 | −0.8253 | −1.3545 ± 0.0003 | 0.5292 ± 0.0002 | 1.0000 | 3.86 ± 0.15 | 5.1e-11 / 5.6e-06 | 7.1e-12 / 1.6e-06 | 0.432 | 30×20 | 370 | 0.616 |
| hub12x12_cb | intel-DEFAULT | −118.779 ± 0.023 | −0.8249 | −1.3547 ± 0.0003 | 0.5299 ± 0.0003 | 1.0000 | 3.42 ± 0.11 | 4.3e-11 / 2.9e-06 | 6.8e-12 / 5.9e-06 | 0.432 | 30×20 | 470 | 0.784 |

Checkerboard and dense differ by ≈0.006–0.007 in E/N: the extra Trotter error of the checkerboard
split at Δτ=0.1. **GPU results must always be compared with the same hopping variant.**

## Cross-checks and what they taught us

1. **GNU `STABLOG` vs GNU default (same seeds):** all observables agree to ~1e-13 relative, i.e. the
   two stabilization schemes produce the *same Monte Carlo trajectory* over 600 sweeps × 16 ranks.
   Round-off differences did not flip a single accept/reject decision that mattered. `compare.py`: PASS.
2. **ifx vs GNU (different RNG streams, statistically independent):** scalars agree
   (E: 1.5σ, Kin: 0.5σ, Pot: 1.9σ). Correlation functions fail `compare.py` (χ²/n up to 4.5 for
   real-space spin correlations). Diagnosis: the per-bin series of S_zz(π,π) drifts slowly over many
   bins in *both* runs (GNU bins 12–18 ≈ 4.5, bins 19–22 ≈ 2.7; ifx drifts from ≈4 to ≈2.7), so with
   28 bins of 20 sweeps the integrated autocorrelation time of the AF mode is several bins and the
   errors are underestimated. In addition, all R/k/τ entries move coherently with that mode, so a
   χ² that treats them as independent is inflated even with correct errors. **This is a flaw of the
   reference design and of `compare.py`, not evidence of a compiler bug.** Rebinning cannot fix it
   with only 30 bins (ALF averages each bin over all MPI ranks).
3. ifx binaries are ~1.3× slower than GNU on this AMD host (12×12: 0.78 vs 0.61 s/sweep). Not a
   statement about Intel CPUs.

### Consequences for validation (feeds PLAN.md / HARD.md)
- **Primary GPU-vs-CPU observable check: same compiler, same seeds.** A CPU build and a GPU build
  from the same compiler (nvfortran for NVIDIA, ifx for Intel) share the RNG stream, so their
  trajectories should coincide until round-off flips an accept/reject; item 1 shows this can hold for
  hundreds of sweeps. That is a far sharper test than error-bar agreement. Needs an nvfortran CPU
  baseline (not built yet).
- **Statistical check (different streams):** compare a curated, small set of nearly independent
  quantities (E, Kin, Pot, n, S_zz and S_den at a few k, G(k, β/2), χ_zz(π,π)) with autocorrelation-aware
  errors, and use parameters with short autocorrelation (e.g. β=2 or U=2) or long runs. Full-array
  χ² is reported for information only. `compare.py` is marked provisional until this is done.
  **Done in ticket 106:** `testsuite/gpu_reference/compare2.py` (curated 13 observables, rebinning plateau,
  Welch t, Šidák α=0.01) and the 250-bin `configs_stat.py` set (8×8, β=2). It passes GNU vs ifx and different
  seeds, and fails cb vs dense and U=2 vs 2.1. The 12×12/20×20 references above (≤ 30 bins) are **not certifiable**:
  the tool returns INCONCLUSIVE for them. Use them for timing and same-seed checks only.
- 20×20 references have only 8 usable bins: their error bars are indicative only.

## Known issues found
- pyALF 2.6 (PyPI) cannot parse ALF's plain-text complex output written by gfortran 13
  (`(a, b)` with a space after the comma; `ana.read_scal` turns it into an unparsable token).
  pyALF builds with HDF5 by default, which is why CI does not see it. We use HDF5 throughout.
- Unprivileged `perf` (`perf_event_paranoid=4`) and VTune (`ptrace_scope=1`) are blocked; profiling uses
  gdb stack sampling and gprof, see [PROFILE.md](PROFILE.md).

## How to reproduce (on nividia)
```sh
cd ~/alf-gpu
./build_ctest.sh gnu ; ./build_ctest.sh intel          # unit tests -> logs/ctest_*.log
./build_prod.sh gnu  ; ./build_prod.sh intel           # bin/ALF-{gnu,intel}-{DEFAULT,LOG}.out
./run_ed.sh gnu-DEFAULT hubbard_finT_chain
./ref_batch.sh                                         # all of the above reference runs (~75 min)
```

# Delayed and submatrix sequential updates in ALF: final report

Branch `feature/delayed-submatrix`, baseline `f81b5cbc` (master, 2026-08-20), report written on top of `93826bc0`
(2026-10-31). Sources: `CLAUDE.md`, the design spec `docs/superpowers/specs/2026-10-30-delayed-submatrix-design.md`,
`docs/algorithm-notes.md` (cited as "notes §N"), `docs/benchmarks/*`, `testsuite/lowrank/README.md` and the review ledger
`.superpowers/sdd/2026-10-30-delayed-submatrix/progress.md`. Facts that are not in these sources are marked "not recorded".

## 1. Summary

Two variants of ALF's sequential single-field update were added: `delayed` and `submatrix`.

- `delayed` keeps accepted rank-1 corrections as G = G0 − X Yᵀ and applies them with one ZGEMM.
- `submatrix` is the DQMC submatrix update of Sun and Xu [SX-S]. Between flushes it keeps only Γ̃⁻¹, of size k×k.

The scheme is chosen at runtime in `&VAR_QMC`:

```
Sequential_scheme = "standard"   ! "standard" | "delayed" | "submatrix", case-insensitive
N_delay           = 32           ! max pending accepted vertices before a flush
```

Both keys are optional. The defaults are `standard` and 32.

- **Bit-identity.** With `standard`, the branch is bit-identical to the baseline on 9 regression sets.
- **Speed.** For the Mz Hubbard model at L=20 (N=400, β=2, 1 thread, Apple M1 Pro, N_delay=32), the time per sweep drops
  by a factor of 2.2 (delayed) and 2.3 (submatrix).
- **Correctness.**
  - Both schemes make the same Metropolis decisions as `standard` in every trajectory test.
  - They pass the statistical comparison on three sets.
  - With the same seed, they reproduce the standard trajectory over the full statistical runs.

## 2. Method

The kernels, wraps and flushes follow notes §1–§5.

These files changed:
- `Prog/Lowrank_update_mod.F90` (new);
- `Prog/Wrapgr_mod.F90` (hooks);
- `Prog/QMC_runtime_var_mod.F90` (namelist, getters, MPI broadcast);
- `Prog/Makefile`;
- `Prog/Trace_updates_mod.F90` (new);
- `Prog/upgrade_mod.F90`, which only gains a trace call under `#ifdef ALF_TRACE_UPDATES`.

`main.F90` and the Hamiltonians are unchanged.

**Delayed** (notes §3). For each computed flavour, G_j = G0 − X_j Y_jᵀ.

| Operation | What it does | Cost |
|---|---|---|
| Proposal at site p | needs only G_j(p,p) = G0(p,p) − X(p,:)·Y(p,:) | O(j) |
| Accept | appends one column to X and one to Y (two ZGEMVs) | O(N j) |
| Flush | G0 ← G0 − X Yᵀ, one ZGEMM | BLAS-3 |
| Diagonal vertex wrap | rescales rows of X and Y | O(j) |

Repeated sites need no special treatment.

**Submatrix** (notes §4). For distinct pending sites S with D = diag(Δ):

- G_k = G0 − G0(:,S) Γ̃⁻¹ D (1 − G0)(S,:), with Γ̃ = I + D (I − G0)_SS.

| Operation | What it does | Cost |
|---|---|---|
| Proposal | solve with Γ̃⁻¹ | O(k²) |
| Accept | borders Γ̃⁻¹; the Schur complement equals the fermionic ratio factor d | O(k²) |
| Flush | two ZGEMMs | BLAS-3 |
| Diagonal wrap | rescales one row and one column of Γ̃⁻¹ for every pending site it touches | O(k) per site |

For Hubbard, the wrap rescale happens only in the down-sweep.

**Ratio, weight and RNG order** (notes §1). The only thing either scheme computes differently is G_pp.

- These repeat `Upgrade2`'s FINAL-mode expressions, in the same order:
  - the ratio and `weight_reconstruction`;
  - the `^N_SUN` power, `e^{Z1 α}` and the γ ratio;
  - the Metropolis weight and the phase update;
  - the `Control_upgrade*` calls.
- Random numbers are drawn in the standard order: `nsigma%flip`, then the `T0_proposal` draw, then the Metropolis draw.
  `Prog/Wrapgr_mod.F90` is unchanged around these calls.

**Eligibility and fallback** (design §3.3).

- Vertex n is eligible iff, for every flavour nf = 1..N_FL, `Op_V(n,nf)%N_non_zero ≤ 1` and `Op_V(n,nf)%diag` both
  hold.
- This covers all Hubbard variants: Mz, density channel and continuous HS.
- An ineligible vertex flushes the pending state, then goes through the unchanged `Upgrade2`.
- Eligibility is computed once, at init.
- For a non-standard scheme, `Lowrank_init` writes `Sequential scheme: <NAME>, N_delay = <n>` to stdout and to `info`.
- If ineligible vertices exist, it also writes one line, once:
  `Sequential scheme <NAME>: <m> of <M> vertices are not diagonal rank-1; they use the standard update (first: n=<i>)`.
- An unknown scheme name, or `N_delay < 1`, terminates the run.

**Flush points.** The pending state is applied to G:

1. before an ineligible vertex;
2. at the end of the vertex loop of `WRAPGRUP`, which is before its Global_tau moves and before the return;
3. at the end of the vertex loop of `WRAPGRDO`, which is before the hopping propagation and the return.
   `WRAPGRDO`'s Global_tau moves run at entry, on the flushed G of the previous call;
4. for submatrix, before proposing at a site that is already pending in any flavour (notes §4, restriction 2);
5. when `N_delay` pending vertices, or the column capacity, is reached.

Because of 2 and 3, every reader outside the two routines (stabilisation, `CGR`, measurements, time-displaced, global
moves) sees a flushed G. It also means the effective delay is min(N_delay, accepts per slice).

**Stale-G protection.**

- `Lowrank_assert_flushed` runs at entry and exit of `WRAPGRUP`/`WRAPGRDO`, in every build. It is one integer compare,
  and it terminates the run if anything is pending.
- `Lowrank_propose`, `Lowrank_accept` and `Lowrank_wrap_hook` terminate on an ineligible vertex, or on an accept without a
  matching proposal.
- Submatrix accept also checks that the site is not already pending (an O(1) slot check).
- The debug build (`-DALF_LOWRANK_DEBUG`) carries a shadow G with standard rank-1 updates.
  - After every accept, hook and flush it checks that max|G_lowrank − G_shadow|/max|G_shadow| ≤ 1e-10.
  - It compares the proposal's d with the shadow's d, which catches a stale G_pp read.
  - On a violation it aborts with n, τ and the flavour.
  - It writes the check count and the largest error to `lowrank_debug_summary.txt`.

**Memory cap.**

- The capacity in columns is N_cols_max = min(N_delay·N_P_max, M·N_P_max), with N_P_max = 1 and M the number of vertices
  per slice. For submatrix it is further capped at Ndim.
- Delayed allocates X and Y, each of size Ndim × N_cols_max × N_FL_eff.
- Submatrix allocates Γ̃⁻¹, of size N_cols_max² × N_FL_eff, plus three Ndim × N_cols_max work arrays.
- So a very large N_delay does not over-allocate: the `start` set at N_delay=100000 peaked at 8 MiB RSS.
- The data layout counts columns, not vertices, so a rank-N_P (Woodbury) extension only changes the bookkeeping
  (notes §3).

**Initialisation.** `main.F90` calls `Wrapgr_alloc` only when Global_tau moves are on. The module is therefore
initialised lazily, at the first `WRAPGRUP`/`WRAPGRDO` call (Task 8 ruling). This supersedes the init in `Wrapgr_alloc`
described in design §3.1.

## 3. Validation

The levels and their pass criteria are defined in design §4. The tolerances and their justification are in `CLAUDE.md`,
section "Tolerances". The Fortran under `Prog/` and the unit test are unchanged since `ab81f3fd`: `git diff ab81f3fd..HEAD
-- Prog testsuite/Prog.tests` is empty. Final re-verification on the report commit: see §3.6.

### 3.1 Unit (`testsuite/Prog.tests/36-lowrank-unit.F90`, ctest `36-lowrank-unit`)

Setup:
- a synthetic G = (1+A)⁻¹, with Ndim 32–64;
- three slices up and down, with the real `Op_Wrapup/do`;
- a fixed accept mask;
- the real `Upgrade2` as reference.

The table shows the maximum, over up/down and N_delay ∈ {1, 7, 32, 1000}, of two errors: the relative ratio error, and
max|G_cand − G_ref|/max|G_ref|. The values come from the ctest output of 2026-10-31 07:39. That was the Task 15 run, with
Fortran identical to `ab81f3fd`; the values agree with `CLAUDE.md` Phase 3/4.

| case | what it exercises | cond(1+A) up / down | std drift vs exact (max) | delayed ratio / G | submatrix ratio / G | tol |
|---|---|---|---|---|---|---|
| mz | Mz pair, N_FL=2, type 2 | 4.6e6 / 6.3e6 | 1.6e-13 | 2.9e-13 / 2.3e-14 | 3.0e-13 / 2.0e-14 | 1e-10 |
| density | complex g, N_FL=1, N_SUN=2 | 2.7e2 / 2.0e2 | 1.3e-14 | 1.5e-14 / 7.3e-15 | 1.5e-14 / 7.5e-15 | 1e-10 |
| continuous | type 3 field | 1.5e3 / 1.2e3 | 7.2e-15 | 3.6e-16 / 1.0e-15 | 3.6e-16 / 1.0e-15 | 1e-10 |
| repeat | same site twice per slice (submatrix repeated-site flush) | 2.6e6 / 7.7e6 | 1.5e-12 | 8.6e-13 / 7.4e-14 | 3.8e-13 / 7.5e-14 | 1e-10 |
| mixed | ineligible V_SUN bonds between Mz vertices (mid-slice flush) | 5.7e6 / 1.9e6 | 1.6e-12 | 2.1e-14 / 1.0e-14 | 2.5e-14 / 9.6e-15 | 1e-10 |
| asym | different site permutation per flavour, Δ = 0 in flavour 2 | 2.5e6 / 5.2e6 | 1.1e-13 | 7.5e-14 / 1.7e-14 | 5.3e-14 / 1.9e-14 | 1e-10 |
| repeat_skip | repeat with vertices wrapped but not proposed | 2.6e8 / 5.4e7 | 9.4e-14 | 5.2e-14 / 1.2e-14 | 7.5e-14 / 1.3e-14 | 1e-10 |
| mz_illcond | mz with an ill-conditioned inter-slice similarity | 1.5e7 / 1.9e8 | 3.5e-11 | 2.0e-12 / 1.4e-12 | 2.6e-10 / 8.8e-11 | 1e-8 |

- **Standard.** `standard` vs `Upgrade2` gives exactly 0 in every case. The standard drift against an exact (1+A)⁻¹
  stays between 4e-5 and 0.22 × cond·ε (`CLAUDE.md`, Tolerances).
- **Submatrix in `mz_illcond`.** Its error is the forward error of the explicit Γ̃⁻¹ (cond(Γ̃) up to ~1e6), not an
  ordering error. Two independent reviews (Task 11, Phase 4) reached this conclusion.
- **Kernel equivalence**, `Upgrade2` FINAL vs `Upgrade_lowrank` on case mz, 200 proposals, both schemes, all N_delay:
  0 toggle mismatches, |ΔPhase| ≤ 2e-15, identical next ranf.
- **Eligibility.** Diagonal rank-1 vertices, rank-2 V_SUN bonds and non-diagonal rank-1 vertices are classified
  correctly. The one-time message and the `info` line are checked.
- **Mutations.** Removing the Γ̃⁻¹ rescale, or the repeated-site flush, makes 15–30 rows fail (Tasks 11, 13).
- **Verified on** `3de8a5f9` (delayed) and `ab81f3fd` (submatrix). ctest passed 40/40: 39 baseline tests plus
  `36-lowrank-unit`.

### 3.2 Trajectory (`compare_trajectory.sh`, trace build, NSweep=3 NBin=2, same seeds)

Pass criteria:
- identical accept/reject decisions for every proposal;
- `confout_0` identical field by field;
- all bins within 1e-8 relative.

A decision difference within 1e-10 of ranf would count as ROUNDOFF. None occurred.

| set | model / feature | std vs delayed | std vs submatrix |
|---|---|---|---|
| start | Hubbard Mz 6x6 (Start example) | IDENTICAL | IDENTICAL |
| doped | Mz 4x4, μ=−0.5, β=4 | IDENTICAL | IDENTICAL |
| density | `Mz=.F.` (complex g, N_FL=1) | IDENTICAL | IDENTICAL |
| continuous | `Continuous=.T.` | IDENTICAL | IDENTICAL |
| kondo | Kondo bilayer, 16 of 32 vertices ineligible per slice | IDENTICAL | IDENTICAL |
| pv_adiab | Plain_Vanilla, projector, adiabatic g_t, N_FL_eff=1 < N_FL=2 | IDENTICAL | IDENTICAL |
| z2gt | Z2_Matter, Global_tau moves, N_delay=5 | IDENTICAL | IDENTICAL |
| start, N_delay=3 | frequent capacity flushes | IDENTICAL | IDENTICAL |
| start, `Propose_S0=.T.` | a no-op for type-2 fields | IDENTICAL | IDENTICAL |
| start, N_delay=100000 | memory cap (peak RSS 8 MiB) | not recorded | IDENTICAL |
| tv | all vertices ineligible: fallback | not recorded | byte-identical, message once |

- **Further identical runs.**
  - Delayed vs submatrix is IDENTICAL on start and doped.
  - Delayed with N_delay 1…5 is IDENTICAL on start (Task 8).
  - The Phase-3 reviewer also ran delayed on a non-projector Plain_Vanilla run with flavour reconstruction. It was
    IDENTICAL, but it is not a repeatable set.
- **Proposal counts.** Start has 21600 proposals per run. For the other sets the counts are not recorded.
- **Ratio agreement on start.** The median relative difference is 3e-15. The maximum (2.5e-8 submatrix, 4.7e-8 delayed)
  occurs only at |R| ~ 1e-5…1e-4, from cancellation in d, and is the same for both schemes.
- **Verified on:**
  - `3de8a5f9`: delayed, 8 sets;
  - `ab81f3fd`: 14 comparisons. These are all sets std vs submatrix, N_delay 3 and 1e5, `Propose_S0`, delayed vs submatrix
    on start and doped, and std vs delayed on start;
  - Task 15 re-ran start std vs delayed and std vs submatrix, with the Fortran unchanged.

### 3.3 Statistical (`stat_runs.sh` + `stat_compare.py`, independent seeds)

All sets are 6x6 with U=4:

| set | β | μ | field | sweeps |
|---|---|---|---|---|
| half | 5 | 0 | Mz | 2000 |
| doped | 5 | −1.0 | Mz, average sign ≈ 0.59 | 4000 (NBin=80) |
| density | 4 | 0 | `Mz=.F.` | 2000 |

Pass rule:
- every scalar must have |z| < 3;
- every k-resolved array must pass a binomial 3σ bound on the fraction of |z| > 2 over symmetry-unique entries, plus
  max|z| < 4.5;
- a replicate is required only after an array-only failure. This protocol was pre-registered in
  `testsuite/lowrank/README.md`.

The rule was calibrated on 4 standard chains: 16 of 18 standard-vs-standard comparisons pass, and no scalar failed.

| comparison | half | doped | density | replicate needed |
|---|---|---|---|---|
| delayed s3 vs standard s1 | PASS; Ener z=1.92, Kin z=2.06, Pot z=−0.29; Green_eqJK 3/26 entries with abs(z)>2 | PASS; sign 0.603(11) vs 0.590(10), z=−0.92 | PASS (z not recorded) | no |
| submatrix s4 vs standard s1 | PASS; Ener z=1.18, Kin z=1.37 | PASS; sign 0.601(10) vs 0.590(10), z=−0.78; Part z=−2.00 | PASS (z not recorded) | no |

- **Scalars.** All scalar |z| are < 3 for delayed and ≤ 2.00 for submatrix. The other individual z values are not
  recorded.
- **Commits.** The submatrix runs were verified on `ab81f3fd` (re-run `submatrix_s4_v4`, PASS ×3). The commit of the
  binary used for the delayed runs is not recorded separately; Phase 3 was verified at `3de8a5f9`.

**Same-seed identity.** With the same seed, the statistical runs follow the standard trajectory exactly. For all three
sets, over the full run:
- `confout_0` of `submatrix_s4` equals that of `standard_s4`;
- `delayed_s3` equals `standard_s3`.

So the schemes are not just statistically equivalent: over thousands of sweeps they make the same decisions.

### 3.4 Regression (`regress_vs_baseline.sh`)

- **Result.** With default parameters, all 9 sets are IDENTICAL to baseline `f81b5cbc`. Every output file except `info`
  is compared byte by byte with `cmp`.
- **Sets.** start, explicit (`sequential_scheme='STANDARD' N_delay=7`), doped, density, continuous, tv, kondo, pv_adiab
  and z2gt.
- **Harness sanity check.** With `BRANCH_EXTRA=NSweep=6`, every set reports DIFFERENT, as it should.
- **Verified on** `3de8a5f9`, `ab81f3fd` and again in Task 15, with ctest 40/40 each time.

### 3.5 Stability (`Precision Green/Phase` in `info`)

Paired means same seed, so the same trajectory, and the ratio isolates round-off. The submatrix rows use the final binary
`ab81f3fd` (`submatrix_s4_v4` vs `standard_s4`).

| ratio of Precision ... max | half | doped | density |
|---|---|---|---|
| delayed s3 / standard s1 (different trajectories), Green | 1.11 | 0.25 | 0.61 |
| delayed s3 / standard s3 (paired), Green | 1.00 | 0.94 | 0.89 |
| delayed s3 / standard s3 (paired), Phase | 0.98 | 0.94 | 0.94 |
| submatrix s4 / standard s4 (paired), Green | 0.99 | 1.03 | 1.92 (2.2e-11 vs 1.1e-11) |
| submatrix s4 / standard s4 (paired), Phase | 0.98 | 0.94 | 1.90 (2.0e-12 vs 1.0e-12) |

- **Means.** The paired Green means are 1.03–1.12 for submatrix.
- **Earlier binary.** The first submatrix statistical run used a binary built before the Phase-4 fix round (b3caa12c with
  uncommitted debug-only edits). On the same trajectory it gave a density Phase max ratio of 2.33. The final binary gives
  1.90.
- **Criterion.** All of these meet the criterion in design §4: maxima within a factor of a few, and comparable means.
- **Benchmark runs** (L=12 and 20; N_delay 4–256). The Green max of every run is within 0.5 % of standard. The submatrix
  Green mean rises mildly with N_delay: +2.6 % (L=12) and +4.2 % (L=20) from N_delay=4 to 256 (`ndelay_scan.md`).
- **Debug shadow mode.** The maximum relative error is 2.0e-14 for delayed (~1.35e5 checks; start, doped, N_delay=3) and
  3.0e-14 for submatrix.
- **Self-tests.** Both fire: delayed at 4.1e-2, and submatrix at 2.7e-1 in the down-sweep only. In Hubbard the up-sweep
  hook never touches a pending site.

### 3.6 Final verification

This is filled in after the whole-branch review (Task 16 Step 2); see the `CLAUDE.md` section "Final report".

## 4. Performance

**Setup**
- Model: Mz Hubbard on the square lattice, U=4, β=2, Δτ=0.1 (Ltrot=20), one measurement slice, default Nwrap=10.
- Seeds: the same seed for all schemes. The trajectories are therefore identical, and the acceptance agrees to all
  printed digits.
- Timing: wall time of the sweep phase (`CPU Time` in `info`); setup is excluded.
- Benchmark binary: `26b354c7`, with Fortran identical to `ab81f3fd`.
- Machine: Apple M1 Pro (8P + 2E cores, 32 GB), macOS 27.2.
- Compiler: GNU Fortran 16.2.0 (Homebrew GCC), flags
  `-std=f2008 -cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`.
- BLAS/LAPACK: Apple Accelerate (vecLib), via `-llapack -lblas`.
- Threads: OMP_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1 unless stated.
- Details: `size_scan.md`, `ndelay_scan.md`, `profile.md`, `results.csv`, and the raw data `raw_*.txt`.

### 4.1 Time per sweep vs N (N_delay=32, 1 thread; `size_scan.md`)

| L | N | standard [s] | delayed [s] | submatrix [s] | delayed speedup | submatrix speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 36 | 0.00724 | 0.00661 | 0.00701 | 1.10 | 1.03 |
| 8 | 64 | 0.0236 | 0.0200 | 0.0210 | 1.18 | 1.13 |
| 10 | 100 | 0.0725 | 0.0544 | 0.0546 | 1.33 | 1.33 |
| 12 | 144 | 0.1815 | 0.1227 | 0.1185 | 1.48 | 1.53 |
| 14 | 196 | 0.368 | 0.219 | 0.210 | 1.68 | 1.75 |
| 16 | 256 | 1.098 | 0.728 | 0.700 | 1.51 | 1.57 |
| 18 | 324 | 1.423 | 0.714 | 0.677 | 1.99 | 2.10 |
| 20 | 400 | 2.560 | 1.178 | 1.121 | 2.17 | 2.28 |

- **Second L=20 measurement.** The N_delay scan has an independent L=20 run at N_delay=32: standard 2.592 s, delayed
  2.21x, submatrix 2.33x. This agrees with the table within the noise.
- **L=16 outlier.**
  - All three schemes jump by ×3.0–3.3 from L=14. The N^2.46 trend of standard predicts ×1.9.
  - The low-rank schemes jump relatively more, so their speedup dips; it recovers at L=18.
  - This is consistent with a power-of-two (leading dimension 256) cache effect. That is a guess, not a proof: there is
    one run per point.

### 4.2 N_delay scan (`ndelay_scan.md`)

Time relative to N_delay=32 (positive = slower):

| L | scheme | 4 | 8 | 16 | 32 | 64 | 128 | 256 | speedup at 32 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | delayed | +12.9 % | +1.9 % | −1.6 % | 0 | +6.5 % | +5.3 % | +5.3 % | 1.48 |
| 12 | submatrix | +23.8 % | +15.7 % | −0.8 % | 0 | +10.9 % | +14.2 % | +11.6 % | 1.56 |
| 20 | delayed | +47.6 % | +17.8 % | +4.9 % | 0 | +5.1 % | +18.8 % | +33.0 % | 2.21 |
| 20 | submatrix | +59.7 % | +24.0 % | +8.7 % | 0 | +1.1 % | +16.6 % | +44.4 % | 2.33 |

- **Optimum.** The best value is 16 at L=12 and 32 at L=20, for both schemes. Within the 2.3 % noise of the best are:
  - 16–32 at L=12;
  - 32 for delayed at L=20;
  - 32–64 for submatrix at L=20.
- **Plateau and rise.**
  - At L=12 (about 62 accepts per slice), 64–256 sit on a plateau, because the per-slice flush caps the effective delay.
  - At L=20 (about 172 accepts per slice), the rise sets in below the cap.
  - We attribute the rise to O(N k) per-accept costs (delayed) and O(k²) per-proposal and per-accept costs (submatrix).
    This is an interpretation; no profile was taken as a function of N_delay.
- **Precision.** The Green max is within 0.5 % of standard for every N_delay. At N_delay=32, the Green mean differs from
  standard by +1.2 % / 0.0 % for submatrix (L=12 / L=20), and +0.1 % / −0.7 % for delayed.

### 4.3 Profile at L=20 (`profile.md`, macOS `sample`, ~12.5k samples each)

| category | standard | submatrix |
|---|---:|---:|
| rank-1 ZGERU in `Upgrade2` | 63.1 % (1.64 s) | – |
| flush ZGEMM | – | 19.0 % (0.21 s) |
| ratio / border / propose | 1.4 % | 2.2 % (0.025 s) |
| hopping / interaction wraps | 8.5 % (0.22 s) | 18.5 % (0.21 s) |
| stabilisation (QR/UDV, `CGR`, `WRAPUR/UL`) | 26.7 % (0.69 s) | 59.5 % (0.66 s) |
| measurements | 0.3 % | 0.7 % |

- Submatrix cuts the update cost from 1.64 s to about 0.24 s per sweep, a factor of 7.
- Stabilisation and wraps cost the same absolute time in both schemes, and make up 78 % of the submatrix sweep.
- Even a free update would give at most 2.592 / (0.692 + 0.220 + 0.036 + 0.008) ≈ 2.7x at this size. Submatrix reaches
  2.3x.
- Delayed was not profiled.

### 4.4 Threads

| threads | standard [s] | delayed [s] | submatrix [s] |
|---:|---:|---:|---:|
| 1 | 2.592 | 1.178 | 1.121 |
| 8 | 2.596 | 1.185 | 1.121 |

8 threads give no gain for any scheme, so the speedups (2.19x, 2.32x) are unchanged. The reason was not investigated.

## 5. Recommendation

- **Which scheme.** For production runs of Hubbard-type models with N ≳ 150, use `Sequential_scheme = "submatrix"` or
  `"delayed"` with the default `N_delay = 32`.
  - Submatrix is about 3–6 % faster than delayed for L ≥ 12. That is at or slightly above the 2.3 % noise.
  - Delayed has the smaller round-off, both in the ill-conditioned unit case and in the paired density precision (§3.5).
    It is the safer choice when precision matters more than a few percent of run time.
- **N_delay.** 32 is best, or within noise of best, at N=144 and N=400 for both schemes.
  - It was chosen as a single default for N = 144–400 and is untested outside that range.
  - Values of 16–64 cost at most ~11 %. Avoid N_delay ≤ 8 or ≥ 128 at large N.
- **Expected gain.**
  - It is modest at small N: 1.03–1.18x for N ≤ 64, about 1.3x at N=100 and about 1.5x at N=144.
  - At large N it is capped by stabilisation and wraps, which the update does not touch. At N=400 the bound for a free
    update is about 2.7x, and 2.2–2.3x is reached.
  - More stabilisation work (smaller Nwrap) or more measurement slices lower the speedup further.
- **When to keep `standard`.**
  - For models with mostly ineligible vertices (rank > 1, non-diagonal), where the schemes only add flushes.
  - For MPI and tempering runs, until those are validated (§6).

## 6. Caveats and limitations

**Scope**
- **Global updates and Langevin/HMC** were out of scope. They read G only after `WRAPGRUP`/`WRAPGRDO` return, that is,
  after a flush. They were not tested with the new schemes, except Global_tau through z2gt.
- **Rank > 1 and non-diagonal vertices** are not accelerated. They fall back to `Upgrade2` after a flush (notes §4,
  design §7 Q1). A model whose vertices are mostly ineligible sees no gain.
- **`Sequential_scheme`** is silently ignored when `Sequential = .F.`.
- **`N_delay`** is validated (≥ 1) only for the non-standard schemes.

**Coverage gaps**
- **MPI and tempering** were only compiled; the Phase-3 review built the MPI and Tempering+LOG versions.
  - No MPI or tempering run was made.
  - `Lowrank_init` assumes the tempering `Temp_<g>/info` layout.
  - In the debug build, every MPI rank rewrites the same `lowrank_debug_summary.txt`, which is a race.
  - The trace output is single-rank only.
- **Projector code** is covered only by the `pv_adiab` trajectory set (delayed and submatrix, IDENTICAL decisions). There
  is no statistical test and no unit case.
- **`N_FL_eff < N_FL`** (`weight_reconstruction`):
  - It is covered by the `pv_adiab` trajectory set (both schemes), and by one ad-hoc non-projector Plain_Vanilla run in
    the Phase-3 review (delayed only, not repeatable).
  - It is not covered at unit level or statistically.
- **Other gaps at unit level.**
  - The ledger (Task 7) also lists `g_t`, an `N_non_zero = 0` flavour and type-1 Ising fields as untested at unit level.
  - An eligible type-1 vertex with `Propose_S0` and T0 < 1 is not covered. `Propose_S0` was only tested on Hubbard type-2
    fields, where it is a no-op.
  - Whether kondo or z2gt exercise eligible type-1 vertices is not recorded.
- **The capacity-trigger flush** (capacity reached while `n_pend < N_delay`) is not demonstrably covered by a test
  (Task 13).
- **`LOWRANK_DEBUG_SELFTEST_PROPOSE`** was exercised only for delayed, not for submatrix. The submatrix hook self-test
  fires only in the down-sweep for Hubbard.
- **The `asym` unit case** does not separate the per-flavour rescale branch of the submatrix hook (Task 13 minor).
- **The determinant ratio** is not checked directly against an exact determinant; it is checked against `Upgrade2`
  (Task 3 minor).

**Numerics**
- **Explicit Γ̃⁻¹ in submatrix.**
  - In the ill-conditioned unit case its error is 2.6e-10 / 8.8e-11 (ratio / G), against 2.0e-12 / 1.4e-12 for delayed.
    This is within TOL_ILL = 1e-8 and consistent with cond(Γ̃)·ε.
  - In the paired statistical runs, the submatrix Green and Phase precision maxima reach 1.9x standard on the density set.
  - These are within the approved criterion, but they are a systematic, if small, loss compared with delayed (paired
    ≤ 1.00).
- **Green precision mean.** For submatrix it rises by 2.6–4.2 % from N_delay=4 to 256.
- **Trajectory ratio difference.** It reaches ~5e-8 at |R| ~ 1e-5, from cancellation in d. It is the same for both
  schemes and does not affect the decisions.

**Implementation**
- **Lazy init.** The scheme name is validated only at the first `WRAPGRUP`/`WRAPGRDO`.
  - A typo therefore aborts the run after setup (about 24 s at L=20), not at parameter reading.
  - No finalize runs at program exit, which is harmless.
- **Stack arrays.** The automatic arrays `v`, `b`, `r`, of size N_cols_max, live on the stack in `Lowrank_propose` and
  `Lowrank_accept`. With N_cols_max ≤ M this is harmless, but a large capacity could exceed a small stack limit (Task 11
  minor).
- **Effective delay.** The per-slice flush caps the effective delay at the number of accepts per slice. An N_delay larger
  than that has no effect, except on allocation (up to the cap).
- **NaN test in debug builds.** `-ffast-math` may optimise away the NaN test of the debug check (Task 9 minor).
- **pyALF.** pyALF's `Prog/default_parameters_generic.py` does not know `Sequential_scheme` or `N_delay`. That file is
  outside the allowed set (Task 2 minor).

**Benchmarks**
- **One platform.** One machine (M1 Pro), one compiler, one BLAS (Accelerate). OpenBLAS is installed but was not
  compared.
- **Noise.** Timing noise is about 2 %: a 2.3 % spread over 3 runs of standard at L=12, the only repeated point.
  Differences of a few percent between delayed and submatrix are not significant.
- **L=16 anomaly** (§4.1). There is one run per point, and the cause is not established.
- **Run sizes.** L=20 has 24 sweeps per run; β=2 only; one measurement slice. Runs that measure on every slice, or use a
  different Nwrap, change the time fractions.
- **Setup cost.** The ~24 s setup at L=20 is `validate_ham_variables` → `Op_V_is_symmetric`, using gfortran's `matmul`.
  This is a pre-existing ALF cost, unrelated to this work, but it dominates short runs at large N.

## 7. Reproducing

Environment variables do not persist between shell calls, so always source `configure.sh` in the same command. Run from
`ALF/`.

```bash
# build (normal binary) and the CTest suite incl. 36-lowrank-unit (40 tests)
source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && make all -j8
source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && mkdir -p testsuite/tests && cd testsuite/tests \
  && cmake .. && make -j8 && ctest -j8          # ctest -R 36-lowrank-unit -V  for the unit table
```

Rebuild the normal binary (no `ALF_FLAGS_EXT`) and re-run `cmake ..` before `ctest`, because test 36 links `Prog/*.o`.

```bash
# regression vs baseline f81b5cbc (9 sets, exit 0 iff all IDENTICAL)
testsuite/lowrank/regress_vs_baseline.sh
BRANCH_EXTRA=NSweep=6 testsuite/lowrank/regress_vs_baseline.sh   # harness sanity: all DIFFERENT

# trajectory (sets: start doped density continuous tv kondo pv_adiab z2gt)
testsuite/lowrank/compare_trajectory.sh start standard delayed
testsuite/lowrank/compare_trajectory.sh start standard submatrix
testsuite/lowrank/compare_trajectory.sh z2gt  standard submatrix N_delay=5
testsuite/lowrank/compare_trajectory.sh start standard submatrix N_delay=100000
testsuite/lowrank/compare_trajectory.sh start delayed  submatrix
python3 -I testsuite/lowrank/trace_diff.py --selftest

# statistical (each job 25-105 s; run in background, <= 8 concurrent)
for set in half doped density; do
  testsuite/lowrank/stat_runs.sh $set standard 1; testsuite/lowrank/stat_runs.sh $set delayed 3
  testsuite/lowrank/stat_runs.sh $set submatrix 4
  testsuite/lowrank/stat_runs.sh $set standard 3; testsuite/lowrank/stat_runs.sh $set standard 4   # for the paired checks
done
python3 -I testsuite/lowrank/stat_compare.py runs/stat/half/standard_s1 runs/stat/half/submatrix_s4
cmp runs/stat/half/submatrix_s4/confout_0 runs/stat/half/standard_s4/confout_0     # same-seed identity
python3 -I testsuite/lowrank/stat_compare.py --selftest

# debug shadow mode / self-tests (rebuild the normal binary afterwards)
export ALF_FLAGS_EXT="-DALF_LOWRANK_DEBUG"     # add -DLOWRANK_DEBUG_SELFTEST: the run must abort with "stale G"
source configure.sh GNU noMPI NO-INTERACTIVE >/dev/null 2>&1 && make -C Prog clean && make program -j8
```

**Benchmarks.** Copy the normal binary to `runs/bench/ALF.out`; it stays frozen there, and `BENCH_BIN` overrides it.
Each job prints one raw line:

```bash
testsuite/lowrank/bench.sh SCHEME L N_DELAY [THREADS]   # e.g. bench.sh submatrix 20 32 >> docs/benchmarks/raw_size_scan.txt
# size scan: L in 6..20 (even), all three schemes, N_delay=32; N_delay scan: L=12,20, N_delay in 4..256; threads: THREADS=8 at L=20
python3 -I testsuite/lowrank/bench_report.py docs/benchmarks/raw_size_scan.txt --repeat docs/benchmarks/raw_repeat.txt
python3 -I testsuite/lowrank/bench_report.py --ndelay docs/benchmarks/raw_ndelay_scan.txt \
    --size-scan docs/benchmarks/raw_size_scan.txt --repeat docs/benchmarks/raw_repeat.txt \
    --threads docs/benchmarks/raw_threads.txt
```

- Regenerating from the committed raw files must leave git clean.
- **Profile.** Run `sample <pid> 20 -file <out>` on a running L=20 job after its setup (see profile.md), then
  `python3 -I testsuite/lowrank/profile_samples.py <sample files>`. The exact launch command of the profiled jobs is not
  recorded (profile.md: NSweep=25, sampled from 29 s after launch).

## 8. Commits (`git log --oneline f81b5cbc..93826bc0`)

Design and plan
- `f9bb1707` Add design notes for delayed/submatrix local updates
- `c138d41c` Add implementation plan for delayed/submatrix updates

Phase 2: validation harness (Tasks 1–5)
- `50435365` test: add bit-identity regression harness vs baseline
- `c661c799` test: regress harness compares outputs only and sanity-checks runs
- `dbbfae2d` feat: add Sequential_scheme and N_delay runtime parameters
- `def28365` test: make the explicit regression set exercise parser and N_delay irrelevance
- `92301028` test: add unit-level low-rank update harness (standard vs standard)
- `e42d47d4` test: add update trace and trajectory comparison harness
- `3a91358c` test: restore normal build on trajectory-script exit, harden trace_diff and compare_runs
- `22717b67` test: add statistical comparison harness
- `80ae5322` test: unique-entry binomial array rule, doped set at beta=5 (sign in window)
- `d39e45a7` test: particle-hole merge default in stat_compare, document validation protocol
- `c7c4dddf` test: verify scheme took effect in stat_runs, fail on non-finite z, validate scheme name
- `34f2503f` test: stat_compare fails on non-finite input before the zero-entry skip

Phase 3: delayed (Tasks 6–10)
- `99ce01ba` feat: add Lowrank_update_mod skeleton with eligibility
- `a9d934f3` feat: delayed update kernel, wrap hook and BLAS-3 flush
- `cd9c8f64` feat: route sequential updates through delayed scheme
- `289f7009` feat: shadow-G debug mode for low-rank updates
- `de0c222a` fix: debug summary file instead of main.F90 change
- `e9174f7c` fix: debug shadow uses its own d and checks the proposal's d
- `3de8a5f9` fix: phase 3 review findings (model coverage sets, submatrix init guard, N_P_max, harness reporting)
- `12a27b05` docs: record Phase 3 (delayed) validation results

Phase 4: submatrix (Tasks 11–13)
- `c3444777` feat: submatrix update kernel
- `b3caa12c` docs: fix README wording for submatrix trajectory check
- `ec524b1b` test: submatrix trajectory validation
- `ab81f3fd` fix: phase 4 review findings (capacity cap, slot invariant, unit coverage, hook early exit)
- `26b354c7` docs: record Phase 4 (submatrix) validation results

Phase 5: benchmarks and tuning (Tasks 14–15)
- `40093a44` bench: add benchmark harness and size scan results
- `642f16a2` bench: report provenance from info, commit raw scan data
- `130cc07f` bench: N_delay scan, L=18, profile and threads
- `2bd85893` perf: set default N_delay from scan (its message says "optimum 16-64"; corrected in `93826bc0`)
- `93826bc0` docs: correct N_delay scan claims, single generated source

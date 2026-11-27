# Known issues

Open and recently fixed problems of this branch line (`Code-modularity` and the branches cut from
it). Every entry: where it applies, symptom, cause, status, and the branch or commit of the fix.
New findings of the bug hunt (docs/ai/phase2/ROADMAP.md, P2.8) are added here first.

Classes (owner rule 2026-11-26, AGENTS.md "Working mode"): **A** silently wrong results in the current code path or a
production configuration (or a guard hiding them) — fixed at once as its own feature; **B** real bug elsewhere —
fixed in a batch after each phase merge; **C** cosmetic or theoretical — fixed when the code is touched.

## Index

| ID | Issue | Class | Status |
|---|---|---|---|
| DP | Determinant phase underflow, CPU (CGR / DET_C) | A | fixed in aef50743 |
| DP-GPU | Determinant phase underflow, GPU line | A | fixed on the CM line with feature `gpu-stab` (step 1d, GPU=host: `k_det_phase` + the `LU_PHASE_SAFE` rule, test `52-gpu-stab` mode `phase`); device backends in step 2; not yet in every other GPU branch |
| H1 | NaN guards removed by -ffast-math | A | fixed by feature `nan-guards` (G and scale checks by bit pattern) |
| H2 | Large beta: scale overflow in CGR -> NaN | B | open; since aef50743 the run stops (phase guard) instead of silent NaN |
| H3 | No measurements -> NaN observables | B | fixed by feature `b-batch-1` (start-up stop: NSweep < 1 or LOBS_ST > LOBS_EN) |
| H4 | confin with field value 0 -> NaN | B | fixed by feature `b-batch-1` (two-sided field check, types 1, 2, 4) |
| H5 | Global_tau_moves with partial Nt_sequential range | B | fixed in `cm/feat/b-batch-2` (commit "Wrapgr: …"; test 52) |
| H6 | Tempering with an odd number of replicas hangs | B | fixed by feature `b-batch-1` (start-up stop; upstream issue pending) |
| H7 | Bilayer_honeycomb potential/total energy wrong (also upstream) | B | fixed by feature `b-batch-1` on the CPU line (U/U2 per orbital half); GPU device weights fixed with feature `gpu-obs-tau` (step 1e: `gpu_obs_mod` uses U on the orbitals no <= Norb/2 and U2 above for the bilayer lattices; test `53-gpu-obs` case bil1/bil2 and control 2, e2e `mz-bilayer-honeycomb`); upstream issue pending |
| M1 | Dtau = 0 -> unbounded memory | B | fixed by feature `b-batch-1` (stop if Ltrot overflows) |
| M2 | Projector with Thtrot = 0 and Ltau = 1: out of bounds | B | fixed by feature `b-batch-1` (start-up stop for this combination; tau_p bounds unchanged) |
| M3 | Nwrap <= 0 accepted | B | fixed by feature `b-batch-1` (stop if Nwrap < 1) |
| M4 | Nt_sequential range not checked -> segfault | B | fixed by feature `b-batch-1` (range check when the range is not empty) |
| M5 | Langevin dt = 0 / Max_Force = 0 -> zeros; HMC/Amplitude edge cases | B | fixed in `cm/feat/b-batch-2` (commit "QMC_runtime_var: …"; start-up stop, test 42) |
| M6 | Precision Green up to 10 passes silently | B | fixed in `cm/feat/b-batch-2` (commit "control: …"; warning above 1e-3, test 51) |
| M7 | Phase update overflow for extreme R | C | open (from code reading) |
| M8 | Code-modularity not bitwise identical to master (1e-14); global-tau threshold changed | B | open (to investigate with Fakher's history) |
| L1 | master-only error-format bug | C | fixed on Code-modularity (Fakher) |
| L2 | confin from another L accepted silently | C | open |
| L3 | Default-integer overflow at very large Ndim | C | open |
| L4 | Restart not byte-identical for Langevin Ltau = 1 / HMC | C | open |
| L5 | Tempering messages | C | open |
| L6 | GPU line: confin 0 gives finite wrong numbers | B | open (GPU line) |
| G1 | GPU line: `ALF_GPU_SELFTEST` leaves test data in the device copy of G marked valid | B | fixed by feature `gpu-wrap-update` (step 1c; harmless in 1b, where the sweep did not use the mirror) |
| G2 | GPU line: wrapped G (device wraps and updates, step 1c) less accurate than on the CPU: Precision Green mean 1.4 - 2.5x the CPU value, same chain | B | fixed by feature `gpu-hopd-optin` (cause: dense hopping product `Hd`; now opt-in, `ALF_GPU_HOPD=1`) |
| G3 | GPU line: the device precision check (`Control_PrecisionG_vals`, step 1d) kept the fixed abort threshold 10 and ignored `ALF_UNSTABLE_THRESHOLD` (CPU-line feature `unstable-threshold`, merged in parallel with 1e) | B | fixed by feature `gpu-hopd-optin` (uses `G_thresh` like `Control_PrecisionG`) |
| G4 | `stp/production-mpi-intel-*` tags (Intel production on gpu-port): raw-product LU determinant phase at every CGR | C | open, tags only (production L <= 44); not on the CM line |
| G6 | GPU line: negative control `ALF_GPU_FAULT=upd` acts only on the standard (rank-1) update, not on delayed/submatrix | C | open (test gap; found by the P2.4 step 2 `gpu-intel` gate) |
| MK | Makefile: Prog objects not rebuilt after Libraries/Modules changes | B | fixed in `cm/feat/b-batch-2` (commit "build: …"; test 50) |
| WF | Trial wave function normalisation, large N_part (suspected) | B | open, not reproduced |
| P1 | Op_V_is_symmetric slow (start-up, Symm = .T.) | B | open; fast version on gpu-port |

## Determinant phase underflow for N >~ 3000 (L >~ 54) in CGR / DET_C

- **Applies to:** upstream master (`4c690b04`) and `Code-modularity` before the det-phase merge: the CPU
  code. **Fixed on this line by the det-phase feature** (see Fix below); still open upstream.
- **Status:** **confirmed on the CPU at L = 64, Nt = 16 (acceptance 0.217; control L = 48:
  0.322)** — Helma job 1451287, 2026-11-25 (Hubbard square, U = 5, 1 sweep, 32 threads). That CPU
  binary was built from the `gpu-port` tree, whose `mymats_mod.F90` (`DET_C`) is byte-identical to
  upstream master `4c690b04` and whose `DET_C`/phase lines of `cgr1_mod.F90` are unchanged.
  First found on the GPU port (entry below).
- **Cause:** `DET_C` (`Libraries/Modules/mymats_mod.F90`) multiplies the raw LU diagonal,
  `det_C = det_C*mat(i,i)`. For a unitary matrix |det| = 1, but the running product underflows
  double precision for large N. `CGR` (`Prog/cgr1_mod.F90`) forms the phase as
  `PHASE = CONJG(DET_C(RHS, N_size)); PHASE = PHASE/ABS(PHASE)` with RHS = U_R^dagger U_L (default,
  STAB3 and LOG builds; STAB1/STAB2 call `DET_C` in the same way), so the phase can become
  0/0 = NaN. The Metropolis weight is then NaN and every proposal is rejected, with return code 0,
  a finite G and NaN observables. `Control_PrecisionP` does not notice (comparisons with NaN are
  false). The same function is used by `Op_set` (`Prog/Operator_mod.F90`, normalisation det U = 1),
  the likely cause of the documented `Op_set` failure at Ndim = 4096.
- **Evidence:** the CPU run above; the GPU line (entry below). Numerically, the LU diagonal of a random
  unitary matrix (numpy/scipy, `lu_factor`) has a running product with minimum ln|P_k| = -676
  (N = 2000), -1024 (N = 3000), -1367 (N = 4000); the limit is -708, so the raw product is exactly
  0 for N = 3000 and 4000.
- **Who is affected:** runs with N >~ 2950 (square lattice L >~ 54). L <= 48 is not affected.
- **Fix (det-phase feature, 2026-11-25):** `LU_PHASE_SAFE` (`MyMats`) checks a raw LU product whose
  phase is all that is used: it is kept bit for bit if it is finite, in [1e-250, 1e250] and agrees in
  phase with the product of the normalised factors, otherwise that normalised product (modulus 1) is
  used. `DET_C_PHASE` = `DET_C` + this check; `INV_C` has an optional `DET_PHASE`. Used in `CGR` (all
  STAB variants), `Op_set` and `Compute_Fermion_Det` (`Global_mod`; there the magnitude enters
  `log|det|` and would have been log 0). `Control_Check_Phase` at the end of `CGR` stops the run if
  the phase is NaN, Inf or not of unit modulus. Normal runs are unchanged bit for bit (16 regression
  scenarios); tests `14-det-phase`, `41-check-phase-fail`. `DET_C` itself is unchanged (its
  magnitude matters elsewhere). Large-N check (Helma job 1457891, ifx + MKL, Hubbard N_FL = 2, Mz = .F., U = 5,
  Nt = 16, 1 sweep): L = 64 acceptance 0.2182 before, 0.3197 after; L = 48 identical (0.3193). Merged as aef50743.
  Whether and how this goes to the upstream developers is the owner's decision.
- **Since cpu-s1 (CG-ALF S1 on the CPU line):** the default `CGR` (default, STAB3, LOG builds) takes the phase from
  the det(U) tracked by the UDV states (Householder reflectors of the QR, modulus 1) and calls no LU, so `DET_C`
  is not on this path any more; `ALF_CGR_LU=1` restores the LU above, `ALF_CGR_LU=2` compares both.
- **See:** branch `fix/det-phase-underflow` (GPU line), `docs/gpu/KNOWN_ISSUE_DET_PHASE_UNDERFLOW.md`.

## Build: Prog objects do not depend on the modules of Libraries/Modules

- **Found:** 2026-11-25 (det-phase work). `Prog/.deps.mk` (from `Prog/gen_deps.py`) lists only the
  dependencies between `Prog/` modules. After a change of an interface in `Libraries/Modules`
  (e.g. a new optional argument of `DET_C` in `MyMats`), `make` recompiles the library but not the
  `Prog/` objects that use it; they keep calling the old interface. Seen as a segmentation fault in a
  projector run (`WaveFunction_mod.o` calling `DET_C`), gone after `make clean`.
- **Who is affected:** incremental builds after library changes. CI and fresh builds are not.
- **Workaround:** `make clean` (or `make cleanprog`) after changing anything in `Libraries/Modules`.
- **Fixed** in `cm/feat/b-batch-2` (commit "build: Prog objects rebuilt after Libraries/Modules changes"):
  `Prog/Makefile` makes every object depend on the `.mod` files of `Libraries/Modules` (gfortran rewrites a
  `.mod` only when the interface changes) and `ALF.out` on the two archives. Test `50-make-lib-deps`.

## Trial wave function normalisation for large N_part (suspected, not checked)

- `WaveFunction_mod` normalises the trial wave function with `Det(Mat, N_part)` (= `DET_C`), whose
  raw LU product can under- or overflow for large `N_part` as above. Here the magnitude matters, so
  `DET_C_PHASE` is not a fix; a log-determinant would be. Projector runs with large N_part only.
  Status: suspected from the code, not reproduced (P2.8).

## Bug hunt 2026-11-25 (P2.8, run early): summary

Found by the hunter agents A (input handling, Claude-BugHunter method), B (numerics), C (MPI, restart, differential)
and D (GPU line gpu-port on Helma H100), on upstream master
`4c690b04` and `Code-modularity`; the high-severity ones (H1-H5) re-checked by the main session. "Both" = master
and Code-modularity. Reports: hunter scratch files `bughunt/A/REPORT.md`, `bughunt/B/REPORT.md` (summarised here).
Small runs: 4x4 Hubbard (Mz) from `Scripts_and_Parameters_files/Start` with the changes listed.
None is fixed yet; each fix becomes its own dev/canvas feature.

### H1. NaN guards are removed by -ffast-math (both, high) — re-checked
`configure.sh` builds with `-ffast-math` (GNU, also in Devel) and `-fp-model=fast=2` (Intel). Under these flags
`any(a /= a)` and `ieee_is_nan(a)` are folded to false (10-line test: a NaN read from input is not detected with
`-O3 -ffast-math`, it is without). So "Green function contains NaN" (`Control_PrecisionG`, control_mod.F90) can never
fire, and the silent NaN runs H2, H3, H4 end with return code 0. The phase guard of the det-phase feature tests the
bit pattern (`Control_Is_Finite`) and is not affected. **Fixed** (feature `nan-guards`): `Control_PrecisionG`
and `testscale` (udv_state_mod) use `Control_All_Finite`, a bit-pattern test that also catches Inf.

### H2. Large beta: scale overflow in CGR gives NaN results with return code 0 (both, default build, high) — re-checked
4x4 Hubbard U = 8, Beta = 30, Nwrap = 5 (also U = 4, Beta = 60; U = 20, Beta = 15): Precision Green, Phase and tau are
NaN, every observable file holds NaN, acceptance falls to 0, return code 0. Cause: the default stabilisation multiplies
the two scale vectors directly, `TPUP(:,J) = udvr%D(:)*TPUP(:,J)*udvl%D(J)` (cgr1_mod.F90, default branch), which
overflows once the product passes 1.8e308 (~1e154 per side); the pre-check tests only element (1,1) and `testscale`
warns at 0.1*huge but never checks for Inf. The LOG build is fine on the same input; with the det-phase feature the
run now stops (phase guard), but with a phase message instead of one about the scales.

### H3. No measurements: NaN in every observable, return code 0 (both, high) — re-checked
`NSweep = 0` (also the default when the key is missing, QMC_runtime_var_mod.F90), `LOBS_ST > LOBS_EN`, or
`Sequential = .F.` with no other update. Cause: division by `Obs%N = 0` in observables_mod.F90 (bin averages).

### H4. confin with field value 0: NaN everywhere, return code 0 (both, high) — re-checked
A `confin_0` whose fields are 0. Fields_test checks only `abs(f) - 1 > 1e-8` (one-sided, Fields_mod.F90), so |f| < 1
passes; should be `abs(abs(f) - 1) > eps`.

### H5. Global_tau_moves with a partial Nt_sequential range: G not wrapped (both, high) — re-checked
`Global_tau_moves = .T., Nt_sequential_start = 3, Nt_sequential_end = 16, ham_U = 1`: return code 0, Precision Green Max
1.12 (start = 1: 2e-14). `N_Global_tau = 0` with the default range: 1.92. At U = 4 the same inputs stop with a misleading
"Unstable matrix". Cause: WRAPGRUP/WRAPGRDO wrap only the operators start..end; the others are placed by
`Wrapgr_PlaceGR` only when `N_Global_tau > 0` or `N_Global_tau_MALA > 0`. Same with the delayed/submatrix schemes.
Also with `N_Global_tau > 0`, a range starting above 1 was entered with G at vertex 0.
**Fixed** in `cm/feat/b-batch-2`: WRAPGRUP places G at `Nt_sequential_start - 1` before the range and always at
`Size(Op_V,1)` at the end; WRAPGRDO always places G at the end of the range and back to 0 after it. Test
`52-global-tau-range` (start = 3 / empty range, N_Global_tau = 0, all three schemes): Precision Green Max 5.99 (rc 0)
and abort (rc 3) before, ~1e-13 after. The default range (1..N_op) does not call the new placements.

### H6. Tempering with an odd number of replicas hangs forever (both, high; hunter C)
`mpiexec -n 3` with `mpi_per_parameter_set = 1` (3 replicas): deadlock, no output. One replica: the run finishes, then
MPI_Finalize aborts ("Communicator being freed has 48 unmatched message(s)"). Cause: main.F90 only checks that the rank
count is a multiple of `mpi_per_parameter_set`; in Global_mod `Exchange_Step` the partner loop wraps around for an odd
group count (rank 1 waits for rank 0, which talks to rank 2), `List_masters` (Isize/2 entries) is written out of bounds,
and with one group every rank is paired with itself. Expected: a clear stop "Tempering needs an even number (>= 2) of
parameter sets" (PARALLEL_PARAMS does no exchange and is not affected).

### H7. Bilayer_honeycomb: potential and total energy wrong (upstream master and both lines; hunter D, code checked)
The interaction is set up correctly (`Ham_U_vec`: U on orbitals 1..Norb/2, U2 on the rest), but the measurement of
`Pot_scal` (and so `Ener_scal`) in Hamiltonian_Hubbard_smod uses `no_I == 1 -> U`, `no_I == 2 -> U2` and nothing else.
Bilayer_square (Norb = 2) is right; Bilayer_honeycomb (Norb = 4) weights (U, U2, 0, 0) instead of (U, U, U2, U2).
The GPU device observables (gpu-port, gpu_obs_mod) use yet another weighting (U, U2, U, U). Markov chain unaffected.
Fix: weight with `Ham_U_vec(no_I)` (host and device).

### Medium
- **M1. Dtau = 0** (both): `Ltrot = nint(beta/dtau) = nint(+Inf) = 2147483647` passes the check `Ltrot <= 0`; the run
  allocates ~550 GB for the fields of a 4x4 lattice. Probably what ran the Mac out of memory during the hunt.
- **M2. Projector with Thtrot = 0** (Theta < Dtau/2) and Ltau = 1 (both): `tau_p_mod` loops NT up to Ltrot-Thtrot and
  uses NT+1, reading `nsigma%f(:, Ltrot+1)`. Code-modularity: segfault; master: return code 0 with garbage in the last
  slice. Ham_set only warns that no projection is applied.
- **M3. Nwrap <= 0 or missing** (both): Nwrap = 0 (default when missing) gives no stabilisation (integer division by 0,
  silent on arm64, SIGFPE on x86); Nwrap < 0 segfaults.
- **M4. Nt_sequential_start < 1 or Nt_sequential_end > size(Op_V,1)** with Global_tau_moves (both): segfault (index
  of Op_V out of range in WRAPGRUP); the range is never checked.
- **M5. Langevin with Delta_t_Langevin_HMC = 0 or Max_Force = 0** (both): return code 0, all bins 0 with sign 0. Related:
  HMC `Leapfrog_steps = 0` (acceptance 0, NaN "Phase diff"), `Amplitude = 0` (frozen chain, acceptance 1),
  `N_HMC_sweeps = 0` and `Ltau = 2` silently do nothing. **Fixed** in `cm/feat/b-batch-2`:
  `check_Langevin_HMC_variables` stops at start-up on each of these (Ltau not 0/1, Amplitude = 0; Langevin:
  dt or Max_Force <= 0; HMC: dt <= 0, Leapfrog_steps < 1, N_HMC_sweeps < 1).
- **M6. Precision Green up to 10 passes without a warning** (both): `Control_PrecisionG` aborts only above 10; 0.1 was
  reached silently (continuous fields, U = 20, Nwrap = 10). **Fixed** in `cm/feat/b-batch-2`: `Control_Print`
  writes a warning to info and stderr when Precision Green Max > 1e-3 (results unchanged; abort stays at 10).
- **M7. Phase update `Phase*R/sqrt(R*conjg(R))`** (both, from code reading, not reproduced): overflows for |R| > 1.3e154
  and underflows for |R| < 1.5e-162 (upgrade_mod, Wrapgr_mod); the phase would become 0 or NaN. Use R/abs(R).

- **M8. Code-modularity is not bitwise identical to master** for 7 of 17 inputs (hunter C): Hubbard with Mz = .F.,
  Kondo, LRC, tV, Spin_Peierls. Configurations are byte-identical (same Markov chain); observables and "Precision Phase"
  differ at <= 1e-14 relative. The phase of the first CGR already differs at 1e-16; root cause not found (persists with
  -fno-fast-math). Also a real change in Phase 6: the global-tau threshold in Wrapgr_mod went from `Zero = 10D-8` (1e-7)
  to `Eps_small` (1e-10), so proposals with T0_Proposal_ratio in (1e-10, 1e-7] are now attempted.

### Low
- **L1.** master only (fixed on Code-modularity): format bug `'("A","I0",...)'` makes error messages crash (unknown
  ham_name in Hamiltonian_main_mod, bad confin values in Fields_mod).
- **L2.** A confin from another L or Beta is accepted silently (the first N values are reused); L1 <= 0 gives the
  one-dimensional message.
- **L3.** Default-integer products overflow at sizes not run today: COMPARE (`N*M`, mymats_mod) for Ndim >= 46341;
  tempering MPI count `size(GR)` for Ndim >= 32769 with N_FL = 2; `N_op*Ltrot` above 2.1e9.

- **L4.** Restart is not byte-identical for Langevin with Ltau = 1 (confout <= 4e-16) and HMC (1e-14): the force reused
  from Tau_m is recomputed after a restart (qmc_driver_mod). Langevin with Ltau = 0 and all discrete cases are exact.
- **L5.** Tempering messages: np = 3 with mpi_per_parameter_set = 2 says the reverse of the problem; a missing Temp_n
  directory gives a raw "Cannot open file" error; extra Temp_ directories are ignored silently.
- **L6.** GPU line (gpu-port): a confin with field value 0 gives finite, wrong numbers on the GPU (the device factor of
  field 0 is 1, gpu_state_mod) where the CPU gives NaN (H4); neither build validates confin values.

### Performance (found during CG-ALF, not a correctness bug)
- **P1.** `Op_V_is_symmetric` (start-up check for Symm = .T.) multiplies dense Ndim x Ndim matrices for every operator,
  O(n_op Ndim^3): ~100 s at L = 24 on a CPU, much longer at L >= 40. A fast version exists on the GPU line (`gpu-port`)
  and as fork branch `ai/feat/op-v-symm-fast` (82157fe2).

### Checked and fine
Parameter edges with clear errors (Dtau < 0, Beta <= 0, N_SUN = 0, N_FL = 3, unknown lattice, odd L with checkerboard),
missing/garbled namelists, seeds and confin files (gfortran runtime errors naming the file), all models at small L,
AddressSanitizer/UBSan on Code-modularity for all models (no report), FPE traps in ALF code (with reference LAPACK).

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
- **Status on the CM line:** fixed with feature `gpu-stab` (step 1d): the GPU=host stabilisation takes the phase from the LU
  at every CGR; `k_det_phase` returns the normalised and the raw product and the host selects by the rule of `LU_PHASE_SAFE`
  (bitwise the same rule, tested in `52-gpu-stab phase` at n = 3200 with the raw product underflowing); `Control_Check_Phase`
  after every device CGR. The Helma check at L = 56 .. 96 belongs to the device backends (step 2).
- **To do:** merge the fix into every GPU branch before it is ported (P2.4); a GPU branch that has
  not merged it must not be used for N >~ 2950.

## G4. Intel production tags: raw-product LU phase at every CGR (2026-11-26, found in P2.4 step 2 `gpu-intel` step 0)

- **Applies to:** the archive tags `stp/production-mpi-intel-large-prompt` (4a25ac48), `stp/production-mpi-intel-unstable-guard`
  (200f32c0) and `stp/production-mpi-intel-large-prompt-team-combined` (288b2a43) of the fork (Intel Max 1550 production).
- **Symptom (expected, not run):** the same as "GPU: NaN determinant phase for N >~ 3000" above: for N >~ 2950 the phase becomes NaN,
  almost all proposals are rejected, return code 0.
- **Cause:** these tags have neither the det(U) tracking (oneMKL default on gpu-port since 72a3900f) nor the fix 0b3ffba4: `k_det_phase`
  takes the raw product of the LU diagonal in every CGR (`gpu_stab_mod.F90` :811, :1210 on the tags).
- **Why class C:** production on these tags ran at L <= 44 (N <= 1936), below the threshold (not checked beyond the L = 44 reproducer).
  The CM line is not affected (step 1d `LU_PHASE_SAFE` rule, step 1f tracked det(U); the Intel backend of step 2 uses both).
- **Status:** open on the tags only; do not use them for L >~ 54. No fix planned (archive tags).

## G6. GPU line: `ALF_GPU_FAULT=upd` does not reach the delayed/submatrix update (2026-11-27)

- **Applies to:** the fault injection of `gpu_update_mod` (common GPU code, all backends).
- **Symptom:** with `Sequential_scheme = "delayed"` or `"submatrix"` the negative control `ALF_GPU_DEBUG=1 ALF_GPU_FAULT=upd` is armed
  (info: "GPU fault armed") but changes nothing: the GPU run equals the CPU path to 2e-12 (Max 1550, Hubbard Mz 8x8, gate job 566631).
  With the standard scheme the same fault stops the run (rc 3) as intended.
- **Cause:** the fault scales `Zg` only in the rank-1 branch of the accepted update (`gpu_update_mod.F90` :332); the low-rank branch
  (`gpu_lr_accept`) gets the unscaled `xv`.
- **Why class C:** results are not affected; only a negative control has a gap. A test that relies on `upd` with a low-rank scheme
  would pass vacuously. Workaround: `ALF_GPU_FAULT=obs` or the standard scheme as the control.
- **Status:** open; to be fixed in the common code when `gpu_update_mod` / `gpu_lowrank_mod` are touched next.

## G2. GPU line: Precision Green of the device path (2026-11-26; fixed by feature `gpu-hopd-optin`)

- **Applies to:** GPU=host (device wraps and updates of step 1c; with and without the device stabilisation of step 1d).
- **Symptom:** Hubbard 4x4 (SU(2), Ltrot = 50, 6 sweeps), same chain (confout byte-identical to the CPU run): mean of the
  `Precision Green` line of the info file, GPU path / CPU path: Nwrap = 10: 2.0e-14 / 1.3e-14, Nwrap = 25: 7.6e-11 / 1.7e-11
  (step 1c canvas with host stabilisation: 9.5e-11 / 1.7e-11), Nwrap = 50: 8.9e-5 / 2.2e-5 (1c canvas 9.5e-5 / 2.2e-5). The
  observables agree to 1e-11 .. 1e-5, i.e. within that noise. Not caused by the device stabilisation (the 1c canvas shows it).
- **Cause (found 2026-11-26, branch `rd/g2-precision`):** the dense hopping product `Hd` of `gpu_wrap_mod` (`hopd_build`,
  `gpu_hopd_on`, on for every Ndim <= 800). For full-size matrices (G in the wraps, U in the stabilisation) the device
  multiplies by a precomputed, once-rounded product of all checkerboard factors (one ZGEMM) instead of applying the 2x2
  factors one by one as `Hop_mod` does. It matters only with `Checkerboard = .T.`: with dense hopping there is one element and the device
  path is bitwise the CPU one. Probable mechanism (not proven): the rounding error of `Hd` is the same at every slice, so
  it adds up over the Nwrap slices instead of averaging out.
- **Numbers** (4x4 Hubbard, U = 4, Beta = 5, Ltrot = 50, Symm, 20 sweeps, 12 seeds, same chain in every pair; geometric
  mean of the ratio GPU/CPU of the `Precision Green` mean; "Hd off" = debug switch that makes `gpu_hopd_on` false):

  | case | GPU / CPU | Hd off / CPU | GPU worse in |
  |---|---|---|---|
  | SU(2), checkerboard, Nwrap = 10 | 1.41 | 1.00 | 12/12 |
  | SU(2), checkerboard, Nwrap = 25 | 2.41 | 1.01 | 12/12 |
  | SU(2), checkerboard, Nwrap = 50 | 2.52 | 1.01 | 11/12 |
  | Mz, checkerboard, Nwrap = 10 | 1.15 | 1.00 | 11/12 |
  | SU(2), dense hopping, Nwrap = 25 | 1.01 | 1.01 | 4/12 |
  | Mz, dense hopping, Nwrap = 10 | 1.02 | 1.02 | 7/12 |

  (Mz with Nwrap >= 25 at these parameters stops with "Unstable matrix" on the CPU too.)
- **Ruled out, with numbers:** (1) B B^-1 not inverse: max|B B^-1 - 1| = 2.4e-15 on the device and 2.0e-15 on the CPU (right
  multiplication; 2.4e-15 for both on the left); with Hd off the device and CPU hopping are bitwise equal. (2) Real vs complex: the device has a real
  path for real elements (`el_left_d`/`el_right_d`, ek = 1); the dense run is bitwise the CPU one (B on a random matrix:
  difference 0). (3) Rank-1 loop vs ZGERU: replacing `k_upd_geru` by `la_zgeru` gives 4.085e-11 instead of 4.087e-11 (seed 0,
  6 sweeps, Nwrap = 25). Replacing the device V wraps by `Op_Wrapup/Op_Wrapdo` gives the same result to the last digit.
  Device stabilisation alone (ALF_GPU_NOUPDATE=1, which also moves the stabilisation to the host) = CPU to the last digit.
- **Depends on the BLAS:** the numbers above are from the Mac (gfortran 16, Apple Accelerate). On Linux (melmac,
  gfortran 14, reference BLAS/LAPACK) the same runs give GPU/CPU 0.88 with `Hd` and 0.99 without (SU(2), checkerboard,
  Nwrap = 25, 12 seeds, same chain; manager check 2026-11-26), i.e. `Hd` is not worse there. So the loss probably comes from how
  the library ZGEMM rounds the one dense product (not proven: a Mac build against the reference BLAS failed for other reasons); the per-factor path gives the CPU numbers on both platforms. Check again with
  cuBLAS/oneMKL/hipBLAS in step 2.
- **Fix (feature `gpu-hopd-optin`):** `Hd` is opt-in (`ALF_GPU_HOPD=1`, a debug/benchmark switch like `ALF_GPU_TL/TR/VB`; refused in
  production runs); the default is the per-factor checkerboard path, which gives the CPU numbers. Still to do: measure the
  cost on a real GPU in step 2 (one ZGEMM against one kernel launch per checkerboard group) before deciding whether an opt-in `Hd` is still worth having.
  Debug patch and scripts: side-session notes of the G2 session (`g2-diag.patch`, `run.sh`, `sweep.sh`).

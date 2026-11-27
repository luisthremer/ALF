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
| DP-GPU | Determinant phase underflow, GPU line | A | fixed on `fix/det-phase-underflow` (GPU line); not yet in every GPU branch |
| H1 | NaN guards removed by -ffast-math | A | fixed by feature `nan-guards` (G and scale checks by bit pattern) |
| H2 | Large beta: scale overflow in CGR -> NaN | B | open; since aef50743 the run stops (phase guard) instead of silent NaN |
| H3 | No measurements -> NaN observables | B | fixed by feature `b-batch-1` (start-up stop: NSweep < 1 or LOBS_ST > LOBS_EN) |
| H4 | confin with field value 0 -> NaN | B | fixed by feature `b-batch-1` (two-sided field check, types 1, 2, 4) |
| H5 | Global_tau_moves with partial Nt_sequential range | B | fixed in `cm/feat/b-batch-2` (commit "Wrapgr: …"; test 52) |
| H6 | Tempering with an odd number of replicas hangs | B | fixed by feature `b-batch-1` (start-up stop; upstream issue pending) |
| H7 | Bilayer_honeycomb potential/total energy wrong (also upstream) | B | fixed by feature `b-batch-1` on the CPU line (U/U2 per orbital half); GPU device weights: P2.4 step 1e; upstream issue pending |
| M1 | Dtau = 0 -> unbounded memory | B | fixed by feature `b-batch-1` (stop if Ltrot overflows) |
| M2 | Projector with Thtrot = 0 and Ltau = 1: out of bounds | B | fixed by feature `b-batch-1` (start-up stop for this combination; tau_p bounds unchanged) |
| M3 | Nwrap <= 0 accepted | B | fixed by feature `b-batch-1` (stop if Nwrap < 1) |
| M4 | Nt_sequential range not checked -> segfault | B | fixed by feature `b-batch-1` (range check when the range is not empty) |
| M5 | Langevin dt = 0 / Max_Force = 0 -> zeros; HMC/Amplitude edge cases | B | fixed in `cm/feat/b-batch-2` (commit "QMC_runtime_var: …"; start-up stop, test 42) |
| M6 | Precision Green up to 10 passes silently | B | fixed in `cm/feat/b-batch-2` (commit "control: …"; warning above 1e-3, test 51) |
| M7 | Phase update overflow for extreme R | C | open (from code reading) |
| M8 | Code-modularity not bitwise identical to master (1e-14); global-tau threshold changed | B | open (to investigate with Fakher's history) |
| L1 | master-only error-format bug | C | fixed on Code-modularity (Fakher) |
| L2 | confin from another L accepted silently | C | fixed in `cm/feat/p28-fuzz` (Z6) |
| L3 | Default-integer overflow at very large Ndim | C | open; extended by the P2.8.1 scan (Z9) |
| L4 | Restart not byte-identical for Langevin Ltau = 1 / HMC | C | open |
| L5 | Tempering messages | C | open |
| L6 | GPU line: confin 0 gives finite wrong numbers | B | open (GPU line) |
| MK | Makefile: Prog objects not rebuilt after Libraries/Modules changes | B | fixed in `cm/feat/b-batch-2` (commit "build: …"; test 50) |
| WF | Trial wave function normalisation, large N_part (suspected) | B | open, not reproduced |
| P1 | Op_V_is_symmetric slow (start-up, Symm = .T.) | B | fixed in a9fdbd64 (CPU line) |
| Z1 | NaN/Inf in a real model parameter: NaN results (rc 0) or segfault | B | fixed in `cm/feat/p28-fuzz` |
| Z2 | NaN/Inf in a real `&VAR_QMC` parameter (Amplitude = Inf freezes the chain) | B | fixed in `cm/feat/p28-fuzz` |
| Z3 | NBin < 1 or a measuring interval outside 1..Ltrot: nothing written / NaN, rc 0 | B | fixed in `cm/feat/p28-fuzz` |
| Z4 | NaN field values pass Fields_test, then segfault | B | fixed in `cm/feat/p28-fuzz` |
| Z5 | confin_0.h5 read without error checks: heap overflow, undefined seed | B | fixed in `cm/feat/p28-fuzz` |
| Z6 | Text confin with more values than fields accepted (= L2) | C | fixed in `cm/feat/p28-fuzz` |
| Z7 | data.h5 that HDF5 cannot open: every bin lost with rc 0, or abort at the end | B | fixed in `cm/feat/p28-fuzz` |
| Z8 | "Switch to LOG" message for absurd but finite inputs (U = 1e5, field 1e3) | C | open |
| Z9 | Size limits: Ndim^2 in default integers, O(N^2) lattice tables | C | open (part of L3) |
| Z10 | Namelist errors stop with the bare gfortran message ("End of file") | C | open |

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

## Input fuzzing and size probes 2026-11-27 (P2.8.3, P2.8.1)

Harness `testsuite/fuzz/fuzz_inputs.py` (helper `h5mut.f90` for HDF5 files): each case starts from
`Scripts_and_Parameters_files/Start` (4x4 Hubbard, Beta = 1, 3 sweeps, 2 bins, Ltau = 0), changes one input
(`parameters`, `seeds`, `confin_0`, `confin_0.h5`, `data.h5`) and runs `ALF.out` with 60 s and 2 GB limits. Outcome
classes: clean stop (ALF message), runtime stop (gfortran message with file and line), memory stop (gfortran
"Error allocating N bytes"), crash, hang, NaN (rc 0 with NaN in an observable or the precision lines), frozen
(acceptance 0), runs. 102 text cases and 14 HDF5 cases (GNU, noMPI, melmac): on `Code-modularity-CPU` 9dc3037a
13 text and 6 HDF5 cases fail, with the fixes none. CTest: `54-fuzz-*` to
`61-fuzz-*` (serial builds), unit cases `42-...-lobs_window`, `42-...-nbin`, `43-...-fieldnan`. The
claude-bughunter plugin was not available in the session that ran this (triage by hand, same classes).

- **Z1. NaN/Inf in a real model parameter** (B). `ham_U = NaN`: no interaction operator is built (`abs(U) > eps` is
  false), return code 0, `Pot_scal`/`Ener_scal` NaN, every other observable that of the free system. `Phi_X = NaN`:
  rc 0, `Kin_scal`/`Ener_scal` NaN, Green and spin correlations finite. `ham_T = NaN`: segfault in
  `test_checkerboard_decomposition` (Predefined_Hop_mod). The namelist reader accepts `NaN` and `Inf`. **Fixed:** the
  generated `read_parameters` (`Prog/parse_ham_mod.py`) calls `Check_finite_parameter` (runtime_error_mod, bit-pattern
  test) for every real parameter of every model; `Prog/Makefile` regenerates the readers when the generator changes.
- **Z2. NaN/Inf in a real `&VAR_QMC` parameter** (B). `Amplitude = Inf` (continuous fields): acceptance 0, the chain
  never moves, finite results, rc 0; `Amplitude = NaN` stopped with the wrong message "Amplitude = 0"; `CPU_MAX = NaN`
  accepted; Langevin with `Max_Force = Inf` ran without a force cap. **Fixed:** `read_and_broadcast_QMC_var_and_ham_name` checks CPU_MAX, Amplitude and the Langevin/HMC/MALA
  step sizes and forces.
- **Z3. No bins or no measured time slice** (B, H3 class). `NBin = 0` or `< 0` (without CPU_MAX): rc 0, no bin, no
  output file. `LOBS_ST = -5, LOBS_EN = -1` (also `LOBS_ST > Ltrot`): the window is not empty, so the H3 check passes,
  but no slice is measured: NaN in every equal-time observable, rc 0. **Fixed** in
  `set_default_values_measuring_interval`: stop if the window has no slice in 1..Ltrot or NBin < 1 without CPU_MAX.
- **Z4. NaN field values pass Fields_test** (B). All type tests compare with NaN, which is false, so a NaN field
  (from `confin_0.h5`) passes and the run segfaults in `Op_mmultR` (index from the field value). **Fixed:**
  Fields_test stops on any NaN or Inf field (bit-pattern test).
- **Z5. confin_0.h5 is read without any error check** (B). A file of a larger lattice (L = 6 into L = 4): the
  dataset is read into the smaller buffer, heap corruption, SIGABRT. Missing or short `seed` dataset: HDF5 error
  messages, then the run goes on with an undefined seed, rc 0. An unreadable file (empty, truncated, garbage) left
  seed and fields undefined; only the type test of the discrete fields stopped it (continuous fields would have run on
  garbage). **Fixed:** `Fields_read_conf_h5` checks the open, the existence, rank and shape of both datasets and
  every read.
- **Z6. Text confin with more values than fields** (C, = L2). A confin of L = 6 (or of a larger Beta) in an L = 4
  run: the first values are used, rc 0. **Fixed:** `Fields_read_conf` stops if a non-blank line follows the last
  field. Fewer values already stopped (gfortran "End of file").
- **Z7. data.h5 that HDF5 cannot open** (B). A garbage `data.h5`: every bin's `h5fopen` fails unchecked, 600 KB of
  HDF5 messages, rc 0, no bin saved. A truncated one: the run aborts (SIGABRT) after the last sweep, results lost.
  **Fixed:** `Init_hdf5_file` stops at start-up if an existing data file cannot be opened (or a new one cannot be
  created); Print_bin stops if the open fails during the run.
- **Z8. "Switch to LOG" for absurd but finite inputs** (C, open). `ham_U = 1e5`, `ham_T = 1e200`, `ham_chem = 1e5`,
  a continuous field of 1e3: the run stops with "At least one scale is NaN or Inf. Switch to LOG is required.", which
  points at the build instead of the input.
- **Z9. Size limits** (C, open; P2.8.1, adds to L3). Static scan for default-integer products of sizes:
  `Control_All_Finite(A, Ndim*Ndim)` (control_mod, the nan-guards check) overflows for Ndim >= 46341, the loop is then
  skipped and the guard passes silently; `MPI_Sendrecv_UDV_state` (`n = ndim*ndim`, tempering) and the L3 entries
  likewise. Runs: L = 200 stops in `lattices_v3_mod.F90:350` allocating 6.4 GB for an integer table of Latt%N^2
  entries, before G (25.6 GB) is allocated: the lattice tables are O(N^2) memory. Ltrot = 1e8 / 2e9 (Beta = 1e7 / 2e8)
  stop cleanly with gfortran's "Error allocating ... bytes" (25.6 GB / 512 GB for the fields). `N_SUN = 200`, Ltrot = 1,
  the 1-D chain and L = 1 behave. None of the overflows is reachable at today's sizes (G alone would need >= 32 GB).
- **Z10. Namelist errors** (C, open). A missing namelist, an unknown key or a wrong type stops with the gfortran
  runtime error and the file and line of the READ ("End of file", "Cannot match namelist object name l3"), rc 2.
  Clear enough for a wrong key; "End of file" for a missing `&VAR_lattice` does not name the namelist.
- **Checked and fine:** missing/empty parameters file, CRLF line ends, tabs, a UTF-8 BOM, duplicate keys (last wins),
  Beta/Dtau/Theta NaN, Inf, < 0, 1e-300 (Ltrot/Thtrot checks), Ltau, N_SUN, N_FL, L1/L2 out of range, unknown
  ham_name or lattice, seeds missing/empty/garbage/real/huge/binary (gfortran stops), seeds 0 and negative (run),
  confin empty/truncated/garbage/binary/too short/field 0, 3, 0.5, NaN as text, both confin_0 and confin_0.h5,
  continuous fields NaN/Inf (stopped, see Z8 for the message), `CPU_MAX < 0`, huge Nwrap, N_Phi, NSweep (runs).
- **Not checked:** MPI, tempering and PARALLEL_PARAMS input handling (several ranks reading partly corrupted files),
  models other than Hubbard, the analysis programs, the Intel compiler.

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
  **Fixed in a9fdbd64** (`cm/feat/opv-symm-fast`, CG-ALF T3): all-diagonal operator sets return at once, otherwise each
  factor is multiplied in on its own columns, O(n_op Ndim N^2). Check alone, old -> new (melmac, GNU, 1 thread, loaded
  node): Hubbard L = 24/32: 142 s / 1341 s -> < 0.1 ms; tV L = 24/32/48: 247 s / 2330 s / ~17 h (extrapolated) ->
  0.02 / 0.08 / 0.44 s. Same result as the old check on 100 testsuite scenarios and in unit test 53.

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
- **To do:** merge the fix into every GPU branch before it is ported (P2.4); a GPU branch that has
  not merged it must not be used for N >~ 2950.

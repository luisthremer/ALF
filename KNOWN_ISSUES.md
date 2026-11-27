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
| L2 | confin from another L accepted silently | C | fixed in `cm/feat/p28-fuzz` (Z6) |
| L3 | Default-integer overflow at very large Ndim | C | open; extended by the P2.8.1 scan (Z9) |
| L4 | Restart not byte-identical for Langevin Ltau = 1 / HMC | C | open |
| L5 | Tempering messages | C | open |
| L6 | GPU line: confin 0 gives finite wrong numbers | B | open (GPU line) |
| G1 | GPU line: `ALF_GPU_SELFTEST` leaves test data in the device copy of G marked valid | B | fixed by feature `gpu-wrap-update` (step 1c; harmless in 1b, where the sweep did not use the mirror) |
| G2 | GPU line: wrapped G (device wraps and updates, step 1c) less accurate than on the CPU: Precision Green mean 1.4 - 2.5x the CPU value, same chain | B | fixed by feature `gpu-hopd-optin` (cause: dense hopping product `Hd`; now opt-in, `ALF_GPU_HOPD=1`) |
| G3 | GPU line: the device precision check (`Control_PrecisionG_vals`, step 1d) kept the fixed abort threshold 10 and ignored `ALF_UNSTABLE_THRESHOLD` (CPU-line feature `unstable-threshold`, merged in parallel with 1e) | B | fixed by feature `gpu-hopd-optin` (uses `G_thresh` like `Control_PrecisionG`) |
| G4 | `stp/production-mpi-intel-*` tags (Intel production on gpu-port): raw-product LU determinant phase at every CGR | C | open, tags only (production L <= 44); not on the CM line |
| G5 | GPU=nvidia (nvfortran 25.7): a complex(8) host scalar read in a serial `!$omp target` region loses its imaginary part; hits `k_cgr_phase` (`ph0`) in the LU phase mode | B | fixed in feature `gpu-nvidia` (`map(to: ph0)` in `k_cgr_phase`; gpu-port and the production branches are not affected) |
| G6 | GPU line: negative control `ALF_GPU_FAULT=upd` acts only on the standard (rank-1) update, not on delayed/submatrix | C | open (test gap; found by the P2.4 step 2 `gpu-intel` gate) |
| MK | Makefile: Prog objects not rebuilt after Libraries/Modules changes | B | fixed in `cm/feat/b-batch-2` (commit "build: …"; test 50) |
| WF | Trial wave function normalisation, large N_part (suspected) | B | open, not reproduced |
| P1 | Op_V_is_symmetric slow (start-up, Symm = .T.) | B | fixed in a9fdbd64 (CPU line) |
| R1 | GPU=amd: rocSOLVER's blocked `zgeqrf` (ROCm core-7.14) returns NaN with info = 0 when a 64-column block has tau = 0 (identity-like input, n >= ~200); the own HIP QR covers only m <= 2304 | B | fixed by feature `gpu-amd` (step 2: NaN test after every hipSOLVER `zgeqrf`, recomputed with the unblocked `rocsolver_zgeqr2`; test `55-linalg-hip` identity n = 200, 2400); rocSOLVER report pending |
| R2 | GPU=amd (gfortran 14 amdgcn offload): CM's OpenMP target kernels give wrong results in full runs, depending on layout and timing (OpenMP permutation kernels of gpu_stab_mod at 6x6; device `la_zgetrf` in the CGR of the LU check mode: wrong tau observables (Precision tau 2e-3, silent), NaN G, libgomp mapping fatal, GCN queue error); correct in isolation | A (GPU=amd only; CUDA and GPU=host not affected) | open, investigating, root cause unknown. Mitigation in feature `gpu-amd`: HIP replacements mandatory on GPU=amd (`ALF_HIP_KERNELS=0` debug only), CGR-phase LU of the debug modes on host LAPACK, broad default-path verification |
| R3 | GPU=amd: the HIP slice kernel `hk_vslice_r` (feature `gpu-amd`) read `V_meta` with a fixed stride 4; feature `p26-models` (merged d1e0dce6) widened `V_meta` to (8, Nop, N_FL), so every GPU=amd run of the GPU line from d1e0dce6 on used wrong vertex factors (Hubbard 4x4 Mz: stop "CGR: la_zqrpp (-101: NaN or Inf in the matrix)") | A (GPU=amd only) | fixed by feature `gpu-amd-r3` (stride `size(V_meta,1)` passed to the kernel) |
| G7 | GPU line: the device stabilisation (column-norm pre-pivoted QR, all device paths and GPU=host) moves Mz SpinXY/SpinT equal-time correlations at L >= 24 by ~1e-8 relative against the CPU path (ZGEQP3); CPU with -DQR_PREPIVOT 4.3e-9, GPU=host 6.6e-9, CUDA/AMD/gpu-port of the same order; not a backend bug. The e2e tolerance (100 x Precision Green mean) is too tight there: it needs size scaling for Mz | C | open (property of the algorithm; tolerance rule to be scaled in step 3). Interim rule since feature `p26-g10`: `gpu_e2e.py` accepts max(1e-11, 100 x Precision Green mean, Precision Green max) of the CPU run (only the CPU run's max: with the GPU run's max G10 would have passed). Case `sp-standard` (Spin_Peierls 12x12) has an absolute floor of 5e-10: CI runners give Part_scal 2.2e-10 where melmac gives 2e-11 (run 46033748753, tolerance 1.5e-10; `confout` identical): round-off at the tolerance edge, not a defect |
| G8 | GPU line: the stabilisation self-test `ALF_GPU_SELFTEST=3` reports a CGR phase difference of order 1 (1.8 - 2.0) between device and host for N_SUN = 2 models (Hubbard SU(2) `Mz = .F.`, Kondo, Spin_Peierls) and stops; G agrees to 1e-14, and the production runs agree (phase precision 1e-11, sign columns equal). Found in P2.6.1; reproduced on the unchanged base 5c0af0a4 (Hubbard SU(2)) | C | open (debug self-test only; probably the self-test compares a phase of different convention for N_SUN > 1, not investigated) |
| G9 | GPU line, H100 (GPU=nvidia): Kondo 8x8 bilayer (Ndim 128), every scheme: Dimer_tau of the GPU path differs from the CPU path by 1.26e-11, above the floor 1e-11 of the e2e rule (confout identical, Precision Green mean 7e-14); GPU=host and the 4x4 cases are within it. Round-off of the device tau pass, same family as G7. GPU=amd (MI210, feature `gpu-amd-p26`): the same Kondo 8x8 cases 1.48e-11; tV 12x12 delayed 16: Pot_scal 1.67e-11 against tol 1.66e-11, while the CPU-only null (same binary, submatrix 16 vs delayed 16, `confout` identical) gives 1.88e-11 (round-off at the tolerance edge) | C | open (tolerance rule for tau files to be scaled with size, together with G7) |
| G10 | GPU line common code (P2.6.1, `p26-models`): Spin_Peierls 12x12 (type-4 J vertices), standard update: Precision Green of the device path ~30x the CPU value (GPU=host 2.6e-10 vs 9.0e-12, MI210 likewise), same chain, observables drift up to 1.8e-7. Cause: the device rank-r accept used the Woodbury form G - G(:,P) W R (W = X Tri^T) instead of Upgrade2's grouping (G(:,P) X)(Tri^T R). The H100 evidence (big-sp-std OK) came from nvfortran without -ffast-math; the same tree fails on GPU=host | B | fixed by feature `p26-g10` for the common OpenMP kernels (GPU=host, CUDA, oneMKL; device / CPU 1.10, observables 1.8e-9; test `gpu-state-e2e` case `sp12-precision`); GPU=amd: the HIP accept `alfhk_rk_accept` still takes W, open until it follows |
| G11 | CPU delayed / submatrix updates (`lowrank_update_mod`, Woodbury form with K = Delta Mat^-1): Spin_Peierls 12x12, Precision Green 2.5-2.8e-10 against 9.0e-12 of the standard update (same model, same seeds); the device delayed path uses the same form and gives the same. Same mechanism as G10 | B | open (fix: Upgrade2's grouping in `Lowrank_accept` and in the device `k_gen_accept`; CPU line) |
| CI1 | CI: `41-check-phase-fail` failed once on `bookworm-pgi-24-09 · MPI` (run 45968198665, P2.6.1 canvas): the expected "phase| /= 1: (NaN" line was not in the output before MPI_ABORT (output not flushed); the rerun of the same commit passed, the code path (control_mod) is unchanged | C | open (flaky: flush stdout before the abort, or match the error line instead) |
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
  stop cleanly with gfortran's "Error allocating ... bytes" (25.6 GB / 512 GB for the fields). Dtau = 1e-300 stops with
  the M1 message, except in builds with `-ffpe-trap=invalid` (CI "MALA ... (debug)"): there `Ltrot = nint(Beta/Dtau)`
  (in every model) traps with SIGFPE before the check. `N_SUN = 200`, Ltrot = 1,
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

## G5. GPU=nvidia: complex scalar loses its imaginary part in a serial target region (nvfortran 25.7)

- **Applies to:** `Code-modularity-GPU` built with `GPU=nvidia` (nvfortran 25.7, `-mp=gpu -gpu=cc90,cuda12.9`, Helma H100),
  P2.4 step 2 (`gpu-nvidia`). GPU=none and GPU=host are not affected (no offload).
- **Symptom:** `gpu-state-e2e` on the H100, SU(2) cases (`su2-standard`, `su2-delayed-N4`, `su2-nwrap25`): the run in LU phase
  mode (`ALF_GPU_DEBUG=1 ALF_GPU_CGR_LU=1`) has `Precision Phase, Max` 1.4e-6 (tracked mode: 1.5e-13) and observables 1e-7 .. 1.6e-6
  away from the tracked run and from the CPU path (`ALF_GPU_DISABLE=1`), with the same configurations. The check mode
  (`ALF_GPU_CGR_LU=2`) reports LU and tracked phase equal to 1.4e-15, so the LU phase itself is right.
- **Cause (compiler):** in a serial `!$omp target` region (no `teams`/`parallel do`), an implicitly firstprivate `complex(8)` scalar
  arrives on the device with imaginary part 0. Measured with a 40-line program at -O2, -O3, -O3 -fast, -Mfprelaxed (same result):
  `z = (0.6, 0.8)` read inside `!$omp target` gives `(0.6, 0)` (implicit firstprivate of a dummy, of a local copy, and explicit
  `firstprivate`); `map(to: z)`, a `target teams distribute parallel do` loop and passing `real(z)`, `aimag(z)` as two reals give
  `(0.6, 0.8)`. At -O0 the dummy case stops with CUDA_ERROR_ILLEGAL_ADDRESS. In ALF, `k_cgr_phase` (`Prog/gpu_stab_mod.F90`)
  reads the host phase `ph0` in such a region; it is used only when the phase does not come from the tracked det(U) on the device
  (`dev /= 1`: LU mode, i.e. `ALF_GPU_CGR_LU=1`, or the default after a `git revert` of the S1 commit). The SU(2) phases are nearly
  real, hence the small error; a model with a genuinely complex phase would get a wrong phase. A scan of all serial target regions
  of `gpu_*_mod` and `linalg_cuda_mod` found no other complex scalar read from the host (the others are assigned inside the region).
- **Default path:** not affected (tracked mode reads `cs(7)` on the device; e2e active vs CPU path within 1e-14 for SU(2)).
- **Class:** B (manager, 2026-11-26): only the LU check mode and the S1-revert path are hit, and no production run uses the CM CUDA
  code. gpu-port and its production branches/tags (`gpu-delayed-submatrix-defect-production-*`, `stp/production-mpi-intel-*`,
  `gpu-port-unified`) start `k_cgr_phase` from `cs(3)` on the device and have no such read (scan of all serial target regions).
- **Fix:** first commit of feature `gpu-nvidia`: `!$omp target map(to: ph0)` in `k_cgr_phase` (no-op on GPU=host). Regression test:
  the LU-mode lanes of `gpu-state-e2e` run on the H100 (verified: SU(2) tracked vs LU 8.6e-15 after the fix, job 1465022;
  `docs/ai/phase2/nvidia/RESULTS.md`). When editing device code: pass complex host scalars into a serial target
  region with `map(to:)` or as two reals.

## R1. GPU=amd: rocSOLVER blocked zgeqrf gives NaN on tau = 0 blocks (2026-11-26; fixed by feature `gpu-amd`)

- **Applies to:** GPU=amd, `la_zgeqrf` / `la_zqrpp` on device data whenever hipSOLVER is used: m > 2304 (the own HIP QR
  of `alf_hip_qr.hip.cpp` stops there, so every QR at L >= 49 on the square lattice) or `ALF_HIP_QR=0`. Also gpu-port
  (same hipSOLVER call, no guard). ROCm core-7.14, rocSOLVER behind hipSOLVER, MI210.
- **Symptom:** `hipsolverZgeqrf` returns info = 0, but from column 65 on the factors are NaN (identity n = 200: 9180 NaN below,
  17884 above and 136 on the diagonal; n = 500 likewise; n <= 128 and random matrices are correct; `zunmqr` with tau = 0
  reflectors is correct). Not silent on the GPU path: the NaN reaches the scales (`k_testscale` stops the run) or the tracked
  det(U) of the upload (untracked slot, stop in CGR). Where it can occur: QR of identity-like U (the start-up upload of the
  UDV stacks, `gpu_stab_put_udv` / `k_det_qr`; reset slots) at large N.
- **Cause (probable, not traced into rocSOLVER):** the blocked path (switch at n = 128, block 64): the trailing update with
  the T factor of a block whose reflectors all have tau = 0. The unblocked `rocsolver_zgeqr2` is correct.
- **Fix:** `linalg_hip_mod` keeps a device copy of A before every hipSOLVER `zgeqrf`, tests the result for NaN/Inf with the
  column-norm kernel (n values D->H) and, on a hit, restores A and recomputes with `rocsolver_zgeqr2` (one line in stdout at
  the first occurrence). Cost on the hipSOLVER path: one N^2 device copy and one norm kernel per QR.
- **Test:** `testsuite/GPU.tests/55-linalg-hip.F90`: identity n = 200 and n = 2400 against the host `zgeqrf` (exact),
  run with the own QR and with `ALF_HIP_QR=0`.

## R2. GPU=amd: OpenMP target kernels wrong in full runs (2026-11-26; open, class A for GPU=amd)

- **Applies to:** GPU=amd builds of the GPU line (gfortran 14.2 with amdgcn offload, ROCm core-7.14, MI210). Not CUDA, not
  GPU=host (same OpenMP code passes there).
- **Symptoms (all in runs that end or abort; one silent):** (1) with the HIP replacements switched off (`ALF_HIP_KERNELS=0`),
  6x6 Start input, NSweep 5 x 2 bins: "Green function contains NaN or Inf"; bisected to the OpenMP permutation kernels of
  `gpu_stab_mod` (`k_perm_rows/cols`, `k_iperm_rows/cols`); every other kernel group passes on OpenMP. (2) with
  `ALF_GPU_CGR_LU=1` or `2` (debug: the CGR phase from a device `la_zgetrf`), even with all HIP kernels: wrong time-displaced
  observables with the run ending normally (Precision tau 2e-3 instead of 3e-15, SU(2) 4x4; silent), or unstable/NaN G, or the
  libgomp fatal "Trying to map into device [stack, 16 B) ... already mapped" at the complex `reduction(+:Zkin)` of
  `gpu_obs_mod`, or "GCN fatal error: Asynchronous queue error". The symptom changes with the build (-O3, -O0 -fcheck=all, -g),
  with extra prints, and disappears when an extra target region or a D->H update is inserted before the kernel.
- **Ruled out:** the deferred HIP synchronisation (fails with `ALF_HIP_ASYNC=0`), gpu-port's patched libgomp plugin (fails with
  the stock plugin), SDMA copies (`HSA_ENABLE_SDMA=0`), hipSOLVER out-of-bounds writes (guard-zone test), a stale device copy
  of the pivots (checked after every `la_zqrpp`), the kernels themselves (OpenMP and HIP versions agree on the same data,
  `56-hipk-stab`; device `la_zgetrf` agrees with host LAPACK, `55-linalg-hip`), simple reduction-variable leaks.
- **Mitigation (feature `gpu-amd`):** the HIP replacements are mandatory on GPU=amd (`ALF_HIP_KERNELS=0` only with
  `ALF_GPU_DEBUG=1`, warned in the info file); the CGR-phase LU of the debug modes runs on host LAPACK on GPU=amd; default-path
  verification GPU=amd vs `ALF_GPU_DISABLE=1` over sizes, schemes and builds (see the feature's report). GPU=amd is
  "experimental: R2 open" until the cause is known.
- **P2.6 models (feature `gpu-amd-p26`, 2026-11-27):** a further symptom (3): the libgomp fatal "Trying to map into device
  [a..a) ... already mapped" in the delayed update of Kondo and Spin_Peierls: `k_dly_gpp` / `k_dly_wrap` of `gpu_lowrank_mod`
  declare `X(ld,jcnt)`, `Y(ld,jcnt)` and map them implicitly; with `jcnt = 0` (first update of a block) the zero-length map
  collides with the mapping already present. Mitigation: HIP versions of both kernels with assumed-size dummies (no map).
- **Next:** reproducer below the full code; GCC 15 offload; upstream report.

## R3. GPU=amd: HIP slice kernel read the old V_meta layout (2026-11-27; fixed by feature `gpu-amd-r3`)

- **Applies to:** GPU=amd builds of `Code-modularity-GPU` from d1e0dce6 (merge of `p26-models`) until the fix; every model and
  scheme (the Op_V slice `k_vslice_R` runs in every wrap). CUDA, oneMKL, GPU=host: not affected (they index `V_meta` by name).
- **Cause:** `alfhk_vslice_r` computed the address of `V_meta(:,n,nf)` with a hard-coded leading extent 4; `p26-models`
  added entries 5..8 (N_non_zero, diag flag, uoff, moff). The kernel read entries of other operators.
  Root cause in the process: P2.6.1 changed the shared device data layout of the `gpu_state_mod` mirror without the vendor
  blocks being checked (CI cannot build GPU=amd). From now on the manager tells the vendor sessions (gpu-amd, gpu-intel)
  before merging any change to the mirror layout.
- **Symptom:** Hubbard 4x4 Mz (standard, delayed, submatrix) stops with "CGR: la_zqrpp (-101: NaN or Inf in the matrix)".
  Wrong but finite factors are possible on other inputs, so silent wrong results cannot be ruled out for runs made in that window.
- **Fix:** the leading extent `size(V_meta,1)` is passed to the kernel (feature `gpu-amd-r3`, one commit). Check on the
  MI210, GPU vs `ALF_GPU_DISABLE=1`, same seeds: Hubbard 4x4 Mz standard / delayed / submatrix (`cases_p26.py` hub-*, 2.7e-10,
  before the fix all three stopped) and three cases of the reduced gate of `gpu-amd-hop` (L = 12 Mz Ltau 1 delayed 6.6e-10,
  L = 12 SU(2) submatrix 6.4e-12, L = 16 Mz submatrix 1.6e-10; all PASS); the P2.6 gate of `gpu-amd-p26`.

## G10. GPU line: Precision Green of the rank-r device update (2026-11-27; fixed by feature `p26-g10`)

- **Applies to:** device path of `Code-modularity-GPU` from `p26-models` (d1e0dce6) on, standard update of vertices with rank
  r = N_non_zero > 1 or N > 1 (Kondo, tV, Spin_Peierls); all backends that compile the common kernels. Visible at Spin_Peierls 12x12
  (type-4 J vertices, Beta 2, U 1.5); 6x6 and the Kondo/tV cases are within the e2e rule.
- **Symptom (GPU=host, melmac, gfortran -O3 -ffast-math, same seeds, `confout` identical):** Precision Green mean / max of the GPU
  path 2.6e-10 / 2.0e-7, CPU path 9.0e-12 / 3.3e-9; Part_scal differs by 1.8e-7, Dimer_eq 4.3e-8. MI210 (gpu-amd) the same.
- **Bisection:** the tree of the H100 evidence (cce147cf) gives the same numbers bit for bit on GPU=host, so the H100 pass is the
  platform (nvfortran, no -ffast-math), not the tree. CPU build with -DQR_PREPIVOT: 8.5e-12 (the pre-pivoted QR is not it).
  GPU=host with STABLOG (host ZGEQP3 in CGR): 2.6e-10 (the device stabilisation is not it). `ALF_GPU_NOUPDATE=1`: CPU numbers.
  Debug patch on GPU=host (G_d is host memory): `Op_Wrapup/Op_Wrapdo` in place of the device N > 1 wraps: bit for bit unchanged
  (the wraps and the type-4 factors are exact copies); Upgrade2's accept arithmetic in place of `k_rk_prep/k_rk_upd`: 1.0e-11.
  W from LU (K = Delta Mat^-1, as the CPU delayed path) with the device kernels: 3.7e-10, i.e. the form, not W, loses the digits.
- **Cause:** the device accept applied G <- G - G(:,P) (W R) with W = X Tri^T (r x r, host) and R = E_P^T - G(P,:); Upgrade2
  applies G <- G - (G(:,P) X) (R^T Tri)^T. The two are equal in exact arithmetic; for the J vertices the Woodbury grouping loses
  ~30x. The CPU delayed/submatrix updates use the same Woodbury form and show the same loss (G11).
- **Fix (feature `p26-g10`, `gpu_update_mod`):** the host sends X and Tri (`rk_xtri`, Upgrade2's recursion, unchanged), `k_rk_prep`
  forms C = G(:,P) X and T = Tri^T R, `k_rk_upd` unchanged. Spin_Peierls 12x12 standard: 9.9e-12 vs 9.0e-12 (1.10), observables
  1.8e-9; submatrix 16 (flushes before these vertices, so it uses the standard accept): 1.1e-11; delayed 16, 32 unchanged bit for bit
  (device delayed general path keeps W, like the CPU delayed path: G11). Regression: `gpu_e2e.py` case `sp12-precision` (ratio <= 2;
  before the fix 29.7). GPU=amd: the HIP `alfhk_rk_accept` (merged with `gpu-amd-p26`, f45fdcc0) computes T = W R; the HIP branch keeps receiving W
  (unchanged on GPU=amd, so G10 stays there until the kernel takes X and Tri; not run on the MI210 here).


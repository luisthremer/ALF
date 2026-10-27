# 203 – `gpu_state`: flat device mirror, setup flattening, supported-feature check

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev (review: numerics-reviewer)
- **Depends on:** 201 (done), ADR-002 §2 (accepted). 202 (done; `alf_la` not needed here).
- **Status:** open

## Goal
With `GPU=nvidia`, ALF builds at setup a flat, device-resident copy of everything the sweep needs that does not
change during the run (hopping exponentials in product order, Op_V exponential tables, field types) plus the
device buffers that later tickets fill (G, nsigma, udvst, scratch). A setup check decides whether the run can use
the GPU path and writes the decision and its reason to `info`. A self-test proves that the flat arrays encode
exactly the operators ALF applies. Nothing in the sweep uses the mirror yet; results stay identical.
With `GPU=none` nothing changes (no object code, identical `ALF.out`).

## Context
- ADR-002 §2 (table of arrays), §3 (transfers); HARD.md H4 (flat mirror), H5 (checkerboard product order must
  match the CPU exactly), H10 (guarded early branches only). PLAN §1 transfer table.
- **Setup order** in `Prog/main.F90`: `ham%Ham_set` (≈ 318), Op_V checks incl. time-dependent `g_t` (≈ 321–330),
  `nsigma%make` + types (≈ 394–397), initial field (≈ 404–409), `Hop_mod_init` (≈ 411), `ham%Alloc_obs` (≈ 444),
  `udvst` allocation (≈ 590–603). The hook goes after the `udvst` setup, as one guarded call
  (`#ifdef ALF_GPU` … `call gpu_state_setup(...)`), and one guarded `gpu_state_finalize` before the end.
- **Hopping:** `Prog/Hop_mod.F90` keeps `ExpOpT_vec(N_FL)` (private `DynamicMatrixArray` of
  `class(ContainerElementBase)`), elements `RealExpOpT`/`CmplxExpOpT` (`Prog/OpTTypes_mod.F90`: `mat`, `invmat`,
  `mat_1D2`, `invmat_1D2`, `P(:)`, `Ndim_hop`) and time-dependent variants (`OpT_time_dependent.F90`).
  Application order: `Hop_mod_mmthr/_m1/mmthl/mmthlc/mmthl_m1/Symm` (Hop_mod.F90:143–299). Checkerboard = many
  small elements applied in sequence (families from `Predefined_Hop_mod`, symmetric decomposition
  `Symmetrize_Families` ≈ 1422); dense = one element with full `P`. Flatten **in exactly the order each routine
  applies them**: per flavour an element list with offsets into an index array (`P`) and a value array
  (all four matrices; real stored as real, complex as complex; record which). You may add a small **public
  read-only accessor** to `Hop_mod`/`OpTTypes_mod` (guarded by `#ifdef ALF_GPU`) instead of making
  components public.
- **Op_V:** `Type(Operator)` (`Prog/Operator_mod.F90:56–90`): `P(:)`, `O`, `U`, `E`, `g`, `alpha`, `type`,
  private `M_exp/E_exp`, `diag`, optional `g_t`. MVP subset (H4): `diag` operators, field types 1 and 2 (Ising,
  4-valued), constant coupling (no `g_t`). Build per (n_op, nf) the site list `P` and the exponential table
  `exp(g·φ(s)·E)` for every field value s of that type, computed **the same way** Operator_mod does when it
  applies the operator (find the routine – `Op_mmultR/L`, `Op_Wrapup/do` or similar – and reuse its formula and
  any private table via a guarded accessor so the numbers are bit-identical).
- **Fields:** `nsigma%f(N_op, Ltrot)` (complex), `nsigma%t(N_op)` (`Fields_mod.F90:79`). Device copy:
  `nsigma%f` as is plus `t`. Kept in sync per update later (205); in 203 copy once at setup.
- **Buffers to allocate (not fill):** `G(Ndim,Ndim,N_FL)`, `udvst_d` sized for `NSTM × N_FL_eff` (U, D/logs, V –
  look at `udv_state_mod.F90` for the components of default and STABLOG), τ matrices `GT0/G0T/GTT` only if
  `Ltau = 1`, scratch from a size function. Allocate with `!$omp target enter data map(alloc:)` once.
- **Supported check** `gpu_supported(reason)`: N_FL/N_SUN any; Op_V all `diag` with types 1/2 and no `g_t`;
  hopping elements real or complex, not time-dependent; no Langevin/HMC, no global moves, no tempering/MPI
  (`MPI` builds may run but the GPU path is off), projector code off. Anything else → GPU path inactive, reason
  printed, run proceeds on the CPU unchanged. Environment `ALF_GPU_DISABLE=1` forces inactive.
- **Prog/Makefile:** `GPU=none` must not compile the new module. 201 exports `ALF_GPU=nvidia|intel` only for GPU
  targets; add the source to `SRCS` only when `ALF_GPU` is set (check that `gen_deps.py` handles it).

- Device code rules from 202 / ADR-002 §1: synchronous target regions only (no `nowait`/`depend`), single host
  thread; `alf_la` (once merged) is ordered with them on the same stream.

## Files allowed
- `Prog/gpu_state_mod.F90` (new), `Prog/Makefile` (conditional source), `Prog/main.F90` (two guarded calls + one
  guarded `use`), `Prog/Hop_mod.F90`, `Prog/OpTTypes_mod.F90`, `Prog/Operator_mod.F90` (guarded read-only
  accessors only)
- `docs/gpu/scripts/nividia/t203_*.sh`; `testsuite/gpu_reference/configs.py` only to **add** configs (do not change
  existing ones)
Anything else: stop and report.

## Acceptance
1. **Self-test** (`ALF_GPU_SELFTEST=1`, `GPU=nvidia`, nividia A30) right after setup, then the run continues:
   - hopping: for a random complex N×N matrix, applying B_T, B_T⁻¹ and the half steps from the left and right
     **on the host from the flat arrays** equals `Hop_mod_mmthr/_m1/mmthl/mmthlc/mmthl_m1` (and `Symm` if used):
     max rel. diff reported; expected 0 (same operation order); must be ≤ 1e-14;
   - Op_V: for every operator and every field value, applying the flat table to a random matrix equals ALF's own
     operator application: same bound;
   - device: every mirror array copied back (`target update from`) equals the host-built array bit for bit;
   - any failure terminates with a message (nonzero exit).
   Configs: `hub6x6_cb`, a dense-hopping 6×6 (add to `configs.py` if missing), `hub12x12_cb`, each default and
   `LOG` builds; Mz true and false if both exist for Hubbard.
2. **Unchanged physics with GPU build:** the same three configs with `GPU=nvidia` (mirror built, unused) vs
   `GPU=none` PGI build, same seeds, 20 sweeps: all `data.h5` datasets bit-identical, `info` precision lines equal.
3. **Unsupported model:** one predefined Hamiltonian with a non-diagonal Op_V or type-3 fields (name it) and one
   Hubbard run with `ALF_GPU_DISABLE=1`: `info` says GPU path inactive with the reason; results bit-identical to
   `GPU=none`.
4. **`GPU=none`:** GNU MPI and PGI noMPI `ALF.out` byte-identical to HEAD builds (or only build-metadata
   differences – show what), ctest 52/52 with GNU devel noMPI.
5. **`GPU=intel` compiles** on the `intel` machine (`IntelLLVM noMPI GPU=intel`, build only; run the 6×6 self-test
   on the B580 if it works, report either way – not a blocker).
- Report: diff summary, the flat layout (array, shape, type, bytes for N=144 and N=1024, L=40/200), self-test
  numbers per config, setup time added, unsupported-model output.

## nividia / intel rules
Test only; develop on the laptop and rsync. Work dir `~/alf-gpu/work/t203/` on each machine, `TMPDIR` inside it,
≤ 4 cores on nividia, ≤ 8 on intel after checking load, GPU only while running, runs ≤ 15 min, clean up build trees
at the end. Remote shells are zsh: run scripts with `bash`. On intel use `ONEAPI_DEVICE_SELECTOR=level_zero:0`, never
plain `sycl-ls`. Disk on nividia is tight (≈ 19 GB free): keep the work dir ≤ 4 GB.

## Out of scope
Device kernels for wraps (204), updates (205), stabilization (206); using `alf_la`; any change of results;
MPI/tempering on the GPU.

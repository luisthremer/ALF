# Two GPU branches: `gpu-port` and `gpu-port-unified`

| Branch | What it is | For whom |
|---|---|---|
| `gpu-port` (this branch) | The general GPU port of ALF: upstream ALF (`master` 521ddfca) plus opt-in GPU backends (`GPU=nvidia`, `GPU=intel`, `GPU=amd`), delayed/submatrix updates, MPI/`PARALLEL_PARAMS` on the GPU, `alf_la`, MaxEnt bootstrap, `alf_exit`, `ALF_UNSTABLE_THRESHOLD`. Hubbard model and measurement window as upstream. | Every ALF user. |
| `gpu-port-unified` | `gpu-port` plus the space-time point-defect Hubbard model (`ham_h0` = field on one time slice, observables `SpinZ_local`/`SpinZ_localtau`, `LOBS_ST` default 0, Ltau = 1 switches the field off). | The defect production runs only. |

`gpu-port` was made from `gpu-port-unified` (72a3900f) by removing the defect, so it is a
fast-forward of the old `origin/gpu-port` (43045813). With `GPU=none` it is bitwise identical
to `master` 521ddfca for the upstream features tested (Hubbard `Mz` T/F, `Symm`, `Checkerboard`,
`Ltau` 0/1, projector, static pinning field `ham_h0`, `N_SUN`, Kondo, tV, explicit `LOBS_ST/EN`),
with the stock parameter file and **default** `LOBS_ST = 0, LOBS_EN = 0` (effective start 1 as
upstream).

## What differs between the branches

| File | `gpu-port` | `gpu-port-unified` |
|---|---|---|
| `Prog/Hamiltonians/Hamiltonian_Hubbard_smod.F90` | master's file plus three `ALF_GPU` hooks (`gpu_obs_setup`, `gpu_obs_active`/`gpu_obs_measure`). `ham_h0`: static pinning field at the origin, orbital 1, all slices, observable `SpinZ` local (upstream). | `ham_h0`: field on time slice 1 (`Op_T(2,nf)%g_t(1)`) at site `Latt%N/2 - L/2`; `SpinZ_local`, `SpinZ_localtau` with end-point weights; Ltau = 1 sets `ham_h0 = 0`; refuses projector, Langevin/HMC and a non-default window with `ham_h0 /= 0`. |
| `Prog/QMC_runtime_var_mod.F90` | `LOBS_ST = 0` in the file means 1 (upstream). | `LOBS_ST = 0` means 0 (the down sweep also measures tau = 0). |
| `Prog/main.F90` | down-sweep measurement labelled `NTAU1` exactly as master. | the same label through an alias `nt_obs`. |
| `Prog/Hamiltonian_main_mod.F90` | `Obs_local` private (upstream). | `Obs_local` public (used by `gpu_obs_mod`). |
| `Prog/gpu_obs_mod.F90` | no device accumulators for local observables; with `ham_h0 /= 0` the device measurement declines and the host `Obser` measures. | `local_obs`, `lacc1/lacc2`: device accumulators of `SpinZ_local(tau)`. |
| `Libraries/Modules/alf_hipk_mod.F90`, `alf_hip_kernels.hip.cpp` | no `hk_local` kernel. | `hk_local` / `alfhk_local` (HIP version of `local_obs`). |
| `Prog/observables_mod.F90`, `Analysis/ana.F90`, `ana_hdf5.F90`, `ana_mod.F90` | master's files (local observables with `Nt = 1` only). | tau-resolved local observables (`Nt > 1`, `_localtau` analysis). |
| `testsuite/Run.tests/` | `run-schemes` on the standard Hubbard model (`ham_h0 = 0`). | `run-schemes` on the defect, plus `run-pinning` and `run-local`. |
| `docs/gpu/` | this file, banners. | — |

Kept on `gpu-port` as generic, inert features: the device support of upstream's time-dependent
hopping operator (`Prog/OpT_time_dependent.F90`, `Op_T%g_t`, unchanged from master) in
`Hop_mod`/`gpu_wrap_mod`/`gpu_tau_mod`, and the unit test `36-OpT-time-dependent`. No
Hamiltonian of this branch sets `g_t`, so the code path is the upstream one; a user model can
opt in.

Removed with the defect (back to master): the support for tau-resolved local observables
(`Obser_Latt_Local` with `Nt > 1` in `Prog/observables_mod.F90`, `_localtau` analysis in
`Analysis/ana.F90`, `ana_hdf5.F90`, `ana_mod.F90`); on `gpu-port-unified` it is still there.

Defect-production material that is still in this tree but only makes sense on
`gpu-port-unified`: `docs/gpu/TIMEDEP_HOPPING.md`, `docs/gpu/TUNING_DEFECT.md`,
`docs/gpu/sessions/`, `docs/gpu/scripts/defect/` (its `compare.py` and `run_cmp.py` are generic),
`docs/gpu/scripts/pp/localtau_err.py`, `Scripts_and_Parameters_files/make_PP_Uscan.py`,
`Scripts_and_Parameters_files/Start_PP_GPU/parameters_template` (`ham_h0 = 1` there means the
static upstream field on this branch). Benchmarks in the docs were measured on the defect model;
the GPU code paths are the same, apart from the time-dependent hopping.

## Keeping them in sync

- New general work (GPU, updates, MPI, fixes) goes to `gpu-port`. The defect branch takes it
  with `git merge gpu-port` on `gpu-port-unified`. Never merge `gpu-port-unified` into
  `gpu-port`: that would bring the defect back.
- **One-time step, before the first such merge:** the commits that removed the defect from
  `gpu-port` must not be applied to `gpu-port-unified`. Record them as merged without their
  changes, on `gpu-port-unified`, with the `gpu-port` commit that last removes defect code
  (the commit "docs: two branches ..." or later is fine, it only adds docs; use the tip of
  the cleanup):

  ```bash
  git checkout gpu-port-unified
  git merge -s ours <tip of the gpu-port cleanup> -m "Record the gpu-port cleanup as merged (defect stays here)"
  ```

  After that, `git merge gpu-port` brings only later work. Without this step the merge would
  silently remove the defect from `gpu-port-unified` (it is a fast-forward).
- Defect-only changes go to `gpu-port-unified` only. A fix that touches both (for example in
  `gpu_obs_mod.F90` outside `local_obs`) is made on `gpu-port` and merged.
- Check after a merge: `GPU=none` of `gpu-port` vs `master` bitwise (stock parameters), and
  `gpu-port-unified` vs `gpu-port` bitwise with `ham_h0 = 0` and `LOBS_ST = 1`.

## All branches (state 2026-11-13)

Merged-ness: `git merge-base --is-ancestor <branch> origin/gpu-port-unified` (and `origin/gpu-port`); "+N" = `git rev-list --count origin/gpu-port-unified..<branch>`. `[L]` = also a local branch (all equal to their origin copy; there is no local-only branch), `[WT]` = has a worktree (remove it before deleting the branch). `defect/*` are remote-tracking refs of the sibling repo `../ALF_Defect-2` (as of the last fetch), not branches of `origin`. A local branch named `origin` does not exist: `git branch -a` prints the symbolic ref `refs/remotes/origin/HEAD` as `origin` (likewise `defect`). Deleting means the local branch and the `origin` copy, nothing was deleted.

| Branch | What it is | Merged into | May be deleted? |
|---|---|---|---|
| `gpu-port` [L][WT] | general GPU port, no defect (main branch) | main (contained in `gpu-port-unified`) | no (main branch) |
| `gpu-port-unified` [L][WT] | `gpu-port` plus point-defect model, all GPU lines (main branch) | main (contains `gpu-port`; 10 commits not in `gpu-port`) | no (main branch) |
| `master` [L] | upstream ALF 521ddfca | both | no (upstream mirror) |
| `luis_playground` | your MaxEnt playground (Gaussian bootstrap) | both (all commits; 13 scratch files left out) | no (main branch; scratch files stay) |
| `JAKOB_BACKUP` [L] | Jakob's defect backup, pin on first time slice | both | no (main branch, Jakob's data line) |
| `intel-batched` [L][WT] | Intel bottleneck study, batched rewrite begun | none (+5; +8 vs `gpu-port`) | no (active, not merged) |
| `fix/amd-checkerboard` [L][WT] | AMD checkerboard hang fix (`!GCC$ novector`) | both | no (active, AMD checks still open) |
| `speedup/unified` [WT] | LU-free det phase (S1), `la_sync` removal on CUDA (S2), nsys script (S4) | S1 in both (as 985e1b0e + 72a3900f, enabled for oneMKL only); S2, S4 in none (+3 by ancestry) | no (S2 waits for an NVIDIA run) |
| `speedup/nvidia` [WT] | same three changes, separate commits | as `speedup/unified` | no (S2 waits for an NVIDIA run) |
| `gpu-port-apple-metal` [L][WT] | Apple Metal float-float ZGEMM prototype | none (+1) | no (not merged, no gain on M1) |
| `gpu-port-unified-maxent` [L][WT] | work branch: MaxEnt bootstrap merge, PP tool | `gpu-port-unified` (+2 vs `gpu-port`: `prepare_sim.py`, Helma docs) | yes (merged; remove worktree) |
| `gpu-port-unified-tests` [L][WT] | work branch: Run.tests, `run-local` | both | yes (merged; remove worktree) |
| `gpu-delayed-submatrix-defect-production-old-Ntau` [L] | defect variant with old tau labelling of `SpinZ_localtau` | none (+4) | owner decides (4 commits unmerged; kept for old data) |
| `production-mpi-intel-unstable-guard-asyncblk-test` [L][WT] | merge of unstable-guard and team-combined at older tips | none (+1, a merge commit only) | owner decides (1 commit unmerged, merge only) |
| `gpu-delayed-submatrix-defect` [L] | first defect GPU port, delayed/submatrix | both | yes (merged) |
| `gpu-delayed-submatrix-defect-fullgpu` [L] | defect port with full GPU path | both | yes (merged) |
| `gpu-delayed-submatrix-defect-production` [L] | production line, JAKOB_BACKUP merged | both | yes (merged) |
| `gpu-delayed-submatrix-defect-production-mpi` [L] | MPI with GPU path; base of `gpu-port-unified` | both | yes (merged) |
| `gpu-delayed-submatrix-defect-production-mpi-amd` [L][WT] | `GPU=amd` HIP line | both | yes (merged; remove worktree) |
| `gpu-delayed-submatrix-defect-production-mpi-pp` [L] | parallel tempering on the GPU path | both | yes (merged) |
| `gpu-delayed-submatrix-defect-production-mpi-pp-helma` [L][WT] | Helma H100 PP production chains | both | yes (merged; remove worktree) |
| `gpu-delayed-submatrix-defect-production-mpi-pp-opt` [L][WT] | PP profiling, `acc0` reduction (worktree in `/private/var/...wt_pp`) | both | yes (merged; remove worktree) |
| `production-mpi-intel-large-prompt` [L] | Intel SuperMUC large-prompt production | both | yes (merged) |
| `production-mpi-intel-large-prompt-team` [L][WT] | Intel team issue, staging-buffer check | both | yes (merged; remove worktree) |
| `production-mpi-intel-large-prompt-team-async` [L][WT] | Intel async block path (phase 3b) | both | yes (merged; remove worktree) |
| `production-mpi-intel-large-prompt-team-combined` [L][WT] | Intel device QR, async, explicit-shape kernels | both | yes (merged; remove worktree) |
| `production-mpi-intel-large-prompt-team-sycqr` [L][WT] | Intel SYCL device Householder QR | both | yes (merged; remove worktree) |
| `production-mpi-intel-unstable-guard` [L][WT] | `ALF_UNSTABLE_THRESHOLD`, SuperMUC scripts | both | yes (merged; remove worktree) |
| `intel-supermuc-ng2` [L] | SuperMUC-NG Phase 2 session, benchmarks | both | yes (merged) |
| `gpu-port-207-tau` | ticket 207: quad-precision tau reference | both | yes (merged) |
| `gpu-port-207-2` | ticket 207-2: gate fix | both | yes (merged) |
| `gpu-port-207-r` | merge of 207-2 into 207-r | both | yes (merged) |
| `gpu-port-207-3` | ticket 207-3: device correct-side test | both | yes (merged) |
| `gpu-port-207-final` | ticket 207 closing sections | both | yes (merged) |
| `gpu-port-207-4` | 207-4: beta=10 NVIDIA calibration scripts | none (+1) | owner decides (1 commit unmerged) |
| `gpu-port-207-harness` | 207-5: harness QREF gates | none (+4) | owner decides (4 commits unmerged) |
| `gpu-port-207-wip` | 207 evidence record, scripts (on hold) | none (+5) | owner decides (5 commits unmerged) |
| `gpu-port-208-1` | ticket 208-1 closed | both | yes (merged) |
| `gpu-port-208-2` | ticket 208-2 closed | both | yes (merged) |
| `gpu-port-208-3` | ticket 208-3 closed (old `gpu-port` tip 43045813) | both | yes (merged) |
| `gpu-port-212-hardening` | 212: pre-registration only | none (+2) | owner decides (2 commits unmerged) |
| `gpu-port-212-hardening-b` | 212: T2 pass, closed | both | yes (merged) |
| `gpu-port-302a-onemkl` | 302a oneMKL, T2 pass | both | yes (merged) |
| `gpu-port-302b` | 302b ticket, pre-registration | both | yes (merged) |
| `gpu-port-302b-docs` | 302b review and T1/T2 record | none (+2) | owner decides (2 commits unmerged) |
| `gpu-port-302b-impl` | 302b candidate 1 | both | yes (merged) |
| `gpu-port-302b-impl2` | 207 split into 207-1 ... 207-7, PLAN | both | yes (merged) |
| `gpu-port-302c` | 302c-1: device tau pass, `GPU=intel` | both | yes (merged) |
| `gpu-port-302c-2` | 302c-2: B580 control ladders | both | yes (merged) |
| `gpu-port-delayed-submatrix` | NOT_IMPLEMENTED.md, PARKING.md | both | yes (merged) |
| `gpu-port-harness-302b` | acceptance harness after review (c) | both | yes (merged) |
| `gpu-port-harness-gpunone` | harness exclusion list, `GPU=none` | both | yes (merged) |
| `gpu-port-process-lanes` | harness fixes, fourth review | both | yes (merged) |
| `feature/delayed-submatrix` | early delayed-submatrix work (SuperMUC tests) | none (+265) | owner decides (265 commits unmerged) |
| `feature/delayed-submatrix-clean` | same, large text files removed from history | none (+132) | owner decides (132 commits unmerged) |
| `s-wave-Pairing_correlations` | on-site pairing correlation (2026-09) | none (+1) | owner decides (1 commit unmerged) |
| `defect/master` | sibling repo: master mirror | both | yes (merged; ref only) |
| `defect/Point_Pin` | sibling repo: first point-pin version | both | yes (merged; ref only) |
| `defect/Point_Pin_updated` | sibling repo: Point_Pin plus master (its HEAD) | both | yes (merged; ref only) |
| `defect/localtau` | sibling repo: local SpinZ correlation observable | none (+6) | owner decides (6 commits unmerged) |
| `defect/localtau_Obs` | sibling repo: two bug fixes | none (+1) | owner decides (1 commit unmerged) |
| `defect/Predefined_Hoppings_set_OPT-pinning` | sibling repo: predefined hoppings, pinning runs (2025-11) | none (+37) | owner decides (37 commits unmerged) |

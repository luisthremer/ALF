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
| `testsuite/Run.tests/` | `run-schemes` on the standard Hubbard model (`ham_h0 = 0`). | `run-schemes` on the defect, plus `run-pinning` and `run-local`. |
| `docs/gpu/` | this file, banners. | — |

Kept on `gpu-port` as generic, inert features: `Prog/OpT_time_dependent.F90` (time-dependent
hopping operator `Op_T%g_t`, with its device support in `gpu_wrap_mod`/`gpu_tau_mod` and unit
test `36-OpT-time-dependent`). No Hamiltonian of this branch sets `g_t`, so the code path is the
upstream one; a user model can opt in.

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

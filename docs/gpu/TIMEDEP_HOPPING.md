# Time-dependent hopping on the GPU path (space-time point defect)

> **Point-defect material.** The space-time point-defect Hubbard model described here lives on branch
> `gpu-port-unified`; on `gpu-port` `ham_h0` is upstream's static pinning field. See [BRANCHES.md](BRANCHES.md).


Branch `gpu-delayed-submatrix-defect-fullgpu` (from `gpu-delayed-submatrix-defect`, which is `gpu-delayed-submatrix` plus the
merged `Point_Pin_updated`, NTAU1 fix at 58fc5372). Work done 2026-11-07. Scripts: `docs/gpu/scripts/defect/`.
Tuning guide and parameter files: [TUNING_DEFECT.md](TUNING_DEFECT.md).

> **Update 2026-11-09 (merge of `JAKOB_BACKUP` 3852a5f6..7e40b650 into production):** the pin is now one slice, `g_t(1) = -1`
> (tau = 0), not `-0.5` on `Ltrot/2, Ltrot/2+1`. The down sweep measures and labels `NTAU1` (0..Ltrot-1), default `LOBS_ST = 0`;
> the slice-0-as-`Ltrot` remapping described below is removed. `SpinZ_localtau` has `Ltrot+1-2*Thtrot` columns (tau = 0..beta),
> index `Ntau+1-Thtrot`, endpoints weighted 2 (host `Hamiltonian_Hubbard_smod.F90` and device `gpu_obs_mod.F90:local_obs`).
> Statements below about a pulse at beta/2 and the `Ltrot`-column grid describe the state before this merge.

## 1. Starting state

- `gpu-delayed-submatrix-defect` (58fc5372): defect merged, down sweep passes `NTAU1`. A GPU run with `ham_h0 /= 0` was **refused**
  (`gpu_supported`: hopping element is time dependent) and ran entirely on the CPU. `plan/phase1/REPORT.md`: the CPU path of this
  branch is bit-identical to `Point_Pin_updated` with the NTAU1 fix, `standard` vs `submatrix` agree to 1e-11.
- Not contained in any branch although the task assumed it: the second part of the down-sweep fix (store slice 0 as slice `Ltrot`).
  Added here (section 4, `main.F90`).
- Device observables were off for `ham_h0 /= 0` (host `Obser`, one D->H of G per measured slice).

## 2. Gap table (defect run, production parameters: `Checkerboard = T` is ignored on the CPU when `h0 /= 0`, `N_FL = 2`, `Symm = T`)

| Operation in the sweep | Before | After |
|---|---|---|
| Whole GPU path | refused (all on host) | active (`gpu_supported` accepts time-dependent elements that follow all static ones) |
| Static hopping `Op_T(1,nf)` (full dense, `Checkerboard` flag ignored) | dense GEMM (`Hd` for Ndim <= 800, `dense_left/right` above) | unchanged, verified at Ndim = 1024 (L = 32) |
| Time-dependent `Op_T(2,nf)` in the 5 propagation routines | host only | device: per-slice elements, support-reduced (1 site for the point defect), tiny kernel |
| Same in the stabilisation products (`wrapur/wrapul`), `WRAPGRUP/DO`, tau pass | host | device (time argument `t` through `gpu_mmth*`) |
| Update schemes standard / delayed / submatrix | host | device (they only touch the diagonal `Op_V`; G is flushed before each hopping step) |
| `Mz = .F.` (`N_FL = 2`, density coupling, complex couplings) | host | device, verified |
| `Symm` for the measurement | host `Hop_mod_Symm` | device: static `A M B` (once) + time-dependent half steps at the measured slice |
| Scalars, `Green_eq`, `Den_eq`, `SpinZ/XY/T_eq` | host (D->H of G per slice) | device accumulators (`gpu_obs_mod`) |
| `SpinZ_local`, `SpinZ_localtau` | host | device accumulators (diagonal of G only), D->H once per bin |
| `Ltau = 1` pass (`gpu_tau_mod`) | host | device propagation; the tau observables with the defect still use the host `ObserT` (D->H of the four G matrices per tau slice), see section 6 |
| RNG, accept/reject, field bookkeeping, bin I/O | host | host (allowed) |

## 3. Design

Host semantics (`OpT_time_dependent.F90`): with `A(g) = U exp(g E) U^dagger` the element is applied as `A(g)` (`lmult`, `rmult`),
`A(-g)` (`lmultinv`, `rmultinv`) with `g = g_t(t)`, or with `Symm` `g = (g_t(t) + g_t(t-1))/2` (`t-1 = Ltrot` for `t = 1`); the
element is skipped when `|g|^2 <= Zero`. `adjointaction` (Hop_mod_Symm) applies `A(g_t(t1)/2)` from the left and `A(-g_t(t2)/2)` from
the right (no averaging; `t = 0` means `Ltrot`). In `Hop_mod_mm*` the time-dependent elements are the last operators of the flavour
(`nc` is applied in descending order for `mmthr/mmthl_m1/Symm`, ascending for the others).

- `Hop_mod_gpu_td_element` exports `U, E, P, g_t, symm`. `gpu_supported` accepts time-dependent elements only when no static element
  follows them (otherwise: refused with a reason); `Op_V` with `g_t` stays refused.
- **Per-slice elements instead of per-slice dense matrices.** For every active slice and for each of the two coupling kinds (propagation:
  matrices `A(g)`, `A(-g)`; Symm: `A(g_t/2)`, `A(-g_t/2)`) one ordinary complex hopping element is appended to `hop_meta/hop_P/hop_Cval`
  and one single-element group to `hop_gptr`; the host table `td_grp(t,kind,j,nf)` selects it (0 = inactive = skipped exactly as on the
  host). The existing kernels (`el_left_d/el_right_d`, `dense_left/right` for `nh > 8`) apply it. No change to the static element list,
  `Hd`, the groups of the static elements or their kernels: static runs are bit-identical to the previous binary (section 5).
- **Support reduction (`td_reduce`).** For the point defect `Op_T(2,nf)` is an Ndim x Ndim operator with one nonzero eigenvalue
  (`U` is Ndim x Ndim). `A(g) - 1 = sum_k (exp(g E_k) - 1) u_k u_k^dagger` only involves the eigenvectors with `E_k /= 0` and only
  the rows where they are nonzero (threshold 1e-13 relative). The element is stored on that support (1 site, `nh = 1`) instead of an
  Ndim x Ndim matrix; without this a set of matrices at L = 32 needs more than 2 GB. Generic dense time-dependent operators fall back
  to their full support, and `gpu_supported` refuses a set larger than 2 GB.
- **Time argument.** `gpu_mmthr/_m1/lc/l/l_m1(M, ld, n, nf, t)`, `gpu_symm(M, ld, n, nf, t1[, t2])`, `gpu_symm_td(M, ld, n, nf, t)`. The
  argument is required (compile errors found every caller). Slices used: `wrapur/wrapul` the slice `t`, `WRAPGRUP` `ntau+1`, `WRAPGRDO`
  `ntau`, tau pass `nt1`, measurement `Ntau` (label of `Obser`; `0 -> Ltrot`).
- **Device measurement.** `gpu_obs_measure`: with `Symm` and an active time-dependent half step at the label `Ntau`, `G <- A_td G B_td`
  (device, in place on a work copy) before the static `A (.) B` GEMMs; the saved G is restored afterwards. `local_obs` adds
  `n_up - n_dn = G(i,i,2) - G(i,i,1)` (times sign/phase) into `lacc1(I,no)` and `lacc2(I,Ntau-Thtrot,no)` (device), counters on the host
  as in `Obser`; `gpu_obs_flush` adds them to `Obs_local` once per bin. `Obs_local` is now `public` in `Hamiltonian_main` (needed by `main.F90`).
- **Down-sweep label (`main.F90`).** After `WRAPGRDO(GR,NTAU)` G is at `NTAU1 = NTAU-1`; `nt_obs = NTAU1`, and for Hamiltonians that
  have `Obs_local` (the defect) at finite temperature `nt_obs = Ltrot` for `NTAU1 = 0` (G(0) = G(beta)); window check and `Obser` use
  `nt_obs`. Without it slice `Ltrot` gets half the weight in `SpinZ_localtau`. **This changes the CPU result of defect runs** (one more
  measurement per down sweep, in every observable); runs without `Obs_local` are untouched (checked bit-identical, section 5).

## 4. Changed files per work package (commits on the branch)

| WP | Content | Files | Commit |
|---|---|---|---|
| WP1/WP2/WP3/WP5 | time-dependent hopping on the device (equal time, tau pass, Symm), support reduction, `t` argument | `Hop_mod.F90`, `gpu_state_mod.F90`, `gpu_wrap_mod.F90`, `gpu_update_mod.F90`, `gpu_tau_mod.F90` | b676b1da |
| WP4 | device observables for the defect, N_SUN handling, `nt_obs` label fix | `gpu_obs_mod.F90`, `main.F90`, `Hamiltonian_main_mod.F90`, `gpu_state_mod.F90` | 975b7b21 |
| tests | `td_selftest` (+ negative control) in `ALF_GPU_SELFTEST=2` | `gpu_wrap_mod.F90` | 52c4d172 |
| docs/scripts | this file, `TUNING_DEFECT.md`, `FEATURES.md`, `scripts/defect/` (drivers, parameter files, benchmark data) | `docs/gpu/**` | a8e67977, 8f88bc2f |

## 5. Verification (A30 `nividia`, nvfortran 26.9; CPU reference = same source, `GPU=none`, nvfortran, default stabilisation)

The CPU reference has to be built with the same compiler: the random number stream of the GNU build differs (initial configurations
differ), so GNU-vs-nvfortran trajectories are not comparable (known, CLAUDE.md). Fixed seed (`Start/seeds`), `CPU_MAX = 0`, no `confin_0`,
pi flux (`N_Phi = L^2/2`), `U = 5.0`, `h0 = 1`, `beta = L`, `Dtau = 0.1`, `Symm = T` unless stated, `N_SUN = 2` (internally 1 per flavour).
Deviation = max over a file of |A-B| divided by the largest entry of that file (`compare.py`); acceptance and `confout_0` identical in all
passing rows ("id" below).

| Test | Result |
|---|---|
| **U = 0** (G independent of the fields), L = 4 (4 sweeps) and L = 8 (3 sweeps), all files incl. `SpinZ_localtau` (2634 / 20746 numbers) | 3.0e-14 / 8.5e-14 (earlier L = 4 run: 1.1e-15) |
| L = 4, 10 sweeps, Mz = T: standard / delayed(4) / submatrix(4) | id; 3.5e-11 / 4.0e-11 / 3.8e-11 |
| L = 4, Mz = F: standard / delayed / submatrix | id; 1.7e-15 / 1.5e-15 / 2.2e-15 |
| L = 4, Ltau = 1, submatrix: Mz = T / Mz = F | id; 9.4e-11 / 3.0e-15 |
| L = 8, 10 sweeps, Mz = T: standard / delayed / submatrix | id; 5.1e-10 / 5.1e-10 / 5.1e-10 (largest in `SpinZ_localtau`, see below) |
| L = 8, Mz = F: standard / delayed / submatrix | id; 8.4e-15 / 8.4e-15 / 8.4e-15 |
| L = 8, Ltau = 1, submatrix: Mz = T / Mz = F | id; 1.3e-10 / 1.4e-14 |
| L = 4, `Checkerboard = T` (ignored with the defect), Mz = T | id; 5.9e-11 |
| L = 4, `Symm = F`, Mz = F, delayed | id; 5.2e-15 |
| **L = 32** (Ndim = 1024 > `HOPD_NMAX`), beta = 4, 1 sweep: Mz = T submatrix(8) / Mz = F standard / Mz = F delayed(8) | id (confout 40994 entries); 2.5e-12 / 6.0e-15 / 1.4e-14 |
| L = 32, beta = 4, Mz = F, submatrix(8), 2 sweeps (earlier run, same code up to logging) | id; 1.2e-14 |
| Noise level of the Mz = T comparison without the defect: L = 8, `h0 = 0`, Mz = T, 10 sweeps | 7.5e-11 (Ener_scal) ; `Precision Green` mean 1.3e-10 (CPU) / 2.2e-10 (GPU) on the defect run: the Mz = T deviations are the intrinsic accuracy of G (not the defect, not a code difference: Mz = F with the same code is at 1e-14) |
| Physics (L = 6, Mz = T, 60 bins x 40 sweeps, CPU): `SpinZ_localtau` on the pin site | symmetric about `nt = Ltrot/2 - 1/2` (partner `nt <-> Ltrot-1-nt`): max/rms z = 2.3 / 1.3 over 29 pairs; symmetric about `Ltrot/2`: max z = 12.9 (see section 6) |
| Device-vs-host measurement (`ALF_GPU_NOOBS=1`) | equal to 1e-14 |

Self-tests (`ALF_GPU_DEBUG=1 ALF_GPU_SELFTEST=2 ALF_GPU_EXIT=1`, L = 4 and 8, Mz = T/F, with the defect): all five hop routines against
`Hop_mod_mm*` at every active slice, slice 1 and `Ltrot`, and `gpu_symm` against `Hop_mod_Symm` for `t1 = t2`, `t2 = Ltrot`, `t2 = 1`
(80 cases per flavour set): **max rel. diff 0 / 2.4e-16**. Negative control built in (device at an inactive slice against the host at an active
one): diff 6.6e-2, the test fails if it were below 1e-9. `ALF_GPU_SELFTEST=1` passes with and without the defect.
`SELFTEST=3/4` need a `-DQR_PREPIVOT` build (existing requirement) and were not run.

**No regression.**
- GPU, static models (h0 = 0): new binary against the binary of the starting commit, same seed: N_FL = 1 and 2, Mz = T/F, standard/delayed/submatrix,
  Ltau = 0/1, `Checkerboard` T/F, `Symm` T/F, L = 4..8: **bit-identical output files** (dev 0.0) in 6 of 7 cases. The exception: `N_FL = 2`
  with `N_SUN = 4` (internal 2) now takes the device measurement (previously host): 1.5e-15.
- CPU (gfortran build of the new tree against the old tree), h0 = 0: N_FL = 1/2, Ltau 0/1: **bit-identical**.
- The CPU result of defect runs changes by the intended `nt_obs` fix (section 3).

## 6. Findings that are not bugs of this work

1. **The pulse sits at the start of a slot.** `Op_T(2,nf)` is applied before `Op_T(1,nf)` (operator `nc` descending in `mmthr`). I checked
   the CPU Green function at U = 0 against an independent matrix-product model built from the printed `Op_T(1)`: exact agreement (8 digits,
   all slices 17..24, `Symm` T and F) when the pulse `g_t(Ltrot/2), g_t(Ltrot/2+1)` is applied *before* the hopping in slots 20 and 21 (L = 4,
   Ltrot = 40). Consequence: the labelling of G is correct (no off-by-one in the code), but the response is symmetric about
   `nt = Ltrot/2 - 1/2`, i.e. about `tau = beta/2 - Dtau/2`, not about `Ltrot/2` (U = 0, exact for `Symm = F` and `T`; Monte Carlo
   z-test above). The maximum on the pin site sits at slice `Ltrot/2 - 1` or `Ltrot/2` (0.576(15) / 0.569(16) at L = 6). This is a
   definition of the defect on the grid (O(Dtau)), identical on CPU and GPU, and relevant for the finite-size scaling only as a constant
   half-step offset. If a pulse centred on `Ltrot/2` is wanted, `Op_T(2)` has to be placed after `Op_T(1)` or `g_t` shifted; not changed here.
2. **`Mz = .T.` at beta = L = 24, U = 5, `Nwrap = 10`: Precision Green mean 8e-8, max 2e-3** (GPU), while `Mz = .F.` gives 4e-14 / 1e-11. Use
   `Nwrap = 4..5` for `Mz = .T.` (not needed for production, `Mz = .F.`).
3. In `Hamiltonian_Hubbard_smod` `N_SUN` is halved for `N_FL = 2` (file value 2 means 1 per flavour); the device SpinMz formula of `gpu_obs_mod` is for
   `N_SUN = 1` and the host routine measures nothing for `N_SUN /= 1`; `gpu_obs_measure` mirrors this (counts only).
4. **`Op_set` fails for Ndim = 4096 (L = 64)**: `Det_C` of the eigenvector matrix underflows (`|det| = 0`, "Cannot normalize to SU(N)"), on CPU and GPU, for the dense `Op_T(1)`
   of the defect model. Not fixed (outside this task); see TUNING_DEFECT.md.
5. `Obs_local` was private in `Hamiltonian_main` (public only for PGI/GPU builds in the list of the observables); now public for all builds.

## 7. Still on the host, refused or untested

Host (inside the sweep loop): RNG and accept/reject; the field array `nsigma` (bookkeeping, 1 flip at a time, passed to the device in
kernel arguments); in `submatrix` mode the principal block `G0(B,B)` (k x k, `k = N_delay`) and the matrix `M` (k x k) per block, one upload of
`S, M, factors` per block (size O(N k)); per proposal nothing of size Ndim^2. Per bin: D->H of the observable accumulators and G, `confout`.
Initial uploads. The nsys profile (L = 32, beta = 4, 1 sweep, submatrix N_delay 8): H->D 361 MB in 82 673 copies (of which 134 MB are the
initial maps), **D->H 11.6 MB in 20 938 copies, largest 0.66 MB**: no O(Ndim^2) transfer in the sweep; no host BLAS in the sweep (the
timeline contains only device kernels and cuBLAS/cuSOLVER calls). 65% of the host time in the CUDA API is `cuStreamSynchronize` (latency
bound dispatch of 175 000 small kernel launches per sweep, see TUNING_DEFECT.md).

Refused (with a message in `gpu_supported`): time-dependent `Op_V` (`g_t`); a static hopping element after a time-dependent one; a
time-dependent element whose support needs more than 2 GB; everything that was refused before (MPI, STAB1-3, Langevin/HMC, global moves,
`Propose_S0`, Projector, rank > 1 vertices).

Tau-resolved observables with the defect (`Ltau = 1`): the tau pass and propagation are on the device, but `ObserT` is the host routine
(`o_tau` is switched off for `h0 /= 0` because the tau-displaced Symm needs different left/right slices; `gpu_symm(M,..,t1,t2)` supports them,
the device tau observables do not use them yet). Production uses `Ltau = 0`.

Not tested: MPI, HDF5, intel/oneMKL build, `SELFTEST = 3/4`, L = 32 with Mz = T and `Ltau = 1`, Projector (refused), runs longer than 2 sweeps at L >= 32
against the CPU (the CPU reference costs 3.5 min per sweep at beta = 4 with 4 threads), `Nwrap` /= 10 against the CPU.

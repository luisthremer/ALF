# ALF code map: model/physics layer and libraries

Written at `Code-modularity` af22b976 (Phase 6 + 6G). All paths relative to the repo root. `file:N` = line N.
Status: facts read from the code on 2026-11-25. "Inferred" marks statements not written down in the repo.

---

## 1. Modules

### 1.1 `Prog/Hamiltonian_main_mod.F90` (module `Hamiltonian_main`)
Purpose: owns the whole model state and the abstract model class `ham_base`; every model is a submodule of it.

State (all module variables; `PROTECTED_ARRAY` = `, protected` except for PGI/NVHPC, where it is empty, :139-144):
| Variable | Line | Access | Meaning |
|---|---|---|---|
| `ham` (`class(ham_base), allocatable`) | :204 | public | the active model object |
| `Op_V(:,:)` Operator | :206 | public, protected | interaction ops, `(n_op, nf)` |
| `Op_T(:,:)` Operator | :207 | public, protected | hopping ops, `(n_hop, nf)` |
| `WF_L(:)`, `WF_R(:)` WaveFunction | :208-209 | public, protected | trial wave functions per flavour (projector only) |
| `Calc_Fl(:)` logical | :210 | public, protected | flavours actually computed (optional, set by model) |
| `Calc_Fl_map(:)` | :211 | public, protected | indices of computed flavours, size `N_FL_eff` |
| `nsigma` Fields | :212 | public, **writable** | auxiliary fields `f(n_op, nt)` |
| `Ndim, N_FL, N_FL_eff, N_SUN, Ltrot, Thtrot` | :213-217 | public, protected | sizes |
| `Projector, Group_Comm, Symm, reconstruction_needed` | :218-221 | public, protected | flags / MPI group communicator |
| `leap_frog_bulk` | :222 | public, writable | HMC flag (set .false. in `qmc_setup_mod.F90:109`) |
| `Obs_scal, Obs_local, Obs_eq, Obs_tau` | :226-229 | private to module + submodules | observable containers filled by the model |
| `LOG_T0_REJECTED = -huge(1d0)` | :167 | public | sentinel for rejected global moves (log scale) |
| `HAM_VAR_UNSET_INT = -huge(0)` | :170 | private | sentinel for "not set by Ham_set" |
| `override_symmetric_decomposition_check` | :173 | private | silences the Symm warning |

"protected" binds only code outside the module; Hamiltonian submodules write `Ndim`, `N_FL`, `Op_V`, ... directly.

Public routines:
- `Alloc_Ham(ham_name)` :251 — `select case(str_to_upper(name))` from generated `Hamiltonians_case.h` (:256) → `Ham_Alloc_<Name>` (allocates `ham` as the model type), then `Init_Ham_Variables` (:1077, sets sentinels, `Projector=Symm=.false.`).
- `Validate_Ham_Variables()` :1098 — after `ham%Ham_set`: sizes set and positive, `Thtrot>=0`, `Op_V/Op_T` allocated with 2nd dim `N_FL`, `WF_L/R` size `N_FL` if Projector, all `P(i)` in `1..Ndim`, Symm check via `Op_V_is_symmetric` (warning only, :1211-1226).
- `Setup_flavor_map()` :1241 — if `Calc_Fl` unallocated → all `.true.`; sets `N_FL_eff`, `Calc_Fl_map`, `reconstruction_needed = (N_FL_eff /= N_FL)`.
- `Clear_Ham_state()` :1277 — `WF_clear` (projector), `Op_clear` on all Op_V/Op_T, deallocate `Calc_Fl_map`. Called `main.F90:392`.
- `Set_Group_Comm(comm)` :1306 — only writer of `Group_Comm` (called `main.F90:251` after MPI_COMM_SPLIT).
- `Set_model_for_testing(Ndim,N_FL,N_SUN,Ltrot,Thtrot,Projector,Symm,Op_T_in,Op_V_in)` :1324 — unit tests only; `move_alloc` of the operator arrays; call `Setup_flavor_map` afterwards (used by `testsuite/Prog.tests/40-walker.F90:26`, `testsuite/MALA.tests/measurement_cache.F90:114`).

Type `ham_base` (:175-202): **no deferred procedures**; all are `nopass` with a `_base` default. A model overrides by re-binding the same name in its extended type.
| Binding | Default (line) | Default behaviour | Called from |
|---|---|---|---|
| `Ham_set` | :273 | terminate ("Ham_set not defined") — must override | `qmc_setup_mod.F90:110` |
| `Alloc_obs(Ltau)` | :316 | warning only | `main.F90:302` |
| `Obser(GR,Phase,Ntau,Mc_step_weight)` | :541 | warning once | sweep / `Langevin_HMC_mod.F90:180` |
| `ObserT(NT,GT0,G0T,G00,GTT,Phase,Mc_step_weight)` | :582 | warning once | `tau_m_mod.F90:138` |
| `Init_obs(Ltau)` | :650 | init all allocated `Obs_*` | `main.F90:332` |
| `Pr_obs(Ltau)` | :608 | `Print_bin_*` all allocated `Obs_*` on `Group_Comm` | `main.F90:368` |
| `S0(n,nt,Hs_new)` | :290 | returns 1; terminates if `Op_V(n,1)%type==1` | `upgrade`, `Delta_S0_global_base` |
| `Global_move(T0_ratio,nsigma_old,size_clust)` | :348 | terminate | via `Global_move_log_T0` |
| `Global_move_log_T0(...)` | :387 | log of `Global_move` with safe-range warning | `Global_mod.F90:552` |
| `Delta_S0_global(nsigma_old)` | :443 | product of `1/S0` over all fields (prints notice once) | `Get_Delta_S0_global_base` |
| `Get_Delta_S0_global(nsigma_old)` | :496 | `log(Delta_S0_global)` | `Global_mod.F90:269` |
| `Global_move_tau(...)` | :721 | terminate | `Wrapgr_mod.F90:529` |
| `Global_MALA_move_tau(Flip_list,Flip_length,ntau)` | :758 | all fields, warning once | `Wrapgr_mod.F90:656` |
| `Global_MALA_move(Flip_list)` | :803 | `Flip_list=1`, warning once | `Langevin_HMC_mod.F90:622` |
| `Hamiltonian_set_nsigma(Initial_field)` | :837 | no-op (allocate `Initial_field` to force a start config) | `qmc_setup_mod.F90:180` |
| `Overide_global_tau_sampling_parameters(...)` | :862 | no-op | `qmc_setup_mod.F90:212` |
| `Ham_Langevin_HMC_S0(Forces_0)` | :884 | zero, warning once | Langevin/HMC |
| `Ham_Langevin_HMC_S0_single(Force_0,n,nt)` | :915 | via full version | `Wrapgr_mod.F90:151` |
| `Apply_B_HMC(Forces_HMC,ltrans)` | :975 | no-op | `Langevin_HMC_mod.F90:455` |
| `weight_reconstruction(weight)` | :990 | terminate | `Global_mod.F90:184` (only if `reconstruction_needed`) |
| `GR_reconstruction(GR)` / `GRT_reconstruction(GT0,G0T)` | :1013 / :1039 | terminate | `Langevin_HMC_mod.F90:179`, `tau_m_mod.F90:136` |
| `write_parameters_hdf5(filename)` (HDF5 only) | :1050 | no-op | `qmc_setup_mod.F90:459` |

Call order of setup (`qmc_setup_mod.F90:103-156`, `Setup_model`): `Fields_init(Amplitude)` → `Alloc_Ham` → `ham%Ham_set` → `Validate_Ham_Variables` → (if any `.not.Calc_Fl`) `check_compatibility_reconstruct_greens_function` → check `size(g_t)==Ltrot` for every Op with `g_t` → `Setup_flavor_map` → projector WF check → `set_default_values_measuring_interval` → `Hop_mod_init`. Then `Setup_fields` (:163): `nsigma%make(N_op,Ltrot)`, `nsigma%t(n)=Op_V(n,1)%type`, `flip_protocol` from `Op_V(n,1)`, RNG from file `seeds`, `ham%Hamiltonian_set_nsigma`, `nsigma%in(Group_Comm[,Initial_field])`.

### 1.2 `Prog/Hamiltonians/` (submodules of `Hamiltonian_main`)
| File | Lines | Overrides beyond Ham_set/Alloc_obs/Obser/ObserT/write_parameters_hdf5 | Namelists |
|---|---|---|---|
| `Hamiltonian_##NAME##_smod.F90` | 443 | template | `VAR_##XYZ##` (empty) |
| `Hamiltonian_Hubbard_smod.F90` | 1170 | S0, Ham_Langevin_HMC_S0(_single), Get_Delta_S0_global, Global_move_tau, Overide_..., Global_move, Global_MALA_move | VAR_lattice, VAR_Model_Generic, VAR_Hubbard |
| `Hamiltonian_Hubbard_Plain_Vanilla_smod.F90` | 866 | weight/GR/GRT_reconstruction | VAR_lattice, VAR_Hubbard_Plain_Vanilla (sets N_SUN=1, N_FL=2 in code :234-235) |
| `Hamiltonian_Kondo_smod.F90` | 850 | — | VAR_lattice, VAR_Model_Generic, VAR_Kondo |
| `Hamiltonian_LRC_smod.F90` | 887 | Global_move_tau, Overide_..., S0, Ham_Langevin_HMC_S0(_single) | ..., VAR_LRC |
| `Hamiltonian_Spin_Peierls_smod.F90` | 990 | weight/GR/GRT_reconstruction, S0, Hamiltonian_set_nsigma | ..., VAR_Spin_Peierls |
| `Hamiltonian_Z2_Matter_smod.F90` | 1756 | Global_move_tau, Hamiltonian_set_nsigma, Overide_..., Get_Delta_S0_global, S0 | VAR_lattice, VAR_Z2_Matter |
| `Hamiltonian_tV_smod.F90` | 778 | — | ..., VAR_tV |
| `LRC_mod.F90` | 555 | module (not submodule) for long-range Coulomb: `LRC_Set_VIJ`, `LRC_draw_field`, `LRC_Print`, `LRC_Clear` | — |
| `Bid_N2*.F90, Bid_N3, Bid_N4` | ~20 | generated formula snippets (Doxygen-excluded, `Prog/Doxyfile:828`) | — |
| `Not_supported/` | — | 8 old models, not in the build | — |

Hubbard as reference (`Hamiltonian_Hubbard_smod.F90`):
- :136-154 `type, extends(ham_base) :: ham_Hubbard` with the overrides; :155-187 three `#PARAMETERS` blocks; :190-194 submodule-private state (`Latt`, `Latt_unit`, `Hopping_Matrix(:)`, `List`, `Invlist`).
- :198 `module Subroutine Ham_Alloc_hubbard` → `allocate(ham_Hubbard::ham)`; :204 `#include "Hamiltonian_Hubbard_read_write_parameters.F90"` (generated).
- `Ham_Set` :213: `read_parameters()` → `Ltrot=nint(beta/dtau)`, `Thtrot=nint(theta/dtau)` if Projector, `Ltrot=Ltrot+2*Thtrot` (:242-245) → N_FL/N_SUN checks; if N_FL=2 then `N_SUN=N_SUN/2` (:270-272) → `Ham_Latt` (:360, `Predefined_Latt`) → `Ham_Hop` (:376, `Set_Default_hopping_parameters_<lattice>` + `Predefined_Hoppings_set_OPT`, optional pinning) → `Ham_V` (:492, `Predefined_Int_U_MZ[_continuous_HS]` or `..._U_SUN[...]`) → `Ham_Trial` (:471, `Predefined_TrialWaveFunction`) → append to `info` on `irank_g==0`.
- Optional flavour reconstruction hook: commented `Calc_Fl` lines :288-292.
- `Alloc_obs` :568: `Obs_scal(4)` Kin/Pot/Part/Ener via `Obser_Vec_make`; `Obs_eq` Green/SpinZ/SpinXY/SpinT/Den (N_FL=2) via `Obser_Latt_make`; `Obs_tau` only if `Ltau==1`, `Nt=Ltrot+1-2*Thtrot`, channel `'T0'` if Projector.
- `Obser` :700ff: `ZP=Phase/Re(Phase)`, `ZS=sign(Re Phase)*Mc_step_weight`; builds `GRC = 1 - GR^T`; increments `Obs%N`, `Obs%Ave_sign`; uses `Predefined_Obs_eq_*_measure`.
- `S0` (continuous HS): `exp((-Hs_new^2 + f^2)/2)`; HMC forces = `Re f` for type-3 ops.

### 1.3 Build glue: `Prog/Hamiltonians.list(.d)`, `parse_ham.py`, `parse_ham_mod.py`, `Prog/Makefile`
- `Hamiltonians.list`: one model per line `Name [path]`; default path `Hamiltonians/Hamiltonian_<Name>_smod.F90` (`parse_ham_mod.py:17-42`); `#` lines ignored. Every file in `Hamiltonians.list.d/` (directory currently absent) is appended.
- `parse_ham.py` flags (:20-50): `--print_obj_list`, `--print_src_list`, `--create_hamiltonians_interface` (→ `Hamiltonians_interface.h`, `module subroutine Ham_Alloc_<Name>` interfaces, :70-79), `--create_hamiltonians_case` (→ `Hamiltonians_case.h`, `Case('<NAME UPPER>') call Ham_Alloc_<Name>`, :81-88), `--create_read_write_par` (→ `Hamiltonians/Hamiltonian_<Name>_read_write_parameters.F90`, :90-99).
- `parse(filename)` (`parse_ham_mod.py:44`): text between `!#PARAMETERS START# <NamelistName>` and `!#PARAMETERS END#`; each line with `::` is `type :: name = default ! comment` (:78-103). A line whose type starts with `!` (e.g. `!Integer :: N_FL = 1`) is `defined_in_base`: not declared locally (it is the `Hamiltonian_main` variable) but still in the namelist, defaulted and broadcast.
- Generated `read_parameters()` (template `parse_ham_mod.py:154-213`): declares `NAMELIST /<block>/`, assigns all defaults at run time, rank `irank_g==0` opens `parameters` (or `Temp_<igroup>/parameters` with TEMPERING), `REWIND; READ(NML=<block>)` for each block, then `MPI_BCAST(...,0,Group_Comm)` per variable (raw MPI, :305-315). `write_parameters_hdf5` writes each variable as attribute under `parameters/<block>` in the HDF5 file (HDF5 builds).
- Makefile: `Hamiltonian_main_mod.o` depends on the two `.h` files (:30); all model objects depend on `.parse_parameters` (:32), which is rebuilt only when `Hamiltonians.list` or `.list.d/*` change (:43-46). **Changing a `#PARAMETERS` block in an existing model does not regenerate the reader**: run `make tidy` (removes generated files, :66-68) or `./parse_ham.py --create_read_write_par` in `Prog/`.
- `gen_deps.py` writes `.deps.mk` (OBJS/MODS + USE dependencies) from `SRCS` (:4-10) + model sources.
- `default_parameters_generic.py`: Python dict of generic namelists (`VAR_QMC`, `VAR_errors`, `VAR_TEMP`, `VAR_Max_Stoch`) with defaults/comments for pyALF; not used by the Fortran build.

### 1.4 `Prog/Operator_mod.F90` (module `Operator_mod`) — see section 4 for fields
No `private` statement: everything public. Uses `mpi_shared_memory` (U and M_exp may live in node-shared memory).
- `Op_make(Op,N)` :196 — allocate `O(N,N)`, `P(N)`; zero; `type=0`, `diag=.false.`, `flip_protocol=1`.
- `Op_set(Op)` :262 — (**known issue:** its det U normalisation uses `DET_C`, which underflows for large N; see `KNOWN_ISSUES.md`) validates (`type` in 0..4, `P>=1`, `O` Hermitian), sets `diag`, diagonalises `O` → `U`,`E`, sorts non-zero eigenvalues first (`N_non_zero`), normalises `det U = 1` (column 1 / det); for type 1/2 precomputes `E_exp`, `M_exp`; marks `g_t_alloc` if `g_t` allocated. Must be called after filling an operator (also for Op_T).
- `Op_clear(Op,N)` :227 — deallocate.
- `Op_phase(Phase,Op_V,nsigma,nf)` :163 — multiplies by `exp(i Im(g alpha Phi))` over all `(n,nt)`.
- `Op_exp(g,Op,Mat)` :494 — `Mat = exp(g O)` (full N×N).
- `Op_mmultL(Mat,Op,HS_Field,cop,nt,sign)` :558, `Op_mmultR(Mat,Op,HS_Field,cop,nt)` :651 — multiply `Mat` by `exp(±g Phi P^T O P)` from left/right; `cop` = 'N'/'T'/'C'.
- `Op_Wrapup` :746 / `Op_Wrapdo` :865 `(Mat,Op,HS_Field,Ndim,N_Type,nt)` — Green-function wrapping; `N_Type`=1 full, 2 basis change only (not `Op%type`).
- `Op_equal(a,b)` :966, `Op_V_is_symmetric(Op_V,N_FL,Ndim)` :1037 (all-diagonal shortcut, else column-local products), `Op_is_real(Op)` :1151, type-bound `get_g_t_alloc()` :1170.

### 1.5 `Prog/Fields_mod.F90` (module `Fields_mod`) — see section 4
Public: `Fields`, `Fields_init` (:71-72). Private module state: `Phi_st(-2:2,2)`, `Gama_st(-2:2,2)`, `Flip_st(-2:2,3)`, `Amplitude` (:75-77), set in `Fields_init(Amplitude_in)` :258.

### 1.6 `Prog/WaveFunction_mod.F90`
`Type WaveFunction` (:53): `P(Ndim,N_part)` complex, `Degen` real. `WF_overlap(WF_L,WF_R,Z_norm)` :76 (rescales WF_R so ⟨L|R⟩=1, avoids fake sign), `WF_alloc(WF,Ndim,N_part)`, `WF_clear(WF)` (pure, ~:110-125).

### 1.7 Hopping propagation: `Hop_mod`, `ContainerElementBase_mod`, `DynamicMatrixArray_mod`, `OpTTypes_mod`, `OpT_time_dependent.F90`, `matTypes_mod`
- `ContainerElementBase_mod.F90`: abstract `ContainerElementBase` (:41-50) with deferred `rmult, lmult, rmultinv, lmultinv (this,arg,t)`, `adjointaction(this,arg,t1,t2)`, `dump`, `dealloc`.
- `DynamicMatrixArray_mod.F90`: `DynamicMatrixArray` (:47) = growable vector of `class(ContainerElementBase)` pointers; `init`, `dealloc`, `pushback(itm)`, `at(pos)`, `back()`, `length()`.
- `OpTTypes_mod.F90`: `RealExpOpT` (:50) / `CmplxExpOpT` (:74) store `mat=exp(gO)`, `invmat`, and half-step `mat_1D2`, `invmat_1D2` (symmetrised upper part), plus `P`, `g`, `Ndim_hop`; mult via `ZDSLSYMM`/`ZSLHEMM`.
- `OpT_time_dependent.F90` (module `OpT_time_dependent_mod`): `OpT_time_dependent` (:50) stores `U`, `E`, `g_t(:)`, `symm`; exponent built per time slice `t`.
- `Hop_mod.F90`: private `ExpOpT_vec(N_FL)` (:59). `Hop_mod_init` (:117) → `OpT_postprocess` (:77) per `Op_T(nc,nf)`: time-dependent if `g_t` allocated, else real or complex by `Op_is_real`. Apply routines (`In` modified in place, flavour `nf`, slice `t`): `Hop_mod_mmthr` (:142, `e^{-dtau T} In`), `Hop_mod_mmthr_m1` (:161, `e^{dtau T} In`), `Hop_mod_mmthl` (:184, `In e^{-dtau T}`), `Hop_mod_mmthlc` (:207), `Hop_mod_mmthl_m1` (:230, `In e^{dtau T}`), `Hop_mod_Symm(Out,In,t1[,t2])` (:278, `e^{-dtau T/2} In e^{dtau T/2}` for computed flavours).
- `matTypes_mod.F90` (`RealMat`, `CmplxMat`): **not compiled** (not in `Prog/Makefile:4-10`; excluded in `.github/workflows/tests.yml:37`).

### 1.8 `Predefined_*` (helpers called by models)
- `Predefined_Latt_mod.F90` (`Predefined_Lattices`): `Predefined_Latt(Lattice_type,L1,L2,Ndim,List,Invlist,Latt,Latt_Unit)` :100; types SQUARE, N_LEG_LADDER, BILAYER_SQUARE, TRIANGULAR, KAGOME, HONEYCOMB, BILAYER_HONEYCOMB, PI_FLUX (:114-226).
- `Predefined_Hop_mod.F90` (`Predefined_Hoppings`): `Hopping_Matrix_type` (:57-82: bonds `T`, `T_loc`, `list`, flux `N_Phi, Phi_X, Phi_Y, Bulk`, checkerboard families `N_Fam, L_Fam, List_Fam, Prop_Fam`); `Set_Default_hopping_parameters_{square,triangular,kagome,N_Leg_Ladder,honeycomb,Bilayer_square,Bilayer_honeycomb}` (:255-1226); `Predefined_Hoppings_set_OPT(this,List,Invlist,Latt,Latt_unit,Dtau,Checkerboard,Symm,OP_T[,pinned_vertices,pinning_factor,pinning_offset])` :1625 (sets `Op_T%g=-Dtau` or `-Dtau*Prop_Fam`, calls `Op_set`); `Predefined_Hoppings_Compute_Kin(...,GRC,Z_Kin)` :1851; `Predefined_hoppings_clear` :96; `inquire_hop` :122.
- `Predefined_Int_mod.F90` (`Predefined_Int`): one `Operator` per call, each ends in `Op_set`: `_U_SUN` :61 (type 2), `_U_SUN_continuous_HS` :92, `_U_MZ(Op_up,Op_do,...)` :123, `_U_MZ_continuous_HS` :161, `_V_SUN` :200, `_VJ_SUN` :232, `_Ising_SUN` :264, `_LRC` :296, `_Jz` :328, `_SaSb` :381.
- `Predefined_Obs_mod.F90` (`Predefined_Obs`): equal time `Predefined_Obs_eq_{SpinSUN,SpinMz,Green,Den}_measure(Latt,Latt_unit,List,GR,GRC,N_SUN,ZS,ZP,Obs...)` :78-267; time displaced `Predefined_Obs_tau_{Green,SpinSUN,SpinMz,Den}_measure(...,NT,GT0,G0T,G00,GTT,...)` :338-534; scalars `Predefined_Obs_scal_Renyi_Ent_*`, `..._Mutual_Inf_*` :665-759.
- `Predefined_Trial_mod.F90` (`Predefined_Trial`): `Predefined_TrialWaveFunction(Lattice_type,Ndim,List,Invlist,Latt,Latt_unit,N_part,N_FL,WF_L,WF_R)` :115 (ground state of a hopping problem; Kekulé distortion for HONEYCOMB :189).

### 1.9 `Prog/observables_mod.F90` (module `Observables`)
Types: `Obser_Vec` (:53: `N`, `Ave_Sign`, `Obs_vec(:)`, `analysis_mode`, `description(:)`; bindings `init`, `print_bin`, `measure`), `Obser_Latt_Local` (:70: `Obs_Latt(i,tau,norb)`), `Obser_Latt` (:91: `Obs_Latt(i-j,tau,no,no)`, `Obs_Latt0(no)`, pointer to `Latt`/`Latt_unit` — linked, not copied, `dtau`, `Channel` in {T0,P,P_PH,PH,PP} :107-112).
API: `Obser_Vec_make(Obs,N,Filename[,analysis_mode,description])` :269, `Obser_Latt_make(Obs,Nt,Filename,Latt,Latt_unit,Channel,dtau)` :194, `Obser_Latt_local_make` :121, `*_Init` :184/:258/:319, `Obser_vec_measure(obs,value,Phase)` :329, `Print_bin_Latt(Obs,Group_Comm)` :355 (suffix `_eq`/`_tau`), `Print_bin_Latt_Local` :567 (`_local`/`_localtau`), `Print_bin_Vec` :739 (`_scal`). Print routines reduce with `comm_mod` (`comm_reduce_mean[_array]`, `is_root(Group_Comm)`, :415-419, :624-627, :782-785) and write text or `data.h5`.

### 1.10 `Prog/QMC_runtime_var_mod.F90` (module `QMC_runtime_var`) — see section 3
State (all private, :47-77): `Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX, Propose_S0, Global_moves, N_Global, Global_tau_moves, Nt_sequential_start/end, N_Global_tau, Sequential, Amplitude, Langevin, HMC, Leapfrog_steps, N_HMC_sweeps, Delta_t_Langevin_HMC, Max_Force, *MALA*`, `ham_name`, `mpi_per_parameter_set`, `Tempering_calc_det`, TEMPERING: `N_exchange_steps, N_Tempering_frequency`; `settings_locked` (:60).
Namelists: `VAR_QMC` :80-87, `VAR_HAM_NAME` :89, `VAR_TEMP` :77 (TEMPERING).
Public: `get_X/set_X` pairs (:114-153), `set_QMC_runtime_default_var` :167, `read_and_broadcast_QMC_var_and_ham_name(Group_Comm)` :917, `broadcast_QMC_runtime_var(comm)` :840 (MPI), `read_and_broadcast_TEMPERING_var` :886, `set_default_values_measuring_interval(Thtrot,Ltrot,Projector)` :217, checks `check_langevin_schemes_and_variables` :253, `check_update_schemes_compatibility` :285 (ends with lock, :312), `check_MALA_variables` :370, `check_compatibility_reconstruct_greens_function` :186, policies `apply_CPU_MAX_policy` :322 (`NBin=10^7` if CPU_MAX>0), `apply_Langevin_overrides` :333, `set_default_global_tau_sampling(N_op)` :353, `lock_QMC_runtime_settings` :424.

### 1.11 `Libraries/Modules/` (archive `modules_90.a`, list in `Libraries/Modules/Makefile`)
| File | Module | Purpose |
|---|---|---|
| `comm_mod.F90` | `comm_mod` | MPI-optional communication layer (section 5) |
| `linalg_mod.F90` | `linalg_mod` | 1:1 typed wrappers `la_z*` of BLAS/LAPACK used by the QMC core (single GPU-backend switch point, :37-49): `la_zgemm, la_zhemm, la_ztrsm, la_ztrmm, la_zgemv, la_zgeru, la_zscal, la_zcopy, la_zlacpy, la_zgetrf, la_zgetrs, la_zgetri, la_zgeqrf, la_zgeqp3, la_zungqr, la_zunmqr` (:54-55). Used by cgr1, cgr2_2, udv_state, QDRP_decompose, wrapul, UDV_WRAP, upgrade, Global, WaveFunction, mymats, Mat_subroutines. New BLAS/LAPACK calls should go through it. |
| `mymats_mod.F90` | `MyMats` | generic linear algebra: `MMULT`, `INITD`, `COMPARE`, `DET`(=`DET_C`), `CT`, `INV`, `UDV`, `QR`, `SVD`, `DIAG` (real sym / Hermitian), `DIAG_GEN`, `SECONDS` (:49-85) |
| `Mat_subroutines_mod.F90` | `mat_subroutines` | small×large products with a sparse operator `P^T A P`: `ZSLGEMM(side,op,N,M1,M2,A,P,Mat)` :51, `ZSLHEMM` :519 (Hermitian A), `ZDSLSYMM` :968 (real symmetric A) |
| `matrix_mod.F90` | `Matrix` | `Mat_C`, `Mat_R` (pointer `el(:,:)`, `dim`), `Make_Mat`, `Clear_Mat` |
| `random_wrap_mod.F90` | `Random_Wrap` | Fortran intrinsic RNG wrapper: `Get_seed_Len`, `Ranset(seed_vec)`, `Ranget(seed_vec)`, `ranf_wrap()` U(0,1), `rang_wrap()` N(0,1) (Box-Muller), `nranf(N)` integer in 1..N, `lcg(seed)` |
| `lattices_v3_mod.F90` | `Lattices_v3` | `Unit_cell` (`Norb, N_coord, Orb_pos_p`), `Lattice` (`N, Ns, list, invlist, nnlist, listk, invlistk, imj, a1_p, a2_p, b1_p, b2_p, L1_p, L2_p, ...`, :53-68); `Make_lattice(L1_p,L2_p,a1_p,a2_p,Latt[,nnlist_range])` :106, `Clear_Lattice`, `Inv_R(x,Latt)` :541, `Inv_K` :495, `npbc`, `Iscalar`, `Xnorm`, `Fourier_K_to_R`, `Fourier_R_to_K`, `Print_latt` |
| `lattices_interface_mod.F90` | `lattices_interface` | C interface to lattices for pyALF |
| `mpi_shared_mem_mod_v2.F90` | `mpi_shared_memory` | node-shared allocations: `mpi_shared_memory_init(comm,chunk_size)` :98, generic `allocate_shared_memory` 1-4D real/cmplx :81, `deallocate_all_shared_memory` :592, flag `use_mpi_shm` :67 |
| `runtime_error_mod.F90` | `runtime_error_mod` | error enum `ERROR_NONE..ERROR_RUNNING_FILE_FOUND` (:56-68); `Terminate_on_error(code,__FILE__,__LINE__)` :81 → `error stop` or `MPI_ABORT(MPI_COMM_WORLD)`; `alf_exit(code)`: plain exit with status 0 / 2 / else 1, no message (for statuses job scripts read) |
| `natural_constants_mod.F90` | `Natural_Constants` | `pi, twopi`, tolerances `Eps_machine=1d-12, Eps_small=1d-10, Eps_convergence=1d-6` (:8-10), SI constants |
| `precdef_mod.F90` | `precdef` | kind params `byte,long,single,double`; `rone, rzero, cone, czero` |
| `files_mod.F90` | `Files_mod` | `File_i(file,i)`, `File_add(file,file1)`, `str_to_upper(s)` |
| `alf_hdf5_mod.F90` | `alf_hdf5` | HDF5 helpers: `init_dset`, `append_dat`, `write_latt`, `write_attribute`, `read_attribute`, `test_attribute` |
| `entanglement_mod.F90` | `entanglement_mod` | Rényi entropy / mutual information over replicas: `Init_Entanglement_replicas(Group_Comm)` :61, `Calc_Renyi_Ent`, `Calc_Mutual_Inf` (header @brief is a stale copy) |
| `errors_mod.F90` | `ERRORS` | jackknife/bin error analysis: `ERRCALC`, `ERRCALCJ`, `COV`, `COV_ERR`, `INTERGRATE(_F)` |
| `fourier_mod.F90` | `Fourier` | Matsubara ↔ tau transforms (`Matz_tau*`, `Tau_Matz*`, splines) |
| `maxent_mod.F90`, `maxent_stoch_mod.F90` | `MaxEnt_mod`, `MaxEnt_stoch_mod` | classic and stochastic analytic continuation |
| `histograms_mod.F90`, `histograms_v2_mod.F90` | `Histograms(_v2)` | histograms (v2 on a log mesh) |
| `log_mesh_mod.F90` | `Log_Mesh` | `Make_log_mesh`, `Lookup_log_mesh`, `Inter_log_mesh` |
| `check_omp_num_threads_mod.F90` | `check_omp_num_threads_mod` | `check_omp_num_threads()` |
| `iso8601_datetime_mod.F90` | `iso8601_datetime_mod` | `iso8601_datetime()` timestamp string |
`Libraries/libqrref/`: reference LAPACK QR sources (`zgeqrf`, `zungqr`, `zlarf*`, ...). There is no module named `mat_mod`; linear algebra lives in `linalg_mod`, `MyMats`, `mat_subroutines`.

---

## 2. How to add a new Hamiltonian `Foo`
1. Copy `Prog/Hamiltonians/Hamiltonian_##NAME##_smod.F90` → `Prog/Hamiltonians/Hamiltonian_Foo_smod.F90`; replace every `##NAME##` with `Foo` (submodule name `ham_Foo_smod`, type `ham_Foo`, `Ham_Alloc_Foo`, include file name). The `Ham_Alloc_<Name>` spelling must equal the name in the list (the generated interface uses it verbatim).
2. Add a line `Foo` to `Prog/Hamiltonians.list` (or a file in `Prog/Hamiltonians.list.d/`); add the path as 2nd column only for a non-default location.
3. Parameters: declare them inside `!#PARAMETERS START# VAR_<block>` … `!#PARAMETERS END#` as `type :: name = default ! comment` (one per line; arrays as `[a,b]`). Reuse `VAR_lattice` and `VAR_Model_Generic`; base variables (`N_SUN`, `N_FL`, `Symm`, `Projector`) are listed as commented declarations `!Integer :: N_FL = 1` so they are read into `Hamiltonian_main` (Hubbard :162-175). Name the model block `VAR_Foo`: the smoke test `testsuite/Prog.tests/36-Ham-smoke.F90` writes only `&VAR_lattice/`, `&VAR_Model_Generic/`, `&VAR_<Name>/` and fails if the model reads any other namelist.
4. Generated, do not edit: `Prog/Hamiltonians_interface.h`, `Prog/Hamiltonians_case.h`, `Prog/Hamiltonians/Hamiltonian_Foo_read_write_parameters.F90` (contains `read_parameters()` with MPI broadcast and `write_parameters_hdf5`), `Prog/.deps.mk`. After editing a `#PARAMETERS` block later, run `make tidy` in `Prog/` (the reader is regenerated only when the list changes, `Prog/Makefile:43`).
5. Implement in the type-bound list (`type, extends(ham_base) :: ham_Foo`):
   - required: `Ham_Set` — call `read_parameters()`; set `Ndim`, `N_FL`, `N_SUN`, `Ltrot` (incl. `+2*Thtrot`), `Thtrot`, `Projector`, `Symm`; allocate and fill `Op_T(:,N_FL)` and `Op_V(:,N_FL)` (each op: `Op_make` → fill `P`, `O`, `g`, `alpha`, `type` → `Op_set`, or use `Predefined_*`); if Projector allocate `WF_L(N_FL)`, `WF_R(N_FL)`; optionally allocate `Calc_Fl(N_FL)` for flavour reconstruction; write `info` on `irank_g==0`.
   - usual: `Alloc_obs(Ltau)` (allocate `Obs_scal/Obs_eq/Obs_tau/Obs_local` with `Obser_*_make`), `Obser`, `ObserT`; `write_parameters_hdf5` (inside `#ifdef HDF5`; implementation comes from the generated include).
   - only if needed: `S0` (any type-1 field or a bare action for type 3), `Global_move`/`Delta_S0_global`/`Get_Delta_S0_global` (Global_moves), `Global_move_tau` + `Overide_global_tau_sampling_parameters` (Global_tau_moves), `Ham_Langevin_HMC_S0(_single)`/`Apply_B_HMC` (Langevin/HMC with bare action), `Global_MALA_move(_tau)`, `Hamiltonian_set_nsigma` (start configuration), `weight_reconstruction`, `GR_reconstruction`, `GRT_reconstruction` (if some `Calc_Fl=.false.`; not allowed with Langevin/HMC/MALA, `QMC_runtime_var_mod.F90:186`).
6. Rules checked at run time: `Op_V`/`Op_T` second dim = `N_FL`; all `P(i)` in 1..Ndim; `O` Hermitian; `Op_V(n,nf)%type` and `flip_protocol` are taken from flavour 1 (`qmc_setup_mod.F90:171-173`); `g_t` (if used) has size `Ltrot`.
7. Build `make` in `Prog/`; `testsuite/Prog.tests/CMakeLists.txt:140ff` adds `36-Ham-smoke-Foo` automatically. Run in a directory with a `parameters` file containing `&VAR_ham_name ham_name='Foo' /` plus `VAR_QMC` and the model namelists, and a `seeds` file.

---

## 3. How to add a runtime (QMC) parameter `Xyz` (`Prog/QMC_runtime_var_mod.F90`)
1. Declare the private module variable (:47-77 block).
2. Add it to `NAMELIST /VAR_QMC/` (:80-87) (or `VAR_TEMP` :77 for tempering-only).
3. Default: assign in `set_QMC_runtime_default_var` (:167-183). It is called inside `read_and_broadcast_QMC_var_and_ham_name` on the reading rank before `READ(5,NML=VAR_QMC)` (:955).
4. Broadcast: add `CALL MPI_BCAST(Xyz,1,MPI_<type>,0,MPI_COMM_i,ierr)` in `broadcast_QMC_runtime_var` (:840-882). The communicator is `MPI_COMM_WORLD`, or `Group_Comm` with `PARALLEL_PARAMS` (:935-943). (Raw MPI here; `comm_bcast` in `comm_mod` only has a real8 version.)
5. Accessors: `get_Xyz()` function and `set_Xyz(val)` subroutine that first calls `ensure_settings_unlocked("set_Xyz")` (pattern :442-450); add both to the `public` list (:114-153).
6. Lock: settings are locked by `lock_QMC_runtime_settings` at the end of `check_update_schemes_compatibility` (:312), called in `qmc_setup_mod.F90:247`. After that every `set_*` and policy routine terminates (`ensure_settings_unlocked` :428). Derived values or overrides must be applied before (pattern: `apply_*` routines :322-360, called in `main.F90:287`, `qmc_setup_mod.F90:205,235`). Consistency checks go into `check_*` routines.
7. Consumers `use QMC_runtime_var, only: get_Xyz`; no direct variable access.
8. Also update: `Prog/default_parameters_generic.py` (`VAR_QMC` dict, pyALF defaults), `Scripts_and_Parameters_files/Start/parameters`, `Documentation/*.tex`, and `testsuite/Prog.tests/37-QMC-runtime-var.F90` (defaults / set-get round trip) and `38-QMC-runtime-lock-fail.F90` (lock).

Model-specific parameters are NOT added here: they go into the model's `#PARAMETERS` block (section 2.3).

---

## 4. Operator and Fields

### 4.1 `Type Operator` (`Prog/Operator_mod.F90:59-93`)
Represents `g · Phi(field) · (c† P^T O P c + alpha)` with exponent `exp(g Phi P^T O P)`.
| Field | Type | Meaning |
|---|---|---|
| `N` | int | size of the local block (number of sites in `P`) |
| `P(N)` | int ptr | site indices 1..Ndim; column i of P has its 1 at row `P(i)` |
| `O(N,N)` | cmplx ptr | Hermitian local matrix (checked in `Op_set`) |
| `U(N,N)` | cmplx ptr | eigenvectors of O (det U = 1); may be MPI shared |
| `E(N)` | real ptr | eigenvalues; the first `N_non_zero` are non-zero |
| `N_non_zero` | int | number of non-zero eigenvalues (set by `Op_set`) |
| `diag` | logical | O diagonal (set by `Op_set`) |
| `g` | cmplx | coupling; for Op_T `g=-Dtau` (·`Prop_Fam` for checkerboard), `Predefined_Hop_mod.F90:1714,1829`; Hubbard U: `g=sqrt(-Dtau U/N_SUN)`, `alpha=-1/2`, type 2 (`Predefined_Int_mod.F90:70-74`) |
| `g_t(:)` | cmplx alloc | optional time-dependent coupling, size `Ltrot` (replaces `g` at slice nt) |
| `alpha` | cmplx | constant shift; enters the weight phase (`Op_phase`) |
| `type` | int | 0 hopping/no field (`Op_make` default); 1 Ising ±1; 2 discrete HS ±1,±2; 3 continuous real; 4 two fields (Re discrete ±1,±2, Im continuous). Op_set rejects other values |
| `Flip_protocol` | int | for type 4 (and copied to `nsigma`): 1 flip one of the two at random (default), 2 both, 3 real part, 4 imaginary part |
| private `M_exp(N,N,·)`, `E_exp(N,-type:type)` | | precomputed `exp(±g Phi O)` for type 1 (3 slices) and 2 (5 slices) only |
| private `U_alloc, M_exp_alloc, g_t_alloc, win_U, win_M_exp` | | bookkeeping / MPI windows |

### 4.2 `Type Fields` (`Prog/Fields_mod.F90:79-99`); global instance `nsigma` in `Hamiltonian_main:212`
| Component | Meaning |
|---|---|
| `f(N_op,Ltrot)` complex | field value per operator n (same index as `Op_V(n,:)`) and slice |
| `t(N_op)` | field type per operator (= `Op_V(n,1)%type`) |
| `Flip_protocol(N_op)` | per operator (= `Op_V(n,1)%flip_protocol`) |
Type-bound:
- `make(N_op,N_tau)` :236, `clear()` :249, `set()` :588 (random ±1, type 4: random ±1 + i·Amplitude·U(-½,½)).
- `Phi(n,nt)` :112 → complex: type 1 `f`; type 2 `Phi_st(f,2)` = ±sqrt(2(3∓√6)) for |f|=1,2; type 3 `Re f`; type 4 `Phi_st(Re f,2)·sqrt(1+Im f)`.
- `Gama(n,nt)` :143 → real: 1 for types 1,3; `1±√6/3` for types 2,4.
- `Flip(n,nt)` :173 → proposed new value: type 1 `-f`; type 2 random other of ±1,±2 (`Flip_st`); type 3 `f + Amplitude·(U-½)`; type 4 by `Flip_protocol`. `Amplitude` = `VAR_QMC` parameter, given to `Fields_init`.
- `i(n,nt)` :221 → integer `nint(Re f)` (types 1,2 only, else terminate).
- `in(Group_Comm[,Initial_field])` :337 — reads `confin_<rank>` / `confin_0.h5` (or `Temp_<g>/...`); terminates if both text and h5 exist; otherwise uses `Initial_field` or random `set()` seeded from `seeds`.
- `out(Group_Comm)` :474 — writes `confout_<rank>` (or `.h5`).
Weight per vertex: `Gama · exp(g Phi O)`; the sweep (`upgrade_mod.F90:169,205,290`) uses `Phi`, `Gama` and assigns `nsigma%f`.

---

## 5. `comm_mod` (`Libraries/Modules/comm_mod.F90`)
Thin layer: the caller passes the communicator; no `#ifdef MPI` in callers; without MPI every routine acts as one rank (rank 0, size 1, root) and reductions leave the value unchanged (:27-30).
| Routine | Line | Semantics |
|---|---|---|
| `comm_rank(comm)` | :78 | rank in comm (0 serial) |
| `comm_size(comm)` | :89 | size (1 serial) |
| `is_root(comm)` | :100 | `comm_rank(comm)==0` |
| `comm_world()` | :107 | `MPI_COMM_WORLD` handle (0 serial), for callers without `use mpi` |
| `comm_world_rank()`, `comm_world_size()` | :116, :126 | rank/size in MPI_COMM_WORLD |
| `comm_reduce_sum(x,comm)` | :53 (r8 :135, i8 :147) | MPI_SUM to rank 0, valid only on rank 0 |
| `comm_reduce_max(x,comm)` | :58 (r8 :159) | MPI_MAX to rank 0 |
| `comm_bcast(x,comm)` | :63 (r8 :171) | broadcast from rank 0 (real8 only) |
| `comm_reduce_mean(x,comm)` | :71 (r8 :181) | sum to rank 0 divided by `comm_size` |
| `comm_reduce_mean_array(buf,n,comm)` | :196 | same for a complex array of any rank (pass `size`) |
Users: `observables_mod.F90` (Print_bin_*), `control_mod.F90`, `main.F90`, `qmc_setup_mod.F90`; test `testsuite/Prog.tests/39-comm-mod.F90`.

Rule for MPI use (inferred from the comm_mod header and Phase-5 commits 5.1-5.12, e.g. 7eb97370, 58b419b1, 0708cda5; no written rule found in the repo):
- New communication goes through `comm_mod` with an explicit communicator (`Group_Comm` for the parameter group, `comm_world()` for all ranks), not raw `MPI_*` behind `#ifdef MPI`.
- If a type/operation is missing, add a specific routine/overload to `comm_mod` (same MPI type and op as the old call, no-op serial branch) and extend test 39.
- Raw MPI still remains in: generated `read_parameters` (`parse_ham_mod.py:154-213`), model `Ham_Set` info blocks, `QMC_runtime_var_mod.F90` broadcasts, `Fields_mod.F90` (`in/out`), `Operator_mod.F90` (shared-memory fences), `Langevin_HMC_mod.F90`, `Predefined_Hop_mod.F90`, `Set_random_mod.F90`, `main.F90`, `runtime_error_mod.F90` (MPI_ABORT), `mpi_shared_mem_mod_v2.F90`.
- Termination always via `Terminate_on_error(code,__FILE__,__LINE__)` (MPI-safe abort); `alf_exit(code)` only where the exit status itself is the message (e.g. 2 for an invalid `ALF_UNSTABLE_THRESHOLD`).

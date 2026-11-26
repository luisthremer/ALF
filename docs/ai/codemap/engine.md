# ALF code map: QMC engine (Prog/, written at Code-modularity af22b976)

Scope: `main.F90`, `qmc_*_mod.F90`, `Wrapgr_mod`, `upgrade_mod`, `lowrank_update_mod`, `wrapul/wrapur`, `cgr1/cgr2_2`,
`UDV_WRAP`, `udv_state`, `QDRP_decompose`, `tau_m/tau_p`, `Global_mod`, `Langevin_HMC_mod`,
`control_mod`, `Set_random_mod`. All paths relative to `Prog/`. Line numbers are for this HEAD.
Conventions: `nf` = flavour index 1..N_FL; `nf_eff` = index over computed flavours,
`nf = Calc_Fl_map(nf_eff)`; every flavour loop in the engine runs over `nf_eff = 1..N_FL_eff`.
After every Green-function/weight evaluation the engine calls `ham%weight_reconstruction` /
`ham%GR_reconstruction` when `reconstruction_needed` (fills non-computed flavours).

---------------------------------------------------------------------------------------------
## 1. Modules

### main.F90 — `Program Main` (115–418)
Driver: MPI init, RUNNING lock file (193–227), tempering/group comm split (229–266), setup,
bin/sweep loop, teardown. Owns the single `Type(qmc_walker_t) :: walker` (158). No routines.
USEs: QMC_runtime_var (getters), Hamiltonian_main (ham, Ltrot, nsigma, Group_Comm), Control,
Global_mod, Langevin_HMC_mod, qmc_setup_mod, qmc_driver_mod, qmc_sweep_mod, qmc_walker_mod, Wrapgr_mod.

### qmc_walker_mod.F90 — walker state (32–277)
Purpose: one type holding all large per-walker arrays + init/stabilisation helpers.
- `Type qmc_walker_t` (74–87): `GR(Ndim,Ndim,N_FL)`, `GR_Tilde` (same shape; symmetrised G),
  `Test(Ndim,Ndim)` (precision-check copy), `Phase` (complex scalar), `Phase_array(N_FL)`,
  `udvl(N_FL_eff)`, `udvr(N_FL_eff)` (Class UDV_State), `NSTM`, `Stab_nt(0:NSTM)`,
  `udvst(NSTM,N_FL_eff)` (storage). All ALLOCATABLE because callees declare allocatable dummies (55–57).
- `Walker_layout(walker, Ltrot, Nwrap)` (100): NSTM = ceil(Ltrot/Nwrap); Stab_nt(n)=n*Nwrap, Stab_nt(NSTM)=Ltrot.
- `Walker_alloc(walker)` (125): allocates arrays; calls `CGR_work_alloc(Ndim)`, `Tau_m_work_alloc`
  (Ltau=1, finite T) or `Tau_p_work_alloc` (Ltau=1, projector); udvl/udvr `%init` (identity or WF_L/WF_R%P);
  `udvst(NSTM,:)%reset('l')`.
- `Walker_dealloc(walker)` (157).
- `Walker_init_stack(walker)` (181): left propagation Ltrot→0 via WRAPUL, storing udvst(NST,:) for NST=NSTM-1..1.
- `Walker_init_GR(walker)` (203): CGR (NVAR=1) per flavour + `Op_phase` → Phase_array, Phase = prod^N_SUN.
- `Recompute_GR_phase(walker, NVAR)` (229): stabilisation point: Test=GR, CGR(Z1,NVAR,GR,udvr,udvl),
  `Control_PrecisionG(GR,Test)`, Op_phase, `Control_PrecisionP(Z,Phase)`, Phase=Z. (GPU stabilisation hook.)
- `Reset_udv(udv(:), side)` (258): reset to identity or trial WF (projector).
USEs: UDV_State_mod, Hamiltonian_main, Operator_mod(Op_phase), wrapul_mod, cgr1_mod, Control, Tau_m/Tau_p_mod.

### qmc_sweep_mod.F90 — sequential sweeps (32–211)
- `Sweep_up(walker)` (75): τ = 0→Ltrot; see §2. Equal-time measurement at NTAU1 ∈ [LOBS_ST, LOBS_EN].
- `Sweep_down(walker)` (135): τ = Ltrot→1; stabilisation, Tau_p mid-sweep (projector), final G on slice 0,
  Tau_m (finite T, Ltau=1) at the end (196–198), Tau_p fallback (206–208).
USEs: qmc_walker_mod, Wrapgr_mod(WRAPGRUP/WRAPGRDO), wrapur/wrapul, Tau_m/Tau_p, Hop_mod(Hop_mod_Symm), QMC_runtime_var.

### qmc_driver_mod.F90 — continuous-field global updates (32–188)
- `Langevin_update(walker)` (71): `Langevin_HMC%update(...)`, then Tau_p (NST=0) or Tau_m if Ltau=1; set_L_Forces(.true.).
- `HMC_update(walker)` (96): N_HMC_sweeps × `Langevin_HMC%update`; between trajectories calc_Forces +
  `Langevin_HMC_Reset_storage`; measurements only if `.not. sequential`.
- `MALA_update(walker, NSW)` (138): same with `Metropolis_Langevin` and N_Global_MALA_sweeps;
  invalidates cached forces if any other update ran (145–150).
USEs: Langevin_HMC_mod, Tau_m/Tau_p, QMC_runtime_var.

### qmc_setup_mod.F90 — one-time setup + info file (35–465)
- `Setup_model(N_op)` (103): `Fields_init(Amplitude)`, `Alloc_Ham(ham_name)`, `ham%Ham_set`,
  `Validate_Ham_Variables`, g_t size checks for Op_V/Op_T (116–136), `Setup_flavor_map`,
  projector WF check, `set_default_values_measuring_interval`, `Hop_mod_init`. Returns N_op=size(Op_V,1).
- `Setup_fields(N_op)` (164): `nsigma%make(N_op,Ltrot)`, sets nsigma%t / flip_protocol from Op_V(n,1),
  `Set_Random_number_Generator("seeds")`, `ham%Hamiltonian_set_nsigma`, `nsigma%in` (reads confin or random).
- `Setup_update_schemes(N_op)` (196): global-τ sampling params (default or `ham%Overide_global_tau_sampling_parameters`),
  `Wrapgr_alloc` if N_Global_tau>0 or N_Global_tau_MALA>0 or Sequential_MALA (218–220), Langevin overrides,
  `Langevin_HMC%make` / `Metropolis_Langevin%make`, `check_MALA_variables`,
  `check_update_schemes_compatibility` (locks runtime settings, QMC_runtime_var:312).
- `Write_info_parameters` (268), `Write_info_start` (401), `Write_info_end` (423), `Init_hdf5_file` (444, HDF5).

### Wrapgr_mod.F90 — slice propagation of G + local and global-τ updates (43–909)
State: `GR_ST(Ndim,Ndim,N_FL)` private module array (64), backup for MALA / multi-flip rejection;
allocated by `Wrapgr_alloc` (74) only when needed (qmc_setup:218), freed by `Wrapgr_dealloc` (79, main:388–390).
- `WRAPGRUP(...)` (85): G at slice NTAU (m=0) → NTAU+1, updates fields of slice NTAU+1. See §4.
- `WRAPGRDO(...)` (237): G at slice NTAU (m=N_op), updates fields of NTAU, → NTAU-1. See §4.
- `Wrapgr_PlaceGR(GR,m,m1,ntau)` (407): moves G inside a slice from "after vertex m" to "after vertex m1"
  (Op_Wrapup N_type 1+2 forward, Op_Wrapdo 2+1 backward).
- `Wrapgr_Random_update(GR,m,ntau,PHASE,N_Global_tau)` (477): N_Global_tau moves from
  `ham%Global_move_tau(T0_Proposal_ratio,S0_ratio,Flip_list,Flip_length,Flip_value,ntau)`;
  Upgrade2 "Intermediate" for all but last flip, "Final" for last; on rejection with Flip_length>1 restore GR_ST.
- `Wrapgr_Langevin_update(GR,m,ntau,PHASE,N_Global_tau_MALA,Delta_t,Max_Force)` (598): MALA on a
  `ham%Global_MALA_move_tau` flip list (type-3 fields only); forces via `calculate_Force`; all flips "Intermediate";
  own Metropolis accept (755), Control_upgrade/_eff.
- `wrapgr_sort` generic (59): `Wrapgr_sort_value` (778), `Wrapgr_sort_langevin` (823, also drops non-type-3).
- `Wrapgr_Test` (874): debug only.
USEs: Hamiltonian_main, Operator_mod(Op_Wrapup/Op_Wrapdo), Control, Random_Wrap, Hop_mod, upgrade_mod(Upgrade2),
Langevin_HMC_mod(Calculate_Force).

### upgrade_mod.F90 — single-vertex Metropolis update (33–307)
- `Upgrade2(GR,N_op,NT,PHASE,Hs_new,Prev_Ratiotot,S0_ratio,T0_proposal_ratio,toggle,mode)` (106). See §4.
No module state (all locals allocated per call, 155–158, 293–296). USEs (inside routine 109–115):
Hamiltonian_main, Random_wrap, Control, Fields_mod, MyMats(Det), linalg_mod.

### lowrank_update_mod.F90 — delayed and submatrix sequential update (added after af22b976)
Selected by `Sequential_scheme` = "standard" (default, nothing allocated) | "delayed" | "submatrix" and
`N_delay` in `&VAR_QMC` (QMC_runtime_var). Read the module header ("HOW IT WORKS") first.
- Pending form: delayed `G = G0 - X Y^T`; submatrix `G = G0 - G0(:,S) M (E_S^T - G0(S,:))`; flushed by one
  ZGEMM (`Lowrank_flush`) after N_delay accepts, before any non-eligible vertex and at the end of every slice.
- `Lowrank_init(scheme, N_delay, Nt_start, Nt_end, file_info)` from `Setup_update_schemes`;
  `Lowrank_finalize` from main.
- Called from WRAPGRUP/WRAPGRDO per vertex: `Lowrank_eligible(n)` (not for sequential-MALA vertices), else
  `Lowrank_flush`; `Lowrank_wrap_hook` after each Op_Wrapup/Op_Wrapdo; `Upgrade_lowrank` instead of
  Upgrade2 "Final"; `Lowrank_tune_slice_begin/_end` (timed N_delay = -1).
- Module state (private): the pending X, Y / S, M per computed flavour, eligibility, N_delay and timing.
  For a GPU port this state belongs per walker.
- Tests: `testsuite/Lowrank.tests` (36-lowrank-unit; lowrank-compare-schemes: same configurations as
  standard). Derivations: `Documentation/lowrank.tex`, `Documentation/lowrank/*.md`.

### wrapul_mod.F90 / wrapur_mod.F90 — propagation + UDV (no state)
- `WRAPUL(NTAU1, NTAU, UDVL)` (wrapul:37): left prop. $(B(\beta,\tau_1))^\dagger \to (B(\beta,\tau))^\dagger$, NTAU1>NTAU.
  Default build (112–137): `Op_mmultR(...,'c')` + `Hop_mod_mmthlc` on udvl%U, then `%decompose`, `%testscale`.
  STAB1/STAB2 (58–110): explicit TMP product + `UDV_WRAP_Pivot`.
- `WRAPUR(NTAU, NTAU1, UDVR)` (wrapur:37): right prop. $B(\tau,1)\to B(\tau_1,1)$, NTAU1>NTAU.
  Default (103–128): `Hop_mod_mmthR` + `Op_mmultR('n')` on udvr%U, `%decompose`, `%testscale`. STAB1/2: 56–102.
  Note argument order differs: WRAPUL(high, low), WRAPUR(low, high).

### cgr1_mod.F90 — equal-time G from UDV pair
> **Known issue:** the phase comes from `DET_C` (raw LU-diagonal product), which can underflow to a
> NaN phase for N >~ 3000; see `KNOWN_ISSUES.md` ("Determinant phase underflow").
State: private work arrays `TPUP, RHS (N,N), DUP, TAU (N), IPVT, VISITED (N)` (41–43), allocated once by
`CGR_work_alloc(N)` (48, from Walker_alloc), `CGR_work_dealloc` (54). Used only by the default (QR) branch.
- `CGR(PHASE,NVAR,GRUP,udvr,udvl)` (59): $G=(1+U_R D_R V_R V_L D_L U_L)^{-1}$ (in ALF's storage convention),
  PHASE = phase of det. NVAR=1: large scales in DL; NVAR=2: in DR. Dispatches to CGRP if `.not.allocated(udvl%V)` (projector).
  Branches: STAB1/2 w/o STABLOG (81–197), default/STAB3/STABLOG (199–478, uses `QDRP_decompose`).
- `CGRP(phase,GRUP,udvr,udvl)` (495): projector $G = 1 - U_R (U_L^\dagger U_R)^{-1} U_L^\dagger$.

### cgr2_2_mod.F90 — time-displaced G at a stabilisation slice
- `CGR2_2(GRT0,GR00,GRTT,GR0T,udv2,udv1,LQ)` (199): from right UDV `udv2` (0→τ) and left `udv1` (β→τ)
  builds the $2LQ\times2LQ$ block inverse; returns G00, G0T, GT0, GTT. Branch on larger D(1).
- `solve_extended_System` (97 STAB1/2 | 157 default), `get_blocks(A,B,C,D,INP,LQ)` (56). No state.

### UDV_WRAP_mod.F90 (43–240) — `UDV_Wrap_Pivot(A,U,D,V,NCON,N1,N2)` (54 STAB1 | 125 else; column-norm pivoting
+ phase fix via `Pivot_Phase`), `UDV_Wrap(A,U,D,V,NCON)` (212; QR then SVD). Used by STAB1/2 paths,
CGR2_2 STAB1/2, `Compute_Fermion_Det` (Global_mod:929).

### QDRP_decompose_mod.F90 — `Module QDRP_mod` (40)
- `QDRP_decompose(Ndim,N_part,Mat,D,IPVT,TAU,WORK,LWORK)` (61): `la_zgeqp3` pivoted QR; D = |R_ii|, R rows scaled;
  allocates WORK (caller deallocates).
- `Pivot_phase(Phase,IPVT,N_size)` (104): sign of permutation.

### udv_state_mod.F90 — `Type UDV_State` (86–111)
Fields: `U(ndim,n_part)`, `V(n_part,n_part)` (finite T only; absent in projector), `D(n_part)` or, with STABLOG,
`L(n_part)` (log scales), `ndim`, `n_part`, `side` ('l'/'r').
Methods: `alloc(t[,t_part])` (132), `init(t,side[,P])` (169), `dealloc` (295), `reset(side[,P])` (321),
`assign` / `=` (398; deep copy), `decompose` (449; QR with pivoting of U·D, V updated, STABLOG variant),
`print`, `setscale`, `getscale`, `testscale` (254; NaN abort / overflow warnings), `MPI_Sendrecv` (602, MPI).
"Projector" ⇔ `.not. allocated(%V)`.

### tau_m_mod.F90 — time-displaced measurements, finite T (43–283)
State: private work arrays `GT0,G00,GTT,G0T (Ndim,Ndim,N_FL)`, `_T` copies (Symm only), `HLP4/5/6` (61–63);
`Tau_m_work_alloc` (68) / `Tau_m_work_dealloc` (77) from Walker_alloc/dealloc.
- `TAU_M(udvst,GR,PHASE,NSTM,NWRAP,STAB_NT,LOBS_ST,LOBS_EN)` (84): needs G on slice 0 and udvst full of
  LEFT propagations. Propagates GT0 (PROPR), G0T (PROPRM1), GTT (both, or `Langevin_HMC%Wrap_Forces` for
  Langevin/HMC) slice by slice, calls `ham%OBSERT(NT,GT0,G0T,G00,GTT,PHASE,Mc_step_weight)` for NT=0..Ltrot;
  at Stab_nt recomputes via local udvr + WRAPUR + `CGR2_2(...,udvr,udvst(NST,:))` and `Control_Precision_tau`.
  For Langevin/HMC also does `ham%Obser` on GTT (188, 199). Mc_step_weight = Delta_t_running for Langevin.
- `PROPR(AIN,NT)` (232): $A\leftarrow B(NT)A$; `PROPRM1(AIN,NT)` (257): $A\leftarrow A B^{-1}(NT)$. Public (used by tau_p).

### tau_p_mod.F90 — time-displaced measurements, projector (42–361)
State: private `GRUPB,GRUP,G00UP,G0TUP,GT0UP,GTTUP (LQ,LQ,N_FL)`, `_T` copies, `TEMP` (60–63);
`Tau_p_work_alloc` (68) / `Tau_p_work_dealloc` (79).
- `Tau_p(udvl,udvr,udvst,GR,PHASE,NSTM,STAB_NT,NST_IN,LOBS_ST,LOBS_EN)` (104): input G on Stab_nt(NST_IN) ≤ Thtrot+1,
  storage full of left props for Stab_nt > that. Propagates to Thtrot+1, then measures OBSERT for
  τ ∈ [0, Ltrot-2·Thtrot] with CGRP checks at Stab_nt (217–289); NCHECK fallback (300–317); Langevin/HMC tail (319–352).

### Global_mod.F90 — global (space-time) updates, parallel tempering (45–1044; everything public)
State: private `Tempering_acceptance` (Obser_Vec, 62).
- `Global_Updates(Phase,GR,udvr,udvl,Stab_nt,udvst,N_Global)` (453): precondition/postcondition: storage full of
  left props, G on slice 0. Old weight via `Compute_Fermion_Det(...,"Full")`; N_Global proposals
  `ham%Global_move_log_T0(log_T0,nsigma_old,size_clust)`; new weight `Compute_Fermion_Det(...,"Empty")` (refills udvst);
  `Compute_Ratio_Global`; accept/restore nsigma; at end rebuilds storage if last move rejected (600–621) and recomputes
  G/Phase with CGR (622–633). Counters: `Control_upgrade_Glob`, `Control_PrecisionP_Glob`.
- `Compute_Ratio_Global(...)` (675, function): $T_0\,W(\sigma')/W(\sigma)$ = Ratio(1)·exp(Ratio(2)), incl.
  Gama, $e^{N_{SUN}\, g\,\alpha\,\Delta\phi}$, `ham%Get_Delta_S0_global`.
- `Compute_Fermion_Det(Phase_det,Det_Vec,udvl,udvst,Stab_nt,storage)` (795): log|det| per flavour as Det_Vec(:,nf_eff)
  (+phase); storage="Empty" recomputes udvl/udvst from scratch via WRAPUL.
- `Exchange_Step(Phase,GR,udvr,udvl,Stab_nt,udvst,N_exchange_steps,Tempering_calc_det)` (111, TEMPERING only): replica
  exchange; swaps nsigma via MPI; on accept recomputes storage/G (or MPI-swaps GR/udv*).
- `Global_Tempering_setup` (994), `_init_obs` (1010), `_obser` (1022), `_Pr` (1038); `npbc_tempering` (function ending :984).
USEs: Hamiltonian_main, MyMats(DET_C, DET_C_LU), Operator_mod(Op_phase), Control, Observables, Fields_mod, wrapul, cgr1.

### Langevin_HMC_mod.F90 — Langevin / HMC / global MALA (36–1052)
- `Type Langevin_HMC_type` (69–93, private components): `scheme` (enum none/Langevin/HMC/MALA, 63–68),
  `Update_scheme` (string), `L_Forces` (cached forces valid), `Delta_t_running`, `Delta_t_Langevin_HMC`, `Max_Force`,
  `Leapfrog_Steps`, `Det_vec_old`, `Phase_Det_old`, `Forces(N_op,Ltrot)` (complex), `Forces_0(N_op,Ltrot)` (real).
- Module instances: `Langevin_HMC` (Langevin or HMC), `Metropolis_Langevin` (global MALA) (95).
- Methods: `make`=`Langevin_HMC_setup` (742), `clean` (862), `Wrap_Forces`=`Wrapgrup_Forces(this,Gr,NT1)` (201;
  propagates G NT→NT1 like WRAPGRUP without updates and stores forces of slice NT1), `Update`=`Langevin_HMC_update`
  (345; select case scheme), `calc_Forces`=`Langevin_HMC_Forces(this,Phase,GR,GR_Tilde,Test,udvr,udvl,Stab_nt,udvst,
  LOBS_ST,LOBS_EN,Calc_Obser_eq)` (114; full up-sweep computing forces, stabilising, optionally measuring),
  `set_L_Forces` (919), `set_det_OLD` (934), `get/set_Update_scheme` (952/970), `get_Delta_t_running` (1001).
- Public free routines: `Langevin_HMC_Reset_storage(Phase,GR,udvr,udvl,Stab_nt,udvst)` (267; refill storage with left
  props, G on slice 0), `calculate_Force(n,ntau,gr)` (1016; force of one type-3 vertex; used by Wrapgr_mod MALA).
USEs: Hamiltonian_main, UDV_State_mod, Control, Hop_mod, Operator_mod, wrapur/wrapul, cgr1; Global_mod (inside update).

### control_mod.F90 — `module Control` (46–719): statistics & sanity checks
State: private SAVE counters (53–85): G precision (XMEANG, XMAXG, NCG), tau precision, phase precision
(normal/Glob/HMC/MALA), acceptance counters NC_up/ACC_up, NC_eff_up/ACC_eff_up, Glob, HMC, Temp, MALA, force stats,
CPU clock.
- `control_init(Group_Comm)` (91; called twice: main:301 and main:324), `Control_upgrade(toggle)` (288),
  `Control_upgrade_eff(toggle)` (295), `_Temp` (302), `_Glob(toggle,size_clust)` (309), `_HMC` (321), `_MALA` (330),
  `Control_PrecisionG(A,B,Ndim)` (aborts if max|ΔG| > G_thresh or NaN; G_thresh = 10, or env `ALF_UNSTABLE_THRESHOLD` > 0 read by `control_read_threshold` in `control_init`), `Control_Precision_tau` (433),
  `Control_PrecisionP*(Z,Z1)` (447/456/466/476), `Control_Langevin` (171), `Control_MALA_sequential` (201),
  `Control_MALA_Global` (221), `Control_MALA_Global_tau` (258), `Control_Print(...)` (487; writes info),
  `make_truncation(prog_truncation,cpu_max,count_bin_start,count_bin_end,group_comm)` (680; stop if <1.5 bins left).

### Set_random_mod.F90 — `module set_random` (33)
- `Set_Random_number_Generator(File_seeds,Seed_in)` (48): rank 0 reads `seeds` file, sends one seed per rank, `RANSET`.

---------------------------------------------------------------------------------------------
## 2. Call graph of one run

```
main (main.F90)
├─ read_and_broadcast_QMC_var_and_ham_name            main:281  (QMC_runtime_var:917; reads VAR_QMC, VAR_HAM_NAME)
├─ qmc_setup_mod:Setup_model(N_op)                    main:284  → ham%Ham_set, Setup_flavor_map, Hop_mod_init
├─ apply_CPU_MAX_policy                               main:287
├─ qmc_setup_mod:Setup_fields(N_op)                   main:289  → nsigma%make/%in, Set_Random_number_Generator
├─ Init_hdf5_file (HDF5)                              main:297
├─ Control:control_init; ham%Alloc_obs(Ltau)          main:301–302
├─ qmc_walker_mod:Walker_layout                       main:304
├─ qmc_setup_mod:Setup_update_schemes(N_op)           main:306  → Wrapgr_alloc?, Langevin_HMC%make, lock settings
├─ Walker_alloc → Walker_init_stack (WRAPUL Ltrot→0) → Walker_init_GR (CGR)   main:311–315
├─ Control:control_init (again)                       main:324
└─ DO NBC=1,NBin                                       main:326
   ├─ ham%Init_obs                                    main:332
   ├─ DO NSW=1,NSweep                                 main:337
   │  ├─ [TEMPERING] Global_mod:Exchange_Step          main:342  (every N_Tempering_frequency sweeps)
   │  ├─ [Global_moves] Global_mod:Global_Updates      main:346
   │  ├─ [scheme=LANGEVIN] qmc_driver:Langevin_update  main:350 → Langevin_HMC%update (Langevin_HMC_update case Langevin)
   │  │        → calc_Forces (Langevin_HMC_Forces: Wrap_Forces, WRAPUR, CGR, ham%Obser) → field step → Reset_storage
   │  │        → Tau_m / Tau_p if Ltau=1
   │  ├─ [scheme=HMC]  qmc_driver:HMC_update           main:354 → leapfrog with ham%Apply_B_HMC, Compute_Fermion_Det
   │  ├─ [MALA]        qmc_driver:MALA_update          main:358 → Metropolis_Langevin%update (case MALA, ham%Global_MALA_move)
   │  └─ [sequential]  qmc_sweep:Sweep_up; Sweep_down  main:361–365
   ├─ ham%Pr_obs; [TEMPERING] Global_Tempering_Pr     main:368–371
   ├─ nsigma%out(Group_Comm)   (writes confout)       main:373
   └─ make_truncation (CPU_MAX)                        main:378
teardown: Walker_dealloc, Wrapgr_dealloc?, Clear_Ham_state, Control_Print, Write_info_end   main:387–406
```

Sequential sweep detail:
```
Sweep_up (qmc_sweep:75)
  Reset_udv(udvr,'r')
  DO NTAU=0,Ltrot-1
    WRAPGRUP(GR,NTAU,PHASE,...)                               qmc_sweep:87
      Hop_mod_mmthr + Hop_mod_mmthl_m1 (kinetic B_T(ntau1) wrap)    Wrapgr:127–128
      DO n=Nt_sequential_start..end:  Op_Wrapup(N_type=1) → [MALA prop.] → Upgrade2 → Op_Wrapup(N_type=2) [→ MALA accept]
      [N_Global_tau>0]      Wrapgr_Random_update  → ham%Global_move_tau, Wrapgr_PlaceGR, Upgrade2
      [N_Global_tau_MALA>0] Wrapgr_Langevin_update → ham%Global_MALA_move_tau, calculate_Force, Upgrade2
      Wrapgr_PlaceGR(GR,m,N_op,ntau1)  (bring G to m=N_op)
    IF NTAU1==Stab_nt(NST):  WRAPUR(Stab_nt(NST-1),NTAU1,udvr); udvl←udvst(NST); udvst(NST)←udvr
                             Recompute_GR_phase(walker, NVAR) (NVAR=2 if NTAU1>Ltrot/2) → CGR, Control_PrecisionG/P
    IF NTAU1∈[LOBS_ST,LOBS_EN]: [Symm: Hop_mod_Symm→GR_Tilde]; [ham%GR_reconstruction]; ham%Obser(G,PHASE,NTAU1,1.0)
Sweep_down (qmc_sweep:135)
  Reset_udv(udvl,'l')
  DO NTAU=Ltrot,1,-1
    WRAPGRDO(GR,NTAU,PHASE,...)                               qmc_sweep:145
      [global-τ moves first at m=N_op] → Wrapgr_PlaceGR(GR,N_op,m1)
      DO n=Nt_sequential_end..start,-1: Op_Wrapdo(N_type=2) → Upgrade2 → Op_Wrapdo(N_type=1, new field)
      Hop_mod_mmthl + Hop_mod_mmthr_m1                                  Wrapgr:399–400
    IF NTAU1∈[LOBS_ST,LOBS_EN]: ham%Obser(...)     (measurement BEFORE stabilisation on the way down)
    IF Stab_nt(NST)==NTAU1≠0: WRAPUL(Stab_nt(NST+1),NTAU1,udvl); udvr←udvst(NST); udvst(NST)←udvl
                              Recompute_GR_phase; [projector, Ltau=1, Stab_nt(NST)≤Thtrot+1<Stab_nt(NST+1)] Tau_p
  WRAPUL(Stab_nt(1),0,udvl); Reset_udv(udvr,'r'); Recompute_GR_phase(walker,1); Reset_udv(udvst(NSTM,:),'l')
  [Ltau=1, finite T] TAU_M(udvst,GR,PHASE,NSTM,Nwrap,Stab_nt,LOBS_ST,LOBS_EN) → PROPR/PROPRM1, WRAPUR, CGR2_2, ham%OBSERT
  [Ltau=1, projector, Stab_nt(1)>Thtrot+1] Tau_p(...,NST_IN=0)
```
Invariant between sweeps: G on slice 0, udvst holds LEFT propagations (β→Stab_nt(n)), udvl on slice 0.
After Sweep_up the storage holds RIGHT propagations (0→Stab_nt(n)); Sweep_down swaps back.
Global_Updates, Langevin/HMC/MALA, Exchange_Step all assume and restore the "G on slice 0, storage left" invariant.

---------------------------------------------------------------------------------------------
## 3. Where data lives

| Data | Owner / declaration | Written by | Read by |
|---|---|---|---|
| `GR(Ndim,Ndim,N_FL)` | walker (qmc_walker:76) | Walker_init_GR, Recompute_GR_phase (CGR), WRAPGRUP/WRAPGRDO (Hop_mod, Op_Wrapup/do, Upgrade2), Wrapgr_PlaceGR, Global_Updates (622–630), Exchange_Step, Langevin_HMC_Forces/Reset_storage/Wrap_Forces | ham%Obser (via qmc_sweep), Tau_m/Tau_p (as G00 input), calculate_Force |
| `GR_Tilde`, `Test` | walker (76–77) | Hop_mod_Symm (qmc_sweep:114,154; Langevin_HMC:178); Recompute_GR_phase / Langevin_HMC_Forces | ham%Obser; Control_PrecisionG |
| `GR_ST` | Wrapgr_mod private (64) | WRAPGRUP/DO MALA path (144,316), Wrapgr_Random_update (546), Wrapgr_Langevin_update (671) | rejection restore |
| `Phase`, `Phase_array` | walker (79–80) | Walker_init_GR, Recompute_GR_phase, Upgrade2 (Final+accept, upgrade:227), MALA paths in Wrapgr (188, 339, 720), Global_Updates, Langevin_HMC_* | ham%Obser/OBSERT, Upgrade2 weight, Control_PrecisionP |
| `udvl`, `udvr` (N_FL_eff) | walker (82) | Reset_udv, WRAPUL/WRAPUR, sweeps (copy from udvst), Compute_Fermion_Det, Reset_storage | CGR, Tau_p |
| `udvst(NSTM,N_FL_eff)`, `Stab_nt(0:NSTM)`, `NSTM` | walker (84–86) | Walker_layout (Stab_nt), Walker_init_stack, Sweep_up (right props), Sweep_down (left props), Global_Updates/Compute_Fermion_Det("Empty"), Langevin_HMC_Reset_storage, Exchange_Step | sweeps, TAU_M (left), Tau_p (left), Compute_Fermion_Det("Full") |
| `nsigma` (Type Fields: `f(N_op,Ltrot)` complex field, `t(N_op)` type 1..4, `flip_protocol`; methods Phi, Gama, Flip, in, out) | Hamiltonian_main:212 (public, not protected) | Setup_fields (make/in), Upgrade2 on accept (upgrade:290), MALA rejection restore in Wrapgr (207, 370, 763), Wrapgr_Random_update restore (584), ham%Global_move_tau / Global_MALA_move_tau / Global_move_log_T0 / Global_MALA_move (model code), Global_Updates restore, Langevin/HMC field steps (Langevin_HMC:408), Exchange_Step (MPI swap); written to disk by `nsigma%out` (main:373) | everything that builds B matrices: WRAPUL/WRAPUR, Op_Wrapup/do callers, PROPR/PROPRM1, Op_phase, ham%S0, Compute_Ratio_Global |
| `Op_V(N_op,N_FL)` (Type Operator: P, O, U, E, E_Exp, g, g_t, alpha, type, N_non_zero, diag) | Hamiltonian_main:206 (protected except under nvfortran, :141–145) | ham%Ham_set (model) | everywhere: Upgrade2, WRAPGR*, wraps, forces, Op_phase |
| `Op_T(N_ops,N_FL)` | Hamiltonian_main:207 | ham%Ham_set | Hop_mod_init (kinetic exponentials); engine only checks g_t size (qmc_setup:127) |
| `WF_L/WF_R(N_FL)` (%P: Ndim×N_part) | Hamiltonian_main:208–209 | ham%Ham_set | Walker_alloc, Reset_udv, Global_mod, Langevin_HMC (projector only) |
| `Ndim, N_FL, N_FL_eff, N_SUN, Ltrot, Thtrot, Projector, Symm, reconstruction_needed, Calc_Fl_map, Group_Comm` | Hamiltonian_main:210–221 (protected) | Ham_set / Setup_flavor_map / Set_Group_Comm | all |
| Runtime parameters | `Module QMC_runtime_var` (QMC_runtime_var_mod.F90:32; private vars, `get_*`/`set_*`; namelists VAR_QMC 80, VAR_TEMP 77, VAR_HAM_NAME 89) | read: `read_and_broadcast_QMC_var_and_ham_name` (917); modified by `apply_CPU_MAX_policy` (322, NBin), `set_default_values_measuring_interval`, `set_default_global_tau_sampling` (353), setters in Setup_update_schemes (213–216), `apply_Langevin_overrides` (333), `check_MALA_variables`; frozen by `lock_QMC_runtime_settings` at end of `check_update_schemes_compatibility` (312); later `set_*` abort via `ensure_settings_unlocked` | main, qmc_sweep (all WRAPGR* args), qmc_driver, qmc_setup, qmc_walker (get_Ltau) |
| Kernel work arrays | cgr1 (TPUP…), tau_m (G00…), tau_p (GRUP…) module privates | *_work_alloc from Walker_alloc | own routine only |
| Acceptance / precision stats | Control module SAVE vars | Control_* | Control_Print |
| Langevin/HMC forces | `Langevin_HMC%Forces/Forces_0`, `Metropolis_Langevin%…` | Wrapgrup_Forces, ham%Ham_Langevin_HMC_S0 | update step, tau_m/tau_p via Wrap_Forces |

Note: sequential parameters reach WRAPGRUP/WRAPGRDO only as arguments (qmc_sweep:87–90, 145–148); Wrapgr_mod
itself does not USE QMC_runtime_var.

---------------------------------------------------------------------------------------------
## 4. Hook points for a new local-update scheme (delayed / submatrix updates)

### Signatures (verbatim)
```fortran
! Wrapgr_mod.F90:85
SUBROUTINE WRAPGRUP(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
           &        Sequential_MALA, Delta_t_MALA_sequential, Max_Force_MALA_sequential, &
           &        N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau)
  COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  GR(:,:,:)
  COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
  INTEGER, INTENT(IN) :: NTAU                       ! 0..Ltrot-1; updates slice NTAU1=NTAU+1
  LOGICAL, INTENT(IN) :: Propose_S0, Sequential_MALA
  INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau, N_Global_tau_MALA
  Real    (kind=Kind(0.d0)), intent(in) :: Delta_t_MALA_sequential, Delta_t_MALA_global_tau
  real    (kind=kind(0.d0)), intent(in) :: Max_Force_MALA_sequential, Max_Force_MALA_global_tau

! Wrapgr_mod.F90:237  -- identical argument list; NTAU in Ltrot..1, updates slice NTAU
SUBROUTINE WRAPGRDO(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
           &        Sequential_MALA, Delta_t_MALA_sequential, Max_Force_MALA_sequential, &
           &        N_Global_tau_MALA, Delta_t_MALA_global_tau, Max_Force_MALA_global_tau)

! upgrade_mod.F90:106
Subroutine Upgrade2(GR,N_op,NT,PHASE,Hs_new, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, toggle,  mode)
  Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim, N_FL)   ! explicit shape
  Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Prev_Ratiotot
  Integer                  , INTENT(IN)    :: N_op, Nt
  Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Phase
  Complex (Kind=Kind(0.d0)), INTENT(IN)    :: Hs_new
  Real    (Kind=Kind(0.d0)), INTENT(IN)    :: S0_ratio, T0_proposal_ratio
  Character (Len=64)       , INTENT(IN)    :: Mode                   ! "Final" | "Intermediate" (case-insensitive)
  Logical                  , INTENT(INOUT) :: toggle                 ! .true. if accepted
```
Callers: qmc_sweep:87 / :145 only (both pass all 13 args from QMC_runtime_var getters).
Upgrade2 callers: Wrapgr:176 (UP, both modes), :338 (DO, MALA "Intermediate"), :380 (DO "Final"),
:557/:568 (Random_update), :711 (Langevin_update).

### What Upgrade2 does (upgrade:106–305), per call
1. op_dim = max_nf Op_V(N_op,nf)%N_non_zero; allocates Mat(op_dim,op_dim), Delta(op_dim,N_FL_eff), u,v,x_v,y_v,xp_v(Ndim,op_dim).
2. Per flavour: $Z_1 = g(nt)\,[\phi(\sigma')-\phi(\sigma)]$, $\Delta_m = e^{Z_1 E_m}-1$,
   $M_{nm} = \delta_{nm}e^{Z_1E_m} - \Delta_m\, G_{P_n P_m}$ (indices via Op_V%P), ratio = det M (closed forms for dim 1,2) ·
   $e^{Z_1\alpha}$. Ratiotot = $(\prod_{nf} \text{Ratio})^{N_{SUN}}\cdot\gamma(\sigma')/\gamma(\sigma)$.
3. "Final": Ratiotot *= Prev_Ratiotot, Weight = S0_ratio·T0_proposal_ratio·|Re(Phase·Ratiotot)/Re(Phase)|.
   "Intermediate": Prev_Ratiotot *= Ratiotot, Weight = 1.5 (always accepted).
4. If Weight > ranf_wrap(): Phase updated (Final only), rank-op_dim Sherman–Morrison/Woodbury update of
   GR(:,:,nf) via la_zgemv/la_zgeru/la_zgemm (230–287), `nsigma%f(N_op,NT) = Hs_new` (290).
5. Final only: Control_upgrade(toggle), Control_upgrade_eff(toggle) (298–301).
Immediate full Ndim×Ndim update per accepted flip ⇒ this is what a delayed/submatrix scheme replaces.

### Structure of the per-slice loop (where a new scheme plugs in)
WRAPGRUP, slice NTAU1 = NTAU+1 (Wrapgr:124–231):
```
GR ← B_T(ntau1) GR B_T(ntau1)^{-1}                     (Hop_mod_mmthr, Hop_mod_mmthl_m1; 125–129)
DO n = Nt_sequential_start, Nt_sequential_end           (130)
   L_MALA = Sequential_MALA .and. Op_V(n,1)%type == 3
   Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,1,ntau1)   ! N_type=1: apply vertex n (old field), go to its eigenbasis (136)
   [L_MALA] save nsigma_st,phase_st,Gr_st=Gr; Gr_st rotated back (N_type=2); force_old=calculate_force; Langevin proposal (141–159)
   [else]   Hs_New = nsigma%flip(n,ntau1)                                                         (161)
   S0_ratio = ham%S0(n,ntau1,Hs_New); Propose_S0 (type 1) sets T0_proposal/T0_Proposal_ratio (163–167)
   IF T0_proposal > ranf: Upgrade2(GR,n,ntau1,PHASE,HS_new,Prev_Ratiotot,S0_ratio,T0_Proposal_ratio,Acc,mode)
         mode = "Intermediate" if L_MALA else "Final"                                            (168–176)
   ELSE Control_upgrade_eff(.false.)                                                             (178–179)
   Op_Wrapup(...,N_type=2)                                 ! back to site basis (181–185)
   [L_MALA] phase update, force_new, proposal ratio, own Metropolis; reject ⇒ nsigma,Gr,phase restored (187–215)
ENDDO
m = Nt_sequential_end (or 0)                                (218–222)
[N_Global_tau>0]      Wrapgr_Random_update(GR,m,ntau1,PHASE,N_Global_tau)                          (223–226)
[N_Global_tau_MALA>0] Wrapgr_Langevin_update(GR,m,ntau1,PHASE,N_Global_tau_MALA,Delta_t,Max_Force) (227–230)
[either>0]            Wrapgr_PlaceGR(GR,m,Size(OP_V,1),ntau1)                                      (231)
```
WRAPGRDO, slice NTAU (Wrapgr:277–401): ORDER REVERSED — global-τ moves first:
```
m = Size(OP_V,1)
[N_Global_tau>0] Wrapgr_Random_update(GR,m,ntau,...); [N_Global_tau_MALA>0] Wrapgr_Langevin_update(GR,m,ntau,...)  (278–286)
m1 = Nt_sequential_end (or 0); [either>0] Wrapgr_PlaceGR(GR,m,m1,ntau)                         (289–294)
DO n = Nt_sequential_end, Nt_sequential_start, -1                                               (297)
   [L_MALA] nsigma_st, phase_st, force_old = calculate_force(n,ntau,Gr)  (G here is in site basis, after vertex n)
   Op_Wrapdo(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,2,ntau)    ! N_type=2: to eigenbasis of vertex n (309)
   proposal (MALA 315–325 | flip 327), S0_ratio, Propose_S0
   IF accepted-proposal:
      [L_MALA] Upgrade2(...,"Intermediate") (338); Op_Wrapup N_type=2 to site basis; force_new; MALA Metropolis;
               accept ⇒ Op_Wrapdo N_type=2 again; reject ⇒ restore Gr_st, nsigma, phase (336–376)
      [else]   Upgrade2(...,"Final") (380)
   ELSE Control_upgrade_eff(.false.)
   HS_Field = nsigma%f(n,ntau)  (new field); Op_Wrapdo(...,N_type=1,ntau)   ! remove vertex n (389–395)
ENDDO
GR ← B_T(ntau)^{-1} GR B_T(ntau)                         (Hop_mod_mmthl, Hop_mod_mmthr_m1; 397–401)
```

### Facts relevant to a delayed / submatrix replacement
- Basis: `Op_Wrapup/Op_Wrapdo` (Operator_mod:746/865) with N_type=1 multiply by $e^{\pm g\phi E}$ and $U^{(\dagger)}$ on rows/cols
  `Op%P`; N_type=2 is the bare $U$ rotation and is a NO-OP when `Op%diag` (Operator_mod:808, 833, 924, 949). Upgrade2 reads
  `GR(P(n),P(m),nf)` in the eigenbasis of vertex n. For non-diagonal vertices G must be rotated between vertices
  (O(Ndim·Op%N) each) — a delayed scheme that leaves G stale must account for these rotations; for `Op%diag`
  vertices (e.g. Hubbard density coupling) only diagonal rescalings of rows/cols P occur.
- Inner loop range: `n = Nt_sequential_start..Nt_sequential_end` (default 1..N_op, set by
  `set_default_global_tau_sampling`, QMC_runtime_var:353). Global-τ moves operate on the same G in the same slice
  and need G "after vertex m" (Wrapgr_PlaceGR semantics) ⇒ any pending delayed updates must be flushed before
  Wrapgr_Random_update / Wrapgr_Langevin_update / Wrapgr_PlaceGR and before Hop_mod wraps at slice end.
- Sequential MALA (type-3 fields) needs exact G immediately (calculate_force on GR, 190/348) and a full G backup
  `Gr_st = Gr` (144/316) ⇒ must flush before force evaluation; or restrict the new scheme to `.not. L_MALA` vertices.
- Upgrade2 "Intermediate" mode (multi-vertex moves) updates G and nsigma unconditionally and accumulates
  Prev_Ratiotot; restoration on reject is done by the caller via GR_ST.
- Phase: updated inside Upgrade2 only for "Final"+accept (upgrade:227): $\text{Phase} \mathrel{*}= R/|R|$.
- Bookkeeping to keep: Control_upgrade + Control_upgrade_eff on every Final decision; Control_upgrade_eff(.false.)
  when proposal is skipped (Wrapgr:178–179, 383–384).
- RNG order: proposal (`nsigma%flip` / `rang_wrap`), then `ranf_wrap` for T0_proposal (168/334), then `ranf_wrap`
  inside Upgrade2 (225). Changing call order changes the Markov chain (bitwise comparisons with reference runs break).
- G on entry/exit of WRAPGRUP: exact equal-time G at slice NTAU / NTAU1 in site basis; callers immediately use it for
  stabilisation compare (`Recompute_GR_phase` → `Control_PrecisionG`, abort if >10) and `ham%Obser` ⇒ G must be fully
  up to date when WRAPGRUP/WRAPGRDO return.
- Flavours: loops over `nf_eff`; Upgrade2 calls `ham%weight_reconstruction(Ratio)` (upgrade:200) for non-computed flavours.
- Upgrade2 allocates/deallocates 7 small arrays per call (155–158, 256/267 inside loop, 274/280) — not per-run work arrays.
- Natural insertion point: replace the `DO n` loop bodies (Wrapgr:130–216 and 297–396) or the `Call Upgrade2(...)`
  sites at 176 / 380 with a new routine taking `GR, n, ntau, PHASE, HS_new, S0_ratio, T0_Proposal_ratio, Acc`, plus a
  flush call before line 218 (UP) / after 396 (DO) and before global-τ calls; selector would be a new VAR_QMC flag in
  QMC_runtime_var passed through qmc_sweep:87/145 (Wrapgr_mod has no QMC_runtime_var USE).

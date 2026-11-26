# ALF code map — index

Written once from a full read of the code at `Code-modularity` `af22b976` (2026-11-25), so that later
sessions need not read the whole code. Line numbers (`file:N`) are for that commit; after later
changes, search for the routine name instead.

| Part | File | Covers |
|---|---|---|
| QMC engine | [`codemap/engine.md`](codemap/engine.md) | `main`, `qmc_walker/sweep/driver/setup_mod`, `Wrapgr_mod`, `upgrade_mod`, wraps, CGR, UDV, `tau_m/tau_p`, global moves, Langevin/HMC/MALA, `control_mod`; call graph of one run; who owns which data; where a local-update scheme plugs in |
| Models and libraries | [`codemap/model.md`](codemap/model.md) | `Hamiltonian_main` (protected state, `ham_base` bindings), the Hamiltonian submodules and their build glue, `Operator`, `Fields`, hopping types, `Predefined_*`, observables, `QMC_runtime_var`, `Libraries/Modules` (`comm_mod`, `linalg_mod`, …); how to add a model or a run-time parameter |
| Build, tests, CI | [`codemap/build.md`](codemap/build.md) | `configure.sh` switches, make targets, ctest, the "bitwise identical" procedure, every CI job, the coding conventions with their commit evidence, commit style |

## The run in ten lines

```
main.F90
 ├─ Setup_model / Setup_fields / Setup_update_schemes  (qmc_setup_mod)   model, fields, update flags
 ├─ Walker_alloc, Walker_init_stack, Walker_init_GR     (qmc_walker_mod) G and UDV stacks at slice 0
 └─ loop bins → loop sweeps
     ├─ Sweep_up   (qmc_sweep_mod): for ntau = 0..Ltrot-1
     │    WRAPGRUP (Wrapgr_mod): hopping wrap, then per vertex n: Op_Wrapup, propose,
     │             Upgrade2 (or Upgrade_lowrank: delayed/submatrix), Op_Wrapup; then global-τ moves
     │    at stabilisation slices: WRAPUR + Recompute_GR_phase (CGR, precision checks)
     │    equal-time measurements (ham%Obser)
     ├─ Sweep_down (qmc_sweep_mod): ntau = Ltrot..1, WRAPGRDO, WRAPUL, time-displaced (TAU_M/Tau_p)
     └─ Langevin / HMC / MALA (qmc_driver_mod), global moves and tempering (Global_mod)
```

## Not covered here

- `Analysis/` (ana.out, maxent) — read its sources when needed.
- `Prog/Hamiltonians/Not_supported/` — not built.
- GPU code — lives on the GPU branches (see [`phase2/ROADMAP.md`](phase2/ROADMAP.md)).

## How to refresh

When a phase changes many modules: let three readers (engine / models+libraries / build+CI) each
read their part in full and rewrite their file; keep the same section headings; update the commit
in the first line of each file and in this index.

# Hooks of the AI features in baseline files

*`master_ai` only. Written with AI assistance; not part of the ALF baseline.*

A hook is the only permitted edit to a file that exists on `master`. Each hook is listed
below, under the file it belongs to, as the exact block of added lines. `testsuite/hook_lint.py`
checks that the diff of every baseline file against `master` consists of these blocks and nothing
else: no removed or changed baseline lines, no unlisted additions.

Rules (draft `AI_DEVELOPMENT_RULES.md`, R2):
- A Fortran hook is at most three lines: `#ifdef ALF_<FEAT>`, one statement, `#endif`. The statement is a
  `use` of the feature's module or a call into it, typically "dispatch and return".
- Build files (`configure.sh`, Makefiles) and `Documentation/doc.tex` may have longer blocks; they are
  listed here in the same way.

## Feature LOWRANK (delayed / submatrix sequential update)

Code: `Prog/Lowrank/`. Tests: `testsuite/Lowrank.tests/`. Switch: `configure.sh ... LOWRANK`.

### Prog/Wrapgr_mod.F90
Module `Wrapgr_mod`: makes the dispatch functions visible.
```
#ifdef ALF_LOWRANK
  use Lowrank_wrapgr_mod, only: Lowrank_wrapgrup, Lowrank_wrapgrdo
#endif
```

### Prog/Wrapgr_mod.F90
`WRAPGRUP`, first executable statement: the low-rank scheme handles the whole slice; with
`Sequential_scheme = "standard"` the function returns `.false.` and the baseline code runs.
```
#ifdef ALF_LOWRANK
    if (Lowrank_wrapgrup(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start,Nt_sequential_end,N_Global_tau,Sequential_MALA,Delta_t_MALA_sequential,Max_Force_MALA_sequential,N_Global_tau_MALA,Delta_t_MALA_global_tau,Max_Force_MALA_global_tau,Wrapgr_Random_update,Wrapgr_Langevin_update,Wrapgr_PlaceGR)) return
#endif
```

### Prog/Wrapgr_mod.F90
`WRAPGRDO`, first executable statement: as for `WRAPGRUP`.
```
#ifdef ALF_LOWRANK
    if (Lowrank_wrapgrdo(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start,Nt_sequential_end,N_Global_tau,Sequential_MALA,Delta_t_MALA_sequential,Max_Force_MALA_sequential,N_Global_tau_MALA,Delta_t_MALA_global_tau,Max_Force_MALA_global_tau,Wrapgr_Random_update,Wrapgr_Langevin_update,Wrapgr_PlaceGR)) return
#endif
```

### Prog/Makefile
Adds the feature's sources when the switch is set.
```
ifneq (,$(findstring -DALF_LOWRANK,$(ALF_FLAGS_PROG)))
SRCS += Lowrank/Lowrank_update_mod.F90 Lowrank/Lowrank_wrapgr_mod.F90
endif
```

### configure.sh
Help text.
```
  LOWRANK: Compile the delayed/submatrix sequential update (master_ai, AI-written)
```

### configure.sh
Reset of the switch (the script is sourced, so a value from an earlier call must not survive).
```
LOWRANK_ENABLED=""
```

### configure.sh
Keyword.
```
    LOWRANK)
      LOWRANK_ENABLED="1"
    ;;
```

### configure.sh
Flag.
```
if [ "${LOWRANK_ENABLED}" = "1" ]; then
  ALF_FLAGS_PROG="${ALF_FLAGS_PROG} -DALF_LOWRANK"
fi
```

### Documentation/doc.tex
Appendix page of the feature.
```
\input{lowrank}
```

### testsuite/CMakeLists.txt
Tests of the feature, only for builds with the switch.
```
if("${FLAG_LIST}" MATCHES "-DALF_LOWRANK")
  add_subdirectory(Lowrank.tests)
endif()
```

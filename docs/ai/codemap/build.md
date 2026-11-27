# ALF codemap: build, tests, CI, conventions (written at Code-modularity af22b976)

Scope: Code-modularity refactor (origin ALF-QMC/ALF PR #656 = Phases 0–5, PR #657 = Phase 6;
Phase 6G commits on top). Base master commit `4c690b04`; 149 commits `4c690b04..HEAD`.
Verified 2026-11-25 by reading files; "unverified" marks anything not run.

## 1. Build and run locally (macOS, Homebrew gfortran)

### Environment (checked on a Mac, 2026-11-25)
- gfortran 16.2.0 (Homebrew GCC), `mpifort`/`mpiexec` = MPICH (Hydra 5.0.1), cmake 4.4.3, python3.
- BLAS/LAPACK: `-llapack -lblas` (macOS system libs). `configure.sh` compiles and runs
  `check_libs.f90` and fails if they don't link/run.
- `. ./configure.sh GNU MPI NO-INTERACTIVE` returned 0 here (ALF_FC=mpifort).
- Set `TMPDIR` to a writable dir if /tmp is restricted (configure.sh uses `mktemp -d`).
- pyALF (needed for `alf_test_branch`, `compare_dirs.py`, test_vs_ed) is not on PATH by default; install it in its own environment.
- pyALF cannot parse gfortran's plain-text complex output → use HDF5 builds for anything pyALF
  analyses.

### configure.sh (must be *sourced*; `make` checks `ALF_FC`)
`source configure.sh MACHINE [MODE] [STAB] [HDF5] [Devel] [NO-INTERACTIVE] [NO-FALLBACK]`
(args case-insensitive, any order; later duplicates overwrite with a warning).
- MACHINE: `GNU`, `Intel` (ifort/mpiifort, MKL), `IntelLLVM`|`IntelX` (ifx; MPI = `mpiifort -fc=ifx`),
  `PGI` (pgfortran / mpifort), `SuperMUC-NG`|`NG`, `JUWELS`, `FRITZ-ifort`, `FRITZ-ifx`, `HELMA`
  (cluster entries do `module load`). Unknown/none → fallback serial gfortran, no MPI, unless `NO-FALLBACK`.
- MODE → macros: `noMPI`|`serial` → none; `MPI` (default) → `-DMPI`; `Tempering` → `-DMPI -DTEMPERING`;
  `PARALLEL_PARAMS`|`PP` → `-DMPI -DTEMPERING -DPARALLEL_PARAMS`. GNU compiler: gfortran / mpifort.
- STAB → macros: none (default stabilisation), `STAB1`, `STAB2` (old UDV), `STAB3` (newest,
  separates scales), `LOG` → `-DSTABLOG` (log storage of scales). `-DQRREF` exists only as a
  commented line.
- `HDF5`: adds `-DHDF5 -DHDF5_ZLIB`; self-builds HDF5 1.14.6 into `HDF5/<compiler_version>/`
  via `HDF5/install_hdf5.sh` (asks unless `NO-INTERACTIVE`); override location with `ALF_HDF5_DIR`.
  Never uses a system HDF5.
- Debug macro `ALF_TRACE_UPDATES` (not set by configure.sh): after `source configure.sh ...` run
  `export ALF_FLAGS_PROG="$ALF_FLAGS_PROG -DALF_TRACE_UPDATES"`, `make clean`, `make`. Every Metropolis decision of
  `Upgrade2` (mode Final) and `Upgrade_lowrank` appends `nt n_op Re(R) Im(R) Weight ranf toggle` to
  `trace_updates.txt` (`Prog/Trace_updates_mod.F90`); diff the files of two runs that differ only in
  `Sequential_scheme`. Without the macro the build is unchanged.
- Debug macro `ALF_TIMERS` (CG-ALF T8 base; not set by configure.sh): `export ALF_FLAGS_PROG="$ALF_FLAGS_PROG -DALF_TIMERS"`,
  `make clean`, `make`. Exclusive stack timers (`Prog/Timers_mod.F90`) split every sweep into update, flush, vwrap,
  hop, stab, meas, global, langevin_hmc and other, with counters (proposals, accepts, flushes by reason); one block
  per bin in `timers.txt` next to the info file (max/mean over the group's ranks). The sum of the parts equals the
  sweep time exactly (checked every bin). Without the macro the build is unchanged.
- `Devel`: GNU adds `-Wconversion -fcheck=all -g -fbacktrace -fmax-errors=10 -pedantic -Werror -Wno-error=cpp`.
- `ALF_FLAGS_EXT`: extra flags appended everywhere (export before sourcing).
- GNU opt flags: `-cpp -O3 -ffree-line-length-none -ffast-math -fopenmp`; `-std=f2008` for Prog/Analysis.
- `configure.sh ... GPU=none|host` (GPU line `Code-modularity-GPU`, P2.4 step 1a): `GPU=host` builds the OpenMP device path (`gpu_*_mod`) running on the host; `GPU=none` (default) gives the CPU code unchanged. Vendor backends (`GPU=nvidia|intel|amd`) come with P2.4 step 2, see docs/ai/phase2/P2_4_GPU_PLAN.md.
- `GPU=amd` (P2.4 step 2), **experimental: KNOWN_ISSUES R2 open** (OpenMP target kernels unreliable on gfortran-14 amdgcn;
  the HIP replacements are mandatory). On the amd machine (MI210): `source configure.sh AMDGPU noMPI GPU=amd NO-INTERACTIVE && make`
  (MACHINE `AMDGPU` = GNU plus the toolchain under `ALF_AMD_HOME`, default `~/amd`: `env.sh` with the gcc 14 offload compiler,
  the patched libgomp plugin as `ALF_GOMP`, ROCm `/opt/rocm/core-7.14` with dev headers; it lists missing pieces and stops).
  Elsewhere MACHINE GNU with gfortran >= 14 built with amdgcn offload, ROCm (hipcc, hipBLAS, hipSOLVER,
  rocSOLVER); env `ALF_AMD_ARCH` (default gfx90a), `ROCM_PATH`, `ROCM_INC`. Flags `-foffload=amdgcn-amdhsa -DALF_GPU
  -DALF_LA_HIP`, no `-static`, no `-ffast-math`. The HIP code is `Libraries/Modules/linalg_hip_mod.F90` (device path of
  `linalg_mod`, presence dispatch), `alf_hipk_mod.F90` + `alf_hip_kernels.hip.cpp` (HIP kernels replacing hot OpenMP kernels
  of `gpu_*_mod` inside `#if defined(ALF_LA_HIP)`), `alf_hip_qr.hip.cpp` (own QR, m <= 2304), `alf_la_hip_shim.cpp`.
  Run-time switches: `ALF_HIP_QR=0` (hipSOLVER QR only), `ALF_HIP_KERNELS=0` (the OpenMP kernels instead of the HIP ones;
  debug only, needs `ALF_GPU_DEBUG=1`, warned in the info file; R2), `ALF_HIP_ASYNC=0` (synchronise after every HIP call). Recommended: `GCN_NUM_THREADS=4`,
  `OMP_TARGET_OFFLOAD=MANDATORY`, and gpu-port's patched libgomp GCN plugin (`docs/gpu/platforms/amd_mi210/*.patch` on
  `fork/gpu-port`: pooled launches, presync hook for the deferred synchronisation; without it every HIP call synchronises).
  Device tests: `testsuite/GPU.tests/55-linalg-hip` (wrappers vs host), `56-hipk-stab` (HIP vs OpenMP stab kernels),
  registered only for `-DALF_LA_HIP`. CI does not build GPU=amd (no ROCm/amdgcn there): check on the MI210 (`ssh amd`).
  Known: KNOWN_ISSUES R1 (rocSOLVER zgeqrf NaN, guarded), R2 (open; the CGR-phase LU of the debug modes
  `ALF_GPU_CGR_LU=1,2` runs on host LAPACK on GPU=amd).
- Exports: `ALF_DIR ALF_FC ALF_LIB ALF_FLAGS_QRREF ALF_FLAGS_MODULES ALF_FLAGS_ANA ALF_FLAGS_PROG`.
- Preprocessor macros actually used in Prog/ + Libraries/Modules: MPI (178 #if lines), HDF5, TEMPERING,
  PARALLEL_PARAMS, HDF5_ZLIB, STAB1/2/3, STABLOG (cgr1, wrapur, udv_state…), OLDNAG, `__PGI`/`__NVCOMPILER`.

### Makefiles
- Top `Makefile`: `all` = `lib ana program`; `lib` → `Libraries/` (Modules → `modules_90.a`,
  libqrref → `libqrref.a`); `ana` → `Analysis/` (`ana.out`, `ana_hdf5.out`); `program` → `Prog/ALF.out`.
  Clean: `clean`=`cleanall`, `cleanprog`, `cleanlib`, `cleanana`, `tidy`.
- `Prog/Makefile`: `SRCS` is the single list of core sources (add new modules there, in any order).
  `Prog/gen_deps.py` scans `module`/`use` statements and writes `.deps.mk` (OBJS, MODS, `x.o: y.o`
  rules) → `make -j` safe. `Libraries/Modules/Makefile` uses the same script (`../../Prog/gen_deps.py`).
  Every Prog object also depends on the `.mod` files of `Libraries/Modules` and `ALF.out` on the two archives,
  so a library change rebuilds or relinks the program (no `make clean` needed; test `50-make-lib-deps`).
- Hamiltonians: registered in `Prog/Hamiltonians.list` (+ `Prog/Hamiltonians.list.d/*`):
  Kondo, Hubbard, Hubbard_Plain_Vanilla, tV, LRC, Z2_Matter, Spin_Peierls. `parse_ham.py` generates
  `Hamiltonians_interface.h`, `Hamiltonians_case.h`, `Hamiltonian_*_read_write_parameters.F90`.
- `git.h` regenerated by `Prog/Git_config.sh` (commit hash, `-dirty`, branch, `git status` into info file).

### Commands
```sh
cd <ALF>                                   # a clean checkout; build products are untracked
. ./configure.sh GNU noMPI NO-INTERACTIVE  # or MPI | Tempering | PP ; add STAB1|STAB2|STAB3|LOG ; HDF5 ; Devel
make -j8 all                               # 'all' (not 'program'): Analysis.tests need Analysis/*.o
# switching variant: make clean (objects carry the macros), re-source configure.sh, make all

# quick run (as CI)
cp -r Scripts_and_Parameters_files/Start run && cd run && "$ALF_DIR/Prog/ALF.out"          # serial
mpiexec -n 4 "$ALF_DIR/Prog/ALF.out"                                                       # MPI
# Tempering/PP: copy Start_Tempering instead; analysis then in Temp_0
"$ALF_DIR/Analysis/ana.out" *              # or ana_hdf5.out for HDF5 builds

# unit tests (uses $ALF_FC, $ALF_FLAGS_PROG, $ALF_LIB from configure.sh)
mkdir -p testsuite/tests && cd testsuite/tests && cmake .. && cmake --build . -j8 && ctest -j8 --output-on-failure

# MALA regressions: serial, non-HDF5 build only (CMake FATAL_ERROR otherwise)
. ./configure.sh GNU noMPI && make -j2 all
cmake -S testsuite -B testsuite/tests -DALF_MALA_TESTS=ON && cmake --build testsuite/tests -j2
OMP_NUM_THREADS=1 ctest --test-dir testsuite/tests --output-on-failure
# CI also runs it as "debug": ALF_FLAGS_EXT='-O0 -fno-fast-math -g -fcheck=all -ffpe-trap=invalid,zero,overflow'
```
- Variant set used per Phase-6 commit (PR #657): clean builds noMPI, MPI, TEMPERING, PARALLEL_PARAMS,
  and MPI with STAB1, STAB2, STAB3, LOG, each followed by ctest; HDF5 builds for HDF5 steps.
- Valgrind (CI): `. ./configure.sh GNU devel noMPI; make program; cd Scripts_and_Parameters_files/test_valgrind;
  valgrind --leak-check=full --errors-for-leak-kinds=definite,indirect,possible --error-exitcode=1 ../../Prog/ALF.out`.

### testsuite/ contents
- `CMakeLists.txt`: compiler/flags/libs from env (`ALF_FC`, `ALF_FLAGS_PROG`, `ALF_LIB`), fallback
  gfortran defaults. Subdirs: matmod.tests, Prog.tests, Analysis.tests; MALA.tests if `ALF_MALA_TESTS=ON`.
  Tests link the `.o` files of the normal make build → build first.
- `matmod.tests/` (Libraries/Modules): 1 diag real; 2–9 matrix inverses (real/complex, variable sizes);
  10 initd; 11/12 UDV complex/real; 13-linalg: every `la_*` wrapper of `linalg_mod` vs direct BLAS/LAPACK, bitwise.
- `Prog.tests/`:
  - Operator_mod: 9 Op-Phase, 10/11 Op-mmultL/R, 13/14 Op-Wrapup/Wrapdo, 21 Op-make, 22 (file 21-sl-mat-mults),
    27 Op-mmultL-m1, 28 Op-validate, 33 Op-equal, 34 V-symm-pass; expected failures (WILL_FAIL):
    29 hermitian, 30 type, 31 projector, 32 1x1, 35 V-symm.
  - 26 polymorphic OpT types (OpTTypes/DynamicMatrixArray/ContainerElementBase).
  - cgr: 15-cgr, 20-qdrp, 23-cgrp; cgr2_2: 16-get-blocks, 17-solve-extended-system;
    udv_state: 18/19 ul/ur-update-matrices, 24-udv, 25-assign-UDV-state.
  - 36-Ham-smoke-<model>: Alloc_Ham + Ham_set + Validate with default parameters, one ctest per model in
    Hamiltonians.list (auto-discovered), own run dir, 60 s timeout. Skipped in TEMPERING builds.
  - 37-QMC-runtime-var (defaults, set/get, CPU_MAX policy); 38-QMC-runtime-lock-fail (setter after
    `lock_QMC_runtime_settings` must abort, WILL_FAIL).
  - 39-comm-mod (comm_mod rank/size/reduce/bcast); 40-walker (qmc_walker_mod layout/alloc/dealloc,
    skipped in TEMPERING builds). MPI builds run 36–40 via `mpiexec -n 1`.
  - Files 1–4, 6–8 exist but are commented out in CMake; `16-UDV/`, `18-UDV/` are old standalone Makefile tests.
- `Analysis.tests/1-maxent-kernel-overflow` (links `Analysis/maxent_wrapper_mod.o`).
- `MALA.tests/`: `mala-measurement-cache` (Fortran, links all Prog objects except main.o via `Prog/.deps.mk`),
  `mala-measurements`, `mala-restart` (`integration.py` runs ALF.out on a 4-site Hubbard chain). Label MALA.
- `test_vs_ed/`: `test_specs.yaml` (5 specs: hubbard_finT_chain, hubbard_T0_chain, 3× MALA variants);
  each runs 5 Δτ values, fits E(Δτ)=y0+aΔτ², passes if |y0−E_ED| ≤ `max_delta` (0.001) and ≤ `max_sigma` (3)σ.
  `prep_test.py` → dirs + CI matrix, `simulate.py`, `analysis.py`.
- `test_branch_parameters.json`: 17 runs for alf_test_branch: Hubbard_N_Leg_Ladder(_temper), Langevin,
  Langevin_1, HMC_1, Hubbard_PAM(_1,_2), Hubbard_Bilayer, Hubbard_Plain_Vanilla(_1), Kondo, LRC, tV,
  Z2_Matter, Spin_Peierls, Hubbard_Global_moves (added in 0.11 `c2e52844`).
- `test_branch_generate_config.py`: generates a *GitLab* CI pipeline (`generated-config.yml`, runner tag `k8s`,
  images on git.physik.uni-wuerzburg.de, `macGNU` with tag `macos`); `--github` raises NotImplementedError.
- `compare_dirs.py DIR_R DIR_T [--results f]`: pyALF analysis of both dirs, `np.allclose` per observable.

## 2. Regression procedure "bitwise identical to master"

Tool facts (read from pyALF source, `py_alf/cli/alf_test_branch.py`):
- `alf_test_branch` runs `git checkout <branch>` inside `--alfdir`, compiles R and T, runs each parameter
  set, analyses with pyALF and compares with `np.allclose` (rtol 1e-5, atol 1e-8). So it is NOT a
  bitwise test by itself; it writes `test.txt` (`<name>: True/False`, `Total`). Exit 1 on failure.
- Because it checks out branches, use a separate clone (git refuses to check out a branch that is
  checked out in another worktree; it also needs a clean tree).

Procedure used for every Code-modularity commit (PR #656/#657 descriptions):
```sh
git clone <repo> /path/alf_tb && cd /path/alf_tb   # both master and the test branch available
export ALF_DIR=$PWD
alf_test_branch --alfdir $PWD --sim_pars $PWD/testsuite/test_branch_parameters.json \
  --machine GNU --mpi --n_mpi 4 --branch_T <branch> --branch_R master --devel
# (CI adds --mpiexec_args=--use-hwthread-cpus for Open MPI only)
cat test.txt                                         # all 17 must be True
```
Then, for "bitwise":
1. compare the 18 `info` files of R vs T runs (dirs `<name>` vs `<name>_test`, Tempering: `Temp_*`);
   ignore lines that always differ: `START TIME`, `FIN TIME`, `CPU Time`, git hash/branch/status block
   (`qmc_setup_mod` writes them; `git.h`/`git_status.h`). Everything else must be identical.
2. a `Scripts_and_Parameters_files/Start` run with master and branch binaries (same configure line):
   output files identical (plain-text builds: `diff`/`cmp`; exact file list compared: unverified).
3. for steps touching stabilisation/sweeps/updates: also Start runs with STAB3 and LOG builds, and a
   CPU_MAX truncation run on 3 ranks (clean stop).
- Never compare across compilers/builds: RNG streams differ between compilers.
- Commit messages record the outcome ("Bitwise identical to master." / "Results unchanged: bitwise
  identical to master."); 54 of 149 commit bodies say "bitwise".
- CI equivalent: workflow `test-branch.yaml` (manual dispatch), branch_T = the dispatched ref.

## 3. CI (.github/workflows, all on GitHub-hosted runners)

No `self-hosted` runner in any workflow. Runners: `ubuntu-latest` (mostly in containers
`ghcr.io/alf-qmc/alf-container/{alf,pyalf}-requirements/<image>`), `ubuntu-slim`, `macos-latest`.
The GitLab mirror (git.physik.uni-wuerzburg.de) has its own runners (tags `k8s`, `macos`) used by
the generated GitLab pipeline.

| Workflow / job | Trigger | What it does / enforces |
|---|---|---|
| tests.yml `lint_stopcommands` | every push, dispatch | fails on `stop`/`error stop` at line start in `*.F90` (excl. Analysis, HDF5, testsuite, Not_supported, runtime_error_mod.F90) |
| tests.yml `report_use_only` step 1 | push | **fails** on `USE` without `ONLY:` in `Prog/*.F90`, excl. `Prog/Hamiltonians/*`, `Prog/matTypes_mod.F90`; allowed bare: `mpi`, `hdf5`, `h5lt`, `alf_hdf5` (since 1.34 `af6be0f1`) |
| tests.yml `report_use_only` step 2 | push | **fails** on direct `call z(gemm|hemm|trsm|trmm|gemv|geru|scal|copy|lacpy|getrf|getrs|getri|geqrf|geqp3|ungqr|unmqr)(` in `Prog/*.F90` + `Libraries/Modules/Mat_subroutines_mod.F90`, excl. Hamiltonians, `Operator_mod.F90`, `matTypes_mod.F90` (since 6G.C10 `618674ac`) |
| tests.yml `report_use_only` step 3 | push | report only: count of bare USE in Prog/ (step 0.14) |
| tests.yml `report_save_state` | push | awk scan of `Prog/*.F90` + `Libraries/Modules/*.F90` (excl. Not_supported, **incl. Hamiltonians**): reports `local_save`, `module_save`; **fails** on `local_implicit` (local initialised in declaration, non-PARAMETER) (since 2.13 `bc88e1a3`) |
| tests.yml `test_linux` | push | matrix {trixie GNU, bookworm-pgi-24-09 PGI, bookworm-intel-2024.2 INTEL, bookworm-intel INTELLLVM} × {noMPI, MPI, TEMPERING, PARALLEL_PARAMS}, always `Devel HDF5`; action `alf_test` |
| tests.yml `test_macos` | push | macos-latest, GNU × 4 modes, `Devel HDF5`; brew mpich + gcc |
| tests.yml `mala_regressions` | push | trixie, GNU noMPI, `optimized` and `debug` builds, `-DALF_MALA_TESTS=ON`, ctest |
| tests.yml `valgrind_gnu_trixie` | push | GNU devel noMPI, valgrind on `test_valgrind` run, fails on leaks |
| tests-broad.yml | manual dispatch | 9 images (bullseye/bookworm/trixie/bullseye-openblas GNU, 2 PGI, 2 INTEL, INTELLLVM) + macOS × 4 modes × {Devel,""} × {HDF5,""} |
| test-branch.yaml | manual dispatch | alf_test_branch vs master, 8 container images + macOS, `--mpi --n_mpi 4 --devel`; uploads `test.txt` |
| test-vs-ed.yaml | manual dispatch | prepare matrix → simulate (Linux containers + macOS) → analyze; uploads `results.png`, `results.json` |
| doc.yml | every push, dispatch | `latexmk -pdf doc.tex` in container `aergus/latex`; doxygen `doxygen/Doxyfile`; deploys Pages only on master |
| mirror.yml | push to main, master, `ALF-**`, tags | force-push to GitLab mirror (only in repo ALF-QMC/ALF) via mirror-template.yml |

Composite action `.github/actions/alf_test/action.yml` (args = configure.sh args):
macOS toolchain (brew mpich, gfortran symlink, HDF5 cache) → `. ./configure.sh $ARGS; make all` →
Start (or Start_Tempering) run, serial or `mpiexec -n 4` → analysis → unit tests (`cmake ..; ctest -VV`).
Facts to know:
- Unit tests are skipped when `ALF_FC` is `mpiifort -fc=ifx` (IntelLLVM MPI modes).
- The analysis step tests `$CONFIG_ARGS` for HDF5, but nothing sets `CONFIG_ARGS`, so it always runs
  `ana.out *` (also in HDF5 builds).
- Concurrency: tests.yml/tests-broad/doc cancel in-progress runs of the same ref.
- Not enforced by CI: "MPI only via comm_mod" (`use mpi` still in 11 Prog files incl. main, control_mod,
  qmc_setup_mod), "no implicit SAVE at module level", PROTECTED state (enforced by the compiler instead).

## 4. Coding conventions (Code-modularity rules) with evidence

1. Every `USE` in `Prog/*.F90` outside `Prog/Hamiltonians/` has an `ONLY:` list; exceptions `mpi`, `hdf5`,
   `h5lt`, `alf_hdf5`. — steps 1.1–1.33, CI since 1.34 `af6be0f1` (tests.yml `report_use_only`).
   Hidden re-exports were removed; modules must import what they use themselves (1.11 `28f45e1e`, 1.21 `ccb3265f`).
2. No implicit SAVE: never initialise a local in its declaration; use `parameter` for constants or an
   executable assignment. Applies to Hamiltonians and Libraries/Modules too. — 2.1–2.12, CI fails since
   2.13 `bc88e1a3`.
3. No `STOP`/`ERROR STOP`; terminate via `Terminate_on_error(ERROR_*, __FILE__, __LINE__)` (or `alf_exit(code)` for a bare exit status) from
   `Libraries/Modules/runtime_error_mod.F90`. — tests.yml `lint_stopcommands`; 5.13 `46c84bab`.
4. Runtime (QMC control) parameters live in `Prog/QMC_runtime_var_mod.F90` with get_/set_ accessors and
   `lock_QMC_runtime_settings` (setter after lock aborts). Policy logic belongs there, not in main:
   `apply_CPU_MAX_policy`, `apply_Langevin_overrides`, `set_default_global_tau_sampling`. — 3.1–3.5
   (`d2337a3f`, `80d25e4c`, `f366ceae`, `408618d4`, `b403b024`); tests 37, 38.
5. Model state in `Hamiltonian_main` is `PUBLIC, PROTECTED` (Ndim, N_FL, N_FL_eff, N_SUN, Ltrot, Thtrot,
   Projector, Symm, Group_Comm, Op_T, Op_V, WF_L, WF_R, Calc_Fl, Calc_Fl_map, reconstruction_needed);
   only Hamiltonian submodules write it. Outside: `Setup_flavor_map`, `Set_Group_Comm`, `Clear_Ham_state`,
   `Set_model_for_testing`. Writable: `nsigma`, `leap_frog_bulk`, `ham`. Allocatable arrays use macro
   `PROTECTED_ARRAY` (empty under `__PGI`/`__NVCOMPILER`, compiler bug). — 4.1–4.18, CHANGELOG 2026-11-21.
6. MPI communication goes through `Libraries/Modules/comm_mod.F90` (`comm_rank/size`, `is_root`,
   `comm_world*`, `comm_reduce_sum/max/mean(_array)`, `comm_bcast`; serial no-ops without MPI).
   Root guards: `if (is_root(Group_Comm))` instead of `#ifdef MPI / Irank_g == 0`. — 5.1–5.13, 6.8 `d7411a7a`.
   (Convention, not CI-checked.)
7. Walker state in `qmc_walker_t` (`Prog/qmc_walker_mod.F90`); work arrays allocated once per run by
   `*_work_alloc` called from `Walker_alloc`, not per call. — 6.1a `11e647cb`, 6G.D1–D3 (`6c65e82a`…).
8. QMC core calls BLAS/LAPACK only through `la_*` wrappers in `Libraries/Modules/linalg_mod.F90`
   (same argument list as the reference routine), with whole arrays, not element addresses. — 6G.A*, 6G.B1
   `d2a3bb64`, 6G.C1–C9, CI since 6G.C10 `618674ac`.
9. main.F90 only orchestrates; steps live in `qmc_setup_mod`, `qmc_walker_mod`, `qmc_sweep_mod`,
   `qmc_driver_mod`. — 6.5–6.7, documented in sampling.tex (6.9 `0d5844c9`).
10. Every refactor step is behaviour-neutral and bitwise identical to master (Section 2); when not
    obviously so, the message argues why (e.g. "None of the moved calls draws random numbers", 6.7h `37ea48ea`).
11. New/moved functionality gets a unit test in `testsuite/` (37–40, 13-linalg); tests must compile in all
    STAB/LOG variants (`5f1c3855`); core routines must not need MPI_Init in unit tests (`af22b976`).
12. New core source file → add to `SRCS` in `Prog/Makefile` (or `Libraries/Modules/Makefile`); deps are automatic.
13. Interfaces of existing routines, Hamiltonians, input and output stay unchanged; CHANGELOG only for
    user-visible/backward-incompatible changes (PR #657: "the CHANGELOG, because nothing changes for users").
14. Exempt directories: `Prog/Hamiltonians/` (USE/BLAS rules), `Prog/Hamiltonians/Not_supported/`
    (all scans), `Prog/matTypes_mod.F90` (not compiled), `Prog/Operator_mod.F90` (BLAS rule).

## 5. Documentation, CHANGELOG, CONTRIBUTING

- `Documentation/doc.tex` (SciPost class) `\input`s: intro; model, updating, langevin, hmc,
  stabilization (Sec. finite T); pqmc; sampling (Sec. Monte Carlo sampling, Alg. 1 `alg:1`
  "Basic structure of the QMC implementation in Prog/main.F90 and the qmc_*_mod modules");
  implementation (+time_dependence), files; running; Hubbard_Plain; predefined (+predefined_lattices,
  generic_hopping, interaction_vertices, predefined_observables, trial_wave_function); Model_classes
  (+Kondo_SUN, LRC, Z2_Matter, Spin_Peierls); maxent; conclusion; appendices wick, performance, license.
  MALA is in updating/sampling/files.tex; protected model state/flavour map in implementation.tex.
  Build: `cd Documentation && latexmk -pdf doc.tex` (CI does this on every push).
- Doxygen: `doxygen doxygen/Doxyfile` → `doxygen/html`.
- `CHANGELOG.md` = "Log of backward compatibility changes and critical bugs", newest first. Entry:
  `### YYYY-MM-DD <title>`, then `Author: ... <br>` and links (GitHub PR/issue, GitLab MR) or `Branch ...`,
  then what breaks and a numbered list of what to change (with Fortran snippets), optional
  `#### Optional changes`.
- `CONTRIBUTING.md` (27 lines): community process only (core-dev voting on GitLab, release cycle:
  minor releases April/October, major ~3 years). No coding rules; those live in CI + commit messages.

## 6. Commit-message style (Code-modularity branch)

- Subject: `<step>: <what changed>`; step = `<phase>.<n>[letter]` from the refactor plan
  (`0.11`, `1.34`, `4.5`, `6.1a`, `6G.D1`). Non-plan commits: `<area>: <subject>` (`testsuite: ...`,
  `UDV_Wrap: ...`, `main: ...`, `fix: ...`). Subject is a statement of the new state, no period,
  e.g. "6G.D1: CGR: work arrays allocated once per run", "4.5: Hamiltonian_main: Ndim PROTECTED".
- Body: short plain-English paragraphs wrapped at ~72: what moved where (routine and module names),
  why (often GPU preparation), why results cannot change, what is NOT changed, which test covers it,
  numbers where relevant ("before 6G.C1 it would have found 103 calls").
- Verification line: "Bitwise identical to master." or "Results unchanged: bitwise identical to master.";
  test-only commits say "Only test sources change; ALF itself is unchanged."
- Trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` (147/149 commits).
- One logical step per commit; fixes of earlier steps name them ("fixes 1.24", "(fixes 1.31)").
- Author on this branch: "Fakher Assaad Home" (a few web-UI commits by `fassaad`).

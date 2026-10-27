# 201 – `configure.sh` GPU switch (`GPU=none|nvidia|intel`)

- **Phase:** 2
- **Owner agent:** fortran-gpu-dev
- **Depends on:** 105 (nvfortran via the PGI profile), ADR-001/002 (accepted)
- **Status:** open

## Goal
`source configure.sh MACHINE MODE [STAB] [HDF5] [Devel] [GPU=none|nvidia|intel]` selects the GPU target at
configure time. `GPU=none` (and leaving the argument out) changes **nothing**: every exported `ALF_*` variable is
identical to today's for every profile. `GPU=nvidia` builds unmodified ALF with nvfortran plus OpenMP offload and the
CUDA libraries linked. `GPU=intel` does the same with ifx (JIT `spir64`), for compile checks.

## Context
- ADR-001 §3, PLAN.md §1 "Configure", HARD.md H10. `configure.sh` argument loop at lines ≈ 233–275: each argument
  is upper-cased (`tr`), unknown tokens fall through to `*)` and become `MACHINE`. A new `GPU=*` case must come
  before `*)`.
- Flags flow: `F90OPTFLAGS` / `PROGRAMMCONFIGURATION` → `ALF_FLAGS_MODULES`, `ALF_FLAGS_PROG`, `ALF_FLAGS_ANA`,
  `ALF_FLAGS_QRREF` (lines ≈ 535–555). **The link line** (`Prog/Makefile:26`, `$(ALF_FC) -o ALF.out $(OBJS)
  $(ALF_LIB)`) gets no compile flags, so offload link flags must go into `ALF_LIB`. The CMake testsuite splits
  `ALF_FLAGS_PROG` and `ALF_LIB` on spaces (`testsuite/CMakeLists.txt`): **no flag may contain a space**.
- nvfortran is used through the existing `PGI` profile (105: works unchanged with NVHPC 26.9; `PGIOPTFLAGS` already
  has `-mp`). `noMPI` uses `pgfortran`, `MPI` uses `mpifort`.
- ALF has no `!$omp` directives today, so offload flags alone must not change numerics.

### Specification
| `GPU=` | allowed MACHINE | added compile flags (modules + prog + ana) | added to `ALF_LIB` | defines |
|---|---|---|---|---|
| `none` / absent | any | – | – | – |
| `nvidia` | `PGI` only | replace `-mp` by `-mp=gpu -gpu=${ALF_GPU_CC:-cc80}` | `-mp=gpu -gpu=${ALF_GPU_CC:-cc80} -cudalib=cublas,cusolver` | `-DALF_GPU -DALF_LA_CUDA` |
| `intel` | `IntelLLVM`/`IntelX` only | `-fiopenmp -fopenmp-targets=spir64`; **remove `-static`** (offload needs dynamic libs) | `-fiopenmp -fopenmp-targets=spir64 -qmkl` | `-DALF_GPU -DALF_LA_ONEMKL` |

- Any other combination: print an error naming the allowed MACHINE and `return 1` (no fallback).
- `ALF_GPU_CC` (e.g. `cc80,cc90`) is an optional environment override, documented in the usage text. AOT
  (`spir64_gen`, PVC) is **not** part of this ticket (301).
- Print one line saying which GPU target is active. Export `ALF_GPU=nvidia|intel` for the GPU targets; for `none`
  export nothing new (identity).
- Add `GPU=` to the usage text in `configure.sh`. The user manual (`Documentation/`) is out of scope; document the
  switch in `docs/gpu/BUILD.md` (new).

## Files allowed
- `configure.sh`
- `docs/gpu/BUILD.md` (new): how to configure each target, flags table, known limits
- `docs/gpu/scripts/nividia/` (new helper scripts, e.g. `env_diff.sh`, `build_gpu.sh`)
Anything else: stop and report.

## Acceptance
On nividia, work dir `~/alf-gpu/work/t201/` (rsync the tree there; `TMPDIR=~/alf-gpu/work/t201/tmp`):
1. **`GPU=none` identity:** for MACHINE ∈ {GNU, IntelLLVM, PGI} × MODE ∈ {noMPI, MPI} × STAB ∈ {default, LOG} ×
   {with/without HDF5, with/without Devel} (HDF5 only where it already builds; skip combos that fail on the
   unmodified tree and list them), dump `env | grep '^ALF_' | sort` after sourcing the **unmodified** `configure.sh`
   and after sourcing the new one with no `GPU=` and with `GPU=none`. All three dumps identical. Script committed in
   `docs/gpu/scripts/nividia/`. (HDF5 self-build: reuse the existing `HDF5/<compiler>` build in the rsync'd tree,
   don't rebuild it per combo.)
2. **`GPU=nvidia`:** `PGI noMPI GPU=nvidia` and `PGI noMPI HDF5 GPU=nvidia`: `make -j8 all` succeeds; `ldd ALF.out`
   shows libcublas/libcusolver; the CMake testsuite (`ctest`) passes (39/39 as in 105). A 6×6 Hubbard run
   (`hub6x6_cb` from `testsuite/gpu_reference/configs.py`, fixed seeds, ≤ 2 min) with the `GPU=nvidia` binary and
   with a `GPU=none` PGI binary: compare `info` precision lines and `Ener_scal` bins. Expected bit-identical (no
   directives in the code); if not, report the differences, don't hide them.
3. **Offload smoke:** a 20-line standalone program (in the work dir, not committed) compiled with the exported
   `ALF_FC`/`ALF_FLAGS_PROG`/`ALF_LIB` that runs one `!$omp target teams distribute parallel do` and checks
   `omp_is_initial_device()` is false inside it, plus one `cublasZgemm` call: works for `GPU=nvidia`.
4. **`GPU=intel` – on the `intel` machine** (Arc B580; owner's rule: Intel targets are tested on `intel`, NVIDIA
   targets on `nividia`). Work dir `~/alf-gpu/work/t201/` there, `TMPDIR` inside it, ≤ 8 cores, check `uptime` and
   `ps` first. `IntelLLVM noMPI GPU=intel`: `make -j8 all` succeeds with ifx 2025.3 (`. ~/alf-gpu/env/intel.sh`);
   the smoke program's target region runs **on the B580** with `OMP_TARGET_OFFLOAD=MANDATORY` and
   `ONEAPI_DEVICE_SELECTOR=level_zero:0` (the B580; `level_zero:1` is the iGPU); `ctest` passes. Don't run plain
   `sycl-ls` there (it aborts in the driver and leaves a 110 MB crash report in `/var/crash`).
5. Error cases: `GNU GPU=nvidia`, `PGI GPU=intel`, `GPU=foo` all return 1 with a clear message.
- Report: diff summary, the identity table, build times, `ldd` excerpts, the 6×6 comparison.

## nividia rules
Test only; develop on the laptop. Write only in `~/alf-gpu/work/t201/`. ≤ 4 cores for builds (`make -j4` is fine if
others run), GPU only for the smoke test and the 6×6 run; check `nvidia-smi` first. Each run ≤ 15 min. Clean up
build trees larger than 1 GB at the end (keep one `GPU=nvidia` binary as `~/alf-gpu/bin/ALF-nvhpc-gpu-none-code.out`).
Ticket 202 runs concurrently (different work dir).

## Out of scope
Any Fortran source change; `alf_la` (202); AOT/PVC flags and `SUPERMUC-NG-PVC` profile (301); CMake testsuite
changes; MPI + GPU.

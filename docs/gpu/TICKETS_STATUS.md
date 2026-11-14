# Status of the old tickets (original ticket-driven GPU development)

State 2026-11-14, checked against `gpu-port` at `b7366be2`. Read-only audit; nothing was implemented.

The ticket list of the original development (2026-10-26 to 2026-11-07, checkout `CLAUDE/ALF`) is the set of tables in
[`PLAN.md`](PLAN.md) plus the files in [`tickets/`](tickets/). The ticket process was dropped on 2026-11-06. The
branches of that time are kept as tags `archive/tickets/*` (see [`BRANCHING.md`](BRANCHING.md)).

Closed at the time and not listed below: 101-108, 201-206, 207 (207-1 to 207-5), 208-1 to 208-3, 212, 302a, 302b,
302c-1. Their code is in the trunk (`archive/tickets/gpu-port-208-3` and `archive/tickets/gpu-port-delayed-submatrix`
are ancestors of `gpu-port`).

Verdicts: **done** (in the trunk, possibly in another form), **obsolete** (waived or replaced by a decision),
**open** (not in the trunk), **partly** (the useful part exists, the rest is named).

## Summary

| Verdict | Tickets |
|---|---|
| done | 207-7, 209-2, 301, 304-1, 304-2, 401, 402 |
| obsolete | 207-6, 302c-2, 303-2 |
| partly | 209-1, 209-3, 210, 211, 303-1, 305, 306, 403, 404, 405 |
| open | 213, 406, 407 (with the unmerged work "401-A") |

The only open item with existing unmerged work is **407 / "401-A"**: general vertices on the device, tag
`archive/tickets/gpu-port-delayed-submatrix-rank-nd` (24c85e9e), not an ancestor of the trunk.

## Table

| Ticket | Content | Verdict | Evidence in the trunk | Missing |
|---|---|---|---|---|
| 207-6 | pre-registration of `207.json` on both GPUs | obsolete | `DECISIONS.md` row 41 (2026-11-02): pre-registration, T2 run and acceptance waived for 207 | nothing; the ticket file still says "open" |
| 207-7 | candidate, T1, T2, verification, merge of 207 | done (functional checks) | `tickets/207-closing-state.md`; commit 27b215a7 "Close 207 (decision 41)"; `archive/tickets/gpu-port-207-final` is an ancestor of the trunk | Intel "check 4" (Szz(pi,pi) 6.3 sigma on the B580) was accepted, no later follow-up found |
| 209-1 | `test_vs_ed` on the GPU | partly | `DELAYED_SUBMATRIX.md` "Minimal test round" (ED, 4-site chain, A30 and B580); `scripts/minitests/` | the upstream `testsuite/test_vs_ed` was never run on a GPU build |
| 209-2 | same-compiler trajectory test | done (other form) | `VERIFIED_2026-11-13.md` (GPU vs same-compiler CPU, configurations identical); `scripts/defect/run_cmp.py`; `README.md` "How correctness is checked" | no first-divergence sweep in the current docs (`scripts/nividia/t206_firstdiv.py` exists) |
| 209-3 | one runnable validation suite | partly | pieces: `ALF_GPU_SELFTEST`, `testsuite/alf_la.tests/`, `testsuite/Run.tests/run_checks.py`, `scripts/defect/run_cmp.py` | no single driver; CI builds `GPU=none` only |
| 210 | `BENCHMARKS.md`: s/sweep CPU vs GPU, N = 64..1024 | partly | `TUNING_DEFECT.md` "Time per sweep"; `README.md` speed table; `SPEEDUPS.md` | no full-sweep table in `BENCHMARKS.md`; CPU column is 1 thread; no crossover |
| 211 | CPU regression gate | partly | `.github/workflows/tests.yml` (GNU, PGI, Intel, IntelLLVM; ctest on push); `CONSOLIDATION.md` (`GPU=none` bitwise identical to upstream) | no CI job `GPU=none` vs `master`; no s/sweep timing gate |
| 213 | STABLOG same-seed excess 1.214e-9 > 1e-9 | open (recommended: drop) | `tickets/213-log-same-seed-excess.md` "draft, not started" | cause not found; not re-measured after 212. STABLOG is not used on the GPU (owner, 2026-11-11) |
| 301 | `SUPERMUC-NG-PVC` profile in `configure.sh` | done (other form) | `platforms/intel_pvc_supermuc_ng2.md`; `scripts/supermuc_ng2/build.sh` and `ifx_aot`; `BUILD.md` | no named machine profile; the build goes through the script |
| 302c-2 | accuracy gates of 207 on the B580 | obsolete | `DECISIONS.md` row 41: QREF accuracy work becomes a report, not a gate | round 4 never recorded; ticket file still says "open" |
| 303-1 | AOT build for the B580 | partly | PVC AOT: `BUILD.md`, `scripts/supermuc_ng2/ifx_aot`; B580 small test: `MACHINES.md` | no AOT build of the full GPU path on the B580 (the B580 is no longer a target) |
| 303-2 | the 209 suite on the B580 | obsolete | the 209 suite does not exist as such (see 209-3); the B580 was replaced by the Max 1550 | - |
| 304-1 | Slurm scripts for PVC | done (other form) | `Scripts_and_Parameters_files/JobfileSuperMUC-NG-Phase2-GPU.sh` and `...-GPU-PP.sh`; `scripts/supermuc_ng2/*.sbatch` | binding by `ZE_AFFINITY_MASK`, not `I_MPI_OFFLOAD_CELL` |
| 304-2 | validation recipe with tolerances | done (other form) | `scripts/supermuc_ng2/verify_vs_cpu.sbatch`; `platforms/intel_pvc_supermuc_ng2.md` "Numerical verification"; `DECISIONS.md` row 43 | tolerances are in three places, no single page |
| 305 | evaluate PVC results | partly | `README.md` speed table; `platforms/intel_pvc_supermuc_ng2.md`; `sessions/supermuc_large_prompt/ISSUE_INTEL_VS_A30.md` | no PVC column in `BENCHMARKS.md` |
| 306 | small dense LA on PVC, host vs device | partly | `sessions/supermuc_large_prompt/FINDINGS.md`; own device QR `ALFSH_DEVQR` (`RUNTIME_SWITCHES.md`) | no routine-by-size table; device `getrf` |
| 401 | delayed updates on the CPU | done | commit 80c1c80d; `Prog/Lowrank_update_mod.F90`; `DELAYED_SUBMATRIX.md` "CPU" | no CPU test with op_dim > 1 in `testsuite/` |
| 402 | submatrix updates on the CPU | done | same | same |
| 403 | delayed/submatrix on the GPU | partly | `Prog/gpu_lowrank_mod.F90`; `DELAYED_SUBMATRIX.md`; `FEATURES.md` (NVIDIA, Intel, AMD) | rank-1 diagonal vertices only (see 407) |
| 404 | performance work | partly | several chains per GPU (`MPI_GPU.md`), device observables (`Prog/gpu_obs_mod.F90`), dense hopping (`gpu_wrap_mod.F90`) | no flavour batching (tried, no gain: `NOT_IMPLEMENTED.md`); no CUDA graphs |
| 405 | MPI and several GPUs | partly | `MPI_GPU.md`, `PP_GPU.md`; `gpu_state_mod.F90` (device per rank); production on SuperMUC and Helma | tempering with exchange is refused on the GPU path (`NOT_IMPLEMENTED.md` row 2); device uses the global rank |
| 406 | accept decision on the device | open (idea only) | decision is on the host by design: `gpu_update_mod.F90:428`; `NOT_IMPLEMENTED.md` section 3 | the feature; not a priority since the update is block based |
| 407 ("401-A") | vertices with `N_non_zero > 1` / non-diagonal on the device, `Calc_Fl` on the device | open | refusal in `Prog/gpu_state_mod.F90:573`; `NOT_IMPLEMENTED.md` rows 10 and 15 | everything on the device side. Unmerged work: tag `archive/tickets/gpu-port-delayed-submatrix-rank-nd` (24c85e9e WIP, 34630588 `Calc_Fl` and ED references for tV, Kondo, Spin_Peierls, f72cfe9f ADR-003); K1 passes, K2/K3 untested |

## Notes

- Why several tickets are "partly": the verification standard changed. The statistical acceptance gates were
  abandoned; the standard is the same-compiler, fixed-seed comparison (`README.md`, `VERIFIED_2026-11-13.md`).
- The ticket files of 207-6, 302c-2 and 213 still show "open" or "draft"; they were not edited by this audit.
- Background and the history of the ticket process: `plan/EVAL_ORIGINAL_GPU_BRANCH.md` (outside the repo).

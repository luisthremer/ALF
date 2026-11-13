# Possible speedups (state 2026-11-13)

Where the time goes, from the measurements in the platform documents:

- **One chain per GPU is latency-bound.** A30, L=32: 44 % of kernel time in about 32 000 small
  flush GEMMs, 65 % of host API time in stream synchronisation. Max 1550: device busy 70 to 78 %
  of the time but only about a third of the tile is used; the QR panel and the update phase
  wait on launches. MI210: the QR panel (one work-group) is 36 % of kernel time.
- **Production layouts are throughput-bound.** A30 with 9 ranks under MPS: flush GEMMs 47 %,
  stabilisation 25 %, measurements 17 %; a kernel is running 99.97 % of the time. Max 1550 with
  4 chains per tile: observables 2.7 s and updates 2.7 s of 9.3 s per sweep.

So a change that removes launches helps single chains; only a change that removes device work
helps production.

## 1. Available now, no code change

| Item | Gain | Numerics | State |
|---|---|---|---|
| Intel: `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576` | L=24 single chain 10.6 to 7.1 s per sweep | bitwise identical | Set in the team jobs. Not found in `sessions/supermuc_large_prompt/prod/job_preempt.sh`. Owner decision for new production sets (a running set keeps its settings). |
| Intel: `ALFSH_ASYNCBLK=1` | One chain on one card, L=24: about 0.5 s per sweep, 15 to 20 % per card (`ISSUE_INTEL_VS_A30.md`, V1 to V2). Full node, 32 ranks on 8 tiles (several ranks per tile): 2.57 to 2.65 sweeps/s per node, +3 % (`phase2/team/p3/STATUS.md`) | bitwise identical (measured) | Measured, not yet cleared for production by the owner (same-seed check of the production binary first); waits for the M12 timing test. Both gains are right for their test: the single-chain number hides latency that other ranks on the same tile would fill anyway. |
| All: `ALF_GPU_OBS_SKIP` (no `Green_eq`, `Den_eq`) | 7 % at L=24 under MPS | the skipped outputs stay zero | Needs `ALF_GPU_DEBUG=1` today. Owner decision whether the analysis needs these outputs. |
| Nwrap 10 instead of 5 at L=32 | Nwrap 5 costs about 45 % more | precision rule of `/jakob_params` | Owner decision 44/45 used Nwrap 5; the abort that motivated it was traced to a faulty GPU tile. |
| Intel: `acc0` parallel reduction in `corr` (arrived with the merge from pp-opt) | the two serial `corr` kernels were 0.47 s of 4.95 s at L=24 | summation order of the `acc0` terms changes (1e-13) | In `gpu-port-unified`. Measured 2026-11-13, L=24, one chain, V7 switches: 4.92 to 4.98 s per sweep before, 4.42 to 4.48 after (observables 1.18 to 0.74 s), about 10 %. |

## 2. Implemented on `speedup/nvidia`, not yet measured

The A30 was occupied by the PP production during this work, so both changes compile
(nvfortran 26.9) but have no timing and no numerics run. `~/alf-gpu/work/spd1/ab.sh` on
`nividia` runs the checks and the timing once the GPU is free.

| Commit | Change | Expected gain | Numerics claim to verify |
|---|---|---|---|
| a3f257c6 (in `gpu-port-unified` for oneMKL) | **Phase of the determinant without an LU.** In `dev_cgr` the `la_getrf` of `U_R^H U_L` is used only by `k_det_phase`; nothing is solved with it. The phase follows from the Householder scalars: det(H_i) = -tau_i / conj(tau_i). It is tracked per stack slot in `dev_decompose`. Debug switch `ALF_GPU_CGR_LU=1` restores the LU, `=2` runs both and reports the difference. | **Measured on the Max 1550 (2026-11-13, L=24, one chain): 4.45 to 4.04 s per sweep with the V7 switches (CGR 1.23 to 0.82 s), 10.06 to 9.59 with default switches.** A30 estimate from profiles: 1.2 % at L=32 single chain, about 3 % in the L=24 MPS production; not measured. HIP keeps the LU. | Verified on Intel: phase equal to the LU phase to 5.8e-15 (default and V7 switches), configurations and acceptance bitwise identical, observables differ at 1e-14 relative; with `ALF_GPU_CGR_LU=1` every file is bitwise identical to the old code. CUDA: enable after one run of `spd1/ab.sh`. |
| 50f837f7 | **`la_sync` in `blk_end` removed on CUDA.** `blk_begin` has already waited for the previous flush, and nothing asynchronous touches `cbuf` / `bS_d` in between. | 10 to 20 ms per sweep. | bitwise identical |
| e58302c0 | `scripts/defect/nsys_solver_calls.py`: per-call times of cuSOLVER from an nsys profile. | tool | |

## 3. Measured and decided

| Question | Result |
|---|---|
| Port the AMD device QR to CUDA? | No, for now. cuSOLVER `geqrf` takes 7.8 ms at n=1024 on the A30 (single chain), `ormqr` 2.2 ms, `getrf` 4.1 ms; together 8 % of a sweep. The AMD gain was against hipSOLVER, which is much slower. Best case about 3 % per sweep single chain, at most 8 % per rank in production. Worth another look only if single-chain `geqrf` at n=576 is well above 4 ms. |
| Fused block-end kernel (AMD b3923afe) | No gain (L=24 4.41 s per sweep unchanged). Not ported. |
| N_delay 128, Nwrap >= 20, oneMKL batched routines, `nowait` on ifx, CholeskyQR2, STABLOG | Measured or rejected earlier (TUNING_DEFECT.md, sessions/supermuc_large_prompt/SYNC.md). Not re-opened. |

**Compute-unit partitioning on the MI210 (measured 2026-11-13, L=24, steady state):**
`HSA_CU_MASK` is honoured (2 of 104 CUs: 28.4 s per sweep) and results are bitwise identical
with and without a mask (25 of 25 output files, also with 4 concurrent processes), but it does
not raise throughput.

| Configuration | s per sweep per process |
|---|---|
| 1 process, no mask | 4.46, 4.52 |
| 1 process, quarter of the CUs | about 5.7 |
| 4 processes, no mask | 9.16, 9.27, 9.21 |
| 4 processes, disjoint quarters | 9.41, 9.77, 9.42 |
| 4 processes, two halves shared by two each | 8.74, 8.63 (within noise) |

Four processes still deliver about 1.95 times one process. Recommendation: no masking.
Scripts and kept outputs: `~/amd/work/cu1/` on `amd`.

## 4. Larger candidates (design notes, no code)

Ranked by expected effect on production throughput.

| # | Idea | Where | Mechanism | Numerics | Effort |
|---|---|---|---|---|---|
| 1 | Two-level submatrix update | `Prog/gpu_lowrank_mod.F90` `blk_end` / `blk_flush` | Keep the pending low-rank update X Y^H over all blocks of a time slice and correct each block's gather by X(B,:) Y(B,:)^H. G is flushed once per slice with rank about 190 instead of once per block with rank 10 to 21. The k x k algebra has to move to the device. Targets the largest item in production (flush GEMMs, 47 % on the A30). | rounding changes | large (weeks) |
| 2 | Static hopping by FFT or band | `Prog/gpu_wrap_mod.F90` (`Hd`), `symm_apply` | The static `Op_T` is a dense N x N exponential of a translation-invariant matrix. An FFT or a truncated band removes about 5800 dense GEMMs per sweep at L=24. | rounding changes; truncation must be shown to be below one ulp | medium to large; needs an FFT per backend |
| 3 | `corr` with coalesced memory access | `Prog/gpu_obs_mod.F90` `corr` | Threads read G along shifted diagonals. A pre-pass that forms the products with a tiled transpose and then bins by distance removes the uncoalesced loads (Intel: the large `corr` kernels are 3.5 times slower than their data volume). | summation order | medium |
| 4 | Equal-time measurements on fewer slices | main loop, `gpu_obs_measure` | Measuring on every second slice halves 12 to 25 % of a sweep. | no bias, larger variance per sweep. Owner decision: comparability with existing data. | small |
| 5 | QR panel look-ahead | `alf_hip_qr.hip.cpp`, `alf_la_sycl_qr2.hpp` | Factor panel j+1 while the trailing update of panel j runs. Helps single chains on AMD and Intel. | same arithmetic | medium; needs the node race test again |
| 6 | Device `getrf` on Intel | `alf_la_sycl_shim.cpp` | oneMKL `getrf` is host-assisted (13 ms at N=1024). Owner noted it as the next candidate after the QR. Item a3f257c6 above removes half of the calls. | rounding changes | medium |
| 7 | Skip the G download before the device tau pass (small) | `Prog/main.F90:1298` `call gpu_g_host(GR)` before `TAU_M`; consumer `tau_begin`, `Prog/gpu_tau_mod.F90:416-424` | The call copies `G_d` to the host `GR` (D->H, 16 N^2 B per flavour, about 5.3 MB at L=24) once per sweep when the host copy is stale. Evidence that the device pass does not need it: `TAU_M` calls `gpu_tau_pass(GR,...)` (`tau_m_mod.F90:82`); `tau_begin` reads `GR` only (a) in the host loop that fills the host shadows of `G00_d`, `GT0_d`, `GTT_d`, `G0T_d` (lines 416-424), which are read only by the host observable fallback (`tau_obs`, `tau_download`) and by `ALF_GPU_TAUDIAG`, and (b) in `upd_to_c(G00_d)` when `g_dev_ok(nf)` is false, a case in which the host copy is current anyway. With `tobs_dev` (default: device observables, no diag/qref, `prec_report`) the pass takes `G00` from `G_d` by `k_copy_c` (lines 426-428) and everything else from the device. It IS needed when `tobs_dev` is false (`ALF_GPU_NOOBS=1`, `ham_h0 /= 0`, Projector, self-tests, `ALF_GPU_TAUDIAG`) and when `gpu_tau_on` is false (host `TAU_M`). So the saving is conditional: guard the call with a query of the same condition as `tobs_dev` and skip the four host loops in `tau_begin` then. Not measured; the expected gain is the D->H copy plus the host loops, a few ms per sweep (a fraction of one percent at L=24, 4.5 s per sweep): only worth it together with other small items. | none (no arithmetic changes) | small |

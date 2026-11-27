# P2.5 — GPU feature table (GPU=host column filled with run evidence)

Draft of 2026-11-26 (read-only), evidence added 2026-11-26 on `Code-modularity-GPU` tip b0effdbd (P2.4 steps 1a-1f merged, GPU=host only).
Spec: ROADMAP.md section P2.5. It drives P2.6 (first: delayed/submatrix for every model).

## How to read it

Cell values:
- **works** — the device path (host backend, OpenMP host fallback) runs it; a run agrees with the CPU path of the same binary.
- **refused** — `gpu_supported` refuses; `info` says `GPU path (gpu_state): INACTIVE, reason: ...` and the whole run uses the CPU path
  (a run with the refusal gives output byte-identical to `ALF_GPU_DISABLE=1`; never a silent wrong answer).
- **stops** — the model or the input reader itself stops the run (CPU path and GPU path alike).
- **untested** — nothing run or read settles it.
- **n/a** — no code for it on that line.

Provenance (kept apart on purpose):
- **gpu-port (doc)** columns: copied from `fork/gpu-port:docs/gpu/{FEATURES,NOT_IMPLEMENTED,VERIFICATION}.md` and `docs/gpu/platforms/*`
  (state 2026-11-13). *Not re-checked.* `works*` = the doc says the code path exists, no run on that vendor documented. **Vendor columns of the new line
  do not exist yet** (`GPU=nvidia|intel|amd` is P2.4 step 2; step-2 sessions add them).
- **New line, GPU=host (run)**: result of the runs in [`p25_runs/`](p25_runs/) (one TSV row per run, `p25_runs/results.tsv`) or of the named ctest
  (`p25_runs/ctest_gpu_host.txt`: 51/51 GPU ctests pass on this build). Evidence key: `[build/case]` = row of `results.tsv`.
- **GPU=amd** (P2.4 step 2, feature `gpu-amd`): **experimental, KNOWN_ISSUES R2 open** (OpenMP target kernels unreliable on
  gfortran-14 amdgcn; HIP replacements mandatory); its column is added at step 3 with the measured results of the feature.

## How the evidence was made

Build: `Code-modularity-GPU` b0effdbd, `source configure.sh GNU noMPI GPU=host && make -j16` (gfortran, melmac, `OMP_NUM_THREADS=1`); further builds of the
same commit for LOG, STAB1/2/3, QRREF, MPI, TEMPERING, PARALLEL_PARAMS (PP), HDF5 (`p25_runs/build_variants.sh`).
Each case (`p25_runs/p25_run.py`, cases in `cases_*.py`, everything in `run_all.sh`): Hubbard 4x4 (some 6x6, 20 sweeps) from `Scripts_and_Parameters_files/Start/parameters`
with a few edits, 2-4 sweeps, same seeds, run twice by the same binary: GPU path and `ALF_GPU_DISABLE=1`.
- *works*: info `GPU path (gpu_state): ACTIVE`, `confout` byte-identical, every other output file equal to max(1e-11, 100 x Precision Green mean of the CPU run)
  (the rule of `testsuite/GPU.tests/gpu_e2e.py`; the device QR differs from ZGEQP3 in the last bits). Measured max |x-y|/max(1,|y|): 9e-14 .. 3e-13 (SU(2)),
  7e-10 .. 2e-9 (Mz-Hubbard, whose Precision Green mean is 6e-11: wrapping amplifies the rounding), 6e-9 with HDF5/Mz.
- *refused*: info INACTIVE with the quoted reason and every file byte-identical to the CPU run.
Result: **107 runs, 54 works, 49 refused, 4 stops, 0 wrong results, 0 mismatches with the draft's predictions that were run** (the draft had two wrong
predictions, corrected here: LRC is refused with field type 3, not 4; discrete-field Hubbard `Global_moves` stops in the model on every path).
Sensitivity: the comparison cannot see a wrong device update by itself; its sensitivity is the negative controls of ctest `gpu-state-e2e` (`ALF_GPU_FAULT=upd|obs`, tau faults, `ALF_GPU_FAULT_STAB=perm|det`), which pass.

## 1. Feature table

Backends of gpu-port: NVIDIA = nvfortran OpenMP offload + cuBLAS/cuSOLVER; Intel = ifx + oneMKL + SYCL shim; AMD = gfortran amdgcn + hipBLAS + HIP kernels.

| Feature | NVIDIA (gpu-port doc) | Intel (doc) | AMD (doc) | New line, GPU=host (run) | Evidence (run `[build/case]` or test) |
|---|---|---|---|---|---|
| **Models** | | | | | |
| Hubbard, `N_FL` 1 or 2, field types 1, 2, rank-1 diagonal vertices | works | works | works | **works**: Mz (N_FL 2) and SU(2) (N_FL 1), 4x4 and 6x6 (20 sweeps), N_leg ladder, bilayer square | `[serial/hub-mz-4x4]`, `hub-su2-4x4`, `hub-mz-6x6-20sw`, `hub-su2-6x6-20sw`, `hub-ladder`, `hub-bilayer-square`; ctest `gpu-state-e2e` |
| Kondo | host-assisted only if structure tests pass; not run | same | same | **works** (P2.6.1): standard, delayed N_delay 4 / 16 (rank-2 delayed device path), submatrix 4; max scaled diff 2.5e-12; GPU=amd works (`gpu-amd-p26`, MI210 2.5e-12; 8x8 Dimer_tau G9) | P2.6.1 `p26_runs/` (GPU path ACTIVE vs ALF_GPU_DISABLE=1, same seeds: confout identical, observables within the e2e rule); ctest `gpu-state-e2e` cases `kondo-*`, `tV-*`, `sp-*` |
| tV | same | same | same | **works** (P2.6.1): standard, delayed 4 / 16, submatrix 4, and Ltau = 1 (tau pass with operator segments); max scaled diff 1.9e-13; GPU=amd works (`gpu-amd-p26`, MI210 1.2e-13) | P2.6.1 `p26_runs/` (GPU path ACTIVE vs ALF_GPU_DISABLE=1, same seeds: confout identical, observables within the e2e rule); ctest `gpu-state-e2e` cases `kondo-*`, `tV-*`, `sp-*` |
| Spin_Peierls | same | same | same | **works** (P2.6.1, field type 4 on the device): standard, delayed 4 / 16, submatrix 4; max scaled diff 8.9e-11 (tol 2.5e-10); GPU=amd works at 6x6 (`gpu-amd-p26`, MI210 2.1e-10); 12x12 KNOWN_ISSUES G10 (all vendors' device algorithm, also GPU=host) | P2.6.1 `p26_runs/` (GPU path ACTIVE vs ALF_GPU_DISABLE=1, same seeds: confout identical, observables within the e2e rule); ctest `gpu-state-e2e` cases `kondo-*`, `tV-*`, `sp-*` |
| Z2_Matter | same | same | same | **refused**: `Hamiltonian Z2_Matter: Op_V(1,1) is not diagonal` | `[serial/mod-Z2-matter]` |
| LRC | same | same | same | **refused**: `Hamiltonian LRC: Op_V(1,1) has field type 3 (only 1 and 2 supported)` | `[serial/mod-LRC]` |
| Plain_Vanilla | same | same | same | **works** (P2.6 `p26-calcfl`): device updates (standard, delayed, submatrix) and stabilisation with `ham%weight_reconstruction` in the ratio and the CGR phase; observables and tau pass on the host (`GR(T)_reconstruction`; info `GPU observables: HOST`, `GPU tau pass: OFF, reason: Projector or reconstruction`). Same seeds vs `ALF_GPU_DISABLE=1`: confout byte-identical, scaled diff <= 1.3e-11 (U > 0 half filling, 4x4 and 6x6) and <= 7e-16 (U < 0 off half filling); CUDA H100 the same (<= 2.5e-12, also 8x8). GPU=amd works (`gpu-amd-calcfl`, MI210: 4x4 and 8x8, <= 3.2e-10 at half filling (tol >= 2.5e-9), <= 2.9e-14 for U < 0). With `Projector = T`: `projector code is not supported` | ctest `gpu-state-e2e` cases `pv-*` (GPU=host, melmac); Helma job 1468063 (GPU=nvidia); `[serial/mod-plain-vanilla-proj]` |
| LR_Jz_Test | not run | same | same | n/a (not built on this line) | `Prog/Hamiltonians/` |
| Non-diagonal `Op_V` | refused | refused | refused | **works** (P2.6.1; N <= 8): dense-block wraps (N_Type 1 and 2), operator segments in the slice kernels; device kernels against Operator_mod 0 .. 5e-17 (`ALF_GPU_SELFTEST=2`); GPU=amd works (`gpu-amd-p26`: HIP versions of the wrap and rank-r kernels) | `p26_runs/`; ctest `gpu-state-e2e` |
| HS field type 3 (continuous) | refused | refused | refused | **refused**: `Hamiltonian Hubbard: Op_V(1,1) has field type 3 (only 1 and 2 supported)` (Hubbard `Continuous = T`, sequential updates) | `[serial/mv-continuous-seq]` |
| HS field type 4 | refused | refused | refused | **works** for non-diagonal operators (Spin_Peierls; exp(Phi g E) formed on the device); diagonal type-4 operators refused (no stock model) | `p26_runs/` sp-* |
| Rank > 1 vertices (`N_non_zero` > 1) | refused | refused | refused | **works** (P2.6.1): r x r block D->H, Upgrade2's ratio on the host, accept G <- G - G(:,P) W R on the device; delayed: r new columns of X, Y per accept | `p26_runs/` (Kondo, tV, Spin_Peierls: r = 2) |
| Field type differs between flavours | refused | refused | refused | untested (no stock model produces it) | code only |
| Time-dependent coupling of `Op_V` (`g_t`) | refused | refused | refused | untested (no model with `g_t` on this line) | code only |
| `Calc_Fl` / weight reconstruction | refused | refused | refused | **works** (Plain_Vanilla, row above; device observables stay refused with Calc_Fl, the host measures); GPU=amd works (`gpu-amd-calcfl`) | ctest `gpu-state-e2e` `pv-*` |
| **Hopping, Trotter** | | | | | |
| Checkerboard hopping | works | works | works | **works** (default) | all Hubbard runs; ctest `gpu-wrap-selftest` |
| Dense hopping | works | works | works | **works** (`Checkerboard = F`) | `[serial/hub-checkerboard-F]`; ctest `gpu-state-tables-dense`, `gpu-wrap-selftest-dense` |
| Symmetric Trotter (`Symm`) | works | works | works | **works** with `Symm = T` (default) and `Symm = F` | `[serial/hub-symm-F]` |
| Time-dependent hopping `Op_T` with `g_t` (defect) | works | works | works | untested on this line (defect model not on this line, no test for the `td_*` tables) | none; P2.6 "defect model" |
| Static element after a time-dependent one | refused | refused | refused | untested (no model) | code only |
| N_SUN | untested | untested | untested | **works**: `N_SUN = 4` with Mz = F (N_FL 1) and Mz = T (N_FL 2) | `[serial/hub-nsun4-su2]`, `hub-nsun4-mz` |
| Complex hopping | works | works | works | **works**: twist `Phi_X = 0.25` and flux `N_Phi = 1` (complex `Op_T`) | `[serial/hub-complex-phiX]`, `hub-complex-flux` |
| **Update schemes** | | | | | |
| Sequential local updates, rank 1 | works | works | works | **works** | `[serial/hub-mz-4x4]`; `gpu-state-e2e` |
| Delayed updates | works | works | works | **works** (`N_delay` 4, 16) | `[serial/sch-delayed-N4]`, `sch-delayed-N16-su2` |
| Submatrix updates | works | works | works | **works** (`N_delay` 4) | `[serial/sch-submatrix-N4]` |
| Delayed/submatrix for non-Hubbard models, rank > 1, non-diagonal vertices | refused | refused | refused | **works** (P2.6.1): delayed keeps G = G0 - X Y^T across general vertices (rank-r columns, X/Y follow the dense wraps); submatrix flushes before them as the CPU does (K_ND = 1); schemes agree bin for bin within 6.3e-11 (gate K4 1e-9), confout identical | `p26_runs/`, K4 check |
| `Sequential_scheme = "auto"` | n/a | n/a | n/a | **refused**: `Sequential_scheme = "auto" is not supported on the GPU path` | `[serial/sch-auto]` |
| Timed `N_delay = -1` with delayed/submatrix | untested | untested | untested | **refused**: `timed automatic N_delay (N_delay = -1) is not supported on the GPU path` | `[serial/sch-timed-N-1]` |
| `Lowrank_cross_slice` | refused | refused | refused | n/a: the input reader stops (unknown namelist variable, exit 2) | `[serial/sch-lowrank-cross-slice]` |
| Global moves (`Global_moves`, `Global_tau_moves`, `N_Global_tau` > 0) | refused | refused | refused | **refused**: `global moves are not supported on the GPU path` (continuous Hubbard ladder; discrete Hubbard `Global_moves = T` stops in the model on every path: `Global_move not implemented`) | `[serial/mv-global]`, `mv-global-tau`, `mv-global-discrete` |
| Langevin, HMC | refused | refused | refused | **refused**: `Langevin/HMC is not supported on the GPU path` | `[serial/mv-langevin]`, `mv-hmc` |
| MALA (sequential, global, global-tau) | refused | refused | refused | **refused**: `MALA updates are not supported on the GPU path` (all three) | `[serial/mv-mala-seq]`, `mv-mala-global`, `mv-mala-global-tau` |
| `Propose_S0` | refused | refused | refused | **refused**: `Propose_S0 is not supported on the GPU path` | `[serial/mv-propose-S0]` |
| **Stabilisation** | | | | | |
| Default (pre-pivoted QR) | works | works | works | **works**: info `GPU stabilization (step 1d): ACTIVE`, CGR phase from the tracked reflector determinants (1f) | all `serial` Hubbard runs; ctest `gpu-stab-selftest`, `gpu-stab-selftest-nwrap7` |
| `LOG` | host-assisted | same | same | **works** (host `ZGEQP3` in CGR, `STABLOG`): the whole case table (41 cases: 20 works, 19 refused, 2 stops) gives the same verdicts as the serial build | `[LOG/*]`; info `GPU tau pass ... CGR2_2 QR: host ZGEQP3 (STABLOG)` |
| STAB3 | refused | refused | refused | **works** (P2.6 `p26-stab`): info `... STAB3: D+/D- split and host ZGEQP3 in CGR`; CGR and CGR2_2 split the scales as `cgr1_mod` / `cgr2_2_mod`, QR by host ZGEQP3 as STABLOG; tracked det(U) phase (S1) applies. GPU=amd: **refused** (`STAB3 is not supported on GPU=amd (KNOWN_ISSUES R2)`) | [§ P2.6 stab](#p26-stab-stab1-stab2-stab3-2026-11-27); ctest 144/144 in a STAB3 GPU=host build (gpu-stab-selftest: device vs host STAB3 CGR 3.1e-15, gpu-state-e2e); CI job GPU=host STAB3 |
| STAB1, STAB2 | refused | refused | refused | **refused** (compile-time, whole build): `STAB1/STAB2 (legacy stabilization, reference only) are not supported on the GPU path: use the default, STAB3 or LOG`; output byte-identical to `ALF_GPU_DISABLE=1` | [§ P2.6 stab](#p26-stab-stab1-stab2-stab3-2026-11-27) |
| `QRREF` | untested | untested | untested | **works** (3 cases, `QRREF is defined` in info, GPU path ACTIVE) | `[QRREF/hub-mz-4x4]`, `hub-su2-4x4`, `sch-delayed-N4` |
| Determinant-phase underflow, N >~ 2950 | fixed only on `fix/det-phase-underflow` | same | AMD own branch | **works as a unit test** (tracked det(U) formula at n = 3200 gives modulus 1 and the analytic phase; LU phase rule bitwise); no full-model run at that N here (needs a lattice of ~3000 sites) | ctest `gpu-stab-phase`, `gpu-stab-refdet`, `gpu-stab-refdet-3200`; KNOWN_ISSUES DP-GPU |
| **Observables** | | | | | |
| Equal-time, Hubbard | works | works | works | **works on the device**: `GPU observables (step 1e): DEVICE` (also with `ham_h0 /= 0` pinning, `LOBS_ST/LOBS_EN = 10/40` window, `Ltau = 0`) | `[serial/hub-h0-pinning]`, `hub-lobs-window`, `hub-ltau0`; ctest `gpu-obs-*` |
| Equal-time, other models | host-assisted | same | same | refused with the model (whole run on the CPU path) | `[serial/mod-*]` |
| Custom user observables | host-assisted | same | same | n/a on the device (observables of a model are on the device only for Hubbard; other models never get there) | as above |
| Time-displaced, `Ltau = 1`, Hubbard | works | code only | code only | **works**: `GPU tau pass (step 1e): ACTIVE`, `Green_tau`, `SpinZ_tau`, ... compared in every Ltau = 1 run | `[serial/hub-mz-4x4]`; ctest `gpu-tau-selftest-*` |
| Defect local observables | works | works | works | n/a (defect = P2.6) | — |
| Entanglement (Renyi) | CPU only | CPU only | CPU only | CPU only (not run) | — |
| **Projector, MPI, I/O** | | | | | |
| Projector (`Projector = .T.`) | refused | refused | refused | **works** (feature p26-projector): device decompose of the Ndim x N_part states and CGRP, host `Tau_p` for Ltau = 1 (udvr and G exported once per sweep); Mz and SU(2), Ltau 0 and 1, standard / delayed / submatrix, Nwrap > Thtrot+1, Theta 2 and 5. Refused on GPU=amd (R2) and GPU=intel (until checked): `projector code is not supported on the GPU path of this backend` | e2e `proj-*` cases (`gpu_e2e.py`); ctest `gpu-stab-selftest-proj`, `gpu-stab-fault-proj-*` |
| MPI independent walkers | works | works | untested | **works** with `mpirun -np 2` (MPI build, GPU=host): Mz, SU(2), delayed, submatrix ACTIVE and equal to the CPU path; Kondo refused as serial | `[MPI-np2/hub-mz-4x4]`, `hub-su2-4x4`, `sch-delayed-N4`, `sch-submatrix-N4`, `mod-kondo`, `mv-global` |
| `PARALLEL_PARAMS` (PP) | works | works | works | **works**: two groups of one rank (Temp_0: U = 4, Temp_1: U = 3.5), Mz and SU(2) delayed ACTIVE, every `Temp_k` file equal to the CPU path | `[PP-np2/grp-mz-4x4]`, `grp-su2-4x4-delayed` |
| Tempering with state exchange | refused | refused | refused | **refused**: `tempering is not supported on the GPU path` (output byte-identical to the CPU path) | `[TEMPERING-np2/temper-mz-4x4]`, `temper-su2-4x4-delayed` |
| HDF5 output | works | works | untested | **works**: `data.h5` equal to the CPU path (`h5diff -d tol`), `confout_0.h5` exact; Kondo refused with exact h5 equality | `[HDF5/hub-mz-4x4]`, `hub-su2-4x4`, `sch-delayed-N4`, `mod-kondo` |
| Restart (`confin`/`confout`) | works | untested | untested | **works**: second run from `confout` -> `confin`, 2 + 2 bins, GPU path and CPU path equal (Mz; SU(2) delayed) | `[serial/restart-mz]`, `restart-su2-delayed` |
| Rank-to-GPU binding | wrapper | `ZE_AFFINITY_MASK` | untested | n/a (no device) | — |
| `ALF_GPU_DISABLE=1` | refused-by-request | same | same | **refused**: `GPU path (gpu_state): INACTIVE, reason: disabled by ALF_GPU_DISABLE=1` (the CPU side of every comparison) | every `disabled` run; ctest `gpu-state-disable` |
| No offload device | refused | refused | refused | allowed in GPU=host: `ACTIVE ..., OpenMP host fallback (no offload device)` | every active run |

## 2. What decides P2.6 (reading of the table)

1. Everything that is not a Hubbard-like rank-1 diagonal vertex falls back to the CPU. Measured: **Kondo, tV, Spin_Peierls, Z2_Matter** stop at "Op_V not diagonal",
   **LRC** at field type 3; **Plain_Vanilla** (weight reconstruction) runs since P2.6 `p26-calcfl`. First: delayed/submatrix with rank > 1 and non-diagonal vertices
   (ROADMAP P2.6, ticket 407); LRC needs field type 3 (continuous), Plain_Vanilla the `Calc_Fl` weight reconstruction.
2. Features that fall back by design and are covered by a test now: auto/timed N_delay, global moves, Langevin/HMC, MALA, `Propose_S0`, projector, tempering,
   STAB1/2 (refused and identical to the CPU path; STAB3 works since P2.6 `p26-stab`).
3. Still untested after this round: HS field type 4, vertices of rank > 1 / differing field types / `g_t` (no stock model reaches those checks),
   time-dependent hopping and the defect model (not on this line), entanglement, a full-model run at N ~ 3000.
4. Vendor columns cannot be filled for the new line until P2.4 step 2 (they re-run `p25_runs/run_all.sh` with the vendor binary: the script takes the binary
   and variants as arguments).

## P2.6 stab: STAB1, STAB2, STAB3 (2026-11-27)

Feature `p26-stab` (design `.handoff/P2_6_stab_design.md`). **STAB3 ported, STAB1/STAB2 kept refused.**

Why STAB1/STAB2 stay refused: (1) nobody runs them: no production run (none of 9185 `info` files on the owner's Mac is from a STAB1/2
build), no CI job, and the upstream manual (`Documentation/running.tex`) calls STAB1 "for reference only"; (2) a port is a second device
stabilisation, not a variant: their own wrap order (`wrapur/wrapul`: B-product from the identity, then x U x D), two decompositions
(`UDV_WRAP_Pivot`), their own CGR (UDV of TPUP, `INV`, `ZGETRF`, `ZGETRS`) and CGR2_2, about the size of P2.4 steps 1b-1d; (3) the tracked
phase (S1) does not apply (their CGR needs the LU for G itself; `Det_C` in STAB1's pivot is the raw product of the underflow class), so the
device would need the LU at every CGR; (4) the refusal already gives output byte-identical to `ALF_GPU_DISABLE=1`: slower, never wrong.

STAB3 on the device: decompose, wraps, scale test and det tracking are the default ones (STAB3 changes only CGR and CGR2_2 on the CPU).
CGR: `k_tpup_stab3` (the four cases of `cgr1_mod`), RHS rows / columns scaled by 1/D where D > 1 (`k_rowscale_d`, `k_colscale_d`); CGR2_2:
D1m/D2m split and rows of V1^-1 / U2^dagger scaled (`k_dm_d`, `k_rowscale_d`). The QR of CGR and CGR2_2 is host ZGEQP3, as for STABLOG
(16 N^2 B each way per CGR): with the device pre-pivoted QR (`la_zqrpp`) the split matrix, whose columns are not graded, gave a Precision Green
twice as noisy as the default path (table). What remains: at L = 16 the STAB3 device Precision Green is 3x the CPU (default device path
1.5x), at L = 8 the same as the default path; the device decompose (pre-pivoted, shared with the default path) is the suspected source, not
proven. The CGR itself is exact against the host STAB3 CGR on the same states (self-test 3e-15).

Evidence (melmac, gfortran, GPU=host, `OMP_NUM_THREADS=1`, Hubbard from `Scripts_and_Parameters_files/Start`, 6 sweeps x 2 bins, GPU path vs
`ALF_GPU_DISABLE=1` of the same binary, same seeds; scripts `~/work/ALF_P2/p26stab/` on melmac):

| Check | Result |
|---|---|
| GPU=none default and STAB3, GPU=host default and LOG: branch vs base ad11a3b4, 16 serial scenarios of `test_branch_parameters.json` | all byte-identical (4 x 16) |
| ctest, STAB3 GPU=host build | 144/144 (gpu-stab-selftest: device vs host CGR max rel diff G 3.1e-15, phase 0; tracked vs LU phase 0; gpu-state-e2e incl. tau pass, faults, phase modes) |
| ctest, STAB1 / STAB2 GPU=host builds | 105/144: the 39 GPU self-tests need an active GPU path; same failing list on the base build (pre-existing, no CI job builds these) |
| STAB1, STAB2: Start run, GPU path vs `ALF_GPU_DISABLE=1` | refusal in `info`; all 33 output files byte-identical |
| Precision Green mean GPU/CPU, Mz Hubbard, geometric mean of the ratio over seeds | L = 8 (8 seeds): default 1.76, STAB3 **1.80** (with the device pre-pivoted QR in CGR: 3.74, not used); L = 16 (4 seeds): default 1.45, STAB3 **2.99** |
| Precision Green mean, Mz L = 16, Start seed (GPU / CPU) | default 8.4e-11 / 6.2e-11; STAB3 1.7e-10 / 6.5e-11; STAB3 with the pre-pivoted QR in CGR 8.7e-10 / 6.5e-11 |
| **CUDA, Helma H100** (nvhpc 25.7, STAB3 GPU=nvidia build of cb921d86, preempt job 1468768, ~0.5 GPU-h in all) | ctest 147/147; GPU vs `ALF_GPU_DISABLE=1`: confout identical and observables within the rule at Mz L = 8 (2.7e-10, tol 2.3e-9), with submatrix, SU(2) L = 8, Mz L = 16 (1.7e-9, tol 4.2e-9); Precision Green GPU/CPU at L = 8 over 8 seeds 1.39 |
| STAB3 GPU vs CPU, `confout` | byte-identical in every case (4x4 e2e set, L = 8 Mz / SU(2) / submatrix, L = 16, L = 8 Beta = 20 U = 8) |
| STAB3 GPU vs CPU, observables (rule of `gpu_e2e.py`: max(1e-11, 100 x Precision Green mean of the CPU run)) | within the rule: 4x4 e2e set, SU(2) L = 8 (5e-13), Mz L = 16 (1.9e-10, tol 6.5e-9), Beta = 20 U = 8 (3.1e-4, tol 4.3e-4); **just above it** on the Start seed at Mz L = 8, also with submatrix: SpinXY_eq 7.5e-9 against 6.5e-9 (default build, same seed: 3.9e-9), because the GPU run's Precision Green there is 3.3x the CPU one |

## 3. Run evidence in the gpu-port docs (not repeated on the new line)

| Backend | Where | What | Scope |
|---|---|---|---|
| AMD MI210 | `CONSOLIDATION.md` (gpu-port) | GPU vs same-compiler CPU 1.15e-14 (L=6), 5.0e-14 (L=12); L=6..24 | Hubbard defect, `Symm = T`, `N_FL = 2` |
| Intel Max 1550 | same, `platforms/intel_pvc_supermuc_ng2.md` | 2.2e-15 (L6F) .. 5.5e-11 (L6T); L6T with V7 switches 1.36e-10 against tolerance 1e-10 (noise-limited, open) | Hubbard defect, L=6..24; 8..32 rank runs |
| NVIDIA A30 | `BENCHMARKS.md`, `MPI_GPU.md`, `PP_GPU.md` | production and acceptance; unified-commit comparison **not run** | Hubbard, defect, MPI, PP |
| NVIDIA H100 (Helma) | `platforms/helma_h100.md` | L=8 5.2e-13, restart, benchmarks to L=32; nvhpc 24.11 aborts at L=8 (`la_qrpp` NaN), use 25.7 | Hubbard defect, PP, MPS |
| New line, GPU=host | `p25_runs/` | 107 runs listed above, 51/51 GPU ctests | Hubbard 4x4/6x6 and the stock models (refusals), no vendor |

## 4. Open items for the owner / manager

- Whether the vendor columns for the new line are filled by re-running `p25_runs/run_all.sh` at step 2 (proposed: yes, same cases, vendor binary), plus the acceptance suite.
- Cases that need new test inputs to leave "untested": field type 4, rank > 1 vertices, differing field types, `g_t`. They would need either a model on this line or a unit test that builds `Op_V` directly (like `gpu-state-refusals`).
- Status of the cells must be refreshed at each merge (ROADMAP P2.5/P2.6).
- Findings on the way: none wrong (class A/B: none). Note for P2.6: LRC is blocked by field type 3, not 4.

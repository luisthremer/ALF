# p8 / obs: parallel observables kernels (`ALF_GPU_OBSPAR=1`)

2026-11-13, SuperMUC-NG Phase 2, ifx/oneMKL 2025.3.0, worktree `ALF_wt_lp_team_combined` (HEAD 7e8b6032 via `git archive`, plus the two files below).
C = measured (job id given), E = estimate. Nothing committed.

## What changed (source)
- `Prog/gpu_obs_mod.F90` (+233 lines, 0 removed except one `if` wrapped around an old kernel): run-time switch `ALF_GPU_OBSPAR` (read in `gpu_obs_setup`; default off = the old
  kernels untouched), new `corr_par` (+ `pq_tree`, `pq_brev`), device work arrays (`o_jct`, `pq`, `zdn`, `zsp`), dispatch in `gpu_obs_measure`
  (`if (o_par) call corr_par else <old loop>`), allocation in `alloc_dev`, release in `gpu_obs_finalize`. Debug-only `ALF_GPU_OBSPAR_PLANT=1` (needs `ALF_GPU_DEBUG=1`):
  one Green output x (1+1e-12) in a separate tiny kernel. `o_skip_acc0`: the old Den branch of `corr` has its serial acc0 kernel wrapped in `if (.not. o_skip_acc0)`.
  Only for N <= 1024 (else the switch is ignored, old kernels). tau-resolved `tcorr` is unchanged (Ltau = 0 in production).
- `Prog/gpu_state_mod.F90` (1 line): `ALF_GPU_OBSPAR` added to the allow-list of `gpu_env_check` (production runs refuse every other `ALF_GPU_*` variable). `ALF_GPU_OBSPAR_PLANT` stays refused without `ALF_GPU_DEBUG=1`.
- Both files are inside `#ifdef ALF_GPU`: the CPU build (GPU=none) is untouched. No Intel-only clause was added (plain `!$omp target teams distribute parallel do`, `declare target` functions, `target`).
  nvfortran was not run here (no machine); on NVIDIA the old reductions use another work-group order, so the new path is a correct sum there but not bitwise to the old one.
- Patch: `obspar.patch`. Scripts: `jobs/`, `bench/obsbench.F90` (micro benchmark), `summ.sh`, `nodeh.sh`. Logs: `data/`.

## Why a bitwise parallel version is possible (C, micro benchmark `bench/obsbench.F90`, jobs 560851, 560960 ... 560977)
The old big kernels are `teams distribute` over the output m with `parallel do reduction(+:s)` over I = 1..N (N <= 1024 work-items, one I each). Measured: with thread_limit 1024/640/576/512 the
result is bitwise the same, with 256 or fewer it is not, and it is run-to-run reproducible. That is a **halving tree on the work-item index** t = I-1 (zero padded):
A_1024(t) = leaf_t, A_s(t) = A_2s(t) + A_2s(t+s), s = 512..1, result A_1(0). New code evaluates exactly this tree in two stages: 32 partial trees over the leaves t = tp + 32 kk
(halving over kk), one work-item per (partial, output) with the outputs on consecutive lanes (coalesced), then the tree over tp (one work-item per output).
Bitwise equal to the old kernels at N = 36, 144, 576, 1024 and Norb = 1, 2 (synthetic data, all kinds, bench). Two things the bench showed matter:
1. The old leaf is `s1 = s1 + term` with the compiler-generated reduction identity. A literal `s1 = 0; s1 = s1 + term` is folded differently (28 of 576 Den outputs differ); the leaf is therefore formed as
   `gpu_obs_zq + term`, `gpu_obs_zq` a run-time zero (bench: opaque zero 0 differences, literal zero 28). It is PUBLIC so that the compiler cannot fold it as a never-written private module variable (a precaution; making it public did not change the Den result below).
2. Even so, the Den kernel (only Den) differed from the old kernel in ALF (L6F: imag parts of Den_eq at 2e-16, Den entry imj = 21 off by 1 ulp in the first call; bench and ALF compile the old kernel differently; five
   different leaf/update expression forms and a separate update kernel did not change it). Rather than chase the compiler, **the Den correlation keeps the old kernel** (called through `corr(kd, 3, ...)`, only its serial acc0 kernel is skipped). 

## Which outputs are bitwise (C, gate G1)
| output | new path | why |
|---|---|---|
| Green (all) | new, parallel | additions only: independent of how the compiler fuses; tree order as old |
| Den correlation | old kernel | by construction |
| Den acc0, SpinMz acc0 (Obs_Latt0) | new: terms in parallel, serial sum in the old order I = 1..N | same expressions; equal in every gate run |
| SpinZ / SpinXY / SpinT correlations (N_FL = 2, N_SUN = 1) | new, parallel | same expressions as old; equal in every gate run. **Compiler dependent in principle** (complex multiply-add contraction): the Den case above shows that a different fused form can differ by 1 ulp in rare entries; every SpinZ/SpinXY/SpinT entry of all 12 gate configs (hundreds of entries x tens of calls each) matched exactly |
| SpinSUN (N_FL = 1) | new, parallel | same remark, matched in L6 nfl1 |
| scalars (kin, pot, rho), local observables, tau-resolved | unchanged old code | not touched |
So everything checked is bitwise; the only unproven-in-principle parts are the compiler-dependent fusing of the spin kernels (G1 is the guard; re-run it after any compiler change).

## Gates
Same-binary pairs `ALF_GPU_OBSPAR=0` vs `=1`, same seed, fresh start, `run_cmp8.py` (copy of run_cmp_gg.py, key `lat`), `compare.py` with `CMP_TOL=0` on every observable file (E25.17 text = full double) and `confout_0` and the acceptance.
Text build `build_txt` (JIT), job **561036** (build + all gates; 9.6 min), final source.
- **G1 PASS (all, worst rel dev 0.000e+00)**: L6F (10x2 sweeps), L6F long (60x3), L12F (4), L12F long (12x2), L12F with V7 switches, **L24F (1x2 sweeps)**, L6T (Mz=T), L6 N_FL=1 (SpinSUN), L6 N_SUN=4 (no SpinMz), L6 bilayer (Norb = 2).
  Default build (no switch) vs the old p5 binary (`p5/build_txtnew`): L6F and L12F **bitwise PASS** (G1d).
- **G2 (planted 1e-12 in one Green output) FAIL as required**: L6F 6.8e-13, L12F 1.05e-15, L24F 4.97e-16 (rel. to file scale, tol 0). Also shown: with the earlier ALF builds the Den 1-ulp difference was caught by the same check (G1 FAIL 2.3e-17), i.e. the gate sees 1 ulp.
- **G3 PASS**: node layout 32 ranks, 4/tile, CCS 4, OMP 3, 3 sweeps, job **561062**: Precision Green identical to the last digit with and without the switch: V7 mean 1.960957926459005E-14 max 6.536e-10, V2 mean 1.681e-14 max 3.68e-10 (8 ranks: mean 1.81e-14 max 3.33e-11). All PASS (< 1e-8).
- **G4 timing** (job **561062**, 14.9 min, alternating 0/1, 2 repetitions; binary `build_aottt`, L=24 reference input, V7 = staging 1 MB + ASYNCBLK=1 + DEVQR=2 + DEVUNMQR=1, V2 = without the two DEVQR switches):

| layout | OBSPAR 0 | OBSPAR 1 | change |
|---|---|---|---|
| one chain (1 rank, tile 0, OMP 8), V7, 3 runs each, sweep / obs [s] | 4.92 / 1.18 | **4.38 / 0.617** | -11 % / -0.57 s |
| 8 ranks, 1 per tile, OMP 14, V7, sweeps/s per card (CPU-time diff 1 vs 3 sweeps) | 0.3590, 0.3547 | **0.4005, 0.3955** | +11.6 % |
| same, wall (TEAM_TIMERS mean over ranks) sweep / obs | 5.59 / 1.36 | 5.00 / 0.735 | -10.6 % |
| node 32 ranks 4/tile CCS 4 OMP 3, V7, sweeps/s per card | 0.8506, 0.8529 | **0.9201, 0.9179** | **+8.1 %** |
| same, wall sweep / obs [s] | 9.31 / 2.72 | 8.70 / 1.88 | -6.6 % (other phases +0.1 s each: shared tile) |
| node, V2, sweeps/s per card | 0.6755, 0.6851 | 0.7106, 0.5979 (outlier) | rep 1 +5.2 %; rep 2 CPU-time metric polluted (host QR spin) |
| same, wall sweep / obs [s] (both reps) | 11.97 / 2.59 | 11.30 / 1.75 | -5.6 % (-> +6 %) |
Phase means: `data/g4_phase_means.txt`. Expected (E) was -0.4 s one chain and +5-8 % throughput: **measured -0.57 s (-11 %) and +8 % (V7) / +6 % (V2)**.
Micro benchmark (one tile, isolation, L24 Norb 1, N_FL 2): all kinds incl. acc0 old 1680 us -> new fused 449 us (3.8x; with Den old kernel kept about 550 us E); L32 2974 -> 817 us; Norb 2 3479 -> 1754 us.
For small systems (L6, N = 36) the new path is slightly slower (240 vs 214 us per call, C bench): irrelevant for production sizes.

## Budget
Node-h used (sum of job elapsed from the log start/end lines, `data/nodeh.txt`): **1.10** (builds 0.21, benches 0.04, gates G1/G2 0.43, debugging the Den difference 0.22, G4 0.25), limit 3. 30 jobs, all single-node `-p preempt --qos=preempt -x i20r02c06s02`, at most 2 at a time.

## Jobs
bench 560849-51, 560959-60, 560964, 560967, 560971, 560975-77 (+ g1 variants 561009, 561012, 561020, 561021, debug 561024, 561029 that located the Den difference); builds 560981-85, 560994, 561000, 561011, 561018;
G1/G2: 560984, 560987, 560995 (intermediate, failing Den), **561036** (final, all pass); G3/G4: **561062**.

## Files
Remote `$SCRATCH/large_prompt/team/p8/obs/`: `src_head` (git archive HEAD + the two files), `build_txt`, `build_aottt` (final), `build_txt9` (debug variants, superseded), `work_g1`, `work_g4`, `jobs`, `bench`.
Local `docs/supermuc_large_prompt/phase2/team/p8/obs/`: `NOTES.md`, `obspar.patch`, `bench/obsbench.F90`, `jobs/` (g1.sh gates, g4.sh timing, build8.sh, b*.sh benches, run_cmp8.py, wrap_*.sh, p8lib.sh, g1b.sh/dbg.sh debug experiments), `sync.sh`, `rr.sh`, `waitjob.sh`, `summ.sh`, `nodeh.sh`, `data/` (job logs, G1 outputs, phase means).

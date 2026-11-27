# Leverages: what made ALF faster or more robust, measured

Owner (2026-11-26): "we want to report what did leverage", for a report to the ALF project (ALF-QMC/ALF). One row per change
that gains speed, size or robustness, with the measured effect, the input, the machine and where it is. A row is added when
the change is merged (status "merged") or when a measurement decides it (status "measured, not merged").
Numbers are copied from the merge messages, the plans and the session reports; every row names its source.

## Speed and size

| # | Leverage | Effect (measured) | Input / machine | Where | Status | For ALF upstream |
|---|---|---|---|---|---|---|
| 1 | **Fast `Op_V_is_symmetric`** (start-up check for `Symm = .T.`; CG-ALF T3, KNOWN_ISSUES P1) | check alone, old -> new: Hubbard L = 24/32/48: 142 s / 1341 s / ~9.5 h (extrapolated, n_op Ndim^3) -> < 0.1 ms (all operators diagonal); tV L = 24/32/48: 247 s / 2330 s / ~17 h (extrapolated) -> 0.02 / 0.08 / 0.44 s. Same result (symmetric or not) as the old check on 100 testsuite scenarios (all 7 models + 3 test models; 79 symmetric, 21 not) and on 270 random sets (test 53) | tV and Hubbard square, `Symm = .T.`, melmac, GNU, 1 thread (loaded node) | `cm/feat/opv-symm-fast` (CPU line), a9fdbd64 | in review | **yes**: pure CPU, every user with Symm = .T. at large N |
| 2 | **Delayed and submatrix local updates**, native in `lowrank_update_mod`, `Sequential_scheme = auto` picks the faster | Hubbard 16x16 beta 6: standard 12.6 s, delayed 10.1 s, submatrix 9.8 s (1.3x); Kondo 12x12 beta 6: 4.5 / 2.7 / 3.8 s (1.7x); tV 16x16: 12.3 / 10.3 / 12.3 s | Mac, GNU, 1 core | Code-modularity 2bb0258a, eaf94e27 (P2.2) | merged | **yes**: CPU, all models with eligible vertices |
| 3 | **Determinant phase without LU** (S1): phase from the tracked det(U) of the Householder reflectors, det H_i = -tau_i / conj(tau_i) | H100 L = 32: 8.57 -> 7.49 s per sweep (-12.6 %); L = 24 within noise. Tracked vs LU phase within 8.1e-15 on 62 cases incl. L = 56 (T5) | Hubbard, Helma H100 | GPU line a3e1220f (step 1f); CPU line 1ab8ad7d (CPU: sweeps 1-5 % faster at L = 8/12/16, standard and delayed; Intel Max 1550 L = 24: 21.3 -> 20.2 s/sweep, -5 %) | merged (GPU and CPU) | **yes**: also removes one LU per stabilisation on the CPU |
| 4 | **Larger Nwrap** (fewer stabilisations) | Nwrap 10 -> 20: -20 % per sweep (L = 24: 2.83 -> 2.18 s, L = 32: 7.82 -> 6.26 s); thermalised Precision Green max at Nwrap 20: 4.0e-11 (L = 24), 5.2e-10 (L = 32), ~20x below the 1e-8 limit; Nwrap 30 only 2x margin, Nwrap 40 breaks the limit at L = 32; PG grows ~10x from L = 24 to 32, re-check needed at L >= 40 | Hubbard, beta = L, H100 (jobs 1463455/1463834/1463941) | docs/ai/phase2/T1_stab_results.md (GPU line); rd/cgalf-stab 7ad95cf6 | measured (thermalised): Nwrap 20 safe at L <= 32 | **yes**: a recommendation per model and size, no code change |
| 5 | **CholeskyQR2 after scale separation** instead of Householder QR in the decomposition | decompose -16 % (L = 24) / -26 % (L = 32) at the same Precision Green; only 1-3 % of a sweep | Hubbard, H100 | rd/cgalf-stab 7ad95cf6 (T1) | measured, not merged (small gain) | maybe, low priority |
| 6 | **Intel batched switches** (`ALFSH_OBSCORR=1`, `ALFSH_BLKFUSE=5`, `ALFSH_BLKPF=1`) | Max 1550 L = 24 one chain: 4.55 -> 3.75 s per sweep; 4 chains per tile 0.96x of the H100 | Hubbard, SuperMUC-NG Phase 2 | gpu-port `intel-batched`; default on Intel after its own check (owner Q4b) | measured; part of `gpu-intel-speed` (planned) | GPU only |
| 7 | **Many chains per GPU** (MPS, 8 processes) | H100 L = 24: 1.18 sweeps/s per GPU, 2.7x the A30 at the same layout | Hubbard, Helma H100 | gpu-port docs `helma_h100.md` | measured (layout, no code) | GPU only; next: walker batching (CG-ALF T12/I1) |
| 8 | **Device observables with h0 /= 0** (GPU step 1e; gpu-port measured on the host when h0 /= 0) | Intel Max 1550 L = 24, h0 = 1: 40.7 -> 20.2 s per sweep against gpu-port (2x) | Hubbard, SuperMUC-NG Phase 2, 8 runs per node | GPU line b75f9d91 (step 1e) | merged | GPU only |
| 9 | **GPU=amd backend on the new line** (HIP device path of `linalg_mod`, HIP kernels for the stabilisation, wraps, updates and observables, own HIP QR) | MI210, s per sweep, against gpu-port: L = 24 with the dense hopping product (`ALF_GPU_HOPD=1`) 1.85 -> 1.10 (1.7x); L = 32 18.2 -> 14.2 .. 16.9 (1.1-1.3x). With the per-factor default (G2) L = 24 is 2.9x slower (5.36 s; one libgomp launch per checkerboard group); large N L = 56/60/64 run (acceptance 0.41) | Hubbard Mz, beta 5, Nwrap 10, submatrix 32, `Symm = F`, MI210 | GPU line, feature `gpu-amd` (b418419b); docs/ai/phase2/P2_4_step2_amd_verification.md | merged; experimental while KNOWN_ISSUES R2 is open | GPU only |

## Precision and robustness

| # | Leverage | Effect (measured) | Where | Status | For ALF upstream |
|---|---|---|---|---|---|
| R1 | **Determinant phase for N >~ 2950** (L >~ 54): normalised LU phase instead of the raw product (which underflows) | before: NaN phase, silently wrong results at large N; after: runs at L = 56/60/64 (acceptance 0.41) | Code-modularity aef50743 (det-phase), GPU step 1d | merged | **yes**: posted upstream as an issue |
| R2 | **NaN guards** in G and the phase (class-A bug H1) | a NaN stops the run instead of giving silently wrong numbers | Code-modularity (nan-guards) | merged | **yes** |
| R3 | **Per-factor hopping product on the GPU** (dense `Hd` opt-in, G2) | Precision Green GPU/CPU 1.4-2.5x -> 1.00-1.01 (12 seeds) | GPU line a576c1b6 | merged | GPU only |
| R4 | **Configurable precision abort** (`ALF_UNSTABLE_THRESHOLD`) | production can set the abort threshold per run | CPU line 73fb16b8 | merged | yes |
| R5 | **G wrapped outside the sequential range** (H5, `Global_tau_moves` with partial `Nt_sequential`) | before: Precision Green 1.12 with return code 0; after: 2e-14 | CPU line bce3cfbe (B batch 2) | merged | **yes**: real bug in upstream |

## Rule
Every merge that gains speed, size or robustness adds or updates a row here (the manager checks it at the merge).
Measuring tools that will not stay (e.g. `Timers_mod`, owner: temporary) are not leverages.

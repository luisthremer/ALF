# QR_PREPIVOT variants 1–3 on the CPU line (CPU-line queue item 6, qr-prepivot, 2026-11-27)

Test-only compile flag, default off: the CGR QR (`cgr1_mod::CGR`, both sites of `cgr2_2_mod::CGR2_2`) uses a
pre-pivoted QR instead of ZGEQP3. From gpu-port tickets 104 (variant 1) and 108 (variants 2, 3; report
`docs/gpu/reports/108-stablog-pivot.md` on gpu-port). Owner decisions (2026-11-27): one shared helper, CGR sites only
(`udv_state_mod::decompose` unchanged on the CPU line), precision reported against the report-108 rule but not gated,
the gate is the bitwise-unchanged default build.

## Use
`ALF_FLAGS_EXT="-DQR_PREPIVOT"` (= 1), `"-DQR_PREPIVOT=2"`, `"-DQR_PREPIVOT=3 -DQR_PPNB=n"` (n default 32) before
`configure.sh`. Any other value is a compile error.

| variant | what it does | flops |
|---|---|---|
| 1 static | column 2-norms, stable descending sort, unpivoted QR (CPU form of the device `la_zqrpp` on the GPU line) | 1 geqrf |
| 2 two-pass | variant-1 order, unpivoted QR of a copy, stable re-sort by abs(R_jj), unpivoted QR of the re-sorted original | 2 geqrf |
| 3 blocked re-ranking | per panel of nb columns: re-sort the remaining columns by their exact residual norm, ZGEQRF of the panel, ZUNMQR on the trailing columns, recompute their norms (DZNRM2, no downdating); nb = 1 gives ZGEQP3's pivots | 1 geqrf + O(N^3/nb) |

## Code
- `Prog/QDRP_decompose_mod.F90`: `QDRP_prepivot(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK, variant, nb)`, same
  contract as `QDRP_decompose` (D separated, one ZGEQRF-format reflector set, WORK/LWORK for the later ZUNMQR, IPVT
  with ZGEQP3 semantics); variant 0 calls `QDRP_decompose`. `QR_PP_VARIANT` / `QR_PP_NB` are module parameters set
  from the macros. The variant is a runtime argument, so the default build compiles and tests all variants.
- Call sites: `#if defined(QR_PREPIVOT)` around the three CGR QRs; `#else` is the unchanged `QDRP_decompose` call.
- On the GPU line variant 1 is three inline copies (cgr1, `prepivot_qdrp` in cgr2_2, udv_state); the CPU -> GPU merge
  replaces the cgr1/cgr2_2 copies by the helper (same arithmetic: DZNRM2 norms, the same stable sort, IPVT = 1 then
  `QDRP_decompose`) and keeps the udv_state block.
- Fixes against gpu-port's study code (108 review nits): partial last panel, nb > N_part, rectangular N_part < Ndim
  (helper takes the leading Ndim x N_part block), no copies of array sections into LAPACK calls.

## Checks (melmac, GNU -O3, reference BLAS)
- Default build: 16 scenarios of `test_branch_parameters.json` byte-identical to `Code-modularity-CPU` 4853a472
  (`bitwise_vs_ref.py`); ctest 94/94.
- Test `54-qdrp-prepivot` (default build, also Devel -fcheck=all -Werror): seeded graded / tied-norm (column = conjg of its neighbour) / Hilbert /
  rectangular (40 x 25) / odd-size matrices, variants 1, 2 and 3 with nb = 1, 4, 7, 64: ZUNMQR with the returned
  WORK/LWORK gives Q (D R) = A P to 1e-13, IPVT a permutation with identity tail; variant 0 bitwise equal to
  `QDRP_decompose`; variant 3 nb = 1 makes exactly ZGEQP3's pivots (3 graded matrices).
- `ALF_CGR_LU=2` CGR phase check, Hubbard L = 8 beta 10, standard and delayed, every variant: 402 CGR calls, max
  difference 0.0; Devel build with variant 3 nb = 7 (partial panels) L = 6: 122 calls, 0.0, no bounds error.

## Timing
Decomposition alone, complex N x N graded, 1 thread, best of 5, ratio to ZGEQP3:

| N | ZGEQP3 ms | v1 | v2 | v3 nb 32 | v3 nb 64 |
|---|---|---|---|---|---|
| 64 | 0.22 | 0.91 | 1.77 | 0.96 | 0.94 |
| 128 | 1.66 | 0.95 | 1.83 | 0.97 | 1.09 |
| 144 | 2.42 | 0.95 | 1.89 | 0.90 | 1.04 |
| 256 | 16.4 | 0.84 | 1.66 | 0.76 | 0.86 |
| 288 | 22.3 | 0.83 | 1.67 | 0.76 | 0.86 |
| 512 | 122.9 | 0.85 | 1.71 | 0.79 | 0.87 |

Per sweep (Hubbard Mz, beta 10, Nwrap 10, Ltau 1, delayed N_delay 16; 1 thread, 10 runs at a time, not pinned, one run
each, wall time including start-up; noise about ±5 %), s/sweep standard / delayed:

| L | default | v1 | v2 | v3 nb 32 | v3 nb 64 |
|---|---|---|---|---|---|
| 8 | 0.448 / 0.427 | 0.418 / 0.428 | 0.472 / 0.467 | 0.420 / 0.426 | 0.416 / 0.443 |
| 12 | 3.74 / 3.97 | 3.68 / 3.64 | 4.28 / 4.19 | 3.71 / 3.62 | 3.75 / 3.94 |
| 16 | 28.6 / 29.0 | 28.3 / 28.8 | 30.9 / 30.8 | 27.7 / 28.8 | 28.0 / 28.4 |

The CGR QR is a small part of a sweep: variants 1 and 3 are within noise of the default, variant 2 is 5–15 % slower.
No speed leverage (no LEVERAGES row).

## Precision (same runs, same seeds; ratio variant / default)
Report-108 rule as information, not a gate: Green and tau max <= 3x, mean <= 1.5x of the default.

| L (standard) | default Green mean, max | v1 mean, max | v2 | v3 nb 32 | v3 nb 64 |
|---|---|---|---|---|---|
| 8 | 6.8e-11, 4.4e-7 | 1.02, 1.00 | 1.00, 1.05 | 1.03, 1.01 | 1.02, 1.00 |
| 12 | 3.0e-10, 2.1e-5 | 1.09, 1.34 | 1.08, 1.48 | 1.08, 1.17 | 0.97, 0.67 |
| 16 | 2.4e-10, 1.2e-6 | 0.87, 1.31 | 0.84, 0.53 | 0.88, 1.89 | 0.84, 1.29 |

tau max ratio: 1.00 (L = 8, 16) and 0.47 (L = 12) for every variant; tau mean within 1 %. All variants meet the
report-108 rule here (short runs, one seed; report 108 found cases at L = 6..12 where they do not).
(v3 nb 64 at L = 8: N = 64 fits one panel, so it is variant 1.)

Mz sensitivity (max |variant − default| / max |default| over all bins, same seeds, 6 sweeps): SpinXY_eq at L = 16
1.3e-11 (v2) .. 4.6e-11 (v3 nb 32) standard, up to 6.1e-11 delayed; SpinXY_tau, SpinZ_eq, energy <= 4e-12. At
L = 8/12: 1.2e-11 .. 2.6e-11. (GPU line, device pre-pivoted QR at L = 24: ~1e-8; CPU with variant 1 at L = 24: 4.3e-9.)

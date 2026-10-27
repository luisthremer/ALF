# 104 - pre-pivoted unpivoted QR vs ZGEQP3 in CGR and CGR2_2 (CPU study)

Machine: nividia, GNU 13.3 + OpenBLAS, MPI+HDF5, 1 rank per run, Hubbard square lattice, Mz, Ltau=1,
Symm=F, Checkerboard=T (one dense case), dtau=0.1, random start, 20 bins, NSweep 4..200 (sized for
about 4 min per run), identical seeds for up/pp. Total 133 runs, 16056 CPU s.
Scripts: `docs/gpu/scripts/nividia/qr_study/`.

## Implementation
`-DQR_PREPIVOT` (default off), only inside `#if defined(QR_PREPIVOT)` blocks: in CGR before
`QDRP_decompose(N_size, N_part, TPUP, ...)` and in CGR2_2 before both `QDRP_decompose(LQ2,...)` calls
(the `#else` branch holds the unchanged upstream call). The column 2-norms (DZNRM2) are computed, a
stable insertion sort (descending) gives the permutation perm, the columns are permuted, `IPVT=1`
(ZGEQP3 treats all columns as fixed, so unpivoted QR), and afterwards `IPVT = perm`, i.e. ZGEQP3
semantics (column j of A*P is column IPVT(j) of A). Downstream code (permutation-sign loop, ZLAPMR/ZLAPMT,
`solve_extended_System`) is untouched. The sign of the permutation is obtained from IPVT by the existing
cycle-counting loop, so the determinant phase needs no extra bookkeeping.

## Checks
- Unit tests (devel, noMPI, GNU, `-fcheck=all`): macro off 39/39 pass, macro on 39/39 pass
  (includes `15-cgr`; there is no dedicated CGR2_2 unit test, it is covered by the physics runs).
- Macro-off unchanged: the diff only adds `#if/#else/#endif` lines, original lines are the `#else` part
  (verified: `diff` original vs modified shows no removed line). Run L6 b4 U4 w10 default with the build
  from the modified source (macro off) and from the unmodified source: all of E, Kin, Pot, SpinZ_eq,
  SpinZ_tau, Green_tau bin data bitwise equal; `info` precision numbers identical.
- Macro on, same case: bin-1 relative difference of E, Kin, Pot about 1e-13, Green_tau 5e-11
  (rounding-level, as expected from a different but valid factorisation).

## Table (up = upstream pivoted QR, pp = pre-pivot; "x / y" = up / pp)
Columns: Precision Green/tau mean and max from `info`, Precision Phase max, first bin where E, Kin, Pot or
S_zz(R) differs by more than 1e-8 (relative to the largest entry of the dataset in that bin), E and
S(pi,pi) (maximum of S_zz(K)) with errors. "abort" = the program stopped with "Unstable matrix"
(ALF's own precision check) or NaN; w20 and U=8 runs abort or fail frequently also for upstream.

| case | stab | Green mean up/pp | Green max up/pp | tau mean up/pp | tau max up/pp | Phase max up/pp | first diff bin | E up / pp | S(pi,pi) up / pp |
|---|---|---|---|---|---|---|---|---|---|
| L6_b10_U4_w10 | def | 3.50e-11 / 3.52e-11 | 9.50e-07 / 9.32e-07 | 5.82e-12 / 5.83e-12 | 5.11e-07 / 5.44e-07 | 0.00e+00 / 0.00e+00 | identical | -31.7105(0.0436) / -31.7105(0.0436) | 5.138(0.446) / 5.138(0.446) |
| L6_b10_U4_w10 | log | 3.60e-11 / 5.95e-11 | 1.16e-06 / 1.62e-05 | 8.69e-12 / 1.12e-11 | 8.83e-07 / 1.05e-06 | 0.00e+00 / 0.00e+00 | identical | -31.7105(0.0436) / -31.7105(0.0436) | 5.138(0.446) / 5.138(0.446) |
| L6_b10_U4_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b10_U4_w20 | log | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L6_b10_U8_w10 | def | 1.28e-06 / 1.28e-06 | 2.63e-01 / 2.69e-01 | 2.48e-07 / 2.49e-07 | 5.74e-01 / 4.96e-01 | 0.00e+00 / 0.00e+00 | 1 | -21.3774(0.0358) / -21.3774(0.0358) | 13.716(1.381) / 13.716(1.381) |
| L6_b10_U8_w10 | log | 1.29e-06 / 1.31e-06 | 5.76e-01 / 5.76e-01 | 3.22e-07 / 4.58e-07 | 1.53e-01 / 2.72e-01 | 0.00e+00 / 0.00e+00 | 1 | -21.4087(0.0754) / -21.4087(0.0754) | 12.519(1.564) / 12.519(1.564) |
| L6_b10_U8_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b10_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b20_U4_w10 | def | 3.01e-11 / 3.02e-11 | 1.32e-06 / 1.32e-06 | 5.03e-12 / 5.01e-12 | 5.11e-07 / 5.25e-07 | 0.00e+00 / 0.00e+00 | identical | -31.6834(0.0283) / -31.6834(0.0283) | 6.661(0.662) / 6.661(0.662) |
| L6_b20_U4_w10 | log | 3.36e-11 / 3.82e-11 | 1.72e-06 / 2.22e-05 | 1.11e-11 / 1.35e-11 | 3.52e-06 / 3.53e-06 | 0.00e+00 / 0.00e+00 | identical | -31.6834(0.0283) / -31.6834(0.0283) | 6.661(0.662) / 6.661(0.662) |
| L6_b20_U4_w20 | def | 1.26e-04 / 1.21e-04 | 4.84e+00 / 6.92e+00 | 2.00e-05 / 1.95e-05 | 2.64e+00 / 2.70e+00 | 0.00e+00 / 0.00e+00 | 1 | -31.8242(0.1007) / -31.7317(0.0232) | 4.888(0.482) / 5.871(0.757) |
| L6_b20_U4_w20 | log | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L6_b20_U8_w10 | def | 2.62e-06 / 2.17e-06 | 3.84e+00 / 6.07e+00 | 3.97e-07 / 3.73e-07 | 2.10e-01 / 6.86e-01 | 0.00e+00 / 0.00e+00 | 1 | -21.2982(0.0481) / -21.2483(0.0522) | 3.301(0.323) / 5.015(0.479) |
| L6_b20_U8_w10 | log | 2.23e-06 / 2.07e-06 | 2.29e+00 / 2.33e+00 | 5.16e-07 / 6.89e-07 | 4.46e-01 / 1.68e+00 | 0.00e+00 / 0.00e+00 | 1 | -21.3168(0.0477) / -21.2483(0.0522) | 5.356(0.465) / 5.015(0.479) |
| L6_b20_U8_w10_dense | def | 2.32e-06 / 2.32e-06 | 3.22e-01 / 2.95e-01 | 3.99e-07 / 3.98e-07 | 7.30e-02 / 5.76e-02 | 0.00e+00 / 0.00e+00 | 1 | -21.1322(0.0927) / -21.1322(0.0927) | 4.702(0.503) / 4.702(0.503) |
| L6_b20_U8_w10_dense | log | 2.25e-06 / 2.28e-06 | 2.44e-01 / 3.50e-01 | 5.09e-07 / 6.84e-07 | 6.86e-02 / 1.80e-01 | 0.00e+00 / 0.00e+00 | 1 | -21.1322(0.0927) / -21.1322(0.0927) | 4.702(0.503) / 4.702(0.503) |
| L6_b20_U8_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b20_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b40_U4_w10 | def | 3.60e-11 / 3.59e-11 | 5.84e-07 / 5.69e-07 | 5.40e-12 / 5.39e-12 | 4.34e-07 / 4.28e-07 | 0.00e+00 / 0.00e+00 | identical | -31.7623(0.0234) / -31.7623(0.0234) | 4.957(0.235) / 4.957(0.235) |
| L6_b40_U4_w10 | log | 3.84e-11 / 3.87e-11 | 9.03e-07 / 8.98e-07 | 2.90e-11 / 3.16e-11 | 4.39e-05 / 4.39e-05 | 0.00e+00 / 0.00e+00 | identical | -31.7623(0.0234) / -31.7623(0.0234) | 4.957(0.235) / 4.957(0.235) |
| L6_b40_U4_w20 | def | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L6_b40_U4_w20 | log | up: abort: Unstable matrix; pp: ok | | | | | | | |
| L6_b40_U8_w10 | def | nan / nan | 4.67e-05 / 6.04e-05 | nan / nan | 2.96e-04 / 3.32e-04 | 0.00e+00 / 0.00e+00 | identical | nan(nan) / nan(nan) | nan(nan) / nan(nan) |
| L6_b40_U8_w10 | log | 2.29e-06 / 2.31e-06 | 9.19e-01 / 9.18e-01 | 7.33e-07 / 9.40e-07 | 1.48e+00 / 1.29e+00 | 0.00e+00 / 0.00e+00 | 4 | -21.2400(0.0507) / -21.2400(0.0507) | 3.873(0.286) / 3.873(0.286) |
| L6_b40_U8_w20 | def | up: abort: NaN; pp: abort: NaN | | | | | | | |
| L6_b40_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b4_U4_w10 | def | 2.39e-11 / 2.38e-11 | 6.31e-07 / 6.20e-07 | 4.81e-12 / 4.83e-12 | 4.44e-07 / 4.50e-07 | 0.00e+00 / 0.00e+00 | identical | -29.6570(0.0278) / -29.6570(0.0278) | 3.123(0.124) / 3.123(0.124) |
| L6_b4_U4_w10 | log | 2.46e-11 / 3.31e-11 | 8.84e-07 / 1.67e-06 | 6.17e-12 / 8.19e-12 | 4.32e-07 / 6.48e-07 | 0.00e+00 / 0.00e+00 | identical | -29.6570(0.0278) / -29.6570(0.0278) | 3.123(0.124) / 3.123(0.124) |
| L6_b4_U4_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b4_U4_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b4_U8_w10 | def | 1.65e-06 / 1.65e-06 | 1.08e+00 / 9.83e-01 | 3.19e-07 / 3.15e-07 | 4.56e-01 / 3.11e-01 | 0.00e+00 / 0.00e+00 | 6 | -19.8511(0.1520) / -19.8511(0.1520) | 7.372(1.353) / 7.372(1.353) |
| L6_b4_U8_w10 | log | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L6_b4_U8_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L6_b4_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b10_U4_w10 | def | 3.84e-11 / 3.74e-11 | 1.03e-05 / 1.04e-05 | 5.03e-12 / 5.01e-12 | 1.14e-06 / 5.68e-07 | 0.00e+00 / 0.00e+00 | identical | -56.4971(0.0520) / -56.4971(0.0520) | 9.173(1.438) / 9.173(1.438) |
| L8_b10_U4_w10 | log | 5.92e-11 / 2.77e-10 | 5.55e-06 / 5.93e-05 | 8.95e-12 / 1.15e-11 | 2.80e-06 / 1.98e-06 | 0.00e+00 / 0.00e+00 | identical | -56.4971(0.0520) / -56.4971(0.0520) | 9.173(1.438) / 9.173(1.438) |
| L8_b10_U4_w20 | def | up: abort: Unstable matrix; pp: ok | | | | | | | |
| L8_b10_U4_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b10_U8_w10 | def | 1.97e-06 / 2.08e-06 | 1.27e+00 / 1.69e+00 | 3.00e-07 / 3.04e-07 | 3.80e-01 / 4.18e-01 | 0.00e+00 / 0.00e+00 | 1 | -37.7775(0.1320) / -37.7775(0.1320) | 6.417(0.907) / 6.417(0.907) |
| L8_b10_U8_w10 | log | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L8_b10_U8_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b10_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b20_U4_w10 | def | 5.05e-11 / 5.16e-11 | 1.31e-06 / 1.34e-06 | 6.12e-12 / 6.14e-12 | 7.53e-07 / 5.89e-07 | 0.00e+00 / 0.00e+00 | identical | -56.6188(0.0978) / -56.6188(0.0978) | 4.047(0.464) / 4.047(0.464) |
| L8_b20_U4_w10 | log | 1.19e-10 / 2.24e-10 | 4.41e-06 / 6.18e-05 | 2.18e-11 / 2.54e-11 | 1.30e-05 / 1.30e-05 | 0.00e+00 / 0.00e+00 | identical | -56.6188(0.0978) / -56.6188(0.0978) | 4.047(0.464) / 4.047(0.464) |
| L8_b20_U4_w20 | def | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L8_b20_U4_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b20_U8_w10 | def | 2.38e-06 / 2.41e-06 | 1.65e+00 / 2.20e+00 | 4.24e-07 / 4.22e-07 | 4.63e+00 / 4.63e+00 | 0.00e+00 / 0.00e+00 | 1 | -37.9032(0.1214) / -37.9032(0.1214) | 7.827(0.785) / 7.827(0.785) |
| L8_b20_U8_w10 | log | 2.49e-06 / 2.67e-06 | 9.04e-01 / 7.84e+00 | 6.02e-07 / 8.67e-07 | 5.40e-01 / 1.42e+00 | 0.00e+00 / 0.00e+00 | 1 | -37.8856(0.0940) / -37.9032(0.1214) | 3.908(0.402) / 7.827(0.785) |
| L8_b20_U8_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b20_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b40_U4_w10 | def | 5.18e-11 / 5.26e-11 | 2.89e-06 / 2.90e-06 | 6.41e-12 / 6.41e-12 | 6.39e-07 / 6.38e-07 | 0.00e+00 / 0.00e+00 | identical | -56.5763(0.0691) / -56.5763(0.0691) | 3.559(0.210) / 3.559(0.210) |
| L8_b40_U4_w10 | log | 1.43e-10 / 1.65e-10 | 1.14e-05 / 3.17e-05 | 6.04e-11 / 6.58e-11 | 4.04e-04 / 4.03e-04 | 0.00e+00 / 0.00e+00 | identical | -56.5763(0.0691) / -56.5763(0.0691) | 3.559(0.210) / 3.559(0.210) |
| L8_b40_U4_w20 | def | up: abort: Unstable matrix; pp: ok | | | | | | | |
| L8_b40_U4_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b40_U8_w10 | def | nan / nan | 4.14e-09 / 4.18e-09 | nan / nan | 1.33e-04 / 2.03e-04 | 0.00e+00 / 0.00e+00 | identical | nan(nan) / nan(nan) | nan(nan) / nan(nan) |
| L8_b40_U8_w10 | log | 2.66e-06 / 2.64e-06 | 9.06e-01 / 9.16e-01 | 1.30e-06 / 1.52e-06 | 5.78e+00 / 5.78e+00 | 2.00e+00 / 2.00e+00 | 1 | -38.0134(0.0944) / -38.0134(0.0944) | 3.461(0.199) / 3.461(0.199) |
| L8_b40_U8_w20 | def | nan / nan | 5.10e-04 / 5.11e-04 | nan / nan | 3.58e+07 / 6.75e+07 | 0.00e+00 / 0.00e+00 | identical | nan(nan) / nan(nan) | nan(nan) / nan(nan) |
| L8_b40_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b4_U4_w10 | def | 2.39e-11 / 2.38e-11 | 6.73e-07 / 6.49e-07 | 4.29e-12 / 4.30e-12 | 2.39e-07 / 3.74e-07 | 0.00e+00 / 0.00e+00 | identical | -52.7953(0.0726) / -52.7953(0.0726) | 3.448(0.333) / 3.448(0.333) |
| L8_b4_U4_w10 | log | 2.62e-11 / 6.86e-11 | 4.53e-07 / 5.20e-06 | 5.78e-12 / 7.94e-12 | 4.81e-07 / 2.45e-06 | 0.00e+00 / 0.00e+00 | identical | -52.7953(0.0726) / -52.7953(0.0726) | 3.448(0.333) / 3.448(0.333) |
| L8_b4_U4_w20 | def | 9.02e-05 / 8.22e-05 | 4.39e+00 / 5.05e+00 | 1.45e-05 / 1.35e-05 | 3.11e+00 / 3.22e+00 | 2.00e+00 / 0.00e+00 | 1 | -52.7265(0.0503) / -52.9499(0.1157) | 3.121(0.229) / 3.536(0.359) |
| L8_b4_U4_w20 | log | up: abort: Unstable matrix; pp: ok | | | | | | | |
| L8_b4_U8_w10 | def | 1.46e-06 / 1.46e-06 | 4.24e-01 / 4.48e-01 | 3.07e-07 / 3.13e-07 | 8.04e-01 / 8.04e-01 | 0.00e+00 / 0.00e+00 | 1 | -34.3526(0.1693) / -34.3526(0.1693) | 5.140(0.949) / 5.140(0.949) |
| L8_b4_U8_w10 | log | up: ok; pp: abort: Unstable matrix | | | | | | | |
| L8_b4_U8_w20 | def | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |
| L8_b4_U8_w20 | log | up: abort: Unstable matrix; pp: abort: Unstable matrix | | | | | | | |

## Findings
Limitations: single seed per case, short runs from a random start (no thermalisation), and the upstream
code itself is on the edge for U=8 and for Nwrap=20 at dtau=0.1 (upstream Green max precision is O(1) for
U=8 and many upstream runs abort with "Unstable matrix"; default stabilization gives NaN for U=8 beta=40).
Aborts are stochastic (the trajectories differ), so up/pp abort patterns are only indicative.

Default stabilization (STAB default), Nwrap=10, U=4: mean precision 2.4e-11..5.2e-11 for both, max
5.7e-7..1e-5 for both, ratios 0.5..1.0; trajectories identical to 1e-8 over all bins. U=8 (Nwrap=10):
mean ratio 0.8..1.05, max ratios 0.9..1.6 of values that are O(1) in both; E and S(pi,pi) agree within
errors (E identical to the digits shown in most cases, differing by about 1 sigma in L6 b20 U8). Nwrap=20:
both builds are unstable (L6 b20 U4: mean 1.26e-4 / 1.21e-4, max 4.8 / 6.9). Verdict: **same order**
(mean within 3x, max within 10x) in every completed pair. No degradation with beta (4..40), L or U found.

STABLOG, Nwrap=10, U=4: Green mean precision pp/up = 1.0..1.7 on L6 (b4 1.35, b10 1.65, b20 1.14, b40 1.0),
but 4.7x (L8 b10: 5.9e-11 -> 2.8e-10) and 1.9x (L8 b20), 2.6x (L8 b4), max up to 14x (L8 b20: 4.4e-6 -> 6.2e-5;
L6 b20: 1.7e-6 -> 2.2e-5; L8 b10: 10.7x; L8 b4: 11.5x; L8 b40: 2.8x). U=8: tau mean up to 1.4x, max up to 8.7x (L8 b20 Green
max 0.9 -> 7.8), and pp aborted with "Unstable matrix" in 3 of 8 U=8 Nwrap=10 STABLOG cases where upstream
completed (L6 b4, L8 b4, L8 b10), none the other way. Observables (E, S) still agree within errors where both
complete (one case, L8 b20 U8, S 3.9 vs 7.8 differs because the trajectories diverged, both are
poorly thermalised/unstable there). Verdict: **worse** for STABLOG: max exceeds the 10x criterion at
U=4 for beta=10..20 (L6 and L8) and the mean exceeds 3x at L8 beta=10; degradation starts at moderate
beta and grows with L (L8) and with U=8 (stability aborts). A plausible cause (not tested) is that
sorting by initial column norm cannot replace the norm downdating of true column pivoting when the
columns span many scales (STABLOG scales enter through exp(L) terms).

Dense (non-checkerboard) L6 b20 U8 w10: default and log both same order (mean 2.3e-6 / 2.3e-6; max 0.32 / 0.30,
0.24 / 0.35).

Consequence for the port: the pre-pivot is acceptable for the default stabilization (CPU numerics
unchanged within noise) but not without precision loss for STABLOG; for STABLOG a GPU pivoted-QR substitute
(or a different strategy, e.g. pre-pivot with re-sorting after a first QR pass, or CPU fallback for that
step) needs a separate decision.

## Manager review (2026-10-27)
- **Reproduced:** L8 b20 U4 w10 STABLOG rerun from the same parameters and seeds gives bit-identical precision
  numbers (Green max 4.41e-6 up vs 6.18e-5 pp). The study sources on nividia are byte-identical to the laptop diff.
- **numerics-reviewer: approve.** The macro-off code is unchanged. IPVT semantics match ZGEQP3's JPVT. The phase/parity
  logic stays correct. The plain 2-norm is the right norm: in CGR/CGR2_2 the STABLOG scales are already folded into the
  matrix, unlike `decompose`.
- **Wording of the conclusion:** a *static* one-shot pre-sort (no norm downdating) is sufficient for the default
  stabilization and **insufficient for STABLOG**. Pre-pivoting schemes with re-ranking (two-pass sort, block or
  tournament pivoting, QRCP of a small panel on the host) were **not** tested.
- **Study-design caveat (ticket author):** U=8 at Δτ=0.1 with Nwrap=10/20 is outside upstream's own stable regime
  (upstream Green max O(1), aborts). Rows with U=8 or Nwrap=20 are not evidence either way. A follow-up should use
  Nwrap chosen so that upstream reaches max ≲ 1e-6 (e.g. Nwrap=5 at U=8) and several seeds.
- Scope: the macro covers CGR and CGR2_2 only. `decompose` keeps ZGEQP3 in the default scheme and pre-pivots under
  STABLOG (upstream).

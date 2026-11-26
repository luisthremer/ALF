# T1 on H100 (Helma jobs 1463455, 1463834; 2026-11-26; 0.22 GPU-h)

Production Hubbard (N_FL 2, Mz .F., flux, U 5.3, beta = L, dtau 0.1, submatrix N_delay 32, ham_h0 0), cold start, same seeds,
3 sweeps per Precision Green run. ALF_GPU_DEBUG=1 on all runs of job 1463834 (no effect: pgd_qr = pg_qr of job 1463455 to the digit).
Code: rd/cgalf-stab (gpu-port + CholeskyQR2 in dev_decompose, ALF_GPU_CHOLQR=1). Raw logs here; stabtime.txt from ALF_GPU_STABTIME=1.

Decompose cost = intercept a of the measured WRAPUR+decompose block time t(Nwrap) = a + b Nwrap (both flavours, median of 5).

| L (N) | method | a: decompose ms | b: wrap ms/slice | PG max Nwrap 5 / 10 / 20 / 40 | s per 3 sweeps Nwrap 10 / 20 / 40 |
|---|---|---|---|---|---|
| 24 (576) | pre-pivoted QR (now) | 6.9 | 0.13 | 1.2e-11 / 9.0e-12 / 8.3e-11 / 8.6e-9 | 8.38 / 6.67 / 5.79 |
| 24 (576) | CholeskyQR2 | 5.8 (-16 %) | 0.13 | 8.0e-12 / 8.1e-12 / 6.8e-11 / 9.7e-9 | 8.16 / 6.67 / 5.74 |
| 32 (1024) | pre-pivoted QR (now) | 14.3 | 0.49 | 3.2e-11 / 7.0e-11 / 9.6e-11 / 1.0e-8 | 23.57 / 18.92 / 16.54 |
| 32 (1024) | CholeskyQR2 | 10.6 (-26 %) | 0.49 | 6.7e-11 / 4.6e-11 / 9.9e-11 / 1.0e-8 | 22.86 / 18.44 / 16.24 |

PG mean equal for both methods within 10 % at every point. Cholesky failures in the sweeps: 0 (of 84-896 per run);
in the STABTIME hook 2 per run (its synthetic start states), handled by the QR fallback.

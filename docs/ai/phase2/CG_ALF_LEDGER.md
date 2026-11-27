# CG-ALF GPU-hour ledger

| Date | Machine | Job | GPUs x time | GPU-h | Purpose |
|---|---|---|---|---|---|
| 2026-11-25 | Helma H100 | 1452169 (E1 GEMM ceilings) | 1 x 6 s | 0.002 | library ceilings |
| 2026-11-25 | Helma H100 | 1452222 (E0 + E2) | 1 x 6 min | 0.10 | FLOP count, nsys timeline |
| 2026-11-25 | Helma H100 | 1452309 (E3) | 1 x <= 40 min | <= 0.67 | setup vs sweep, host profile |
| 2026-11-26 | Helma H100 | 1453558 (E4 S1 A/B) | 1 x ~25 min | 0.42 | LU-free phase on CUDA |
| 2026-11-26 | Helma (CPU on GPU node) | 1452907, 1457891 (det-phase large-N, CPU) | 1 x 1 h + <= 2 h | <= 3.0 | CPU check of the det-phase fix (GPU idle) |
| 2026-11-26 | Helma H100 | 1463455, 1463834, 1463941 (T1 stabilisation scan, cold + thermalised) | 1 GPU each | 0.40 | QR vs CholeskyQR2, Precision Green vs Nwrap at L = 24/32 |

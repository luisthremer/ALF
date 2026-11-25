# CG-ALF GPU-hour ledger

| Date | Machine | Job | GPUs x time | GPU-h | Purpose |
|---|---|---|---|---|---|
| 2026-11-25 | Helma H100 | 1452169 (E1 GEMM ceilings) | 1 x 6 s | 0.002 | library ceilings |
| 2026-11-25 | Helma H100 | 1452222 (E0 + E2) | 1 x 6 min | 0.10 | FLOP count, nsys timeline |
| 2026-11-25 | Helma H100 | 1452309 (E3) | 1 x <= 40 min | <= 0.67 | setup vs sweep, host profile |

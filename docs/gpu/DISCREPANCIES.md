# Discrepancies: research claims vs. code / measurement

Leads come from `ALF_GPU_research/R1–R3_report.md`. Each row records a claim that turned out wrong or imprecise.

| Date | Source (report §) | Claim | What we found (evidence) | Consequence |
|---|---|---|---|---|
| 2026-10-27 | R1 §2 table, calls/sweep | `Upgrade2` called 2·L·N_op·F times per sweep | 2·L·N_op calls; the flavour loop is inside Upgrade2 (12×12: 11 520 calls/sweep, instrumented gprof). Work per flavour as stated | none for flops; GPU kernel granularity is per site, both flavours |
| 2026-10-27 | R1 §2 table, calls/sweep | `UDV_State%decompose` 2·(L/Nw)·F per sweep | 3·(L/Nw)·F: the WRAPUR inside TAU_M adds (L/Nw)·F (measured 24 vs 16) | τ pass also needs device QRP |
| 2026-10-27 | R1 §6 | Ltau=1 "roughly doubles sweep cost" | τ pass = 34–41 % of the sweep, i.e. ×1.5–1.7 ([PROFILE.md](PROFILE.md)) | benchmark with and without Ltau |
| 2026-10-27 | R1 §2 scaling | local updates ≈5 of ≈30 F·L·N³ flops, stabilization dominates | by time the rank-1 ZGERU update is the largest kernel at N=1024 (33 % cb, 20 % dense) and its share grows with N; stabilization 26–55 %. Flop counts understate memory-bound BLAS-2 | delayed/submatrix updates matter more than the flop estimate suggests |

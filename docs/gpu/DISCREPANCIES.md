# Discrepancies: research claims vs. code / measurement

Leads come from `ALF_GPU_research/R1–R3_report.md`. Each row records a claim that turned out wrong or imprecise.

| Date | Source (report §) | Claim | What we found (evidence) | Consequence |
|---|---|---|---|---|
| 2026-10-27 | R1 §2 table, calls/sweep | `Upgrade2` called 2·L·N_op·F times per sweep | 2·L·N_op calls; the flavour loop is inside Upgrade2 (12×12: 11 520 calls/sweep, instrumented gprof). Work per flavour as stated | none for flops; GPU kernel granularity is per site, both flavours |
| 2026-10-27 | R1 §2 table, calls/sweep | `UDV_State%decompose` 2·(L/Nw)·F per sweep | 3·(L/Nw)·F: the WRAPUR inside TAU_M adds (L/Nw)·F (measured 24 vs 16) | τ pass also needs device QRP |
| 2026-10-27 | R1 §6 | Ltau=1 "roughly doubles sweep cost" | τ pass = 34–41 % of the sweep, i.e. ×1.5–1.7 ([PROFILE.md](PROFILE.md)) | benchmark with and without Ltau |
| 2026-10-27 | R1 §2 scaling | local updates ≈5 of ≈30 F·L·N³ flops, stabilization dominates | by time the rank-1 ZGERU update is the largest kernel at N=1024 (33 % cb, 20 % dense) and its share grows with N; stabilization 26–55 %. Flop counts understate memory-bound BLAS-2 | delayed/submatrix updates matter more than the flop estimate suggests |
| 2026-10-27 | R2 §2, R3 §3, brief §2 item 3 | STABLOG pre-pivots and uses unpivoted QR, "the natural GPU path" | Only `UDV_State%decompose` pre-pivots under STABLOG (`IPVT=1`, udv_state_mod.F90:503-542). CGR (cgr1_mod.F90:305-312) and CGR2_2 (cgr2_2_mod.F90:346, 405, 420) set `IPVT=0`, i.e. full ZGEQP3 QRCP, in every STAB variant. Together 33–41 % of CPU time | GPU needs pre-pivoted QR in CGR/CGR2_2 as well: new numerics, CPU accuracy study first (PLAN ticket 104) |
| 2026-10-27 | R1 §4 | Hubbard B_T is applied through polymorphic `RealExpOpT` → `ZDSLSYMM` | confirmed (profile stacks: `realexpopt_lmult` → `zdslsymm`; dense: `zlarcm`/`zlacrm` → DGEMM) | – |

# Phase 2 for the "smuc2 long prompt" session: cut the sync time, then fill the GPU

## Goal (in this order)
1. Bring down the sync/launch overhead of a sweep on the PVC (Intel Max 1550). This comes first.
2. Then utilise the GPU as fully as possible: every tile busy ~100 % of the time, high compute-unit (XVE) use, maximum node throughput.
3. Collect the numbers a GCS regular-project application needs (all 4 GPUs used, performance, scaling, GPU-h budget).

You may implement every option; do the small, low-risk changes first.

## Your work plan
Work through `docs/supermuc_large_prompt/wbs.json`, Parts 4 to 8 (new; Parts 1-3 and Deliverables are your earlier work).
Mark a task `done: true` only after its criterion is met. After each edit run
`python3 "<wbs-json skill>/scripts/wbs.py" check docs/supermuc_large_prompt/wbs.json` (the wbs-json skill is in this app; use the `/anthropic-skills:wbs-json` skill to find the script) and keep the existing ids.

## What we already know (from both sessions' docs, read them: FINDINGS.md here, `ALF_GPU_SUBMATRIX/docs/gpu/SUPERMUC_NG2.md` of the short-prompt session)
- Only ~10 % of a sweep is GPU kernel time; ~55-60 % is waiting at sync points; ~26-28 % oneMKL (QR partly on the host); ~9 % host MKL.
- ~76 000 synchronous operations per sweep at L=24; estimated 100-170 us overhead each on PVC vs ~12 us empty sync kernel on the A30.
- No `nowait` in any `Prog/gpu_*_mod.F90`; `!$omp dispatch` blocks; the SYCL shim calls `wait_and_throw()` after every LAPACK call.
- Per submatrix block (`gpu_lowrank_mod.F90` `blk_begin`/`blk_end`/`blk_flush`, lines ~1200-1290): ~15 serial operations, ~4 300 blocks per sweep at N_delay 32.
- Settings that cut the number of operations help the PVC a lot (N_delay 32->128: -29 %, Nwrap 10->20: -21..-32 %) but not the A30.
- Short-prompt session: best single chain 8.2 s/sweep (8 ranks, 14 threads, N_delay 128, AOT). You: 1.90 sweeps/s per node with CCS (32 ranks). The two were never combined.
- A30 reference: 5.25 s/sweep single chain, 0.371 sweeps/s with 4 walkers (MPS), 0.433 with 8.

## Rules
- Test jobs only (partition `test`), never submit a production job. Budget for this phase: at most 100 GPU-h (25 node-h); count your node-hours from job elapsed times and report them.
- The short-prompt session uses the same account (`pn93la`, its dirs `~/alf-gpu`, `$SCRATCH/alf`). Do not touch its files or jobs.
- Every code change: same-seed comparison against the ifx CPU build (Mz=F <= 1e-13, identical configs), ctest 53/53, timing table in `SYNC.md`. One commit per change on `production-mpi-intel-large-prompt`.
- Changes in shared files (`Prog/gpu_*_mod.F90`) must keep the CUDA path building. If you change behaviour there, note it for a later regression on the NVIDIA A30 (`ssh nividia`, one GPU process at a time); do not run it without the user's OK if the GPU is busy.
- Measure marginal cost per sweep from wall clock (runs of 1 and 4 sweeps), repeat once (noise +-15 %).
- Do not send anything to LRZ yourself; drafts only.
- Pushing was blocked from the sandbox before: commit locally and tell the user.

## Report
At the end, and after each Part, give the user a short report in German: numbered sentences, one fact per sentence, grouped under paragraph subheaders, global numbering. Include a before/after table of s/sweep per chain and sweeps/s per node at L=24 and L=32.

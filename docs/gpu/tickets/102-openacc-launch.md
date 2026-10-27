# 102 – OpenACC vs OpenMP target: launch overhead of per-site work (narrow)

- **Phase:** 1
- **Owner agent:** fortran-gpu-dev
- **Depends on:** 101 (done)
- **Status:** done (2026-10-27; reviewed and re-run by the manager)

## Goal
Measure on nividia (A30, nvfortran 26.9) whether OpenACC launches small kernels markedly faster than OpenMP target
for the per-site work of the DQMC update, the case where 101 showed that launch latency, not bandwidth, sets the cost.
The result feeds ADR-001 (programming model). PLAN.md rule: OpenACC is considered only if OpenMP is ≥ 1.5× slower on a
kernel that matters.

## Context
- 101 (`gpu_bench/`, results in `gpu_bench/results/nividia_A30_2026-10-27.txt`, summary in `docs/gpu/BENCHMARKS.md`):
  one OpenMP target region ≈ 12–15 µs, 7–9 µs with `nowait`; cuBLAS call ≈ 5 µs; scalar D→H→D round trip 16–23 µs;
  rank-1 + round trip ≈ 28–34 µs for N ≤ 400. See `sec_upd` and `sec_opv` in `gpu_bench/src/bench_sec.F90`
  (rank-1 kernel, `nowait` variant, round-trip variants at lines ≈ 380–405).
- What the DQMC single-site update does per site (sequential over sites, L·N sites per sweep): read G(i,i) (and a
  few entries), compute the acceptance ratio, draw a random number, accept/reject, and on accept apply a rank-1 update
  of G (N×N complex). Today the decision is on the host.
- Two ways to cut the cost, both to be measured here:
  1. a cheaper launch model (OpenACC `async` queues vs OpenMP `nowait`/`depend`);
  2. **decision on the device**: ratio, comparison with a pre-generated random number and the conditional rank-1
     update all in device code, so the host only enqueues kernels and never waits per site.

## Measurements (all complex(8), device-resident G, N = 64, 144, 256, 400, 1024)
For each programming model, **OpenMP target** (`-mp=gpu`) and **OpenACC** (`-acc=gpu`), same nvfortran, `-gpu=cc80`:
1. Empty kernel: synchronous launch, and async (`nowait` / `async(q)`) launches of 1000 kernels followed by one wait,
   per launch.
2. Rank-1 update kernel (own loop, as 101's `r1_ker`): synchronous, async; per update.
3. Rank-1 + host decision (read G(i,i) to host, decide, launch update): per site.
4. **Device decision**: per site one kernel (or two, ratio + update) that reads G(i,i), compares with a
   device-resident random number and conditionally applies the rank-1 update, launched async, one wait per 256 sites.
   Also as one kernel per site that is a single gang/team if that is faster. Accept rate ≈ 50% (fixed random input).
5. Op_V per-site launches (101 `sec_opv` variant "one launch per site") in both models.
Medians over ≥ 5 repetitions; also report the cuBLAS `zgeru` (async on the same stream/queue) for reference.

## Correctness
Every variant is checked against the host reference (`bench_ref`), as in 101: max abs diff ≤ 1e-12 for rank-1
sequences (256 updates) and the device-decision variant must reproduce the host-decision result exactly for the same
random input (same accept/reject pattern, G diff ≤ 1e-12). A FAIL is a bug, not a tolerance to adjust.

## Files allowed
- `gpu_bench/src/bench_launch.F90` (new; standalone program or new section, your choice, keep it simple)
- `gpu_bench/Makefile` (new target(s), e.g. `make launch MODEL=omp|acc`; existing targets unchanged)
- `gpu_bench/README.md` (new section), `gpu_bench/results/nividia_A30_launch_<date>.txt` (new)
Reuse `bench_util`/`bench_ref` modules by `use`; do not change them. Anything else: stop and report.

## nividia rules
Testing only: develop on the laptop and rsync `gpu_bench/` to `~/alf-gpu/work/t102/`. Write only inside
`~/alf-gpu/work/t102/` (`TMPDIR` = `~/alf-gpu/work/t102/tmp`). `. ~/alf-gpu/env/nvhpc.sh`. Budget: 2 CPU cores and the
GPU; check `nvidia-smi` first and wait if someone else uses it. Runs ≤ 15 min, release the GPU afterwards, install
nothing, no leftovers outside the work dir. Ticket 108 runs CPU-only jobs concurrently.

## Acceptance
- Both models build and all checks PASS on nividia.
- Results file + a table in the report: per variant and N, µs per launch / per update / per site, OpenMP vs OpenACC
  ratio, and the projected s/sweep for a single-site sweep (2·L·N sites, L = 40) for N = 144, 400, 1024, compared with the
  16-core node figure from `docs/gpu/BENCHMARKS.md` (0.038 s per chain-sweep at N = 144).
- Statement: is OpenMP ≥ 1.5× slower than OpenACC on items 2–4? Does the device decision remove the round-trip cost?

## Out of scope
ALF source (`Prog/`, `Libraries/`), CUDA Fortran or CUDA Graphs, ifx/OpenACC (no OpenACC on ifx), multi-chain batching.

## Resolution (manager, 2026-10-27)
- Reviewed `gpu_bench/src/bench_launch.F90` and the Makefile target. Re-ran both models on nividia (N = 144, 400, 1024):
  33/33 checks PASS each, all timings within noise of the agent's results file
  (`gpu_bench/results/nividia_A30_launch_2026-10-27.txt`).
- Result: the difference between the models is the **asynchronous launch**: OpenMP `nowait` ≈ 10 µs per target region
  (task creation; with or without `depend`), OpenACC `async` ≈ 2.6 µs. Synchronous launches, cuBLAS and the host
  decision are model-independent (host decision ≈ 15–18 µs per site for N ≤ 400).
- Device-side decision (2 async kernels per site): OpenACC 6.5–8.8 µs per site (≈ 2.2× better than the host decision),
  OpenMP 17 µs (no gain). Even the OpenACC variant projects to 0.075 s/sweep for the update part alone at N = 144,
  i.e. ≈ 2× slower than the 16-core node's whole sweep (0.038 s per chain-sweep). At N = 1024 the update is
  bandwidth-bound and the model does not matter (±20 %).
- Caveat for any later use: the benchmark's rank-1 vectors are fixed; in ALF they are a row and a column of G, so a
  real device decision also has to copy them before the update (one more short kernel, or inside the decision kernel).
- PLAN.md's 1.5× criterion is met for async launches and the device decision at N ≤ 400. The decision for ADR-001 is
  taken there.

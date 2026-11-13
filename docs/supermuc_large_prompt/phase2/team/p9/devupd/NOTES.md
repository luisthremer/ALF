# p9 / devupd: on-device update loop (decision 44): feasibility check, stopped before implementation

2026-11-13, SuperMUC-NG Phase 2, ifx/oneMKL 2025.3.0. C = measured (job id or file), E = estimate. No source under Prog/ or Libraries/ was edited. Nothing was committed.

## Result
The bitwise gate G1 cannot be met by a device evaluation of the decision and the M update. The brief allows a tolerance only after the owner agrees, so I stopped before writing any code.
The timing target G4 (update phase <= 0.5x) is also out of reach by estimate, because the brief's premise does not match the production path (section 1).

## 1. What the production path does today (C, source + p7/p8 data)
- Production (`Sequential_scheme = submatrix`, N_FL_eff = 2, single-site vertices) runs the **block mode** (`gpu_lr_blockmode`, `gpu_lowrank_mod.F90` lines 1290-1912). It does NOT use the per-site path that the brief describes (gather kernel + D->H of 2k+1 values per flip). That path (`gpu_lr_gpp` without `blk_open`) is used only outside block mode.
- Per block of b = N_delay vertices: one gather of G0(B,B) with a D->H copy and **one** wait. Then all b proposals (blk_gpp, ratio, decision, M update) run on the **host**, with no device access. Then the end-of-block device chain: scale, fgather, 2 GEMM per flavour, async.
- Per-site host round trips therefore do not exist today. The 112k async ops / sweep of p7 are per-block ops (13 x 8,640 blocks at N_delay 32).
- L1, N_delay 64, V7 (p8/blk logs/p3_summary.txt, job 561035): update 1.24 s/sweep = blk_begin 0.14 + blk_end 0.505 + slice level 0.276 + remainder ~0.32. The remainder is the host per-site work, ~1.2 us per site visit, cexp/cpow/csqrt included.
- Node layout, N_delay 64 (p8/blk, job 561098): the update phase is 2.46 s and device/tile bound. p8 showed that removing host work does not shorten it.

## 2. RNG order (C, source)
- `ranf_wrap` is the intrinsic `RANDOM_NUMBER`; its state is `RANDOM_SEED(GET/PUT)` (Libraries/Modules/random_wrap_mod.F90).
- Per vertex in the production config (Hubbard, `Predefined_Int_U_SUN`, field type 2, `Propose_S0 = .F.`, S0 = 1) there are exactly **3 draws**:
  - `nsigma%flip` draws `nranf(3)` (one ranf). The new field value can be computed on the host in advance from that draw.
  - The T0 test `1.5 > ranf_wrap()` always draws and always passes.
  - `site_update` draws once for the decision.
- The count is fixed, so pre-drawing 3b numbers per block in this order is exact, and no restore/redraw is needed. The count varies only with `Propose_S0 = .T.` or field types 3/4; those would need the save/restore scheme.

## 3. Why bitwise is impossible (C)
**Host code** (objdump of `p8/obs/build_aottt/Prog/*.o`, flags `-O3 -fp-model=fast=2 -no-prec-div -xHost`; `logs/blk_gpp_host.s`):
- `blk_gpp` is vectorised with AVX2/AVX-512 gathers and FMA (`vfmsub231pd` x63).
- The m-sum runs in 4-lane partial sums plus multi-accumulator unrolled versions, reduced by `vextractf128` + `vaddpd`.
- `gpu_lr_accept`: `vfmadd231pd` x65, `vfnmadd231pd` x18, `vaddsubpd` x14.
- `site_update` calls the host libm `cpow` (for `Ratiotot**dble(N_SUN)`; N_SUN = 2 in ref_L24), `csqrt` (phase) and `cexp`.

**Bench** `bench/dubench.F90`, job **561178** (logs/p9du_b1.561178.log): the same Fortran source is run on the host and in an OpenMP target region, with ALF's flags, AOT for PVC, 200,000 random samples each.

| operation (as in ALF) | differ (k=24) | max rel |
|---|---|---|
| `z**dble(N_SUN)` (cpow) | 93 % | 9e-16 |
| phase `ph*z/sqrt(z*conjg(z))` | 10 % | 7e-16 |
| xv = `d/(1+v*d)` | 41 % | 2e-15 |
| cexp | 87 % | 7e-16 |
| blk_gpp chain (W = -M Gcol, gpp) | 66 % (26 % at k=8) | 2e-14 |
| Ratiotot (ratio code of site_update) | 97 % | 2e-15 |
| bt dot (gpu_lr_accept) | 95 % | 2e-13 |
| single M element `M + z*a*b` | 44 % | 2e-14 |
| Weight | 83 % | 1e-11 |

Even one complex multiply-add differs in 44 % of the cases: host and device contract FMAs differently. The libm functions (cpow, csqrt, cexp) are different implementations on the two sides.

Exact emulation would mean copying, on the device, the host compiler's multi-versioned vectorised reduction order and FMA pattern, plus glibc/libimf `cpow`/`csqrt`. That is compiler-version dependent and is not a sound basis.

A host-verified speculative scheme also does not help: verifying the decisions exactly needs host-exact gpp and M, which is the whole current host work.

## 4. Timing outlook (E, partly C)
- Bench part 2 (job 561178): a straightforward one-team OpenMP kernel for the sequential Metropolis of one block is 15-35x slower than the host (b=64: 4.3-4.9 ms vs 0.22-0.29 ms; b=32: 0.88 ms vs 0.025 ms). This is a pessimistic bound; a hand-written SYCL sub-group kernel would be much faster (E: 30-60 us/block).
- Even with that, the best case at L1 is update ~0.85-0.95 s vs 1.24 s (~0.7x, E). The per-block device chain stays, and the oneMKL GEMM k must be known on the host (a per-block D->H of k, or padding to k=b, which is not bitwise and does 2-3x more gemm2 work).
- At the node layout (4 ranks/tile) the Metropolis work would move from idle host cores onto the shared, already saturated tile, so it would likely be slower (E).
- So G4 (<= 0.5x) is not reachable at either layout.

## 5. Budget, jobs, files
- Job 561178 (bench, 1 node, preempt, 22:39:09-22:39:26 by its log, ~1 min with allocation): **~0.02 node-h**. No other jobs.
- Local files: `bench/dubench.F90`, `jobs/b1.sh`, `logs/` (job log, err, host disassembly of blk_gpp).
- Remote: `$SCRATCH/large_prompt/team/p9/devupd/{bench,jobs}`.

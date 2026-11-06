# Known issues

Branch `gpu-port-delayed-submatrix`, state 2026-11-06. Issues are listed with their evidence and status.

## 1. Dense hop and `Symm` matrices are slow on GPUs with weak FP64 (regression of the speed-up round)

Reported by the fork session "GPU performance bottlenecks (fork 2)" after profiling on the machine `maksim` (RTX 4060).
The statement is quoted verbatim:

> **1. The dense hop and `Symm` matrices I added hurt on this card.**
> - **Evidence:** the `nsys` run on `maksim` showed the GPU 96% busy. 13 of the 15 s of kernel time was ZGEMM, about 2 ms per N=400 product, which is the 4060's fp64 peak (~0.26 TFlop/s). That is not launch overhead.
> - **Cause:** each dense hop is an N³ product, so it is cheap on the A30 but costs about 2 ms here. The sparse checkerboard kernels are O(N²) and take about 10 µs.
> - **Fix, tested locally:** I switched the dense path off in a copy at `~/luis_gpu_test/src_nodense`, which is not in git. Wall time for the 20×20, 6-sweep run dropped from 46.4 s to 21.0 s for submatrix, from 49.1 s to 23.7 s for delayed, and from 53.0 s to 27.6 s for standard. Acceptance is identical.
> - **Per sweep:** that is about 3.5 s per sweep including start-up and per-bin output. The `nsys` run gives about 3.2 s per sweep of pure GPU busy time, which is close to his 3 s.
>
> **2. What remains is double-precision work.** The kernels are now 89% busy. In the profile, the ~15,800 rank-32 update products of the block flush (`Kernel2`, N=400 by 32) take 3.7 s of 6.5 s. The QR in the stabilization takes 0.85 s, at about 1.3 ms per call.
>
> - **Why the delayed and submatrix schemes don't help here:** they make the update level-3, but the total flop count stays about 8·N² per accepted flip. On an fp64-weak card the arithmetic itself dominates, not the memory or latency.
> - **Not tested:** why his code does 3 s on the same machine. I don't know his precision or algorithm, so I can't say whether he uses lower precision somewhere or simply has fewer flops.
>
> I haven't changed git. The dense paths are in the committed branch with a fixed `HOPD_NMAX = 800`, which is too large for consumer GPUs. The real fix is to pick that threshold from the device's measured fp64 rate; on the 4060 the best setting is effectively "never dense".

Status: **open**. The dense paths (`hopd_*` in `gpu_wrap_mod.F90`, the dense `Symm` in `gpu_obs_mod.F90`) are in this branch with the fixed
threshold `HOPD_NMAX = 800`. The same N³ argument applies to the dense `Symm` (item 10 in `PARKING.md`: sparse `Symm`). On the A30 the
dense matrices helped; the 4060 numbers above are from the fork session and were not re-measured here.

## 2. Other known issues (from the 2026-11-06 test round)

- **nvfortran miscompile:** `blk_end` in `gpu_lowrank_mod.F90` gives wrong results at `-O2` and above (`-O1` is correct). It is pinned to `-O1`
  with `!pgi$r opt=1`; the cause is not isolated. ifx is not affected in any run so far.
- **GPU coverage:** only Hubbard-type models with diagonal rank-1 vertices run on the device. Plain_Vanilla, Kondo, tV, Spin_Peierls,
  Z2_Matter, LRC and `LR_Jz_Test` are refused at setup (reason in `info`) and run on the CPU path. Also refused: global moves, Langevin/HMC,
  MALA, projector, `Propose_S0`, time-dependent couplings, `Lowrank_cross_slice`, MPI builds, and STAB1/2/3.
- **`Nwrap`:** larger values lose precision quickly (12×12, β = 5: Green-function error 1e-10 at `Nwrap` 10, 2e-4 at 20, abort at 30).
- **tV with the default parameter file** (β = 5, V = 4, `Nwrap` 10) is numerically marginal on every build; use `Nwrap` 4.
- **Intel B580:** raw FP64 throughput equals the A30, but oneMKL `zgemm` (about 0.9 TFlop/s at N = 1024) and `zgeqrf` (47 ms at N = 1024)
  are slow, so stabilization and dense products are 3–8× slower than on the A30 and the gap grows with N.
- **Memory:** about 0.4 GB + 1650·N² bytes at β = 5; 62×62 fits a 24 GB A30, 63×63 does not.
- **Unvalidated:** Ponte Vecchio (SuperMUC-NG Phase 2), Helma and SuperMUC builds, and every model other than Hubbard on the device.
- **Output:** per-bin output (`Pr_obs`) is slow at large N (38 s per bin at 56×56) and dominates short runs.

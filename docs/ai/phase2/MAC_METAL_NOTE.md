# Metal on Apple GPUs for ALF: feasibility note (1 page)

Date 2026-11-27. Question (owner, Q19): can the Apple GPU (M1 Pro here) do anything useful for DQMC?
Companion to the Accelerate backend (`configure.sh`, `ALF_MAC_BLAS`, results in
[`MAC_ACCELERATE.md`](MAC_ACCELERATE.md)).

## Facts

- **Apple GPUs have no FP64.** Metal Shading Language has no `double`; M1 to M4 execute only FP32, FP16,
  and integer. There is no emulation in the driver. MPS (`MPSMatrixMultiplication`) and
  MPSGraph support FP32/FP16/BF16 only.
- ALF stays correct only in FP64. The stabilisation (QR with column pivoting of `U D V`, scales spanning
  $e^{\pm\beta\,\Delta}$) needs the 53-bit mantissa; the measured `Precision Green` of the FP64 code is
  $10^{-14}$ to $10^{-12}$ on the CPU path. FP32 has $\epsilon\approx 6\times10^{-8}$, so a single product already
  carries an error seven to nine orders of magnitude above what the FP64 run reaches (and the `D` scales can leave the FP32 range $10^{\pm38}$ for large $\beta$).
- Double-float (two FP32 numbers, about 44 bits) or Ozaki-style FP64 emulation on FP32 units gives
  FP64-class GEMM only at roughly 1/10 to 1/20 of FP32 speed. The M1 Pro GPU peaks at about 5 TFLOP/s FP32,
  so emulated FP64 lands near 0.3 TFLOP/s (estimate, not measured), about what AMX already gives the CPU: measured ZGEMM N=1024 on this M1 Pro is 233 GFLOP/s on 1 thread and 368 on 4 threads through Accelerate (see MAC_ACCELERATE.md).
  There is no QR/LAPACK on Metal at all; Householder QR would have to be written by hand.
- Unified memory removes the PCIe copy, but the ALF kernels are small and latency-bound
  ($N\le 10^3$ matrices, hundreds of calls per sweep); launch and synchronisation costs per call
  (about 10 to 50 us) are comparable to the CPU time of a single $N\!\times\!N$ ZGEMM at $N<500$.

## What could still run in FP32 (and why it does not help)

| Candidate | Verdict |
|---|---|
| Green function, `UDV` stabilisation, wrap products | No: FP64 required (above). |
| Metropolis ratio $\det$ and local update $G \mathrel{+}= \dots$ | No: error accumulates in $G$ between stabilisations. |
| Observables (equal-time correlators, FFT) | Possible in FP32 after a FP64 round-trip, but cost is $<5\%$ of a sweep. Not worth a second code path. |
| Preconditioner / initial guess, then FP64 refinement (mixed precision) | Mixed-precision refinement needs the FP64 residual GEMM, which is the cost we try to avoid. Only research-grade for DQMC, no known stable scheme. |

## Conclusion

Do not build a Metal backend. The Apple GPU cannot run any ALF kernel at the required precision, emulated FP64 is
no faster than the CPU's AMX through Accelerate, and no LAPACK exists there. The useful Mac speed-up is the
CPU-side Accelerate (AMX) backend, which is what `configure.sh` now selects explicitly. Revisit only if
a future Apple GPU generation exposes FP64 (none announced) or if a stable FP32 DQMC scheme appears in the literature.

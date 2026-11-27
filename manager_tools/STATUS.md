## Now
- **P2.4 step 2 done:** all three vendor backends are merged into `Code-modularity-GPU`: Intel oneMKL/SYCL (4bd27bc5), CUDA (75b9c881), AMD HIP (323b2801, experimental: R2 open). Each was checked on its device against the CPU path (Max 1550, H100, MI210). Next: step 3, integrate (one vendor-neutral sync hook, one tolerance rule incl. size scaling for Mz, G7).
- Running: gpu-intel (speed switches, S-c), gpu-amd (follow-up: speed vs gpu-port, L = 64), t3-opv-symm (fast `Op_V_is_symmetric`, CI fix), gpu-speed (GPU vs CPU node, Helma queue).
- Live progress page: [ALF GPU Phase 2](https://claude.ai/artifact/PmMk2VMb8iUzXtnEYUcHc3) ("Needs you" is live).

## Waiting for you
1. A prepared session card: "QR_PREPIVOT variants 2 and 3" (CPU-queue item 6). The app does not let the manager start sessions on a timer; click the card or say "start it".
2. Check the overnight decisions on the progress page (A/B classes, AMD tolerance rule, Mz L = 24 as an algorithm property); overrule any you disagree with.

## Done since the last owner message (night 26/27 Nov)
- Merged: B batch 2 (H5, M5, M6, MK), LU-free CGR phase on the CPU line (sweeps 1-5 % faster), Timers (temporary), Intel, CUDA and AMD backends, T1/T4/T5/P2.5 docs, BRANCHING archive table, LEVERAGES.md (measured gains for the report to ALF).
- AMD R2 (gfortran 14 offload, class A for GPU=amd only): mitigated; AMD production checked clean (production code at L = 36 agrees with the CPU better than two CPU builds agree).
- T1 thermalised: Nwrap 20 safe up to L = 32, -20 % per sweep. CUDA: 5 % faster than gpu-port at L = 24/32. Intel: h0 = 1 measurements 2x faster than gpu-port.
- Your real email was removed from the fork (one commit rewritten); all sessions use the placeholder.

## Next 3 steps
1. P2.4 step 3: integrate the three backends (vendor-neutral sync/finalize hook, common tolerance rule, G7), then step 4 consolidate.
2. Merge the follow-ups as they come: t3-opv-symm, gpu-intel-speed, AMD speed/L = 64, gpu-speed status.
3. P2.6 (delayed/submatrix for every model on the GPU) and CG-ALF T12 (finish the walker) can start now that the vendor code is in.

## Open branches and running jobs
- `cm/feat/opv-symm-fast` (CI failed on a FP trap in its test; being fixed). gpu-intel-speed, AMD follow-up, rd/gpu-speed in work.
- Cluster: SuperMUC (gpu-intel), Helma (gpu-speed), amd (gpu-amd). Plan: [P2_4_GPU_PLAN.md](https://github.com/luisthremer/ALF/blob/Code-modularity-GPU/docs/ai/phase2/P2_4_GPU_PLAN.md)

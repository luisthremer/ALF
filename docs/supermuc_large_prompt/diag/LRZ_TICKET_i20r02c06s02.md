**To:** LRZ Servicedesk (SuperMUC-NG Phase 2)
**Subject:** Suspected faulty GPU tile on i20r02c06s02 (GPU 0, tile 1): silent wrong results

*Transparency note: this message was co-authored with Claude Opus 5.5 (Anthropic), an AI assistant, at medium reasoning level. It drafted the text and helped with the analysis and the stress test. All job IDs and numbers come from our own job logs.*

Dear LRZ team,

We believe one GPU tile on node **i20r02c06s02** returns silently wrong floating-point results. In flat device hierarchy it is device 1 (card 0, tile 1). Please check the node and take it out of service if you can confirm the fault.

**Project / user:** pn93la / di38zuw
**Software:** ALF (finite-temperature quantum Monte Carlo), OpenMP offload with ifx/oneMKL 2025.3 (modules intel/2025.3.0, intel-mkl/2025.3.0, intel-mpi/2021.15.0). The Level Zero backend is used with OMP_TARGET_OFFLOAD=MANDATORY.

**Symptom**
Our code checks its results internally. On this node, device 1 sporadically produces a completely wrong matrix: the mean deviation per entry is of order 1-10, whereas the normal accuracy is about 1e-14. No error, exception or driver message is reported. The program only stops because of its own consistency check.

**Evidence**
1. All 7 aborts in our production runs since 2026-11-11 happened on device 1 (MPI rank mod 8 = 1).
2. The same input (8 ranks, 1 rank per tile, fresh start from fixed seeds) gave these results:

   | Job | Node | Visible devices | Result |
   |---|---|---|---|
   | 560331 (2 runs) | i20r02c06s02 | all | abort at rank 1 |
   | 560620 run A | i20r02c06s02 | all | abort at rank 1 (deviation 15.68) |
   | 560620 run B | i20r02c06s02 | all | abort at rank 1 (deviation 10.05) |
   | 560622 run C | i20r01c01s06 | all | completed, clean |
   | 560623 run D | i20r02c06s02 | ZE_AFFINITY_MASK=0,2,3,4,5,6,7 (tile 1 hidden) | completed, clean; result identical to run C |

   Apart from these events the code is bitwise reproducible. 128 identical replicas on other nodes gave bit-identical results.
3. Dedicated stress test: a fixed sequence of complex double-precision BLAS and LAPACK operations at N=1024 (zgemm, zgetrf, zgeqrf, zunmqr, ztrsm), repeated on each device. Every result is compared bitwise with the first iteration.
   - Positive control (job 560640, one corrupted matrix entry injected): detected on all 8 tiles in exactly the injected iteration. All other iterations were bit-identical.
   - Suspect node (job 560637, 8 ranks, 1 per tile, 20 min): **[TO FILL: bad iterations per device, first occurrence]**
   - Suspect node (job 560638, 32 ranks, 4 per tile, CCS 4, 15 min): **[TO FILL]**
   - Control node **[node]** (job 560639, 20 min): **[TO FILL]**
4. xpu-smi RAS/ECC counters (reset, cache errors, memory errors correctable/uncorrectable) show "N/A" for our user account, so we could not read them.

**Request**
- Please read the RAS/ECC and error counters of GPU 0 (tile 1) on i20r02c06s02 and run your hardware diagnostics.
- Please drain the node if the fault is confirmed.

We have excluded the node from our jobs in the meantime (`-x i20r02c06s02`). Logs and the stress-test source are available on request (`$SCRATCH/large_prompt/diag/`, `$SCRATCH/prod_2026-11/work_node/`).

Kind regards,
**[Name]**

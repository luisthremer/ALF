# P2.6 tempering on the GPU path (TEMPERING without PARALLEL_PARAMS)

Feature `p26-tempering`, 2026-11-27. Before it, `gpu_supported` refused every tempering run with a state exchange; PP (independent
groups, no exchange) already ran.

## What changed

The device owns the Markov state during a sweep (stacks `udw_*`/`udv_*`, `G_d`, the tracked det(U) per slot, the field mirror
`nsigma_d`); `Exchange_Step` (Global_mod) is host code that reads and swaps the host copies. `Exchange_Step` stays unchanged;
main.F90 brackets it on every rank:

- `gpu_temper_pre` (gpu_stab_mod): stops if a low-rank update is pending (the call site is a sweep boundary), D->H of
  udvst/udvl/udvr (`gpu_get_udv`) and of G (`gpu_g_host`).
- `gpu_temper_post`: `nsigma_d = nsigma%f` + H->D; the stacks back through `gpu_stab_put_udv` (re-derives every tracked
  det(U) by QR); G marked host-written (`gpu_g_host_wrote`, uploaded at its next device use).
- `gpu_supported`: tempering refused only on GPU=amd (not yet checked there). New negative control `ALF_GPU_FAULT=temper`
  (debug runs): keeps the old `nsigma_d` after the exchange.

All MPI stays on host buffers; ranks map to devices as before (`rank mod devices`, or one visible device per rank). Cost per
exchange: 2 (NSTM + 2) N^2 x 16 B transfers and (NSTM + 2) device QRs, plus the host `Exchange_Step` itself (with
`Tempering_calc_det = .T.` two host passes of WRAPUL and a host CGR). T12's refusal of N_walkers > 1 with tempering stays.
GPU=none: every change is under `ALF_GPU`; preprocessed main.F90 identical to the base for noMPI / MPI / TEMPERING / PP
(including `__LINE__`).

## Checks

Runner `../../p25_runs/p25_run.py --groups N --mpi N` with `cases_tempering.py`: `Temp_k` with ham_U = 4 - 0.5 k, one rank per
replica, exchange every sweep, GPU path vs `ALF_GPU_DISABLE=1` with the same seeds, every `Temp_k` compared (confout
byte-identical, other files incl. the tempering acceptance within max(1e-11, 100 x Precision Green mean)). `P25_MPIRUN` selects the
launcher (`mpiexec` for Intel MPI).

| where | build | result |
|---|---|---|
| melmac, GPU=host (GNU 14, Open MPI 5) | Tempering | 2 replicas: the 8 "works" cases OK (Mz std / delayed 4 / submatrix 4, SU(2) submatrix, Ltau 0, 3 exchange steps, Spin_Peierls 6x6 submatrix, Hubbard 8x8 submatrix 16); worst 2.5e-10 against tol >= 3.2e-9; 4 replicas: 4/4 OK (std, delayed, submatrix, 3 exchange steps), worst 1.0e-9 against 2.5e-9; `Acceptance Tempering` identical (e.g. 0.125 per replica in tp-mz-std) |
| melmac | Tempering LOG | 3/3 OK |
| melmac | Tempering | P2.5 cases (`../../p25_runs/cases_tempering.py`, refused before): 2/2 works |
| Max 1550 (GPU=intel AOT, ifx 2025.3, Intel MPI, one tile per rank, job 566912) | Tempering | 2 replicas 10/10 OK, 4 replicas 5/5 OK (incl. 8x8 submatrix 16); worst 8.3e-10 against 6.7e-9 |
| H100 (GPU=nvidia, nvfortran 25.7, 2 ranks on 2 GPUs, Helma job 1470493, 2 min) | Tempering | one GPU per rank (`CUDA_VISIBLE_DEVICES` = local rank): 10/10 as expected (8 works OK, nodet stops on both paths, negative control stops); worst 5.6e-10 against 2.7e-9; unmasked (`rank mod devices`, ranks on OpenMP devices 0 and 1): tp-mz-std OK |

- Negative control `tp-neg-temper` (`ALF_GPU_FAULT=temper`): the run stops on every platform (precision check: "Unstable matrix",
  or "G difference 24.4 above 10"), FAIL as required.
- `tp-mz-nodet` (`Tempering_calc_det = .F.`): the S0-only weight is wrong for replicas that differ in U, so the CPU path stops too
  ("calculation unstable"); expectation "stops", both paths stop alike. (The first melmac pass ran it with expectation "works":
  the two MISMATCH lines in `results_gpu_host.tsv` are that run.) A model where S0-only exchange is valid was not run.
- Clean tree on melmac, on GPU tip aa56b093: GPU=host ctest 146/146. After the rebase onto 45cde353: the CI job command OK;
  GPU=none outputs byte-identical to the base in 16 scenarios; Tempering STAB3 (open since `p26-stab`): std/delayed/submatrix OK
  (worst 1.6e-10), negative control stops; GPU=host ctest 108/154: the 46 GPU unit tests stop at the P2.8 guard Z3
  ("NBin = 0"), from the base, not from this feature.

Files: `cases_tempering.py`, `results_gpu_host.tsv`, `results_max1550.tsv`, `results_h100.tsv`, `../../intel/job_temper.sbatch`.

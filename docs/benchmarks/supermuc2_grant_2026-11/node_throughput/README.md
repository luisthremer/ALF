# Node throughput test: GPU node vs CPU node (SuperMUC-NG Phase 2, 2026-11-27)

Layout agreed with the owner: artifact https://claude.ai/artifact/VT1hhGa6pTxgtbZN29CVTg (`../cpu_vs_gpu_layout.html`).
Submitted 2026-11-27 ~10:35 as jobs 566793-566804 (smoke job 566791, GPU L=16, results in `runs/smoke`).

Cluster directory: `/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/` - `jobs/nt/` (these scripts), `runs/` (one directory per run).

| File | Purpose |
|---|---|
| `run_layouts.sh` | Slurm job: all layouts of one side (cpu / gpu = 410c7b04 / gpub = intel-batched d1060305) at one L, one full node |
| `submit_all.sh` | the 12 submissions (GPU jobs with REPEAT=2) |
| `make_params.py`, `template_parameters_RUN_L24.txt` | parameters: production file with L, N_Phi, Beta=4, Ltau=0, h0=0, U=5.25, NSweep=1, NBin=11, CPU_MAX=0 changed |
| `pin.sh` | per-rank NUMA pinning (numactl); GPU: tile mask + the tile's 14 nearby cores |
| `watch.py` | logs mtime (ns) of `Temp_0/data.h5` every 2 ms -> bin write times |
| `analyse.py` | `results.csv`: time per sweep (10 timed sweeps, error = s.e.m.), compute vs write time per bin, node throughput, GPU utilisation/power |
| `plots/plot_all.py` | the 5 figures + `summary.md` (SciencePlots); `plots/make_fake_data.py` + `plots/fake_out/` = FAKE data, plots only |

Timing (from the smoke job): a bin is written as a short burst of HDF5 writes, so the write events form NBin + 1 = 12 bursts
(file creation + 11 bins). Per timed sweep: t_comp = gap between bursts (sweep, measurements, bin reduction), t_io = burst
length. Reported cost per sweep = t_comp + t_io/10, because production writes one bin per 10 sweeps.

Run the analysis and plots:
```
ssh supermuc2 'python3 /hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/jobs/nt/analyse.py /hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27/runs' > results.csv
python plots/plot_all.py results.csv plots/out        # needs matplotlib + scienceplots
```
Smoke job note: at L=16 the production layout (4 walkers/tile, 3 threads) was slower than 4 walkers/tile with 1 thread
(1.18 vs 0.66 s per sweep) and had unusually long bin writes; the GPU jobs therefore run every layout twice.

## Parallel CPU part and preliminary report (2026-11-27, owner: "run all simulations in parallel")
The four single-node CPU jobs were cancelled (566801-566804) and replaced by `multi.sh` jobs (`submit_cpu_parallel.sh`,
566809-566814): one job with 4 nodes runs one layout per node at the same time (`srun` steps; mpiexec with fork bootstrap
on the local host), because the test partition allows only 5 running jobs per user.

`./prelim.sh` builds `report/report_prelim.pdf` from all runs finished so far (`./prelim.sh fake`: fake data).
`report/report.tex` explains the layouts (x axes) in Tables 1-2; results tables come from `plots/make_tables.py`.

Update 11:00: the CPU runs are short (L=16, 1 thread: 2.8 s/sweep), so no CPU job needs `general`. Waiting jobs
566811-566814 (test/general) and 566818-566819 (preempt, estimated start next day) were cancelled and replaced by two
single-node test jobs, 566820 (cpu L=32) and 566821 (cpu L=36), which fit into gaps more easily.
Update 11:05: CPU cost per sweep grows ~L^6 at fixed beta (L=16 -> 24, 1 thread: 2.8 s -> 22 s), so the L=36 job (566821)
was cancelled and the remaining CPU runs split into one single-node job per layout: L=36 t1..t8 = 566823-566827,
L=32 t6, t8 = 566828-566829 (566820 keeps L=32 t1, t2, t4; it may hit its time limit in t6, which is then simply missing there).

Update 11:36: `peek.sh` shows runs in progress (bins so far, sweep times): `srun --jobid=<id> --overlap -N1 -n1 bash jobs/nt/peek.sh`.
L=32 t4 submitted separately (job 566820 would hit its time limit before t4).

Update 11:55: proposal version published (preliminary) to git@git.physik.uni-wuerzburg.de:fassaad/juelich.git, folder
`SuperMUC_Phase2_Conversion_2026/` (commit ee2ab33): main.tex = report/antrag.tex, antrag_data.tex, datastate.tex,
figures/node_performance.{tex,pdf,dat}. Re-copy and push after the last CPU runs.

## beta = L series (2026-11-27, owner: "as fast as possible in wall-clock time")
Same test at production temperature beta = L, only the configurations that were best at beta = 4 (CPU t1, t2, t4 backup at
L=36; ALF GPU w2t7, w4t3; Intel-opt. w4t1, w4t3), NBIN 6, CPU L=32/36 NBIN 3-4. `submit_betaL.sh` (jobs 566844-566848):
three 4-node test jobs (`multi2.sh`, one run per node) + CPU t1 at L=32/36 on general. Runs: `runs_betaL/`.
Reports: `SET=betaL ./prelim.sh` -> `report_betaL/report_prelim.pdf` (internal) and `report_betaL/antrag_prelim.pdf`
(proposal version); the beta = 4 reports in `report/` are unchanged.
Update 12:00: full report also in the juelich repo, SuperMUC_Phase2_Conversion_2026/full_report/ (commit d2abe10).

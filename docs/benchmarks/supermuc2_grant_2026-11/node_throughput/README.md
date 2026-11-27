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

# Data campaign 2026-11-28 (sim-campaign)

Code: `fork/Code-modularity-GPU` at **34a31f2b** (no `Prog/` change; harness scripts and inputs only, this folder). Builds from `git archive` of that sha:
H100 `GPU=nvidia` (nvhpc 25.7, Helma), MI210 `GPU=amd` (`~/p2_amd/camp`, gfx90a), melmac `GPU=host`, `GPU=none` (gfortran 14.2), MPI `GPU=none`, ifx.
Raw outputs stay on the machines: Helma `~/campaign/{tol_work,bench_work,bench_cpu_work,prof_work}`, amd `~/p2_amd/campaign/`, melmac `~/work/ALF_P2/campaign/`.
Jobs: Helma 1476035 (tol), 1476036 (bench + CPU bench), 1476069 (profile), all H100; total well below the 9 GPU-h cap. Scripts: `campaign.py` (tol, bench),
`profile.py`, `ed_stats.py`, `diff_runs.py`, `tol_rule.py`, `bench_table.py`, `*.sbatch`, `amd_run.sh`, `amd_profile.sh`.

## 1. tol-scaling (`tol/`) → P2.4 step 3
`tol_<machine>.tsv`: per model (hub-mz, hub-su2, kondo, tV, sp), scheme (std, dly16, sub16), L, observable class (scal = `*_scal`, eq = `*_eq`, tau = `*_tau`):
max scaled difference |x-y|/max(1,|y|) GPU path vs CPU path of the same binary (`ALF_GPU_DISABLE=1`), 2 sweeps, same seed (confout identical in every cell), the CPU null
(Nwrap 10 vs 5), the CPU scheme null (scheme vs std; 0 to round-off), Precision Green means. `*_files.tsv` per file. 396 cells: H100 L = 8..32 (Hubbard), up to 16 (others),
MI210 Hubbard L = 8..32 (`tol_mi210-refused.tsv` = Kondo/tV/Spin_Peierls L = 8 on MI210: ACTIVE on this tip, the file name is historic), melmac `GPU=host` L <= 16.
Evaluation of rules: `tol/rule_evaluation.txt`.

- No clean size law: the largest GPU-CPU differences are 1e-10 .. 1e-9 (tau, Mz-Hubbard) and 1e-12 .. 1e-11 (SU(2)) at every L; the spread between models is larger than between L.
  The CPU null itself spans 1e-13 (Kondo eq) to 5e-9 (sp L16 eq, melmac). The schemes do not change the difference (dly/sub vs std on the CPU: <= 1e-11 mostly 0).
- The current rule `max(1e-11, 100 x PG mean)` is exceeded in 65 of 396 cells (worst 20x): it is too tight for tau files (small PG mean, e.g. Kondo) and for large Kondo/SP.
- **Proposal** (not applied): `tol = max(1e-11, 100 x PG_mean(CPU), 30 x null)` with `null` = max scaled difference of a CPU-only pair of the same case (Nwrap 10 vs 5), i.e. the rule
  measures the round-off floor of the case itself. 12 of 396 cells remain above (worst 2.6x): H100 hub-su2 L24 scal (7e-11, null 1.3e-11 ... ), MI210 hub-mz L32 eq (5.4e-9, null 1e-10):
  GPU differences 4-50x above the CPU floor at the largest sizes — candidates for G9 follow-up (pre-pivoted QR), not covered by any size-scaled constant.
  Kondo eq/scal: GPU-CPU 1e-11 .. 5e-11 against a null of 1e-14: a null-relative rule alone (c x null) fails there by up to 670x (the QR pre-pivot, G9), so the `1e-11` floor must stay.

## 2. scheme-bench (`bench/`) → P2.6.2
s per sweep = (t(6 sweeps) - t(2 sweeps)) / 4, Ltau = 0, N_delay 16. `fastest_scheme.md`. (`auto` does not exist on this tip: refused.)
- H100 Hubbard: sub16 fastest (L16: 0.33 s vs 1.60 std vs 0.90 dly16); L24 sub16 cell unreliable (t2 > t6, startup noise), dly16 2.0 s vs std 3.6 s.
  Kondo, tV, Spin_Peierls: std ≈ sub16, dly16 slower by 30-60 %; i.e. the device submatrix path only pays off for Hubbard.
- MI210 (Hubbard L = 24): sub16 1.2 s/sweep, dly16 21 s, std 113 s (a factor 94). Kondo/tV/sp L = 8: std fastest.
- CPU (Helma node core, melmac core, single thread): differences between schemes 0-15 %, no consistent winner (sub16/dly16 slightly ahead for Hubbard L = 24: 3-5 %).
  Timings below 0.2 s/sweep are noise-limited (2-sweep resolution, shared nodes).
- Caveat: one Helma node (other jobs on it), melmac loaded (load ~90/96); numbers are ranking-grade.

## 3. ed-stats (`ed/`) → P2.8.7
The repo's ED references exist only for the 4-site Hubbard chain (`testsuite/test_vs_ed/test_specs.yaml`); Kondo, tV, Spin_Peierls have **none** (not covered, no ED written here).
24 independent seeds x 4 dtau per case (finite-T Mz, finite-T SU(2), T = 0 projector), NSweep 1000, NBin 300/200, fit y0 + a dtau^2 as in `analysis.py`. `ed_stats.tsv` (4 dtau), `ed_stats_3dt.tsv` (without dtau = 0.2).

| case | hits 1σ | 2σ | 3σ | mean pull | pull std |
|---|---|---|---|---|---|
| finT-mz (4 dtau) | 12/24 | 18/24 | 24/24 | -0.89 | 1.18 |
| finT-su2 | 12/24 | 21/24 | 24/24 | 0.84 | 1.13 |
| T0-mz | 15/24 | 22/24 | 23/24 | -0.10 | 1.28 |
| finT-mz (3 dtau) | 13/24 | 23/24 | 24/24 | -0.19 | 1.30 |
| finT-su2 (3 dtau) | 17/24 | 22/24 | 24/24 | 0.52 | 0.94 |
| T0-mz (3 dtau) | 14/24 | 22/24 | 24/24 | 0.32 | 1.07 |

Expected 68 / 95 / 99.7 %: 1σ hits 50-71 %, 2σ 75-96 %, 3σ 96-100 %. The -0.89 mean pull of finT-mz with 4 dtau is a fit-model bias (dtau^4 term; gone with 3 dtau), not the code.
Pull std 1.1-1.3 (error bars ~10-30 % small, likely autocorrelation inside bins of 1000 sweeps): not a finding by itself, 24 seeds.

## 4. diff-runs (`diff/`) → P2.8.5
`diff_runs.tsv` (Start seeds, 4x4, NSweep 10, NBin 4; hub-mz, hub-su2, hub-dly4, kondo, tV, sp), melmac CPU.
- MPI np = 1: every file byte-identical to the serial build, all 6 cases. np = 2, 3: rank r gets seed line (np - r) of `seeds` (rank 0 the last), so each rank's `confout_r` is byte-identical to a serial run of
  that seed (all 6 cases, all ranks). The observable files of np >= 2 hold the **average over the ranks** (checked: bin 1 of hub-mz np=2 = mean of the serial seed-1 and seed-2 bins),
  so they differ from a single chain by design ("ROUNDOFF max diff ~1-3" rows are that averaging, not round-off).
- Restart (4 + 4 bins from `confout` -> `confin`) against one run of 8 bins: byte-identical, all 6 cases.
- ifx 2025.3 vs gfortran 14.2: `confout` differs already after 1 sweep. Cause: `Random_Wrap` uses the compiler's `RANDOM_NUMBER`/`RANDOM_SEED` (`Libraries/Modules/random_wrap_mod.F90`),
  so the chains differ by construction; cross-compiler runs are only comparable statistically (not tested here beyond ED above on gfortran). Not a defect; no KNOWN_ISSUES line.
- No real mismatch found, nothing classified A/B/C, `KNOWN_ISSUES.md` untouched.

## 5. profile (`profile/`) → P2.9.5/6
Hubbard Mz, delayed N_delay 32 (L = 24) / 64 (L = 32), single chain, difference of the 4-sweep and the 2-sweep profile (= 2 sweeps after warm-up).

| | wall s/sweep | launches/sweep | kernel s/sweep | GPU busy | host sync calls/sweep | top kernels |
|---|---|---|---|---|---|---|
| H100 L24 | 4.20 | 94 758 | 0.89 | 24 % | 146 454 | dly_down 33 %, dly_up 21 %, dly_accept2 21 %, QR (geqr2) 10 % |
| H100 L32 | 9.38 | 162 616 | 2.52 | 29 % | 255 439 | dly_down 35 %, dly_accept2 26 %, dly_up 19 %, QR 7 % |
| MI210 L24 | 26.1 | 93 852 | 18.5 | 71 % | 966 | dly_accept2 34 %, dly_up 32 %, dly_down 32 % |
| MI210 L32 | 59.6 | 160 100 | 48.9 | 82 % | 872 | dly_accept2 42 %, dly_up 28 %, dly_down 27 % |

H100: the delayed-update kernels (3 per vertex) are 75 % of kernel time and about 3.3 launches and 5 host syncs per vertex (L24: 576 sites x 50 slices = 28 800 vertex updates per sweep); the GPU is idle 70-75 % (launch/sync bound). The MI210 uses far fewer sync calls
(HIP) but its three delayed kernels take 97 % of the time (OpenMP-offload kernels, `gpu_lowrank_mod`): the same code that is 94x slower than submatrix (section 2). Not comparable one to one with E2 of `CG_ALF_PLAN.md`
(that used other inputs and the earlier build; E2's QR share of 39 % does not appear here with delayed N_delay 32: QR is 7-10 %). Standard and submatrix profiles were not taken.

# R_chi^AFM crossing plot on the GPU (pi-flux Hubbard, submatrix)

`rchi_pipeline.py` is one resumable script (needs numpy/scipy/matplotlib: `~/alf-gpu/env/pyplot` on `nividia`) that runs on the
GPU in the background and produces the plot of R_chi^AFM(U, L) with the finite-size crossing, Uc(L), 1/nu(L), eta_phi(L) and
their extrapolations a + b/L^c, and the equilibration plot that fixes Nskip.

```bash
cd ~/alf-gpu/work && export OMP_TARGET_OFFLOAD=MANDATORY
nohup setsid ~/alf-gpu/env/pyplot/bin/python rchi_pipeline.py run --work ~/alf-gpu/work/rchi --deadline-hours 72 > rchi_run.out 2>&1 < /dev/null &
```

Restart = run the same command again (finished runs are skipped, runs continue from `confout_0`).
`analyze` redoes only the analysis and plots from the data on disk. `--quick` runs a few-minute test with L = 4..10.

## Model and observables
Hubbard, square lattice L x L, pi flux (`N_Phi = L^2/2`), `Mz=.F.` (SU(2) symmetric charge decoupling), `N_FL=1`, `N_SUN=2`, `Beta=L`, `Dtau=0.1`,
`Nwrap=10`, `Symm`, `Checkerboard`, `Sequential_scheme="submatrix"`, `N_delay` 32 (L <= 28) / 64, `Ltau=0`, `ham_h0=0` (no defect). L = 12, 16, ..., 36.
S(k) = `SpinZ_eq` (k space, per bin, no 1/N: mean_k S = <(n_up-n_dn)^2>). R = 1 - S(Q+dq)/S(Q), Q = (pi,pi), dq = (+-2pi/L,0),(0,+-2pi/L) averaged.
Uc(L): crossing of R(L,U) and R(L+4,U) (cubic spline in U). 1/nu(L) = ln(R'_{L+4}/R'_L)/ln((L+4)/L) at the crossing, eta_phi(L) = 1 - ln(S_{L+4}(Q)/S_L(Q))/ln((L+4)/L)
(S(Q) ~ L^(1-eta) with C(r) ~ r^-(1+eta) in 2+1 dimensions). x axis 1/L = 1/L_small. Extrapolation a + b/L^c (3 parameters, c in [0.2,6]; c = 1 with < 4 points).
Errors: delete-one-block jackknife (J <= 20 blocks of consecutive bins) of the complete chain bins -> splines -> crossings -> fits.

## Stages and files (in the work directory)
1. `warm/` equilibration study, 1 sweep per bin, (L,U) = (12: 5.0/5.4/6.0, 400 sweeps), (24: same, 300), (36: 5.4/6.0, 150) -> `warmup.png`, `nskip.json`
   (first sweep after which the 10-sweep moving average of S(Q) and E stays within 3 sigma of the second-half mean; x1.5 + 5, min 20; linear in L).
2. `prod/L<L>_U<U>/` production, bins of 5 sweeps; round 1: Nskip + 10 bins, every further round +20 bins; U order 5.4, 5.3, 5.5, ... so that the
   plot is useful early; after every U point (all L) `Rchi_crossing.png/.pdf`, `summary.md`, `summary.json` are refreshed.
3. `logs/pipeline.log` (one line per run), `rchi_run.out`, `DONE` when finished (deadline or `RCHI_MAX_ROUNDS`, default 12).
`Green_eq`/`Den_eq` of every segment are gzipped (`*.segNNNNN.gz`); the pipeline stops if less than 5.5 GB are free.

## Caveats
* L = 4, 6, 8, 10 (quick mode) are not meant to cross: L = 2 mod 4 lattices have a different pi-flux geometry. The production sizes are multiples of 4.
* Statistics grow with the rounds; the first plot (round 1, 10 bins x 5 sweeps) has large error bars, the exponents need several rounds.
* The binary is `ALF-gpu-prod.out` built from branch `gpu-delayed-submatrix-defect-production` (no timing instrumentation).

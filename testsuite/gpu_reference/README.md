# gpu_reference: CPU reference runs and statistical comparison

Provisional Phase 0 tools (see `docs/gpu/BASELINE.md`) plus the Phase 1 comparison (ticket 106).

| File | Purpose |
|---|---|
| `configs.py` | Hubbard reference configs (6x6, 12x12, 20x20; U=4, beta=4) |
| `configs_stat.py` | Short-autocorrelation validation set for `compare2.py` (8x8, beta=2) |
| `prep.py` | Create run directories (parameters + seeds); `--seed-offset N` rotates the seeds file for an independent MC stream |
| `run.sh` | `run.sh ALF.out NRANKS RUNDIR...` run with MPI, record wall time |
| `summarize.py` | pyALF analysis, writes `reference.json` (used by the old `compare.py`) |
| `compare.py` | Old full-array chi2 comparison. Superseded, informational only |
| `compare2.py` | **Observable-level gate**: `compare2.py RUNDIR_A RUNDIR_B [--alpha 0.01] [--json out.json]` |

Needs an HDF5 build (`data.h5` per run directory) and pyALF.

## compare2.py: method

**Curated observables** (scalar per bin, nearly independent): E, Kin, Pot, Part; S_zz(k) at (pi,pi)
and (pi/2,pi/2); Den(k) at (pi,pi); G(k) at (0,0), (pi/2,pi/2), (pi,0); tau-resolved
G(k=(pi,0), beta/2) and S_zz((pi,pi), beta/2); chi_zz(pi,pi) = int_0^beta S_zz((pi,pi),tau) dtau
(trapezoid on the tau grid; needs `Ltau=1`). k-points are located with pyALF's `Lattice.k_to_n` on the
lattice stored in `data.h5` (no assumed index order); points not on the finite lattice (pi/2 on 6x6)
are skipped and printed. At half filling `Part` and `G(k=(pi/2,pi/2))` are constant up to round-off
(particle-hole symmetry): they are compared exactly (relative 1e-9) and not counted in the
test family; with other fillings they become ordinary tested observables. Known weakness of the
list: E, Kin and Pot are not independent (E = Kin + Pot), and S_zz(pi,pi), its tau=beta/2
value and chi_zz are strongly correlated; Sidak over m tests is therefore conservative, not wrong.

**Errors.** Per-bin values come from `py_alf.ana.read_latt/read_scal(bare_bins=True)` (k-space
correlators, background already subtracted, orbital trace taken). The first `n_skip` bins are
dropped, the others are rebinned by factors f = 1, 2, 4, ... while at least 8 bins remain
(`--min-bins`). Each level gives the jackknife error of the sign-aware ratio <O s>/<s>. The plateau
factor f* is the smallest f from which every further doubling raises the error by less than
`1.5/sqrt(2(N-1))` (about 1.5 noise standard deviations of an error estimated from N bins). An
observable is flagged UNDER-RESOLVED if the error is still growing at the largest level, **or** if
fewer than three levels exist (< 32 bins after skipping: a plateau cannot be demonstrated). In that
case sigma is the error at the largest level, which is a lower bound for the true error.
tau_int = 0.5 (sigma_f*/sigma_1)^2 is reported in bins (0.5 = uncorrelated).

**Test.** t = (a - b)/sqrt(sa^2 + sb^2), Welch-Satterthwaite degrees of freedom from N_a(f*)-1 and
N_b(f*)-1, two-sided Student-t p-value. With m tested observables and family-wise false-alarm rate
alpha (default 0.01) the per-test threshold is the Sidak value p < 1 - (1-alpha)^(1/m)
(about 0.01/m). sum t^2 against chi2_m is printed as information. The old per-array chi2/n of all
k and (k,tau) entries is printed as information and never gates.

**Exit code.** 0 PASS: no observable differs *and* every error is resolved (or all means agree to 1e-9
relative, e.g. same-seed runs with identical trajectories). 1 FAIL: a resolved observable differs, or a
constant observable differs. 2 INCONCLUSIVE: nothing resolved differs, but at least one observable is
UNDER-RESOLVED, so its error cannot be trusted (take more or longer bins). Treat any nonzero as "not green".
(Manager decision in the 106 review: a PASS requires trustworthy errors.)

## Validation set (ticket 106)

`configs_stat.py`: 8x8 square Hubbard, half filling, beta=2, dtau=0.1, Ltau=1, Mz, 4 MPI ranks,
250 bins of 40 sweeps per rank. Seeds differ between compared runs (`--seed-offset`), so even
GNU vs ifx is a comparison of independent streams. Results: see ticket report / `docs/gpu/BASELINE.md`.

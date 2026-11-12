# Max_SAC bootstrap test

`2-maxsac-bootstrap.py` is an end-to-end test of the Gaussian (parametric)
bootstrap in `Analysis/Max_SAC.F90`. It writes synthetic `g_dat` data for a known
spectrum (two Gaussian peaks, fermionic `P` channel, β = 10, 41 τ points), runs
`Max_SAC.out` in temporary directories and checks the output.

## Requirements

* A built `Analysis/Max_SAC.out`.
* `python3` (standard library only, no numpy).

## Build

From the ALF root, in **bash** (`configure.sh` does not work when sourced from zsh).
Use the same configuration the libraries were built with, e.g.

```bash
bash -c 'source configure.sh GNU HDF5 && make lib && make ana'
```

(`source configure.sh GNU noMPI` for a serial build; on macOS with Homebrew
OpenBLAS also `source .envrc` first.)

## Run

```bash
python3 testsuite/Analysis.tests/2-maxsac-bootstrap.py
```

Optional arguments:

* a path to another `Max_SAC.out`, e.g.
  `python3 testsuite/Analysis.tests/2-maxsac-bootstrap.py /path/to/Max_SAC.out`
* `--keep`: keep the temporary run directories (they are always kept on failure;
  the path is printed).

It takes about 30 s. The exit code is 0 if all checks pass, 1 otherwise. With
CTest (testsuite build), it runs as `2-maxsac-bootstrap`:

```bash
ctest -R 2-maxsac-bootstrap --output-on-failure
```

## What is checked

| Check | Catches |
|---|---|
| `N_boot = 1`: exactly one MaxEnt run, no `Aom_ps_boot_*` | backward compatibility, data perturbed without bootstrap |
| `N_boot = 8`: 9 runs (8 samples + 1 reference run on the original data) | missing reference run |
| every run starts with `No dump data` | bootstrap samples mixed through dump files |
| the Monte Carlo energy traces of all runs differ pairwise | identical samples (the old fixed-seed reset in `MaxEnt_stoch`) |
| χ² of the samples exceeds the reference χ² by 0.4–1.6 per τ point | noise missing or with the wrong amplitude |
| `Aom_ps_boot_1..N_alpha-10` exist, `Ndis × 3`, finite; error > 0 where A has weight | broken output |
| reference `Aom_ps_*` agrees with the bootstrap mean (median \|ΔA\|/err < 3) | biased mean, outputs from a noisy sample |
| bootstrap error > Monte Carlo error of a single run | error dominated by MC noise only, i.e. no data noise propagated |

Note: the error per frequency bin is **not** expected to scale linearly with the
error bars of the input. Stochastic MaxEnt works in units of σ and smooths more
when the data are noisier, so such a check would be wrong.

## Using the bootstrap

In `&VAR_Max_Stoch` of `parameters`:

* `N_boot` — number of bootstrap samples (default 1 = no bootstrap, original
  behaviour). Only used with `Stochastic = .true.`.
* `N_boot_checkpoint` — if `.true.`, each sample starts from the configurations
  of the previous one (histograms and warmup are reset).

Output: `Aom_ps_boot_<n>` with columns ω, mean A(ω) over the samples, and the
bootstrap standard deviation. All other files (`Aom_ps_*`, `energies`,
`Best_fit`, `Data_out`, `Green`) come from the reference run on the original data.

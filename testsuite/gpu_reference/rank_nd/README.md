# rank_nd: independent ED references for tV, Kondo, Spin_Peierls (ticket 401, gate O1_ed)

Pure numpy/scipy exact diagonalisation of the operators ALF actually samples, validated against the CPU build
(GNU, HDF5, noMPI, `GPU=none`, tree synced from this repo; binary `~/alf-gpu/work/t401ed/bin/ALF_cpu.out`).

| file | content |
|---|---|
| `edcommon.py` | Fock space (Jordan-Wigner), block-diagonal thermal expectation values in (N_up, N_dn) sectors, ALF equal-time lattice observables, ALF `Lattice%list` order (via pyALF) |
| `ed_check_hubbard.py` | convention check: 4-site Hubbard ring, U=4, beta=2 gives -1.47261997 (upstream value) |
| `ed_tV.py`, `ed_Kondo.py`, `ed_SpinPeierls.py` | the models; the docstrings carry the derivation from the Fortran |
| `ed_cases.py` | validated cases: pyALF `sim_dict` + ED with pyALF output names |
| `ed_validate.py` | prep / run / analyse (pyALF jackknife) / compare, dtau scan and weighted `y0 + a dtau^2` extrapolation (as `testsuite/test_vs_ed`) |
| `ed_negctl.py` | planted-defect controls (ED with one ingredient perturbed vs. the same MC data) |
| `gen_models.py` -> `models.py` | ready parameter sets for the GPU harness, ED numbers, validation record |
| `validation/summary_*.json` | raw comparison data (MC values and errors per dtau, ED values) |
| `ed_build_cpu.sh`, `ed_probe.py` | CPU build script, single-run probe |

Run on nividia in `~/alf-gpu/work/t401ed/` (venv `~/alf-gpu/venv`, `OMP_NUM_THREADS=1`).
Example: `python ed_validate.py tV_chain4 --dtaus 0.2,0.1,0.05,0.025 --nsweep 400 --nbin 200 --tag val_tV`.

## Findings that matter for anyone using these numbers

1. **Constants.** ALF's Hubbard-type vertex is `(U/N_SUN)(sum_s (n_s-1/2))^2 = U(n_up-1/2)(n_dn-1/2) + U/4`, so ALF energies contain `+U/4` per site (this is why -1.4726 needs it).
2. **tV, duplicated bonds (L1=2 or L2=2).** `Op_V` exists once per (cell, direction), so a duplicated pair carries two interaction vertices.
   With `Checkerboard=.true.` it also hops with 2t and `Kin = <H_T>`. With `Checkerboard=.false.` the dense `Op_T` is filled with `=`
   (`Predefined_Hoppings_set_OPT`), so the pair hops with t in the propagation while `Compute_Kin` counts it twice: `Kin = 2<H_T,prop>`,
   `Ener` is not `<H>`. Found with this ED (without the mirror, 27 sigma on Ener); upstream behaviour, not changed. Same for Kondo.
3. **Kondo.** `Pot` = Uc term + JK term + constant; the **Uf term is not in Pot/Ener** (separate observable `Constraint = (1/N_SUN) sum <(n_f-1)^2>`).
   `Obser`'s JK estimator is the Wick form of `A^2` **without the pair terms** `sum_{s!=s'}(c+_s f_s c+_s' f_s' + h.c.)`, so `Pot` is not
   `<Uc term> + <H_JK>` (difference 0.04-0.06 on 2x1, 12-17 sigma; with JK=0 no difference). The constant is
   `JK * ((N_cells*Norb)/8)` in INTEGER arithmetic: 0 for 2x1, JK for 2x2 (true value would be N_cells*JK/4). `T_Loc` of the f orbital is 0 (chem acts on c only).
4. **Spin_Peierls.** ALF writes no energy. `Part` is `2<N>` (`ZRho` adds `GRC` twice), `Pot` is the double occupancy `sum <n_up n_dn>`.
   Hamiltonian: no fermion hopping; `H = sum_I (U/2)(n_I-1)^2 + sum_b [P^2/2M + kQ^2/2 + (J/4)Q - (J/4)(1+Q) K_b^2]`, `k = 1/(2 Lambda)`, `M = k/Omega0^2`,
   `K_b = sum_s (c+_i c_j + h.c.)`. Smallest allowed lattice: `Square L1=2 L2=1` (2 sites, 2 bonds = 2 phonons); 2x2 would need 8 oscillators (not done).
5. `Den_eqR`, `Phi_eqR`, `Dimer_eqR` are connected (pyALF subtracts the `Latt0` products), `Green_eqR`, `SpinZ_eqR` are not. `SpinZ = <(n_up-n_dn)_i (n_up-n_dn)_j>`.

## Validation protocol and result

dtau in {0.2, 0.1, 0.05, 0.025} (SP: 0.1, 0.05, 0.025), 200 bins x 400 sweeps (SP 1000), first 5 bins skipped, fresh random seeds, pyALF jackknife errors.
Comparison: weighted fit `y0 + a dtau^2`, `z = (y0 - ED)/err(y0)`; per-dtau z shown in `validation/*.json` (Trotter error is visible and shrinks as dtau^2).
Nothing was tuned to force agreement; the two ED corrections above (duplicated bonds, Kondo estimator) were derived from the Fortran after the CPU run
disagreed, and the corrected ED then agreed (Kondo 2x1: the same MC data re-analysed with the corrected ED; tV 2x2: fresh runs with V=1, Nwrap 2 for both checkerboard modes; the first V=2 Nwrap-4 runs are only in the nividia run dir `runs/val_tV`).
Numbers: see `models.py` (`validation_cpu_vs_ed`) and the report of ticket 401.

Negative controls (`ed_negctl.py`): V, t, JK, Uc, Uf by 2-10 % and the duplicated-bond mode are rejected with |z| of 10-400; for Spin_Peierls J 10 %, lambda 30 %,
Omega0 30 %, U 50 %, beta 20 % are rejected, 2-5 % perturbations are not resolvable at this statistics. beta by 2 % is not resolved for tV energies.

## Limits

* Lag-1 autocorrelation of the phonon bins is up to 0.7 at dtau 0.025 (errors there underestimated); use dtau 0.1 for SP GPU tests.
* tV with V=2 on 2x2 is unstable at Nwrap 4 / dtau 0.2 (G difference > 10); V=1 and Nwrap 2 were used.
* ED of Kondo 2x2 takes about 2.5 min (many correlator operators).
* `Kondo_2x1_dense` and Kondo/SP beyond the listed sizes were not run.

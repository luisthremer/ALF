# Parameter sets shared by regress_vs_baseline.sh and compare_trajectory.sh (sourced, not executed).
# Each entry is "name|common keys|branch-only keys"; the branch-only keys are keys that the baseline
# binary cannot parse (or a deliberate difference) and are only used where a branch binary is run.
SETS=(
  "start||"
  # explicit: parser check for the branch-only keys (baseline cannot parse them): key name in another case,
  # single-quoted upper-case value, non-default N_delay; N_delay must be irrelevant for "standard".
  "explicit||sequential_scheme='STANDARD' N_delay=7"
  "doped|ham_chem=-0.5 L1=4 L2=4 Beta=4|"
  "density|Mz=.F.|"
  "continuous|Continuous=.T.|"
  # tV: the Start default ham_V=4 aborts at baseline (unstable, positive V gives an imaginary
  # coupling), so use the parameters of testsuite/test_branch_parameters.json (ham_V=-1, N_SUN=1).
  "tv|ham_name=\"tV\" Model=\"tV\" ham_V=-1.d0 N_SUN=1 L1=4 L2=4 Beta=4|"
  # kondo: mixes diagonal rank-1 (U) and non-diagonal rank-2 (V_SUN/VJ_SUN, J) vertices in one time slice: all eligible
  # for delayed (support <= 8); the J vertices are ineligible for the rank-1 submatrix (flush + Upgrade2).
  "kondo|ham_name=\"Kondo\" Model=\"Kondo\" Lattice_type=\"Bilayer_square\" L1=4 L2=4 Beta=2.d0|"
  # pv_adiab: projective algorithm with adiabatic switching (time-dependent couplings g_t, N_FL_eff < N_FL
  # via flavour reconstruction). Beta/Theta/Projector/Adiabatic live in VAR_Model_Generic and VAR_Hubbard_Plain_Vanilla.
  "pv_adiab|ham_name=\"Hubbard_Plain_Vanilla\" Model=\"Hubbard_Plain_Vanilla\" L1=4 L2=4 Projector=.T. Adiabatic=.T. Theta=2.d0 Beta=1.d0|"
  # z2gt: Z2 matter with Global_tau moves (global moves while low-rank updates are pending). Run the delayed side
  # with N_delay=5 so that flushes are exercised.
  "z2gt|ham_name=\"Z2_Matter\" Model=\"Z2_Matter\" L1=4 L2=4 Beta=1.d0 ham_U=1.0 Global_tau_moves=.T. N_Global_tau=4 Symm=.F. Checkerboard=.F.|"
  # Rank-N sets (multi-site vertices; the validation targets of the rank-N delayed/submatrix work).
  # tv_flux: tV bonds (N_P = 2, non-diagonal) with a complex hopping matrix (one flux quantum).
  "tv_flux|ham_name=\"tV\" Model=\"tV\" ham_V=-1.d0 N_SUN=1 L1=4 L2=4 Beta=4 N_Phi=1|"
  # tv_proj: tV bonds in the projective algorithm (Theta=2, trial wave function).
  "tv_proj|ham_name=\"tV\" Model=\"tV\" ham_V=-1.d0 N_SUN=1 L1=4 L2=4 Beta=1.d0 Projector=.T. Theta=2.d0|"
  # z2: Z2_Matter with all vertices updated sequentially (no Global_tau moves: then Nt_sequential_start/end default
  # to 1 .. n_op), so the Ising_SUN gauge bonds and matter bonds (type 1, N_P = 2) are proposed.
  "z2|ham_name=\"Z2_Matter\" Model=\"Z2_Matter\" L1=4 L2=4 Beta=1.d0 ham_U=1.0 Symm=.F. Checkerboard=.F.|"
  # spin_peierls: type-4 (continuous, flip_protocol 2) N_P = 2 bond vertices with off-diagonal O.
  "spin_peierls|ham_name=\"Spin_Peierls\" Model=\"Spin_Peierls\" L1=4 L2=4 Beta=2.d0 N_SUN=2 N_FL=1 ham_Jx=0.4d0 ham_Jy=0.4d0 ham_U=1.5d0 ham_Lambda=0.1d0 ham_omega0=0.5d0|"
)

# Trajectory-only sets "name|common keys|": the baseline binary does not know the model, so these are not part of
# the bit-identity regression (regress_vs_baseline.sh), only of compare_trajectory.sh (branch/trace binary).
# jz: LR_Jz_Test (test-only model, Prog/Hamiltonians/Hamiltonian_LR_Jz_Test_smod.F90): Mz-Hubbard sites plus
# one Predefined_Int_Jz bond vertex (N_P = 2) per nearest-neighbour bond.
TRAJ_ONLY_SETS=(
  "jz|ham_name=\"LR_Jz_Test\" Model=\"LR_Jz_Test\" L1=4 L2=4 ham_U=4.d0 ham_Jz=1.d0 Beta=2.d0 Nwrap=5|"
)

# Statistical sets for stat_runs.sh (independent long chains compared with stat_compare.py): "name|keys".
# Common keys (NSweep=50 NBin=40 Ltau=0 n_skip=5) are added by stat_runs.sh; keys given here override them.
# doped: ham_chem=-1.0 with Beta=5 is the smallest beta in {5,6,7} whose average sign lies in [0.3, 0.8]
# (measured 0.58(2), 0.59(2); at Beta=4 the sign is ~0.92, see CLAUDE.md). NBin=80 (4000 sweeps, ~105 s) keeps the
# sign-reweighted error bars usable.
STAT_SETS=(
  "half|L1=6 L2=6 ham_U=4.d0 Beta=5.d0 ham_chem=0.d0"
  "doped|L1=6 L2=6 ham_U=4.d0 Beta=5.d0 ham_chem=-1.0d0 NBin=80 n_skip=10"
  "density|Mz=.F. L1=6 L2=6 ham_U=4.d0 Beta=4.d0 ham_chem=0.d0"
  # Rank-N statistical sets (4x4, one standard run <= ~3 min; measured 12-55 s): multi-site vertices of tV, Kondo-J, Jz
  # and the projective tV. tv: N_SUN=1 ham_V=-1 is nearly sign-free only away from half filling: ham_chem=-1, Beta=2 gives
  # average sign 0.64(2) (at half filling it is ~0.03). tv_proj: the N_SUN=1 sets have average sign ~0.05 in the projective
  # algorithm, so the statistical set uses N_SUN=2, ham_V=0.5, one flux quantum (complex hopping; sign 1).
  # jz: Nwrap=5 (Nwrap=10 gives Precision Green max 1e-5, Nwrap=5 1.5e-9).
  "tv|ham_name=\"tV\" Model=\"tV\" ham_V=-1.d0 N_SUN=1 L1=4 L2=4 Beta=2.d0 ham_chem=-1.d0 NBin=80"
  "kondo|ham_name=\"Kondo\" Model=\"Kondo\" Lattice_type=\"Bilayer_square\" L1=4 L2=4 Beta=2.d0"
  "jz|ham_name=\"LR_Jz_Test\" Model=\"LR_Jz_Test\" L1=4 L2=4 ham_U=4.d0 ham_Jz=1.d0 Beta=4.d0 Nwrap=5"
  "tv_proj|ham_name=\"tV\" Model=\"tV\" ham_V=0.5d0 N_SUN=2 N_Phi=1 L1=4 L2=4 Beta=1.d0 Projector=.T. Theta=2.d0 NBin=80"
  # z2: Z2_Matter, all vertices sequential (no Global_tau): Hubbard, gauge-Ising and matter-Ising vertices. No Ener/Pot/Kin
  # observables (Part, Flux, X, Q): stat_runs.sh analyses the *_scal files the model has. NSweep=200 (bins of 200
  # sweeps): with NSweep=50 the gauge-field autocorrelation underestimated the errors (standard s1 vs s2: Green_eqJK
  # 4 of 16 entries |z|>2, rms z 1.5, FAIL); see README.
  "z2|ham_name=\"Z2_Matter\" Model=\"Z2_Matter\" L1=4 L2=4 Beta=2.d0 ham_U=1.0 Symm=.F. Checkerboard=.F. NSweep=200"
  # pv_adiab: projective Hubbard with adiabatic switching (time-dependent g_t, N_FL_eff < N_FL), as in SETS.
  "pv_adiab|ham_name=\"Hubbard_Plain_Vanilla\" Model=\"Hubbard_Plain_Vanilla\" L1=4 L2=4 Projector=.T. Adiabatic=.T. Theta=2.d0 Beta=1.d0"
)

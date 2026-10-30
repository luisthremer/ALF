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
  # kondo: mixes eligible (diagonal rank-1) and ineligible (non-diagonal / rank>1) vertices in one time slice.
  "kondo|ham_name=\"Kondo\" Model=\"Kondo\" Lattice_type=\"Bilayer_square\" L1=4 L2=4 Beta=2.d0|"
  # pv_adiab: projective algorithm with adiabatic switching (time-dependent couplings g_t, N_FL_eff < N_FL
  # via flavour reconstruction). Beta/Theta/Projector/Adiabatic live in VAR_Model_Generic and VAR_Hubbard_Plain_Vanilla.
  "pv_adiab|ham_name=\"Hubbard_Plain_Vanilla\" Model=\"Hubbard_Plain_Vanilla\" L1=4 L2=4 Projector=.T. Adiabatic=.T. Theta=2.d0 Beta=1.d0|"
  # z2gt: Z2 matter with Global_tau moves (global moves while low-rank updates are pending). Run the delayed side
  # with N_delay=5 so that flushes are exercised.
  "z2gt|ham_name=\"Z2_Matter\" Model=\"Z2_Matter\" L1=4 L2=4 Beta=1.d0 ham_U=1.0 Global_tau_moves=.T. N_Global_tau=4 Symm=.F. Checkerboard=.F.|"
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
)

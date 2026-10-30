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
)

# Statistical sets for stat_runs.sh (independent long chains compared with stat_compare.py): "name|keys".
# Common keys (NSweep, NBin, Ltau, n_skip) are added by stat_runs.sh. doped: ham_chem=-1.0 was chosen
# in the Task 5 pilot (measured sign 0.926 +- 0.011; the requested window 0.3..0.8 is not reachable at
# beta=4, see CLAUDE.md).
STAT_SETS=(
  "half|L1=6 L2=6 ham_U=4.d0 Beta=5.d0 ham_chem=0.d0"
  "doped|L1=6 L2=6 ham_U=4.d0 Beta=4.d0 ham_chem=-1.0d0"
  "density|Mz=.F. L1=6 L2=6 ham_U=4.d0 Beta=4.d0 ham_chem=0.d0"
)

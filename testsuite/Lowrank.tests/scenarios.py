"""Scenarios of compare_schemes.py beyond testsuite/test_branch_parameters.json.

Written with AI assistance on master_ai; not part of the ALF baseline.

Each entry maps a name to a dictionary:
    pars      parameters set in the template file Scripts_and_Parameters_files/Start/parameters
              (key "Model" sets the model name in both namelists that carry it)
    mode      "serial" (default), "mpi" (one parameter set on 2 ranks; MPI or PARALLEL_PARAMS build),
              "tempering" (replicas, mpi_per ranks each, default 1; TEMPERING or PARALLEL_PARAMS build) or
              "params" (one parameter file per replica; PARALLEL_PARAMS build)
    replicas  list of parameter dictionaries, one per replica (tempering, params)
    restart   number of bins of each of two stages; the second continues from the final configuration
    quick     also run with --quick (the MPI builds of the CI)
    force_in_proposal   the field proposal contains the force (MALA), which is computed from the Green
              function: the final configuration is then only required to agree within the tolerance
"""

# Hubbard on a 6x2 ladder: the base of the feature scenarios
HUB = {"Model": "Hubbard", "Lattice_type": "N_leg_ladder", "L1": 6, "L2": 2, "Beta": 4.0,
       "Ham_U": 4.0, "Ham_T": 1.0, "Nsweep": 20, "NBin": 5, "Ltau": 1, "Mz": True}

# Models: every Hamiltonian of Prog/Hamiltonians.list with its operator kinds, lattices, N_SUN / N_FL
# variants and flavour reconstruction (Hubbard_Plain_Vanilla at half filling computes one flavour only)
_SQ = {"Lattice_type": "Square", "L1": 4, "L2": 4, "Beta": 2.0, "Nsweep": 10, "NBin": 4, "Ltau": 1}
_HSQ = dict(_SQ, Model="Hubbard", Ham_U=4.0, Ham_T=1.0)
_KON = dict(_SQ, Model="Kondo", Lattice_type="Bilayer_square", L1=4, L2=2, Ham_T=1.0, Ham_Uf=1.0, Ham_JK=2.0)
_TV = dict(_SQ, Model="tV", Ham_T=1.0, Ham_V=-1.0, N_SUN=1)
_LRC = dict(_SQ, Model="LRC", Ltau=0, Ham_U=4.0, Ham_alpha=0.1, Percent_change=0.1)
_Z2 = dict(_SQ, Model="Z2_Matter", Beta=4.0, Dtau=0.1, Ltau=0, Ham_T=1.0, Ham_TZ2=1.0, Ham_K=1.0, Ham_J=1.0,
           Ham_h=1.0, Ham_g=1.0, Ham_U=0.0, Nwrap=10)

MODELS = {
    # Hubbard: Ising (Mz, type 1, N_FL = 2), discrete HS (type 2), continuous HS (type 3)
    "Hubbard_Square_Mz": dict(quick=True, pars=dict(_HSQ, Mz=True)),
    "Hubbard_Square_N_FL1": dict(pars=dict(_HSQ, Mz=False, N_FL=1)),
    "Hubbard_Square_N_FL2_N_SUN2": dict(pars=dict(_HSQ, Mz=False, N_FL=2, N_SUN=2)),
    "Hubbard_Square_Mz_N_SUN4": dict(pars=dict(_HSQ, Mz=True, N_SUN=4)),
    "Hubbard_Square_Continuous_Mz": dict(pars=dict(_HSQ, Mz=True, Continuous=True)),
    "Hubbard_Square_Continuous_N_FL1": dict(pars=dict(_HSQ, Mz=False, N_FL=1, Continuous=True)),
    "Hubbard_Square_twist_complex": dict(pars=dict(_HSQ, Mz=True, Phi_X=0.2, Phi_Y=0.1)),
    "Hubbard_Square_flux_complex": dict(pars=dict(_HSQ, Mz=False, N_Phi=1)),
    "Hubbard_Square_no_checkerboard": dict(pars=dict(_HSQ, Mz=True, Checkerboard=False, Symm=False)),
    "Hubbard_Square_pinning": dict(pars=dict(_HSQ, Mz=True, Ham_h0=0.5)),
    "Hubbard_Honeycomb": dict(pars=dict(_HSQ, Lattice_type="Honeycomb", L1=3, L2=3, Mz=True)),
    "Hubbard_Honeycomb_flux": dict(pars=dict(_HSQ, Lattice_type="Honeycomb", L1=3, L2=3, Mz=False, N_Phi=1)),
    "Hubbard_Bilayer_honeycomb": dict(pars=dict(_HSQ, Lattice_type="Bilayer_honeycomb", L1=2, L2=2, Mz=True,
                                                Ham_Tperp=1.0)),
    "Hubbard_Triangular": dict(pars=dict(_HSQ, Lattice_type="Triangular", L1=4, L2=4, Mz=True)),
    # Kondo: single-site U_SUN and two-site V_SUN / VJ_SUN vertices (N_SUN > 2: three vertices per bond)
    "Kondo_N_SUN2": dict(quick=True, pars=dict(_KON, N_SUN=2)),
    "Kondo_N_SUN4_symm": dict(pars=dict(_KON, N_SUN=4, Symm=True)),
    "Kondo_N_SUN4_nosymm": dict(pars=dict(_KON, N_SUN=4, Symm=False, Checkerboard=False)),
    "Kondo_Uc": dict(pars=dict(_KON, Ham_Uc=1.0, Ham_Uf=2.0)),
    # Hubbard_Plain_Vanilla: N_FL = 2, flavour reconstruction at half filling
    "Plain_Vanilla_half_filling": dict(pars=dict(_SQ, Model="Hubbard_Plain_Vanilla", Ham_U=4.0, Ham_T=1.0)),
    "Plain_Vanilla_doped": dict(pars=dict(_SQ, Model="Hubbard_Plain_Vanilla", Ham_U=4.0, Ham_T=1.0, Ham_chem=0.3)),
    "Plain_Vanilla_attractive": dict(pars=dict(_SQ, Model="Hubbard_Plain_Vanilla", Ham_U=-4.0, Ham_T=1.0, Ham_chem=0.3)),
    # tV: two-site type 2 vertices
    "tV_N_SUN1": dict(pars=_TV),
    "tV_N_SUN2": dict(pars=dict(_TV, N_SUN=2)),
    "tV_repulsive": dict(pars=dict(_TV, Ham_V=1.0, Nwrap=2)),
    "tV_bilayer": dict(pars=dict(_TV, Lattice_type="Bilayer_square", L1=4, L2=2, Ham_T2=1.0, Ham_Tperp=0.5,
                                 Ham_V2=-1.0, Ham_Vperp=-0.5)),
    # LRC: one-site vertices, long-range Ising-like action
    "LRC_alpha": dict(quick=True, pars=_LRC),
    "LRC_global_tau": dict(pars=dict(_LRC, Global_tau_moves=True, N_Global_tau=2)),
    "LRC_bilayer": dict(pars=dict(_LRC, Lattice_type="Bilayer_square", L1=4, L2=2, Ham_Tperp=0.5, Ham_T2=1.0)),
    "LRC_projector": dict(pars=dict(_LRC, Projector=True, Theta=2.0)),
    # Z2_Matter: Hubbard (type 2), gauge and matter Ising fields (type 1), plaquette and star terms
    "Z2_Matter_full": dict(quick=True, pars=_Z2),
    "Z2_Matter_with_U": dict(pars=dict(_Z2, Ham_U=2.0)),
    "Z2_Matter_global_tau": dict(pars=dict(_Z2, Global_tau_moves=True, N_Global_tau=2)),
    "Z2_Matter_propose_S0": dict(pars=dict(_Z2, Propose_S0=True)),
    "Z2_Matter_gauge_only": dict(pars=dict(_Z2, Ham_T=0.0, Ham_h=0.0, Ham_J=0.0, Ham_K=0.0)),
    "Z2_Matter_projector": dict(pars=dict(_Z2, Projector=True, Theta=2.0)),
    # Spin_Peierls: type 2 and type 4 (three-body) vertices
    "Spin_Peierls_square": dict(pars=dict(_SQ, Model="Spin_Peierls", Beta=2.0, N_SUN=2, N_FL=1, Ham_Jx=1.0,
                                          Ham_Jy=1.0, Ham_U=1.5, Ham_Lambda=0.1, Ham_omega0=0.5)),
    "Spin_Peierls_chain": dict(pars=dict(_SQ, Model="Spin_Peierls", L1=8, L2=1, N_SUN=2, N_FL=1, Ham_Jx=1.0,
                                         Ham_Jy=1.0, Ham_U=1.5, Ham_Lambda=0.1, Ham_omega0=0.5)),
}

# ALF features: space-time global moves, partial sequential ranges, Propose_S0, Checkerboard / Symm, restart,
# MPI ranks, parallel tempering and one parameter file per replica
FEATURES = {
    "Hubbard_Global_moves": dict(pars=dict(_HSQ, Mz=True, Continuous=True, Global_moves=True, N_Global=2)),
    "Hubbard_Global_moves_and_tau": dict(pars=dict(_HSQ, Mz=True, Continuous=True, Global_moves=True, N_Global=2,
                                                   Global_tau_moves=True, N_Global_tau=2,
                                                   Nt_sequential_start=1, Nt_sequential_end=16)),
    # Nt_sequential_start/end select the vertices n of a time slice that get the sequential update; the
    # remaining ones are left to the Global_tau moves. (A range that starts after vertex 1 makes ALF master itself
    # stop with "unstable": the vertices before the range keep their random initial fields.)
    "Hubbard_partial_range_head": dict(pars=dict(_HSQ, Mz=True, Global_tau_moves=True, N_Global_tau=2,
                                                 Nt_sequential_start=1, Nt_sequential_end=8)),
    "Hubbard_partial_range_start": dict(pars=dict(_HSQ, Mz=False, N_FL=1, Nt_sequential_start=5,
                                                  Nt_sequential_end=-1)),
    "Hubbard_Propose_S0": dict(pars=dict(_HSQ, Mz=True, Propose_S0=True)),
    "LRC_Propose_S0": dict(pars=dict(_LRC, Propose_S0=True)),
    "Hubbard_checkerboard_no_symm": dict(pars=dict(_HSQ, Mz=True, Checkerboard=True, Symm=False)),
    "Hubbard_no_checkerboard_symm": dict(pars=dict(_HSQ, Mz=True, Checkerboard=False, Symm=True)),
    "Hubbard_restart": dict(quick=True, restart=2, pars=dict(_HSQ, Mz=True, NBin=2)),
    "Hubbard_restart_Global_tau": dict(restart=2, pars=dict(_HSQ, Mz=False, N_FL=1, NBin=2, Global_tau_moves=True,
                                                            N_Global_tau=2, Nt_sequential_start=1,
                                                            Nt_sequential_end=16)),
    "Z2_Matter_restart": dict(restart=2, pars=dict(_Z2, NBin=2)),
    # 2 MPI ranks, one parameter set
    "Hubbard_mpi_2ranks": dict(mode="mpi", pars=dict(_HSQ, Mz=True)),
    "Kondo_mpi_2ranks": dict(mode="mpi", pars=dict(_KON, N_SUN=2)),
    "Z2_Matter_mpi_2ranks": dict(mode="mpi", pars=dict(_Z2, Global_tau_moves=True, N_Global_tau=1)),
    # parallel tempering: two replicas with their own interaction (TEMPERING builds)
    "Hubbard_tempering": dict(mode="tempering", pars=dict(_HSQ, Mz=True),
                              replicas=[dict(_HSQ, Mz=True, Ham_U=4.0), dict(_HSQ, Mz=True, Ham_U=3.0)]),
    "Hubbard_tempering_Global_tau": dict(mode="tempering", pars=dict(_HSQ, Mz=False, N_FL=1),
                                         replicas=[dict(_HSQ, Mz=False, N_FL=1, Global_tau_moves=True, N_Global_tau=2,
                                                        Nt_sequential_start=1, Nt_sequential_end=16, Ham_U=4.0),
                                                   dict(_HSQ, Mz=False, N_FL=1, Global_tau_moves=True, N_Global_tau=2,
                                                        Nt_sequential_start=1, Nt_sequential_end=16, Ham_U=2.0)]),
    "Hubbard_tempering_2ranks_per_set": dict(mode="tempering", mpi_per=2, pars=dict(_HSQ, Mz=False, N_FL=1),
                                             replicas=[dict(_HSQ, Mz=False, N_FL=1, Ham_U=4.0),
                                                       dict(_HSQ, Mz=False, N_FL=1, Ham_U=3.0)]),
    "Z2_Matter_tempering": dict(mode="tempering", pars=_Z2,
                                replicas=[_Z2, dict(_Z2, Ham_K=0.5)]),
    # PARALLEL_PARAMS: each replica has its own parameter file, also with a different model, lattice and size
    "Hubbard_parallel_params": dict(mode="params", pars=dict(_HSQ, Mz=True),
                                    replicas=[dict(_HSQ, Mz=True), dict(_HSQ, Mz=False, N_FL=1, L1=4, L2=2, Ham_U=2.0)]),
    "Hubbard_Kondo_parallel_params": dict(mode="params", pars=dict(_HSQ, Mz=True),
                                          replicas=[dict(_HSQ, Mz=True), dict(_KON, N_SUN=2)]),
}

# Test-only models (add_test_models.sh): multi-site Jz vertices, and continuous and discrete fields in one run
_MIX_NML = lambda mz, cont=True: {"VAR_Hubbard_Mixed": {"Mz": mz, "Continuous": cont, "ham_T": 1.0, "ham_U": 4.0}}
_MIX = dict(_SQ, Model="Hubbard_Mixed", Beta=2.0, Checkerboard=True, N_SUN=2)
TEST_MODELS = {
    "LR_Jz_Test": dict(testmodel=True, quick=True, pars=dict(_SQ, Model="LR_Jz_Test", Beta=2.0, Nwrap=2, nml={
        "VAR_LR_Jz_Test": {"ham_T": 1.0, "ham_U": 4.0, "ham_Jz": 1.0, "Dtau": 0.1, "Beta": 2.0}})),
    # 10-site vertices exceed the limit of 8 sites of the low-rank kernel: they run the standard update
    "LR_Wide_Test": dict(testmodel=True, quick=True, pars=dict(_SQ, Model="LR_Wide_Test", L1=12, L2=1, Beta=2.0, Nwrap=2, nml={
        "VAR_LR_Wide_Test": {"ham_T": 1.0, "ham_U": 4.0, "ham_Jz": 0.1, "Dtau": 0.1, "Beta": 2.0}})),
    "Mixed_fields_sequential": dict(testmodel=True, pars=dict(_MIX, Mz=False, N_FL=1, nml=_MIX_NML(False))),
    "Mixed_fields_Mz_sequential": dict(testmodel=True, pars=dict(_MIX, Mz=True, nml=_MIX_NML(True))),
    # sequential MALA on the continuous vertices, the discrete ones use the low-rank update in between
    "Mixed_fields_Sequential_MALA": dict(testmodel=True, quick=True, force_in_proposal=True, pars=dict(
        _MIX, Mz=False, N_FL=1, Sequential_MALA=True, nml=_MIX_NML(False))),
    "Mixed_fields_Mz_Sequential_MALA": dict(testmodel=True, force_in_proposal=True, pars=dict(
        _MIX, Mz=True, Sequential_MALA=True, nml=_MIX_NML(True))),
    "Mixed_fields_MALA_Global_tau_MALA": dict(testmodel=True, force_in_proposal=True, pars=dict(
        _MIX, Mz=False, N_FL=1, Sequential_MALA=True, Global_tau_MALA_moves=True, N_Global_tau_MALA=2,
        Nt_sequential_start=1, Nt_sequential_end=16, nml=_MIX_NML(False))),
}

SCENARIOS = {
    # sequential vertices 1..12 (all sites) followed by the Global_tau moves on every slice
    "Hubbard_Global_tau": dict(quick=True, pars=dict(HUB, Global_tau_moves=True, N_Global_tau=2,
                                                     Nt_sequential_start=1, Nt_sequential_end=12)),
    "Hubbard_Global_tau_MALA_seq": dict(force_in_proposal=True, pars=dict(
        HUB, Continuous=True, Global_tau_MALA_moves=True, Nt_sequential_start=1, Nt_sequential_end=12)),
    "Hubbard_Continuous": dict(pars=dict(HUB, Continuous=True)),
    "Hubbard_Sequential_MALA": dict(force_in_proposal=True, pars=dict(HUB, Continuous=True, Sequential_MALA=True)),
    "Hubbard_Global_tau_MALA": dict(force_in_proposal=True, pars=dict(
        HUB, Continuous=True, Global_tau_MALA_moves=True)),
    "Hubbard_Projector": dict(pars=dict(HUB, Projector=True, Theta=4.0)),
}
SCENARIOS.update(MODELS)
SCENARIOS.update(FEATURES)
SCENARIOS.update(TEST_MODELS)

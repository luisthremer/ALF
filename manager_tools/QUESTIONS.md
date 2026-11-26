# Open questions for the owner (2026-11-26)

Short form in `STATUS.md` ("Waiting for you"). Here each question with its background, options and Claude's
recommendation. Answer in the chat with the number, e.g. "Q1 yes".

## A. Block the next step

**Q1. P2.4 decisions D1-D4 (merge of step 1a waits for this).** Plan: `docs/ai/phase2/P2_4_GPU_PLAN.md`, section 4.
Step 1a (`cm/feat/gpu-build`: `GPU=none|host` switch, `alf_la_mod` host backend, 12 tests) is green in CI and
byte-identical in all 16 scenarios.
- D1: the GPU linear-algebra layer `alf_la` sits *beside* Fakher's `linalg_mod`, not inside it (ROADMAP said
  "inside"). Inside would need device pointers, module-owned workspaces and asynchronous error checks in every host
  caller and would still not cover the GPU code that is not BLAS (pivoting, wraps, updates, measurements).
- D2: a `GPU=host` test build that runs the whole GPU code path on the CPU (probe: works; configurations
  byte-identical to the CPU path; slow, so small sizes only). Lets us test steps 1b-1e on the Mac and in CI.
- D3: step 1 covers Hubbard only (the GPU observables exist only for Hubbard); other models keep running on the CPU
  path, as on gpu-port.
- D4: the changes gpu-port made to the CPU path without a guard (`Op_V_is_symmetric` rewrite, `ALF_UNSTABLE_THRESHOLD`,
  `alf_exit`, MaxEnt `N_boot`, `QR_PREPIVOT`, timers, `Lowrank_cross_slice`) come later as separate CPU-line features,
  so `GPU=none` stays byte-identical.
Recommendation: yes to all four. Then: merge 1a, start 1b (GPU state on the walker).

## B. Upstream (ALF-QMC/ALF)

**Q2. Post the two bug reports?** You said yes in principle; the texts are drafted in the manager worktree,
`.handoff/upstream_issues_H6_H7.md`. Each starts with "written by an AI assistant during an automated bug hunt;
Luis Thremer has not checked it". H7: Bilayer_honeycomb potential/total energy uses weights (U, U2, 0, 0) instead of
(U, U, U2, U2) (master 4c690b04, `Hamiltonian_Hubbard_smod.F90` lines 762-763; Markov chain unaffected). H6: Tempering
with an odd number of parameter sets deadlocks, one set crashes in MPI_Finalize (`main.F90` 246, `Global_mod.F90`
164/191). Answer "post" (as is) or corrections.

## C. GPU line

**Q3. Speedup S1 into the GPU line?** The GPU determinant phase from tracked U determinants instead of an LU
(`rd/cgalf-s1`, origin `speedup/nvidia`). Measured on the H100: -12.6 % per sweep at L = 32, configurations
byte-identical. Looks vendor-neutral but was measured on CUDA only. Options: (a) yes, as a feature after step 1d of
P2.4 (stabilisation), with a check on all three vendors; (b) later in CG-ALF; (c) no. Recommendation: (a).
Also decides whether `rd/cgalf-s1` and `speedup/nvidia` stay.

**Q4. Older questions from the GPU-port sessions (whiteboard card 1).** They were asked for the old `gpu-port`
branch. Since P2.4 now ports into `Code-modularity-GPU`, Claude proposes where each goes:
- Q4a. `prepare_sim.py` and the Helma documentation into the GPU line? The script exists in several analysis folders
  (e.g. `DATA_HOLE/DOUBLE_OCC_TEST/`); which version is meant was not checked. Proposal: yes, the newest version into
  `Scripts_and_Parameters_files/` with the Helma doc, as part of the step `gpu-nvidia` (P2.4 step 2).
- Q4b. Intel switches `ALFSH_OBSCORR=1` and `ALFSH_BLKFUSE=5 ALFSH_BLKPF=1` on by default after the production
  gate? On branch `intel-batched`: L = 24 single chain 3.75 s/sweep (H100 3.51), 4 chains per tile 1.128 sweeps/s per
  card (H100 1.176); gates passed at L = 6/12 (`docs/gpu/platforms/INTEL_BATCHED.md`). Proposal: yes, but in the step
  `gpu-intel` of P2.4, after its own gate on Max 1550; until then they stay opt-in.
- Q4c. Investigate the `Mz = T` deviation? On Intel, L = 6 with `Mz = T` differs from the CPU by 1.36e-10 against a
  tolerance of 1e-10 (`docs/gpu/NOT_IMPLEMENTED.md` line 139; "noise-limited"). Proposal: yes, as part of the
  verification in `gpu-intel` (a longer run decides whether it is noise or a bias).
- Q4d. Finish ticket 407 / "401-A" (non-diagonal or rank > 1 vertices on the GPU, `Calc_Fl` on the device) or leave
  it? Unmerged work: tag `archive/tickets/gpu-port-delayed-submatrix-rank-nd` (K1 passes, K2/K3 untested). Needed for
  Kondo, tV, Spin_Peierls on the GPU. Proposal: P2.6 (feature port), not now.
- Info, no question: PR luisthremer/ALF#1 (GPU det-phase fix into `gpu-port`) has your yes and waits for the
  Ltau = 1 check; the fix goes into P2.4 step 1d anyway. Until merged: no GPU runs with N >~ 2950 (L >~ 54).

## D. Things only you can do

- Q5. Delete the online repo `luisthremer/alf-ops` (its content is in `manager_tools/`; the local clone keeps the
  history): `gh auth refresh -h github.com -s delete_repo`, then `gh repo delete luisthremer/alf-ops --yes`.
- Q6. GitHub MCP: set `GITHUB_PERSONAL_ACCESS_TOKEN` (e.g. in `~/.zshenv`) and restart the app; today the plugin
  sends an empty token ("Authorization header is badly formatted").
- Q7. Switch off `research-superpowers` and `research-desk` in claude.ai (account plugins; Claude cannot).

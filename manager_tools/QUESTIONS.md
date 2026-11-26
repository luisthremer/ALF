# Open questions for the owner (2026-11-26)

Short form in `STATUS.md` ("Waiting for you"). Here each question with its background, options and Claude's
recommendation. Answer in the chat with the number, e.g. "Q1 yes".

## A. Block the next step

**Q1. P2.4 decisions — answered 2026-11-26.** D1: dispatch inside `linalg_mod` (no `alf_la_mod`; step 1a is
redone accordingly). D2: yes (`GPU=host`). D3: yes, Hubbard only, now written into the plan. **D4 still open:**
gpu-port changed seven things on the CPU path without a guard; the proposal is to bring each later as its own
CPU-line feature with its own check, so `GPU=none` stays byte-identical to `Code-modularity-CPU` (details in the
manager chat, 2026-11-26). yes/no

## B. Upstream (ALF-QMC/ALF)

**Q2. Upstream issues — done:** posted as ALF-QMC/ALF#662 (H7) and #663 (H6), 2026-11-26.

## C. GPU line

**Q3. Speedup S1 into the GPU line?** The GPU determinant phase from tracked U determinants instead of an LU
(`rd/cgalf-s1`, origin `speedup/nvidia`). Measured on the H100: -12.6 % per sweep at L = 32, configurations
byte-identical. Already default for Intel on gpu-port (Max 1550, L = 24: 4.45 -> 4.04 s, -9 %); AMD kernels do not track
det(U) yet; STABLOG keeps the LU. Decision memo for the professor: https://claude.ai/artifact/12cwQ6G6urYNRZntuZKFzj Options: (a) yes, as a feature after step 1d of
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
- Q6. GitHub MCP token. The plugin sends `Authorization: Bearer ${GITHUB_PERSONAL_ACCESS_TOKEN}`; the variable is
  unset, so the header is empty. The desktop app does not read `~/.zshrc` (checked), so set it in Claude Code's own
  settings: (1) github.com -> Settings -> Developer settings -> Personal access tokens -> Fine-grained tokens ->
  Generate new token; name `claude-github-mcp`, expiry 90 days, owner `luisthremer`, repositories: only
  `luisthremer/ALF`; permissions: Contents read-only, Issues read/write, Pull requests read/write, Actions read-only
  (Metadata is added automatically). (2) In `~/.claude/settings.json` add at top level
  `"env": { "GITHUB_PERSONAL_ACCESS_TOKEN": "github_pat_..." }` (the file stays on the Mac; never commit it).
  (3) Quit and restart the Claude app. (4) Tell Claude "check GitHub MCP"; it reconnects and tests a read.
  Alternative: switch the GitHub plugin off in claude.ai; `gh` already does everything Claude needs.
- Q7. Switch off `research-superpowers` and `research-desk` in claude.ai (account plugins; Claude cannot).
- Q8. Upstream pull requests waiting for human review: ALF-QMC/ALF#660 (CPU determinant-phase fix) and #655 (AI
  development rules). Only you can ask Fakher / the maintainers; Claude does not contact upstream. Until they are
  merged or closed, `fix/det-c-underflow` and `docs/ai-development-rules` must stay.
- Q9. Clean-up only you decide: empty the macOS Trash when you no longer need `understand-anything-traces-2026-11-26`
  (plugin cache + graph clone); keep or remove the local clone `~/Documents/_WORK/Coding/CLAUDE/alf-ops` (history of
  the old status repo); whiteboard card 2 still mentions the understand-anything graph (not edited by Claude).
- Q10. Cluster logins: SuperMUC-NG Phase 2 needs your login (Claude never types credentials); Helma works through a
  shared ssh master connection that expires after 8 h without use. Not needed for the next steps (all tests run on
  the Mac with GPU=host); needed again for P2.4 step 2 (vendor backends).

# Open questions for the owner (2026-11-26)

Short form in `STATUS.md` ("Waiting for you"). Answer in the chat with the number, e.g. "Q11 yes".

## Open

**Q13. `Lowrank_cross_slice`: keep or drop?** CG-ALF T6 result (CG_ALF_PLAN.md section 7): mathematically sound (the
wrap is a similarity transform, so the cross-slice update gives the same Markov chain up to roundoff), but not faster
on the CPU (on/off = 1.015 at L = 16, i.e. slightly slower; within noise at L = 8/12). Delayed/submatrix themselves give
1.5x over standard at L = 16. Options: (a) do not port it to the GPU and do not bring it to the CM lines (it stays on
gpu-port/tags only); (b) bring it to the CPU line as an opt-in feature anyway. Recommendation: (a).

## Decided (interview 2026-11-26; record with your exact words in the session notes)

- **D4 (Q1): option (b).** The 8 unguarded CPU-path changes of gpu-port come later as separate CPU-line features (one
  roadmap line each). Exceptions: (1) `QR_PREPIVOT` may come in P2.4 steps 1d/1e as a test-only flag, default off
  (GPU=none stays byte-identical); (2) `ALF_UNSTABLE_THRESHOLD` goes first in the CPU-line queue as its own small
  feature. `N_boot` rather with the MaxEnt branch; `Op_V_is_symmetric` = CG-ALF T3.
- **S1 (Q3): yes, option (a).** "Okay, just do it." Own feature after step 1d of P2.4. No LU: document the tracked phase
  (why no LU; invariant: U stays the exact product of the tracked reflectors; check mode). CG-ALF task: check that
  tracking is the right choice in general. Also wanted for CG-ALF: the LU-free phase on the CPU line.
- **Q11:** "Yes, that is what I meant": the same tracked-phase method on every vendor where possible, not a per-vendor
  patchwork. And: "the commit for no Lu has to be a single commit! so that it can be undone later!" (one commit that a
  single `git revert` removes).
- **Q4a:** `prepare_sim.py` is a general PARALLEL_PARAMS helper, not GPU- or model-specific: own feature (CPU line,
  newest version from gpu-port `Scripts_and_Parameters_files/Start_PP/`, with the Helma documentation).
- **Q4b:** "Yes, default, but report in info file!" Intel switches `ALFSH_OBSCORR=1`, `ALFSH_BLKFUSE=5`, `ALFSH_BLKPF=1`
  on by default in P2.4 step `gpu-intel` after their own check on Max 1550; active values written to the `info` file.
  CG-ALF: port the Intel ideas to CUDA/AMD; timing flag, fixed benchmark list, start-up auto-selection, mechanism for
  new hardware.
- **Q4c:** investigate the `Mz = T` deviation on Intel (L6 1.36e-10, L12 1.21e-9: grows with L) with claude-bughunter
  (Claude asks you to switch it on then).
- **Q4d:** delayed and submatrix updates for every model first (Kondo, tV, Spin_Peierls: ticket 407 / 401-A), then a
  benchmark of which scheme is best. `Lowrank_cross_slice`: investigate now in a parallel CG-ALF research branch
  (mathematically sound? fastest on CPU and GPU?), within the WIP limit.
- **Q12:** claude-bughunter for Q4c: Claude asks you to switch it on when the `Mz = T` investigation starts (P2.4 step
  `gpu-intel`).
- **Q14 (T1):** "Yes" (button on the progress page, 2026-11-26 20:14 UTC). Measure the stabilisation candidates
  (pre-pivoted unpivoted QR, CholeskyQR2 after scale separation, adaptive Nwrap) on the Helma H100, preempt partition,
  a few GPU-hours: cost per stabilisation and Precision Green against Nwrap at L = 24 and 32. Session cgalf-research.
- **Q6:** GitHub plugin switched off by you (2026-11-26).
- **Q2, Q5, Q7, Q8, Q9, Q10: closed.** Upstream issues posted (#662, #663); `luisthremer/alf-ops` already deleted; research
  plugins off; Fakher asked about ALF-QMC/ALF#660 and #655 (branches `fix/det-c-underflow`, `docs/ai-development-rules`
  stay until merged/closed); Trash emptied, local alf-ops clone moved to the Trash; cluster logins asked for when P2.4
  step 2 needs them.
- Earlier: D1 dispatch in `linalg_mod`, D2 `GPU=host`, D3 Hubbard only (2026-11-26).

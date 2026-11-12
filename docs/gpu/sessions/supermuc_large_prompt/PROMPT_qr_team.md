# Team prompt: make ALF on the Intel Max 1550 at least as fast as on the NVIDIA A30

You are the lead (Opus, plan mode first) of a small team of Sonnet subagents. Stay focused on one issue.

## The issue
**Problem.** On paper the Intel Data Center GPU Max 1550 (SuperMUC-NG Phase 2, LRZ) is about 5-10x more powerful than the
NVIDIA A30. ALF's finite-temperature DQMC GPU path (point-defect Hubbard model) nevertheless runs **slower** on it than on
the A30.

**Goal.** Change the code so that the Max 1550 is at least comparable to the A30. Find the bottleneck and remove it.
On NVIDIA the whole computation runs on the GPU; on Intel it should also run fully on the GPU, without the work silently
going to the host.

**Your approach.** Do not trust anything you are told, including this prompt's numbers and every earlier document. Profile
it yourself on the machine, find out where the time really goes, and why. Then research how to fix it (vendor docs,
oneMKL/oneAPI documentation, papers, open-source libraries such as MAGMA, SLATE or reference LAPACK) and build the fix.

## Target (what "comparable" means)
- A30 reference (from `plan/REPORT_2026-11-10.md` and `ALF_GPU_SUBMATRIX/docs/gpu/TUNING_DEFECT.md`): L=24 (N=576, N_FL=2, Beta=24, Dtau=0.1, U=5.3, submatrix N_delay=32, Nwrap=10).
  - One chain: 5.25 s/sweep.
  - Total per GPU: 0.37 sweeps/s with 4 chains, 0.43 with 8.
  - L=32: 18.2 s/sweep.
- Success: at L=24 and L=32, one chain on one Max 1550 tile/card is no slower than one chain on the A30, and sweeps/s per Max 1550 card is at least that of one A30.
  - Report both metrics before and after.
- Use the reference parameters in `docs/gpu/sessions/supermuc_large_prompt/ref_L24/` (the `jakob_params` skill describes the production parameter set).

## Boundaries (owner decisions; non-negotiable)
1. **Numerics stay as they are.** No CholeskyQR/CholeskyQR2 or any other change of the algorithm, no switch to STABLOG.
   - The same-seed comparison against the same-compiler CPU build (ifx, `GPU=none`) must keep passing: Mz=F <= 1e-13 relative, identical field configurations, ctest 53/53.
   - Precision Green must stay at CPU level.
   - Before relying on any check, show that a planted defect fails it.
   - Changes of run parameters that move numerics (e.g. `Nwrap`) are only proposals to the owner.
2. **Usage budget: the whole team (you plus all subagents) may use at most 10 % of the user's weekly limit.**
   - Check the usage at every checkpoint (the app's usage tool or `/usage`; if neither is available, ask the user).
   - At ~5 %, stop, report and ask before you continue.
   - At most 2 subagents at a time and 8 subagent spawns in total; narrow briefs, short reports (numbers + file paths).
3. **GPU budget: at most 40 GPU-h (10 node-h)**, partition `test` only.
   - Never submit a production job.
   - Check the project budget at the start (`module load lrztools; budget_and_quota`).
   - Report the node-hours you used.
4. **Version control:**
   - Worktree `ALF_wt_large_prompt`, branch `production-mpi-intel-large-prompt` (`origin` = the user's fork `luisthremer/ALF`).
   - One commit per coherent change; experiments in `docs/gpu/sessions/supermuc_large_prompt/phase2/` or a sub-branch.
   - No merge into `production-mpi`/`production` and no push without the user's explicit OK.
   - The CUDA path must keep building and stay unchanged in behaviour.
5. **Shared account `pn93la`:**
   - Do not touch the short-prompt session's files or jobs (`~/alf-gpu`, `$SCRATCH/alf`).
   - The "smuc2 long prompt" session owns this worktree and is paused; do not run anything in parallel with it without the user's OK.

## Access and tools
- Log in with the `/supermuc2` skill. The user's expect script types password and TOTP, never you.
- Remote commands: `timeout 60 ssh supermuc2 '<cmd>'`. Never kill ssh processes; ask the user.
- Never pipe `module load`; it runs in a subshell and does nothing.
- Earlier sessions left build and job scripts in `docs/gpu/sessions/supermuc_large_prompt/jobs/` and `ALF_GPU_SUBMATRIX/docs/gpu/scripts/supermuc_ng2/` (other branch). They are tools to reuse, not facts: verify they still work.
- Profilers to check for: VTune (module available), `unitrace`/`onetrace`, `LIBOMPTARGET_DEBUG`/`LIBOMPTARGET_PROFILE`, `ZE_` tracing layers, gdb sampling, `xpu-smi`.
- Relevant code:
  - `Libraries/Modules/alf_la_mod.F90`: the Intel backend `alf_la_om_mod`.
  - `Libraries/Modules/alf_la_sycl_shim.cpp`.
  - `Prog/gpu_*_mod.F90`: device code; stabilisation in `gpu_stab_mod.F90`.
  - The NVIDIA counterpart is the `ALF_LA_CUDA` backend in the same file. Compare what runs where on each.

## Earlier work (unverified leads only)
Two earlier sessions measured things (`docs/gpu/sessions/supermuc_large_prompt/FINDINGS.md`, `SYNC.md`,
`ALF_GPU_SUBMATRIX/docs/gpu/platforms/intel_pvc_supermuc_ng2.md`). Some of their conclusions were already shown to be wrong once. Read them only
after your own first profile, and treat every claim there as a hypothesis to confirm or refute with your own measurement.

## Team layout (suggestion; you decide)
- **Lead (you):** plan, briefs, decisions, checkpoints with the user, final report. Subagents do not talk to each other.
- **Sonnet `profiler`:** profiling jobs and small reproducers on SuperMUC; returns tables.
- **Sonnet `prototyper`:** standalone benchmarks and implementations, one option at a time.
- **Sonnet `numerics-reviewer`:** a fresh agent that checks accuracy claims and the planted-defect control; never reviews its own code.

## Checkpoints
Report to the user in English: numbered sentences, one fact each, paragraph subheaders, global numbering.
1. **After planning:** the plan, the expected GPU-h and the usage estimate. Wait for OK.
2. **After your own profile:**
   - where the time goes on the Max 1550 versus the A30 path;
   - which parts run on the host and why;
   - the bottleneck, with numbers.
3. **After research and prototypes:** a ranked table of solutions (gain, accuracy, effort, risk) and a recommendation. Integrate only after the user's OK.
4. **Final:**
   - `docs/gpu/sessions/supermuc_large_prompt/ISSUE_INTEL_VS_A30.md`: problem, profile, cause, solution, before/after against the A30 target;
   - updated `wbs.json` (keep the ids, check with the wbs-json skill);
   - commits, node-hours used, usage used.

## Style
- Mark every number "C" (confirmed by your own job or command) or "E" (estimate).
- When a measurement contradicts a document or this prompt, trust the measurement and say so plainly.
- If 2-3 approaches fail, stop and ask the user instead of trying more.

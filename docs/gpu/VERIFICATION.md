# Verification standard (2026-10-29; approved by the owner subject to a second independent review, DECISIONS.md)

Goal: no error reaches a commit. Correctness ranks above speed. Tooling: `testsuite/gpu_reference/acceptance/`
(README there). Owner decisions: [DECISIONS.md](DECISIONS.md). Measurements behind the rules: gate report of
2026-10-29 (determinism, MPS, CPU-alongside) and the research notes R4–R9 (outside the repo).

## 1. Test the artifact we commit
- Acceptance runs only on binaries built from `git archive <candidate-sha>` (`acc_ship.sh`, `acc_build.sh`): the build
  tree is extracted from the tarball after its sha256 is re-checked, and the build record carries that hash, which the
  verdict checks; each party (implementer, manager) has its own candidate root, so T2 never reuses T1's binaries. A working tree
  is not a candidate. The tools that judge a candidate (`acceptance/`, `compare2.py`, `configs*.py`) and its acceptance
  file must be unchanged since the pre-registration commit; `acc_ship.sh` refuses otherwise.
- The committed sha equals the accepted sha, except for `docs/**` (`acc_check_commit.sh`). Any other change after
  acceptance, including a comment in a `.F90`, needs a new acceptance run of the deterministic tier; the reviewer
  decides and records whether the statistical tier must be repeated.

## 2. Provenance stamp on every run
`acc_run.sh` writes `provenance.json` into every run directory: candidate sha, tree hash, tarball sha256, binary
sha256, configure line and flags, compiler/CUDA/driver versions, the complete environment, parameter and seed
hashes, GPU processes before and after, start/end time, exit code, timeout flag, sha256 of every output file and of
every `data.h5` dataset. Runs use `env -i` plus the tool's env file plus an explicit whitelist. Debug and fault
variables (`ALF_GPU_FAULT*`, `_TAUSWAP`, `_TAU_NS`, `_NS`, `_STABCHECK`, `_STABMIX`, `_TAUDIAG*`, `_NOUPDATE`,
`_SELFTEST`, `_UPDTEST`, `_EXIT`, `_TAUTIME`, `_STABTIME`, `_DISABLE`) are refused unless the acceptance file declares
the run as a self-test or control run, or as a diagnostic gate run (the 207 paired test needs `ALF_GPU_TAUDIAG`; the
  file says so in the gate's reason). GPU builds run with `OMP_TARGET_OFFLOAD=MANDATORY`, and every run of a GPU build
  must show the ACTIVE lines of the GPU paths; an INACTIVE line is an input error. A log without a stamp that matches
  the candidate is not evidence.

## 3. Pre-registration
Each ticket has `testsuite/gpu_reference/acceptance/NNN.json`, committed **before** implementation results are seen,
and frozen with the candidate (the driver refuses any other copy). It holds configs, seed sets, bins and sweeps, each
gate's statistic, threshold and reason, α and the multiple-comparison correction, the expected counts, and the
positive controls with their expected outcome. Development uses the committed dev seeds; acceptance seeds are drawn
fresh by the driver (`secrets.randbits`) after the candidate is frozen. Changing the file needs the owner's quoted,
dated approval (in the file's `approvals` and in DECISIONS.md) **and** fresh seeds afterwards. A statistical FAIL goes
to the owner; "re-run until green" is forbidden. When a criterion fails, relaxing it is never put forward as the first
or recommended option.

## 4. Scripted verdicts
`verdict.py` computes PASS/FAIL for every gate from raw outputs and writes `verdict.json`/`verdict.md`; report tables
are generated from it, never assembled by hand. NaN or Inf, a missing file, line or dataset, an empty comparison set,
a count differing from the acceptance file, a nonzero exit, a timeout, a GPU-inactive run and a provenance mismatch
are never a PASS: they give INPUT_ERROR (the gate could not be evaluated), which is also never counted as a
detection. FAIL means the gate's own criterion caught something. BLOCKED items (owner-decided, e.g. 207 acceptance 6)
make the overall verdict BLOCKED, never PASS. `verdict.py --selfcheck` (synthetic data, CPU only) must PASS whenever
the script changes; `verdict.py --lint` must pass for every acceptance file (the driver refuses otherwise).
`compare2.py` is fail-closed since 2026-10-29 (exit 3: NaN, missing observables, m = 0, identical seeds or means
without `--same-seed`). The per-ticket `tNNN_*` scripts are development tools only.

## 5. Positive controls for every gate
Each gate has at least one control, run on the same build profile and judged with the gate's own criterion
(`"like": <gate>`), that must FAIL; `verdict.py --lint` refuses an acceptance file with an uncovered gate. A control
is killed only by the criterion (outcome FAIL): a crash, a timeout or a missing directory is an INPUT_ERROR and never
counts as a detection. Planted NaN or dropped datasets are expected as INPUT_ERROR. A gate whose control is not killed
is not a gate.

| Gate | Control (must FAIL) |
|---|---|
| kernel self-test | the `ALF_GPU_FAULT*` classes; a graded fault near 10× the tolerance (hardening ticket) |
| same-seed vs CPU oracle | `mutate.py scale:<tau dataset>:1e-3` on the smallest entry above 1e-5 of the maximum (FAIL; detection limit ~1e-9 × max\|dataset\|); `nan:`, `drop:` (INPUT_ERROR) |
| paired same-trajectory | `ALF_GPU_TAUSWAP=8` (LOG), `ALF_GPU_TAU_NS=0` (default); `mutate.py taudiag_inflate/taudiag_bias` |
| precision mean | `ALF_GPU_NS=0` (once per toolchain); known blind spot: the 207 ZUNMQR defect (+17 % on the τ mean) |
| compare2 | `mutate.py shift:<obs>:5` (5 standard errors of the mean, FAIL); `nan:` (INPUT_ERROR) |
| GPU=none identity | a planted declaration appended to a compared file (`acc_gpunone.sh --plant`); a copy of the ctest log with one failed test; other seeds (`data.h5` identity) |
| replay / determinism | a 1-ulp-level difference (`mutate.py scale:…:1e-12`) |

## 6. Statistical gates with known power
Calibrate each with a null run (CPU vs CPU, fresh seeds) and a known small bias. Errors are autocorrelation-corrected
from ≥ 30 independent bins. The CPU reference has ≥ ~10× the statistics of the test or errors are combined in
quadrature. Next to every PASS the report states the detectable effect size. Every compare2 gate carries a **null
control** (CPU vs CPU on two fresh seed sets, `"null": true`, expect PASS; enforced by the lint): in the harness test
of 2026-10-29 a short run (40 bins × 5 sweeps, 6×6) gave a false FAIL on `E` (t = −3.5) that the identical CPU-vs-CPU
comparison reproduced digit for digit, while GPU and CPU agreed to 1e-11 on the same seeds. **Approved 2026-10-29:** maximum-error
statistics are information only; acceptance is on the mean (a fixed quantile over a fixed number of stabilizations
only once calibrated), and maxima are compared only between runs of identical length.

## 7. Determinism (from the measurement of 2026-10-29)
- Where bit identity holds — same binary re-run, same commit rebuilt with the same toolchain — it is demanded.
  Measured 2026-10-29 (stamped, `env -i`, 12×12 β = 10, same binary twice): the LOG path is bitwise reproducible
  (all 38 datasets and the full-precision Precision lines); the default path is **not** (20 of 38 datasets differ, first
  in `Den_eq/back` at ~1e-12 relative; Precision Green mean 6.188 vs 6.180e-10). Until the hardening ticket (212) makes
  it reproducible, replay and non-interference gates run on the LOG path only, and this is stated in the file.
- MPS and CPU jobs alongside do not change results: LOG runs under MPS (2 clients), next to a time-sliced second
  process, or next to 2–3 CPU runs gave `data.h5` bit-identical to a solo run (2026-10-29).
- Across source changes: the documented relative bounds (same-seed per-bin 1e-9 overall and in the τ line, no bin
  above 1e-4) and statistical comparisons.
- Precision summary numbers are never compared to 4 digits.

## 8. Who runs what (tiers; tolerances identical in every tier)
| Tier | Who | Content | GPU time (207-sized ticket) |
|---|---|---|---|
| T0 fast | implementer, every iteration | self-tests and fault controls; 6×6 same-seed both paths | ~2.5 min |
| T1 pre-candidate | implementer, once, on its own frozen candidate (`git archive`; `acc_driver.py --who agent` maps every acceptance set accN to devN = dev seeds rotated by 37·N lines) | the full acceptance on dev seeds | ~90 min |
| T2 acceptance | manager, frozen candidate, unattended (`acc_driver.py`) | all deterministic gates in full; every statistical gate once on fresh seeds; all positive controls | ~120 min |

On top of T2 the manager audits the implementer's report against the raw outputs (`audit.py`) and lists every
mismatch, and an independent `acceptance-verifier` recomputes every verdict (section 10). A ~2 h unattended T2
sequence of ≤ 15 min runs needs the owner's standing approval of the budget.

## 9. Evidence retention
Raw outputs (`info`, `data.h5`, `run.log`, `provenance.json`, verdict files) stay until the owner's gate, then go to
`~/alf-gpu/archive/tNNN/` with a sha256 manifest. Budget: ~0.6 GB per acceptance; the driver refuses to start below
5 GB free; archive ≤ 2 GB unless the owner approves more. Agents never delete raw evidence; cleanup removes build
trees only. Run directories are never overwritten (the scripts refuse).

## 10. Independent review and verification
- `numerics-reviewer` (Opus) reviews three things separately: (a) the code diff, (b) the acceptance file, (c) scripts
  and verdict tooling. (b) and (c) are reviewed before T2 starts; a finding on (b) after T2 data exists is a criteria
  change (section 3). Any diff touching tests, references, thresholds or seeds is flagged in the review and in the
  owner's gate report.
- `acceptance-verifier` (Opus, fresh, read-only) gets only the candidate sha, the acceptance file and the manager's
  raw outputs; it recomputes every verdict with its own code, then reads the implementer's report and the manager's
  summary and lists every mismatch. A ticket closes only when the scripted verdict, the verifier and the reviewer agree;
  disagreements go to the owner.
- The implementer never edits the acceptance file after the freeze, and never certifies its own corrections.

## 10b. Agent lifecycle
- Implementers on Sonnet by default; the manager may propose Opus per ticket for numerically delicate work.
- Verification never runs in the implementer's context. Corrections of a report go to a fresh verifier.
- When an implementer passes ~250k tokens, or at a phase boundary (candidate → review fixes), it writes a handoff note
  into the ticket (state, decisions, what was tried, open problems, evidence paths) and a fresh agent continues.
- Resuming the same agent is fine only for interruptions (rate limits, deadline changes). After a change to an agent
  definition, spawn new agents.
- The manager restarts per phase from a written handoff; every owner decision lives in the ticket and DECISIONS.md.

## 11. CPU reference cache
`refcache.py` stores CPU reference runs in `~/alf-gpu/refcache/` keyed by the CPU-path sources, build, toolchain,
parameters and seeds; any mismatch is a miss. Only the manager writes to it; budget ≤ 2 GB (oldest entries removed
first, after the owner's gate of the tickets that used them). Details in [BASELINE.md](BASELINE.md).

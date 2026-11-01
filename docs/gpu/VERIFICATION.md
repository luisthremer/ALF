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
`_SELFTEST`, `_UPDTEST`, `_EXIT`, `_TAUTIME`, `_STABTIME`, `_DISABLE`) and the `alf_la` test variables `ALF_LA_TEST_*`
(since 2026-10-31) are refused unless the acceptance file declares the run as a self-test or control run, or as a
diagnostic gate run (the 207 paired test needs `ALF_GPU_TAUDIAG`; the file says so in the gate's reason). GPU builds run
with `OMP_TARGET_OFFLOAD=MANDATORY`; `GPU=intel` builds also with `ONEAPI_DEVICE_SELECTOR=level_zero:0`. Both are set by
`acc_run.sh` and `acc_ctest.sh` after the env file and asserted, not overridable; an env file that sets
`ZE_AFFINITY_MASK` or another device-ordering variable, or tries to change the harness's own variables, stops the run.
`verdict.py` checks the recorded env of runs, `acc_ctest.sh` the recorded env of each ctest repetition (INVALID
otherwise). Every run of a GPU build must show the ACTIVE lines of the GPU paths; an INACTIVE line, or a line the
acceptance file lists in `forbidden_lines`, is an input error. A log without a stamp that matches the candidate is not
evidence.

## 3. Pre-registration
Each ticket has `testsuite/gpu_reference/acceptance/NNN.json`, committed **before** implementation results are seen,
and frozen with the candidate (the driver refuses any other copy). It holds configs, seed sets, bins and sweeps, each
gate's statistic, threshold and reason, α and the multiple-comparison correction, the expected counts, and the
positive controls with their expected outcome. Development uses the committed dev seeds; acceptance seeds are drawn
fresh by the driver (`secrets.randbits`) after the candidate is frozen. Changing the file needs the owner's quoted,
dated approval (in the file's `approvals` and in DECISIONS.md) **and** fresh seeds afterwards. A statistical FAIL goes
to the owner; "re-run until green" is forbidden. When a sound criterion fails, relaxing it is never put forward as the first
or recommended option. Sound means calibrated: correct code (e.g. CPU variants with another rounding order) passes it
and a planted defect fails it; if measurements show that correct code cannot meet a criterion, the manager tells the
owner plainly, with the numbers, and recommends a calibrated replacement (CLAUDE.md, working principle).

## 4. Scripted verdicts
`verdict.py` computes PASS/FAIL for every gate from raw outputs and writes `verdict.json`/`verdict.md`; report tables
are generated from it, never assembled by hand. NaN or Inf, a missing file, line or dataset, an empty comparison set,
a count differing from the acceptance file, a nonzero exit, a timeout, a GPU-inactive run and a provenance mismatch
are never a PASS: they give INPUT_ERROR (the gate could not be evaluated), which is also never counted as a
detection. FAIL means the gate's own criterion caught something. BLOCKED items (owner-decided, e.g. 207 acceptance 6)
make the overall verdict BLOCKED, never PASS; a ticket with a BLOCKED verdict closes only with the owner's explicit OK. A gate may be `report_only` (run and shown, not deciding) only by an
owner decision cited in its reason (the lint checks it), e.g. 207's LOG 12×12 same-seed pair, split off to ticket 213. `verdict.py --selfcheck` (synthetic data, CPU only) must PASS whenever
the script changes; `verdict.py --lint` must pass for every acceptance file (the driver refuses otherwise).
`compare2.py` is fail-closed since 2026-10-29 (exit 3: NaN, missing observables, m = 0, identical seeds or means
without `--same-seed`). The per-ticket `tNNN_*` scripts are development tools only.

## 4b. QREF accuracy gates (ticket 207-5; decisions 33, 35, 36, 37)
Gate type `qref` judges the accuracy of the device's τ-pass CGR2_2 output against a binary128 reference by calling the approved 207-2 gate
`docs/gpu/scripts/nividia/t207_d5_gate.py --spec --spec-sha --prog --report` on the harness's own evidence (`qref_gate.py`); the harness decides
nothing about accuracy itself. Its exit status maps to the verdict: PASS 0 → PASS, FAIL 1 → FAIL, INPUT_ERROR 2, DEV-ONLY 3 (a gate no
control trips), TEST-ONLY 4 (spec marked `test_only`) and anything else → INPUT_ERROR.
- **Pre-registered** in the acceptance file: the gate's `spec` (the 207-2 spec: path, `host_member`, members per role, frozen `refs`, qs, null
  codes, cover 1, KEEP, HOSTIN, `min_reprocessed`, expected points per measure, `src_sha`, env roles with stride 3 (decision 35), `n_runs`,
  every gate with measure, statistic, pool and bound, every control) with `spec_template_sha256` = sha256 of its canonical form (sorted keys,
  no blanks, newline), `runs` (main), `ctl_runs` (`{role: [run ids]}`), `tools_sha256` of the 207-2 scripts, and per run a `qref` entry
  (path, qs, nul, hostin, keep = the spec's). `acc_ship.sh` also pins `docs/gpu/scripts/nividia/t207_d5_*` with the judging tools.
- **Instance fields** (known only after the candidate is built and the seeds are drawn, section 3): `binary_sha256`, `prog_sha`,
  `parameters_sha256`, `seeds` may be `"@harness"`; the verdict fills them from its records (ALF build record from the verified tarball, offline
  build record `bin/qref_<path>.qbuild.json` from the same tarball, the runs' parameters, the driver's seeds record per role); a concrete value
  must equal the harness value. The filled spec and its sha256 (passed as `--spec-sha`) and the gate report are written next to `verdict.json`.
- **Runs** (`acc_run.sh --qref`): only `ALF_GPU_*` extra variables; the verdict refuses any recorded variable the gate cannot see: anything
  `ALF*` besides the declared `ALF_GPU_*` tokens (`ALFSH_*`, `ALF_LA_TEST_*`, ...), `LD_PRELOAD`/`LD_AUDIT`, `OMP_TARGET_OFFLOAD` ≠ MANDATORY,
  thread counts ≠ 1, and every name outside the harness's own variables and the file's `qref_env_allowed` (what the env file sets; the lint
  refuses debug, fault and steering names there). This holds for control runs too: their faults are the gate's `ALF_GPU_*` control roles.
- **Processing** right after each run, before the next (`acc_qproc.sh` → `t207_d5_proc.sh`): bridge records in the 207-2 layout generated from
  the harness provenance (the reviewed parser is not taught a second format), strict check, manifest, thinning by the gcd keep rule; the
  interpreter (`D5_PYTHON`) and the machine (host, `/etc/machine-id`, the libopenblas the program loads) are recorded.
- **Stamps:** after the processing (and after `mutate.py qref_swap` for a control copy) the driver writes a `qref_stamp` line per run into the
  attempt log outside the run directory: sha256 of every file of the run directory (records, `qref.txt`, `qref_proc.json`, `qref_manifest.txt`,
  `acc_qproc.json`, every kept dump). Before calling the gate the verdict demands exactly one stamp per run, no file added, missing or changed,
  no symbolic link; it checks them again after the gate. It also checks the offline build record (candidate, tarball, source = `src_sha`,
  program hash, scripts), the bridge records against the harness provenance, and that processing and verdict run on the same machine (the
  gate's exact re-processing depends on the system OpenBLAS).
- **The gate** (207-2, approved) re-derives everything from the run directories, re-runs the pinned offline program on the kept dumps
  (exact equality; the dumps are hashed again afterwards), pins binary, parameters and seed sets, and decides device and controls. Its
  limits (only every KEEP-th dump re-computed; only `ALF_GPU_*` visible; binary↔source pairing) are covered by the stamps, the env check and the
  builds from the tarball above.
- **Controls:** inside the gate the stand-in `var8` and the production control role (`TAU_NS=0` / `TAUSWAP=8` / the 207-3 fault
  `ALF_GPU_FAULT_TAU=eps:<x>`, optionally localized by `ALF_GPU_FAULT_TAU_AT=<nst>:<nfe>`, which needs `FAULT_TAU`) must be killed,
  else FAIL; the harness control `{"like": <qref gate>, "replace": {"runs": [...]}}` judges copies of the main runs made by
  `mutate.py qref_swap:dev:var8` (device and var8 exchanged in the kept dumps and in `qref.txt`, so that the re-processing agrees) and must FAIL.
- Gate type `h5_identical`: `data.h5` of a QREF run and of a production run (same build, config, seeds), every dataset byte-identical;
  control: a mutated copy (`mutate.py scale:…:1e-12`).
- `verdict.py --selfcheck --qref-fixture CANDROOT` (on the machine of a finished QREF dry run) runs the cases against the real gate: NaN,
  Infinity, an E-less number, a missing member, a missing point, shifted fields, a changed header, device and var8 swapped (in `qref.txt` only:
  INPUT_ERROR; consistently: FAIL), a control not killed, a defect in one of two statistics, DEV-ONLY, TEST-ONLY, an extra `ALF_GPU_*`,
  `ALF_LA_TEST_*`, `ALFSH_*`, `OMP_TARGET_OFFLOAD=DISABLED`, `LIBOMPTARGET_DEBUG`, `LD_PRELOAD`, an unknown variable, a changed script or
  program hash, altered records, a deleted or symlinked kept dump, another machine, a missing or doubled stamp, a changed spec. Without the
  fixture the self-check ends INCOMPLETE (exit 4), never PASS.

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
| precision mean | `ALF_GPU_NS=0` (once per toolchain; on the B580/oneMKL it cannot degrade precision, so 302b uses the graded device fault `ALF_GPU_FAULT_STAB=eps:<x>` instead, decisions 26–27); known blind spot: the 207 ZUNMQR defect (+17 % on the τ mean) |
| compare2 | `mutate.py shift:<obs>:5` (5 standard errors of the mean, FAIL); `nan:` (INPUT_ERROR) |
| QREF accuracy (207) | the stand-in `var8` and the production control run (`TAU_NS=0` / `TAUSWAP=8` / graded fault) inside the gate; `replace: {runs: [qref_swap copies]}` (`mutate.py qref_swap:dev:var8`, FAIL); a mutated `data.h5` copy for the QREF identity pair |
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
sequence needs the owner's standing approval of the budget. Runs have no fixed time limit (decision 16): the acceptance
file's `timeout_s` is a hang guard, never a criterion, and run times are not acceptance criteria unless the ticket is about
performance.

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

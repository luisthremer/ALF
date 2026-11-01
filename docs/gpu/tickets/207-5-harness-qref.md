# 207-5 – acceptance harness support for the QREF accuracy gates

- **Phase / lane:** 2 (NVIDIA) · part of 207 (decisions 33, 34)
- **Owner agent:** fortran-gpu-dev (Sonnet); manager integrates
- **Review:** numerics-reviewer (scripts, part c) · merged into `gpu-port` before 207-6 (the judging tools must be in
  207's base)
- **Depends on:** 207-2 (output format and parser)
- **Branch:** `gpu-port-207-harness`
- **Status:** first round done (WIP `123f7c89`); second round running on `3466c2f9` (207-2 gate approved)

## Goal
`acc_driver.py` / `verdict.py` can run and judge the QREF gates from a pre-registered acceptance file (reviewer
change 6, retention of decision 33).

## Tasks
1. Run type "QREF diagnostic": the env allow-list of `t207_d5_qparse.py` (`env_role`): main runs exactly `DEBUG=1`,
   `TAUDIAG_QREF`, `_QHOST=1`, `_QVAR=8`, `_FROM`; control runs that set plus `TAU_NS=0` / `TAUSWAP=8` / the 207-3 graded
   fault, in declared control runs only. Stride 3 (decision 35).
2. Build the offline program and the GNU oracle objects from the candidate's `git archive`, with provenance.
3. Process each run's dumps before the next run starts (nividia has ~11 GB free; ~1.6 GB per run): verdict inputs
   computed first, then retention as in 207-2 (manifest, full `qref.txt`, kept dumps by the gcd keep rule, KEEP 39).
   Reuse the 207-2 scripts (`t207_d5_proc.sh`, `t207_d5_qparse.py`, `t207_d5_agg.py`, `t207_d5_gate.py`) rather than
   re-implementing them; the harness calls them and records their hashes.
3a. Both machines (decision 36): QREF runs on `intel` (B580, `GPU=intel`, after 302c-1) with the offline program built
   there (gfortran 15.2 present on `intel`, real(16)); the Intel reference ensemble has an MKL/ifx host member in place
   of `nvh` (name fixed with 302c-2); keep the member list a parameter of the acceptance file.
4. Gate type computing S_ens(pTT) and the pointwise statistic from `qref.txt` with the frozen ensemble list from the
   acceptance file; input checks (qself, qbr ≤ 1e-25; n0 identical to cpu; nvh present; device branch = host branch
   at every point; point count) → INPUT_ERROR.
5. Identity check: `data.h5` with QREF on and off, same seeds, identical.
6. `verdict.py --selfcheck` mutations: NaN, missing member, missing point, shifted fields, device and var8 swapped.

## Files allowed
`testsuite/gpu_reference/acceptance/` (scripts, README), `docs/gpu/VERIFICATION.md` (gate description).
Second round (2026-11-01): also `docs/gpu/scripts/nividia/t207_d5_{gate.py,qparse.py,agg.py,proc.sh,selftest.py,pairs.py}`
(207-2's approved gate; 207-3 does not edit them). Base: `gpu-port-207-harness` `3466c2f9` (= WIP + `gpu-port-207-r`).

## Done when
Self-checks pass, an end-to-end dry run on a dev build passes, review (c) approves, merged into `gpu-port`.

## State (2026-11-01)
- WIP `123f7c89`: `--qref` runs, `acc_qbuild.sh` (offline program from the candidate tarball), `acc_qproc.sh` (processing
  before the next run, keep rule), `qref_gate.py` (gate types `qref`, `h5_identical`), selfcheck with synthetic
  `qref.txt` (PASS). Dry run on nividia (fake candidate, 12×12 β = 10, 3 bins, stride 3, 60 points): retention and the
  identity gate work (PASS, control killed); the QREF gate FAILED as it should because var8 (+0.078) and the TAU_NS=0
  control (+0.044) were not killed at the placeholder 0.15 on this short run (6×6 β = 10: −0.10 / −0.14). Not a
  calibration; passed to 207-3/207-4.
- Second round, after the 207-2 redesign: replace `statistic()` by the spec-driven gate (spec embedded in the acceptance
  file, including the pointwise bound and control sizes); add `docs/gpu/scripts/nividia/t207_d5_*` to `TOOLS` in
  `acc_ship.sh`; test the shipping path with a real commit; Intel member name from 302c-2; a real `~/alf-gpu/env/gnu.sh`
  on `intel` (manager).
- From the 207-2 gate reviews (rounds 3 and 4, 2026-11-01), the gate relies on the harness for:
  - **Stamps:** sha256 of every kept dump, `qref.txt`, `qref_proc.json` and `qref_manifest.txt` in the attempt log,
    outside the run directory; the gate's `--spec-sha` value and the spec itself in the acceptance file.
  - **Env outside `ALF_GPU_*`:** the gate sees only `ALF_GPU_*` (info line). `ALF_LA_TEST_*`, `ALFSH_*`,
    `OMP_TARGET_OFFLOAD` (DISABLED turns "dev" into a host result) and other debug/fault variables must be refused by
    the harness's env allow-list and provenance check, outside declared control runs.
  - **Binary ↔ source:** the gate takes `prog_sha`/`binary_sha256` ↔ `src_sha` from the spec; building both from the
    candidate tarball with provenance (task 2) is what makes the pairing true.
  - **Record format:** `t207_d5_qparse.derive_run` expects the dev-run records (`rc`, env list, `run_log_sha256`,
    `info_sha256`, `stamp.txt`); `acc_run.sh` writes `alf-acc-prov/1` (`exit_code`, env dict, `inputs.seeds_sha256`, no
    `stamp.txt`). Unadapted, every harness run would be INPUT_ERROR. Adapting means changing reviewed gate code:
    re-run the gate self-test and get a fresh review of that diff.
  - **Spec after seeds:** the spec pins seed hashes per role, so it is written after the fresh acceptance seeds are
    drawn (207-6).
  - **Intel host label** (207-2 round 5): gate, qparse and agg take the host member's name from the spec, but
    `t207_d5_qref.F90` prints `nvh` itself. For Intel either the offline program takes the label as an argument (change
    in `t207_d5_qref.F90`, then re-review) or Intel specs use the label `nvh` for the in-run host member of the ifx
    binary; decide in the second round and say which.
  - `t207_d5_proc.sh` calls `~/alf-gpu/venv/bin/python3` directly: needs a venv on `intel` or a configurable interpreter.
  - `t207_d5_pairs.py` (dev5-only tool) breaks on `agg.ref_list`'s new signature: fix or retire.
  - From the round-5 review (approved): correct the comments that say NVIDIA and Intel builds print different host
    labels (`t207_d5_gate.py:164-165`, `t207_d5_qparse.py:51-52`, `t207_d5_agg.py:17`; the program always prints
    `nvh`; reviewer advises against a label argument, so Intel specs use `nvh` as "in-run host member, label fixed by
    the program"); add `ctl_dev` to `HOST_RESERVED`; fix or retire `t207_d5_pairs.py`; run processing and gate on the
    same machine (exact equality depends on the system OpenBLAS); refuse symlinked dumps or re-hash at read time.
- **Round 2 result (2026-11-01, not reviewed):** `gpu-port-207-harness` `f06b35c5` (tree 0350f690). qref_gate.py calls
  the approved gate; stamps in the attempts log; env allow-list refusals shown; bridge to the 207-2 record layout
  (derive_run unchanged); `qref_swap:dev:var8` control; t207_d5_* in TOOLS; FAULT_TAU_AT in the gate's FAULTS table.
  Evidence on tree 1713bb7c: selfcheck 57 QREF cases 0 wrong, gate self-test 493/493, dry run FAIL by placeholder
  bound at 6×6 (var8 −0.146 not killed; not a calibration). **Not yet run on the final tree 0350f690**
  (`nividia:~/alf-gpu/work/t207-5/chain3.sh`, needs ≥ 7.6 GB free on /).
- Manager decision on the spec's instance fields: `binary_sha256` and `prog_sha` (built after the freeze) and `seeds`
  (drawn by the driver after the freeze, VERIFICATION §3) may be `"@harness"`; `parameters_sha256` is known at
  pre-registration and is concrete in the acceptance file.
- Open: a real `acc_ship.sh` + T1 needs the t207_d5 scripts in gpu-port first (BASE vs PREREG); driver QREF disk
  margin (`disk_min_free_GB` + 2) as a per-file parameter.

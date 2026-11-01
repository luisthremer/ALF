# 207-2 – QREF analysis pipeline: no silent data loss, frozen reference list

- **Phase / lane:** 2 (NVIDIA) · part of 207 (decision 34)
- **Owner agent:** fortran-gpu-dev (Sonnet), CPU only
- **Review:** numerics-reviewer (scripts)
- **Depends on:** – (runs in parallel with 207-1; touches no `Prog/` file)
- **Status:** open

## Goal
The offline accuracy analysis (binary128 reference `t207_d5_qref.F90`, aggregation `t207_d5_agg.py`, statistics
`t207_d5_stat.py`, processing `t207_d5_proc.sh`) can never drop bad data silently, and its reference ensemble is an
explicit list (reviewer changes 1, 2 and 5 of the dev5 review, decision 33).

## Context
Review of dev5 (2026-11-01): `ES10.3` prints 3-digit exponents without `E` (`1.798+308`), non-finite errors are mapped
to `huge`, gfortran `MAXVAL` skips NaN, and `NUM.findall` skips `NaN`/`Infinity`/E-less numbers, so a point whose
device result is NaN silently drops out of S_ens. `agg.py` treats every method name not on an exclusion list as a
correct member. `proc.sh` deletes dumps after processing.

## Tasks
1. Offline program: finite checks on every matrix read from a dump (device outputs, B, B⁻¹) and every computed number;
   non-finite → the point is flagged, and the run's verdict becomes INPUT_ERROR/FAIL; format `ES12.4E3` (or wider).
2. Parser: exactly the expected fields per method, every listed member at every point, point count derived from the
   run parameters (not a constant 210); anything else → error with a message.
3. Reference ensemble from an explicit list given as input (default: cpu, nvh, n4, n100004, n1, n100001, n2, n1608,
   n4808, n5, n100005; LOG: cpu, nvh, n4, n100004, n1, n100001, n2, n5, n100005); unknown names → error.
4. Processing keeps the manifest (sha256 of every dump), the full `qref.txt` and every 40th dump; it computes and
   writes the result before deleting anything, and deletes only dumps listed in its own manifest.
5. Mutation self-tests: NaN, Inf, missing member, missing point, shifted field, unknown extra member, device and var8
   swapped; each must give the error or FAIL.
6. Re-run on the kept dev5 dumps: the recorded numbers are reproduced.

## Files allowed
`docs/gpu/scripts/nividia/t207_d5_*`. Work dir `nividia:~/alf-gpu/work/t207-2/` (TMPDIR inside), read dev5 data in
place, never modify it.

## Done when
Self-tests and the dev5 reproduction pass; numerics-reviewer approves the scripts.

## Result (2026-11-01): done
- Gate redesign `305488ec` (spec-driven, recomputed from run directories), fix round 4 `ffab796f` (re-processing of
  kept dumps by the pinned binary128 program, pinned binary/parameters/seeds, main env fixed in code, frozen refs,
  q95/max/by_stratum statistics, `--spec-sha`), fix round 5 `4118a1de` (keep/hostin/min_reprocessed pinned, exact
  equality of re-processed numbers, end-to-end device FAIL, minimal program env, host label from the spec).
- Self-test 485/485 as expected (nividia `~/alf-gpu/work/t207-2/fix5/final/st/results.txt`); dev5 reproduction
  unchanged; re-processing of all 16 dev5 runs: every number equal.
- numerics-reviewer part (c): round 3 changes required (fabricated runs, noisy refs, swapped env, cherry-picking,
  localized defect, spec sha), round 4 changes required (keep unpinned, limits), round 5 **approve** (evidence
  `~/alf-gpu/work/t207-2/review{5,6}/`). Changed test assertions (flagged): one-digit difference now INPUT_ERROR;
  hostin-1 test split; cover-2 test moved to a one-run setup; messages split for control members.
- Stated limits (in every report): only every KEEP-th dump is re-computed; only `ALF_GPU_*` variables are visible;
  the binary↔source pairing and the binding of records to the stamped run are the harness's job (207-5).
- Follow-ups: 207-5 (docs text on host labels, `ctl_dev` reserved label, `t207_d5_pairs.py`, same machine for
  processing and gate, symlinked dumps, record format); 207-3 (first keep > 0 end-to-end run must re-process with all
  numbers equal; multi-stratum by_stratum run); 207-6 (keep and min_reprocessed choice, below).

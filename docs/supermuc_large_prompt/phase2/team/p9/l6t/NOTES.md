# p9/l6t: L6T same-seed test for V7 (decision 43), SuperMUC-NG Phase 2, 2026-11-13

C = measured (job id given), E = estimate. Nothing committed.

## Rule (docs/gpu/DECISIONS.md decision 43, verbatim)
"L6T: 20 new seed sets; V7 passes if the median of its deviation from the CPU is ≤ 2× the median of oneMKL, and its maximum ≤ 3× the maximum of oneMKL.
Mz=F (L6F, L12F) ≤ 1e-13 for every seed. Never on node i20r02c06s02."

## Set-up (written BEFORE any run)
- Source: `git archive HEAD` of this worktree, HEAD 5924614f (contains the QR race fix 7e8b6032), docs/ removed, rsynced to
  `$SCRATCH/large_prompt/team/p9/l6t/src_head` (GIT_HEAD file inside).
- GPU binary: `p9/l6t/build_txt` = recipe `txtnew` of `p5/jobs/build_comb.sh` (configure `IntelLLVM MPI GPU=intel`, JIT, text observables), env `jobs/env_intel_pvc.sh`
  (ifx/oneMKL 2025.3.0), build job 561175.
- CPU reference: existing `$SCRATCH/large_prompt/team/build_cpu-txt/Prog/ALF.out` (ifx 2025.3.0, `IntelLLVM MPI`, 2026-11-11). Checked: `diff -rq` of its Prog/ and
  Libraries/ against src_head differs only in `Prog/gpu_{lowrank,obs,state,update}_mod.F90` (each file entirely inside `#ifdef ALF_GPU`), `Prog/main.F90`
  (all hunks inside `#ifdef ALF_TEAM_TIMERS`, not defined in either build), the SYCL shim / device-QR sources and Libraries/Modules/Makefile (GPU=intel only).
  So the CPU-compiled physics is identical.
- Driver: `run_cmp_sng.py` and `compare.py` (copies, md5 equal to `team/jobs/`), one MPI rank. Each CPU reference is run once per (system, seed set) and every
  GPU variant is compared to it with `compare.py --tol 0`; deviation = its "worst rel dev" (max over all observable files and confout_0 of max|diff| / file scale).
  A missing/shape-mismatched file counts as inf.
- Systems (as p5/g2, g5): L6T = `L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18`; L6F = same with mz=F;
  L12F = `L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72`.
- Variants (all with `LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1`): oneMKL = V2 = `ALFSH_DEVQR=0 ALFSH_DEVUNMQR=0`;
  V7 = `ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1`; Pxx = V7 + `ALFSH_DEVQR_LIB=p9/l6t/plant/Pxx/libalfdevqr.so`.
- Plants (no worktree source edited; copies of src_head's three device-QR files in `plant/Pxx/`, compiled with the Makefile's command for libalfdevqr.so):
  P0 = no edit (control, must equal V7 bitwise); P12 / P11 / P10 = the existing p5 "tau" plant (t_plant.sh) with factor (1 + 1e-12 / 1e-11 / 1e-10):
  the stored tau of the first reflector of every panel, `tau[2*i] = t.r * (i == 0 ? 1.0 + eps : 1.0)` (same for the imaginary part), one line changed.
- Concurrency: the seed sets of a job run in parallel (one tile each, `ZE_AFFINITY_MASK`), `I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0`; runs are single-rank and
  deterministic, so this does not change any number (P0 = V7 and reproduction of the p5/g5 values serve as checks).
- Seed sets: set k = the original seeds file `src/Scripts_and_Parameters_files/Start/seeds` rotated to start at line k (as p5/g5.sh); with one rank only its
  first line is used. Dev sets (calibration): k = 1 (original), 3, 11, 23, 37. Acceptance sets: k = 41, 43, 45, ..., 79 (20 odd values). First seeds:
  dev 790789, 1048, 361162, 472421, 456393; acceptance 46379, 210334, 974820, 274404, 593079, 497340, 472980, 307776, 699773, 708323, 102605, 599781,
  422490, 542867, 417774, 957426, 907106, 759152, 749907, 433810. All 25 distinct (checked on the cluster).
- Median: of 5 values the middle one; of 20 the mean of the 10th and 11th (Python statistics.median). Script `jobs/rule.py`.

## Calibration criterion (declared before running)
On the 5 dev sets, with the L6T part of the rule (median ≤ 2× oneMKL median, max ≤ 3× oneMKL max):
V7 must pass, P0 must equal V7, and P10 (perturbation at the top of the allowed 1e-12..1e-10 range) must FAIL the L6T part (failing only the Mz=F part
does not count: that part is the old 1e-13 check). P11 and P12 are reported as the detection threshold, they do not decide. If P10 does not fail, stop and report.

## Step 1: calibration (job 561177, i20r03c02s04, 22:38:39-22:44:41; build job 561175, i20r03c02s04, 22:36:11-22:37:48) — all C
Deviation (worst rel dev vs the ifx CPU run, same seed); `logs/p9_cal.561177.log`, `logs/res_cal.txt`.

| variant | k=1 | k=3 | k=11 | k=23 | k=37 | median | max | ratio med (≤2) | ratio max (≤3) | L6F max | rule |
|---|---|---|---|---|---|---|---|---|---|---|---|
| oneMKL (V2) L6T | 5.503e-11 | 1.131e-10 | 1.434e-11 | 7.376e-11 | 1.842e-10 | 7.376e-11 | 1.842e-10 | - | - | - | - |
| V7 L6T | 1.358e-10 | 4.510e-11 | 1.926e-11 | 1.008e-10 | 1.480e-10 | 1.008e-10 | 1.480e-10 | 1.37 | 0.80 | 3.11e-15 | **PASS** |
| P0 (unedited copy) | = V7 in every entry (L6T and L6F) | | | | | | | 1.37 | 0.80 | 3.11e-15 | PASS (= V7, libraries md5-equal) |
| P12 (tau 1+1e-12) L6T | 2.197e-10 | 1.136e-10 | 1.313e-10 | 1.274e-09 | 3.512e-10 | 2.197e-10 | 1.274e-09 | 2.98 | 6.92 | 4.35e-13 | **FAIL** (both L6T criteria, and Mz=F) |
| P11 (1+1e-11) | 2.234e-09 | 1.125e-09 | 1.281e-09 | 1.289e-08 | 3.537e-09 | 2.234e-09 | 1.289e-08 | 30.3 | 70.0 | 4.35e-12 | **FAIL** |
| P10 (1+1e-10) | 2.236e-08 | 1.126e-08 | 1.279e-08 | 1.287e-07 | 3.539e-08 | 2.236e-08 | 1.287e-07 | 303 | 699 | 4.35e-11 | **FAIL** |

- Criterion met: V7 passes, P0 = V7, P10 fails the L6T part (so do P11 and P12: the L6T part catches a 1e-12 relative tau error at the dev sets,
  with margin 1.5x on the median and 2.3x on the max). The plant response is linear in eps (x10 per decade).
- Checks: oneMKL and V7 values equal the p5/g5 values (job 560724, other binary, no ASYNCBLK, sequential runs) to all printed digits, and V7 k=1 equals
  the p5 Gate C value 1.358e-10: ASYNCBLK=1, the concurrent runs and the QR race fix change nothing at L6T. Mz=F L6F V7 ≤ 3.1e-15.
- Note for the reader: V7's dev median is 1.37x oneMKL's: V7 is not more accurate than oneMKL at L6T, it is in the same noise band.

## Step 2: acceptance (declared before running, unchanged from above)
Job `jobs/acc.sh`, partition test (not preemptible), 5 workers, seed sets k = 41..79 odd; per set CPU L6T/L6F/L12F, GPU L6T oneMKL and V7, L6F V7, L12F V7.
Judged once from `rule.py` output; no re-run.

## Step 2 result (job 561181, partition test, i20r01c01s06, 22:45:33-22:55:35; 80/80 comparisons, all runs rc=0) — all C
`logs/p9_acc.561181.log`, `logs/res_acc.txt`. Deviation = worst rel dev vs the ifx CPU run, same seed.

| k | oneMKL L6T | V7 L6T | V7 L6F | V7 L12F |
|---|---|---|---|---|
| 41 | 3.822e-11 | 9.314e-11 | 2.788e-15 | 3.022e-14 |
| 43 | 2.861e-11 | 7.802e-11 | 6.107e-15 | 8.888e-15 |
| 45 | 3.233e-11 | 3.476e-11 | 4.974e-15 | 1.392e-14 |
| 47 | 1.224e-10 | 1.178e-10 | 4.863e-15 | 3.770e-15 |
| 49 | 2.031e-10 | 1.026e-10 | 4.176e-15 | 2.158e-14 |
| 51 | 1.163e-09 | 1.411e-09 | 3.202e-15 | 3.737e-15 |
| 53 | 4.193e-11 | 5.300e-11 | **1.771e-13** | 2.356e-15 |
| 55 | 6.860e-11 | 1.081e-10 | 3.312e-15 | 1.108e-14 |
| 57 | 2.609e-11 | 1.803e-11 | 1.957e-15 | 1.424e-14 |
| 59 | 4.477e-11 | 2.866e-11 | 2.501e-15 | 2.432e-14 |
| 61 | 7.914e-11 | 4.460e-11 | 3.158e-15 | 1.135e-14 |
| 63 | 1.003e-10 | 7.507e-11 | 2.494e-15 | 1.855e-14 |
| 65 | 4.374e-11 | 5.090e-11 | 1.019e-14 | 1.134e-14 |
| 67 | 7.069e-11 | 6.948e-11 | 2.994e-15 | 3.989e-14 |
| 69 | 9.052e-11 | 5.761e-11 | 3.422e-15 | 2.464e-14 |
| 71 | 5.396e-11 | 1.389e-10 | 2.059e-15 | 1.218e-14 |
| 73 | 5.379e-11 | 7.827e-11 | 1.039e-14 | 3.657e-15 |
| 75 | 1.478e-10 | 2.952e-10 | 3.246e-15 | 8.022e-15 |
| 77 | 2.303e-10 | 7.682e-11 | 3.423e-15 | 1.839e-14 |
| 79 | 2.205e-11 | 9.942e-11 | 2.205e-15 | 7.783e-15 |
| median | 6.128e-11 | 7.742e-11 | | |
| max | 1.163e-09 | 1.411e-09 | 1.771e-13 | 3.989e-14 |

- L6T part: median ratio 1.26 (limit 2) ok; max ratio 1.21 (limit 3) ok -> **L6T part PASS**.
- Mz=F part: 39 of 40 runs ≤ 1e-13; **L6F k=53 = 1.771e-13 > 1e-13 -> FAIL**. **Verdict under decision 43: V7 FAILS** (judged once, not re-run).
- What the k=53 excess is (read from the existing outputs, `logs/L6F_53_V7_detail.txt`, no new run): it is `Pot_scal` only (bin 2:
  CPU 16.595874328029623, V7 16.595874328024660, |diff| 4.96e-12, i.e. 3.0e-13 of the entry; compare.py divides by the file scale 28.03 -> 1.771e-13).
  Every other file of that run ≤ 1.2e-14, confout_0 identical (0.0), acceptance identical. Precision Green of the runs themselves:
  CPU max 3.59e-13, V7 max 3.91e-13 (Precision Phase 1.63e-13 / 1.89e-13), i.e. the CPU reference's own stabilisation error at this seed is of the
  same size as the excess. The next-largest L6F value over all 25 V7 runs is 1.04e-14 (dev sets ≤ 3.1e-15).
- Not measured: whether oneMKL (V2) also exceeds 1e-13 on L6F k=53 (the rule does not run oneMKL at Mz=F). That would be a diagnostic, not a re-run of the
  acceptance; it needs the owner's/manager's go-ahead.

## Budget
Node-h from the jobs' own start/end lines (single node each): build 561175 97 s, calibration 561177 362 s, acceptance 561181 602 s:
**total 1061 s = 0.29 node-h** (limit 1.5). Never more than 1 job at a time; none on i20r02c06s02.

## Files
Local (this dir): NOTES.md, sync.sh, waitjob.sh, jobs/{build.sh, cal.sh, acc.sh, lib.sh, rule.py, run_cmp_sng.py, compare.py, hdr.txt}, logs/.
Remote `$SCRATCH/large_prompt/team/p9/l6t/`: src_head (git archive HEAD 5924614f), build_txt, plant/{P0,P10,P11,P12}, work_cal, work_acc (all run dirs, cmp.txt per comparison), jobs/.

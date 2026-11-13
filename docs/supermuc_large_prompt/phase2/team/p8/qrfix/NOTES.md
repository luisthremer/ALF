# p8/qrfix: verification of the ex-parity fix (7e8b6032) and delay-sensitised race tests, SuperMUC-NG Phase 2, 2026-11-13

All numbers C (job ID) = measured. Kernel sources untouched: the delay hooks and plants live only in COPIES (`hdr_*/alf_la_sycl_qr2.hpp`, made by `mkhdr.py`
from `git show 7e8b6032^` / `HEAD`), under `#ifdef QR2_TESTDELAY`. Without that macro a copy is the original kernel text (plus the plant edit for hdr_tau/v/ex).
Remote: `$SCRATCH/large_prompt/team/p8/qrfix/` (bin/, out/, build_aottt, build_txtnew, lib_*, src_head = `git archive HEAD` without docs/, GIT_HEAD = 7e8b6032).
Raw logs: `logs/` (qf_*.log job logs, out/ per-case outputs).

## Delay technique
- apply2_t (G>1): sub-group sid==1 spins QR2_DELAY iterations of a dependent FMA chain right after the exchange barrier, before it reads `ex`
  (the other sub-group of its column group runs ahead into the next executed step and writes `ex` before the barrier).
- panel2_t: the sub-group of column jb-1 spins QR2_PDELAY iterations right after the step barrier, before it reads tau / v (G=1: no barrier follows within the step).
- Counts are read from the environment on the host at every launch (no recompile). Delay 0 = same binary without delay.

## Test 1: unmqr2 with tau=0 patterns (qr2_test race_apply; 144 cases per run: nq 129/200/288 (G=2, chunk 7), 289/400/576 (G=4, chunk 3), Left/Right x N/C, k=64, 48 other columns)
Reference: host MKL zunmqr on the same V/tau (a naive host applier agrees with MKL to <= 2.5e-15 in every case). Case fails if any rep deviates > 1e-13 (observed: 0.5..1.4, i.e. O(1)).
Patterns: none (no zeros), alt (zero tau at odd j: runs of 1), z1 (j%3==1: runs of 1), z2 (runs of 2: even count, no collision in either version = control), z3 (runs of 3), rnd (50 %).

| test | header | delay | expected | result (failing cases of 144; by pattern) | job |
|---|---|---|---|---|---|
| race_apply, 20 reps | old (7e8b6032^) plain | no | latent, undetected | 0 | 560961 |
| race_apply, 20 reps | old, hooks, QR2_DELAY=0 | no | same | 0 | 560961 |
| race_apply, 3 reps | old | 300 / 3000 / 30000 | FAIL | 80 / 80 / 80: alt 24, z1 20, z3 12, rnd 24, none 0, z2 0 | 560961 |
| race_apply, 20 reps | new (HEAD) plain | no | PASS | 0 | 560961 |
| race_apply, 20 reps | new, hooks, QR2_DELAY=0 | no | PASS | 0 | 560961 |
| race_apply, 3 reps | new | 300 / 3000 / 30000 | PASS | 0 / 0 / 0 | 560961 |
| race_apply (sensitivity control) | plant race_ex (single slot) | 3000 | FAIL | 120 failing (all patterns incl. none) | 560961 |
| race_apply | plant race_ex | no (20 reps) | (earlier stress: caught after 337k it.) | 0 | 560961 |

The old header fails exactly where theory says: odd runs of zeros inside one chunk; z2 and none pass. 4 of the 24 z1 cases and all G=4 z3 cases pass: with the descending order or the
chunk alignment of G=4 (chunk = 3 steps) no executed step falls on the parity of the previous executed step inside a chunk (consistent with the analysis, not a hole in the test; alt and rnd cover G=4).
Delay 300 iterations (about 2 us) is already enough; the result does not change up to 30000: the window is deterministic. Without delay the race is invisible (0/144 in 20 reps each, old and new).

Gate A of qr2_test (unchanged `correct` mode) with the NEW header (job 560961): clean 0 failed (104 PASS lines); DEVQR_PLANT=1 31 failed, =2 60 failed, =3 56 failed.

## Test 2: panel plants with the delay (qr2_test race_panel: geqrf2 vs host zgeqrf, M = 64..768, nb 16, R/V/tau/orthogonality as Gate A)
| plant / header | ALFDEVQR_G | delay QR2_PDELAY | expected | result (failing of cases) | job |
|---|---|---|---|---|---|
| real (new) | 1 | 0 (70 cases) / 500 / 5000 (21 each) | PASS | 0 / 0 / 0 | 560961 |
| real (new) | default | 5000 | PASS | 0 of 21 | 560961 |
| old header | 1 | 5000 | PASS (panel unchanged) | 0 of 21 | 560961 |
| race_tau (single tau slot) | 1 | 0 (70 cases) | invisible | 0 of 70 | 560961 |
| race_tau | 1 | 500 / 5000 | FAIL | 21 of 21 / 21 of 21 (dR, dV, dtau O(1)) | 560961 |
| race_tau | default | 5000 | FAIL (see caveat) | 21 of 21 | 560961 |
| race_v (single v buffer) | 1 | 0 (70 cases) | invisible | 0 of 70 | 560961 |
| race_v | 1 | 500 / 5000 | FAIL | 21 of 21 / 21 of 21 | 560961 |
| race_v | default | 5000 | FAIL (see caveat) | 21 of 21 | 560961 |
Caveat: with the default panel shape the trailing panels of every factorisation (M-j <= 256) use G=1, so even M=768 contains G=1 panels; "default" does not isolate G=2 panels.

## Stress test with the delay (node, 32 ranks, 4/tile, CCS 4, existing p5/race/t_stress binary read-only; libalfdevqr.so from hdr_tau / hdr_v / hdr_new with -DQR2_TESTDELAY; job 560972, 150 s per case)
| case | BAD lines (self-comparison with iteration 1) | STRESS_HASH (result of iteration 1) vs clean hash | verdict |
|---|---|---|---|
| S1 race_tau N=576 G=1 PDELAY 2000 | 0 in 37053 it. | 7E74..., differs from clean CC9E... (all 32 ranks) | corruption deterministic: invisible to BAD, caught by the hash |
| S2 race_v N=576 G=1 PDELAY 2000 | 9 on 8 ranks (37371 it.) | 459B..., differs | caught by both |
| S3 race_tau N=576 default G PDELAY 2000 | 0 in 37683 it. | 7AF2..., differs from clean C678... | hash only |
| S4 race_tau N=576 G=1, no delay | 0 in 138698 it. | = clean CC9E... | not detected (as in the p5 race test) |
| S5 real N=576 G=1, PDELAY 2000 + DELAY 3000 | 0 in 6081 it. | = clean CC9E... | PASS |
| S6 real N=1024 default, same delays | 0 in 1440 it. | = clean BA95... (identical to p5 real_N1024) | PASS |
Clean hashes are those of the p5 race runs (job 560785/560803 N=1024, 560786 N=576 default G, 560799 N=576 G=1).
FINDING: `t_stress` compares every iteration only with iteration 1, so a race that corrupts the result the SAME way every time (what the delay makes of race_tau) gives 0 BAD. The per-rank
STRESS_HASH, compared with the clean value, does detect it. The stress test needs "hash == known clean hash" as a pass criterion (not only BAD == 0); that criterion plus the delay catches all three plants.

## Test 3: ALF level (job 560965 on i20r02c06s10; new build = git archive HEAD (7e8b6032) built in job 560962; old = p5 build_txtnew / build_aottt)
Only difference of the source trees: alf_la_sycl_qr2.hpp (diff -rq p5/src_alf vs src_head; GIT_HEAD and an unused HDF5/ dir aside); libalfdevqr.so md5 differs, as expected.
V7 = LIBOMPTARGET_LEVEL_ZERO_STAGING_BUFFER_SIZE=1048576 ALFSH_ASYNCBLK=1 ALFSH_DEVQR=2 ALFSH_DEVUNMQR=1, same seeds, CMP_TOL=0 (compare.py) plus byte comparison of every observable file:
| config | compare.py | byte compare | PG (mean, max) old = new |
|---|---|---|---|
| L6F (ndelay 8, 10x2 sweeps) | PASS worst 0.0, acceptance equal | 12 files, 0 differ | 4.7e-15, 6.0e-13 |
| L12F (ndelay 32, 4 sweeps) | PASS 0.0 | 12 files, 0 differ | 7.2e-15, 8.3e-13 |
| L24F (ndelay 32, 2 sweeps, beta 24, nphi 288) | PASS 0.0 | 12 files, 0 differ | 2.29e-14, 3.51e-11 |
Controls: old vs old V7 L12F: identical (the run is deterministic); new V7 vs new with the device QR off (wrapper): worst 5.4e-15, 11 of 12 files differ (the check sees rounding-level changes).
Node Precision Green (32 ranks, 4/tile, CCS 4, OMP 3, V7, 3 sweeps): old build mean 1.960957926459005E-014 max 6.536447432744344E-010; new build identical to all printed digits (and identical to the p5 value of job 560724/g2). PASS (< 1e-8).

## Files
mkhdr.py, hdr_{old,new,tau,v,ex}/ (header copies), qr2_test.cpp (copy of p5 with race_apply / race_panel), q1.sh (tests 1+2 job), q2_build.sh (ALF build), q3_alf.sh (ALF level), q4_stress.sh (stress),
sync.sh, chk.sh, waitjob.sh, logs/. Node-h: about 0.6 (560961 4 min, 560962 2 min, 560965 10 min, 560972 17 min, one node each).

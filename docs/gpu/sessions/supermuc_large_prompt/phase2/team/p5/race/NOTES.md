# Race test of the device QR kernels (ALFSH_DEVQR=2 + ALFSH_DEVUNMQR=1), SuperMUC-NG Phase 2, 2026-11-12

Detector: t_stress.F90 (copy of diag/t_stress.F90 + flush of BAD lines, heartbeat, and a reference-result checksum STRESS_HASH per rank),
linked against build/ = copy of p5/build_aottt (AOT for PVC, sources identical to worktree head 2d5cc615, md5 of qr2.hpp/devqr.cpp/shim.cpp checked).
32 ranks, 4 per tile, ZEX_NUMBER_OF_CCS 4 on every tile, OMP 3, ALFSH_QRNB=16, 1 MB L0 staging buffer (= production variant V7 settings).
Each rank repeats gemm, getrf, qrpp, ormqr, trsm, gemm (NVAR=1 and 2) at N, every iteration compared on the device bitwise with iteration 1.
Scripts: race_build.sh (job), race_lib.sh (runcase), race_plant.sh, race_plant2.sh, race_real.sh, race_amp.sh. Logs in logs/ (job logs with per-case summary).
Remote: $SCRATCH/large_prompt/team/p5/race/ (out/*.out = raw stress output).
Nothing committed. Kernel sources untouched (plants = the existing p5/plant/*/libalfdevqr.so; amplifier = existing runtime knob ALFDEVQR_G).

Results: see the final report of the session (copied below by the lead if wanted); raw numbers are in logs/rc_*.log.

## Results (all numbers C = measured; iterations = sum over the 32 ranks of one node-job case)
Plants (N=1024 unless stated, ALFSH_DEVQR=2 + DEVUNMQR=1, 4 ranks/tile CCS 4):
- onebuf (gross):      CAUGHT at iteration 2 (3.8 s), 17988 BAD lines on 28 of 32 ranks (job 560782, i20r07c02s02); then 31 ranks hung (job cancelled)
- race_ex:             CAUGHT, 4 BAD (all G_NVAR2) in 337443 iterations, first at 127 s (job 560765, i20r08c01s04)
- race_v:              N=1024 NOT caught (342547 it, job 560766); N=576 caught 1 BAD in 851110 it (job 560806); N=576 + ALFDEVQR_G=1 amplifier: 530 BAD in 712908 it, all 32 ranks, first at 2 s
- race_tau:            NOT caught: N=1024 341069 it (560764), N=576 1019146 it (560799), N=576 + G=1 852116 it (560799)
Real kernels (0 BAD everywhere; reference hash identical on all ranks and on both nodes):
- N=1024: 506032 it (i20r08c01s04, 560785) + 343586 it (i20r08c01s06, 560803)
- N=576: 2038862 it (i20r08c01s06, 560786) + 1355161 it (i20r08c01s04, 560804)
- N=576 amplifier G=1: 570331 it (i20r02c02s08, 560807)
Controls in the same jobs (0 BAD): oneMKL N=1024 96664 + 48061, N=576 309267 + 154134; DEVQR=1 N=1024 213179 + 106240, N=576 583074 + 290243.
Node-h used: about 5.4 (from job start/end times).

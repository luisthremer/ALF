# Comparison with the other session (short-prompt session)

Sources read AFTER FINDINGS.md was written: `plan/SUPERMUC_NG_PHASE2.md` (63 lines; it points to `docs/supermuc/` for details, which I did not open), and four
cross-session messages that the other session pushed to me unsolicited while I was working (L=24 timings, VTune, ccs result, budget estimate, a relayed "new
instruction"). **The independence the task wanted is therefore only partial for the speed numbers**: I had seen its 14.6 s / 16.9 s / 1.90 sweeps/s before my own
CCS scan results for 24/32 ranks were read (my scan jobs had already been submitted and ran independently; the 8/16-rank numbers and every L=28/32, tuning, PP and
numerical result were measured by me without its input). I did not read its branch/commits or its remote directories. I modified none of its files.

| item | other session | this session | probably right / how to check |
|---|---|---|---|
| partitions, limits | test <=4 nodes 1 h, general <=96 nodes 2 d, large 97-140, preempt <=40 nodes 7 d | same (`scontrol show partition`) | agree |
| gres / GPU request | none, 1 rank per tile = `--ntasks-per-node=8` | none; whole node; 4 ranks per tile best with CCS (32 per node) | agree on gres; layout differs because of the CCS finding (below) |
| storage | `$SCRATCH=/hppfs/scratch/0C/di38zuw`, `$HOME` 100 GB, `/hppfs/work/pn93la` not provisioned yet | same paths; no `/hppfs/work/pn93la`; I could not read a home quota (`mmlsquota` fails) | agree; its 100 GB comes from the login banner, which I did not look at. Check: `dssusrinfo`/banner |
| budget | pn93la = 0.90 k GPU-h (login banner), expires 2028-11-10 | not seen by me; `sacct` unreachable | its figure is likely right; check the login banner |
| toolchain | `intel-toolkit/2025.1.0` | `intel/2025.3.0 intel-mpi/2021.15.0 intel-mkl/2025.3.0 hdf5/1.14.5-intel25` | with ifx 2025.1/2025.2 the oneMKL backend does not compile (`omp_get_interop_*` missing in `omp_lib`) unless the bind(C) fix is present; its working tree has an uncommitted change in `alf_la_mod.F90` (probably the same fix). Check: `ifx` 2025.1 + `use omp_lib, only: omp_get_interop_int` |
| `ZE_FLAT_DEVICE_HIERARCHY` | "set to FLAT" (task text) | **unset** in my jobs, flat is the default (8 devices from `sycl-ls`) | harmless; set it explicitly in job scripts to be safe |
| `OMP_TARGET_OFFLOAD` | docs: DEFAULT | MANDATORY (a missing device aborts instead of silently running on CPU, tested) | MANDATORY is safer |
| GPU-aware MPI vars (`I_MPI_OFFLOAD*`) | listed from the docs | not needed: no GPU buffers are passed to MPI (only host reductions) | confirm by code (`control_mod.F90` reductions are host arrays) |
| rank -> tile | rank r on tile r (card r//2) verified with xpu-smi | same, verified per tile with `xpu-smi dump -t`; default `rank mod 8` puts 4 ranks on only 2 cards | agree; one rank per card is 1.8x faster for 4 ranks |
| L=24 timings | 14.6 s/sweep/chain (8 ranks), 32 ranks+CCS 16.9 s, 1.90 sweeps/s/node, without CCS 0.53-0.63 | 14.6 s, 16.9 s, 1.896 sweeps/s; 0.53-0.63 | agree to the digit (the exact agreement is expected: same binary, same method, but see the independence note) |
| CCS | `ZEX_NUMBER_OF_CCS=0:4,...,7:4` found | found; scanned 8/16/24/32/48/64 ranks: optimum 32 (4 per tile), 48 and 64 are worse (1.22, 0.90) | agree; I add that 4 per tile = number of CCS is the optimum and that numerics were verified with it on (3.0e-10 Mz=T, 2.2e-15 Mz=F) |
| L=28 / L=32 cost | estimated from A30 ratios 1:1.9:3.5 | measured: 37.5 s (L=28) and 66.8 s (L=32) per sweep per chain at 32 ranks, i.e. ratios 1:2.2:4.0 | mine are measured; its L=28/32 budget (900 node-h) is therefore slightly optimistic before tuning, and ~35-40 % too high after tuning |
| numerical verification | same-seed vs ifx CPU 2e-15 (Mz=F), ctest 53/53 | 2.2e-15..4e-16 (Mz=F), 5.5e-11 / 3e-10 (Mz=T), L=24: 7.3e-14, PP 1.4e-13, ctest 53/53 on GPU and CPU builds | agree |
| VTune | GPU busy 47 %, XVE 93 % idle | 47.4 %, 93.1 % | agree |
| additions only in this session | | host-assisted QR is the bottleneck (t_timing), no rank pinning (+35 % otherwise), Nwrap=20 / N_delay=64 (-35..-40 %), AOT optional, merge of `-mpi-pp`, PP verification, HDF5 shim for the missing internet, production jobscript and generator | |
| push | blocked in its sandbox as in mine | blocked | needs the user |

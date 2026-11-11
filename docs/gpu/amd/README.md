# AMD (MI210, gfx90a) notes

Toolchain: gfortran 14 + GCC amdgcn offload, hipBLAS/hipSOLVER via `alf_la_hip_*` (see configure.sh `GPU=AMD`).

## Faster libgomp GCN plugin (optional, ~3x on the real run)
`plugin-gcn-pool-prof.patch` applies to GCC 14.2.0 `libgomp/plugin/plugin-gcn.c` (after changing the wait state at
`run_kernel` to `HSA_WAIT_STATE_ACTIVE`, which alone gave nothing):
- pools the per-launch kernarg block and completion signal (empty target region 253 us -> 37 us),
- pools small device allocations,
- `GCN_PROF=1` prints per-kernel wall time at exit.
Build only `libgomp-plugin-gcn.la`, put its `.libs` first in `LD_LIBRARY_PATH`.

## Runtime settings
`GCN_NUM_THREADS=4` (default 16): parallel-region/reduction startup on GCN scales with thread count; 4 was best at L=16.

## Measured (L=16, 2 sweeps, 1 process): 66 s CPU stock -> 21 s (plugin + threads + obs kernels)

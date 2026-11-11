# Phase 4 (device-only ZGEQRF in the SYCL shim) - status

## Gate A (2026-11-11): PASSED
- Code: `Libraries/Modules/alf_la_sycl_devqr.hpp` (header, used by `alf_la_sycl_shim.cpp` when `ALFSH_DEVQR=1`, block size `ALFSH_QRNB`, default 32).
  Panel = ONE work-group, one 16-lane sub-group per panel column (ZGEQR2: ZLARFG with dlapy3/safmin/zladiv semantics, right-looking ZLARF,
  one barrier per column, reflector published in SLM), V^H V + ZLARFT recurrence in the same kernel; trailing update ZLARFB = 3 oneMKL device GEMMs on an explicit V.
- Test `devqr_test.cpp` (vs host MKL zgeqrf, N=36/144/576/1024 + non-square + tiny, random / column scales 1e-12..1e12 sorted and shuffled / 1e-150..1e150):
  all strict cases PASS (dR, dV, dtau <= 3e-14, worst 4.5e-14; orth ||Q^H Q-I|| of the same size as host LAPACK, <= 1.6e-15); planted tau(0)*=1+1e-12 -> FAIL (dtau 1e-12, orth 2e-12).
- Type 3 (cond 1e24 > 1/eps): R equal to <= 1e-12 (column-wise), reflectors/tau differ O(1) (rounding noise in a numerically rank-deficient trailing part, as between any two QRs): informational, orth gated.
- Host MKL ZGEQRF leaves the last 1-row reflector as identity (R(k,k) complex, tau(k)=0); reference LAPACK and device oneMKL make it real; devqr does the latter (that row is excluded from the host comparison).
- Harness pitfalls (SYCL_UR_USE_LEVEL_ZERO_V2=1, needed by standalone programs): a device->host copy is complete only ~ms AFTER wait() returns (sanity.cpp: ~1 % of a 1 MB readback still wrong at once, 0 after 20 ms): readbacks sleep 30 ms; hybrid oneMKL zgeqrf gives wrong results there, so no standalone comparison with device oneMKL.
- Standalone timing (V2, not ALF): devqr N=576 / 1024: nb=16 1.1 / 2.0 ms, nb=32 0.54 / 1.0 ms, nb=48 0.36 / 0.69 ms (time ~ 1/nb: per-panel launch overhead). To be confirmed inside ALF (gate B).

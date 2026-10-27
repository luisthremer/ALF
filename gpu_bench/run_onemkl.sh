#!/bin/sh
# Run an onemkl-backend binary (or shim_test) on a machine without Intel GPU: host fallback.
# usage: . ~/alf-gpu/env/intel.sh; ./run_onemkl.sh build/onemkl_dispatch/bench [N-list]
#
# Why the LD_LIBRARY_PATH line: on nividia ifx is 2025.3 but the oneAPI 'latest' environment (MKL 2026.2,
# icpx 2026.2) puts compiler/2026.2/lib first.  Code compiled by ifx 2025.3 calls __tgt_create_interop
# (OpenMP interop) which only the 2025.3 libomptarget exports ("undefined symbol __tgt_create_interop"
# with the 2026.2 one).  libsycl.so.9 (needed by libmkl_sycl_*) exists only in 2026.2 and is still found.
# On a machine with one consistent oneAPI release (SuperMUC) this line is not needed.
IFXLIB=$(dirname $(dirname $(command -v ifx)))/lib
export LD_LIBRARY_PATH=$IFXLIB:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
# Do NOT set ONEAPI_DEVICE_SELECTOR here: the OpenMP offload runtime would then see the OpenCL CPU device as an
# offload device and abort ("getExtensionsInfoForDevice failed ..."). The shim picks sycl::cpu_selector_v itself.
unset ONEAPI_DEVICE_SELECTOR
exec "$@"

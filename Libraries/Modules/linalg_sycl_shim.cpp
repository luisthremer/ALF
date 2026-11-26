// Copyright (C) 2026 The ALF project. GPL v3 or later (see linalg_mod.F90 for the full notice).
//
// SYCL shim of the oneMKL device backend of linalg_mod (GPU=intel, P2.4 step 2): oneMKL LAPACK (USM API) on the
// sycl::queue of the OpenMP offload runtime (handed over by linalg_onemkl_mod through OpenMP interop). Every call waits
// for its queue, so it is ordered with the synchronous OpenMP target regions and dispatched BLAS calls of the Fortran
// side, and its LAPACK info is known when it returns (no info slots). Return value: 0 = success, > 0 = the LAPACK info
// (e.g. a zero pivot of getrf), < 0 = an illegal argument (-i) or -1 for any other error (message on stderr).
// Ported from fork/gpu-port Libraries/Modules/alf_la_sycl_shim.cpp without its device QR, timing table and test hooks;
// getri is new (gpu-port inverted with getrf + getrs).
#include <sycl/sycl.hpp>
#include <oneapi/mkl/lapack.hpp>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <stdexcept>

namespace {
using cd = std::complex<double>;
namespace lp = oneapi::mkl::lapack;

sycl::queue *g_q = nullptr;   // the queue of the OpenMP runtime (not owned)

oneapi::mkl::transpose tr(int t) {  // 0 = N, 1 = T, 2 = C
  return t == 0 ? oneapi::mkl::transpose::nontrans
                : (t == 1 ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans);
}
oneapi::mkl::side sd(int s) { return s == 0 ? oneapi::mkl::side::left : oneapi::mkl::side::right; }

template <class F> int guarded(const char *what, F &&f) {
  try {
    if (!g_q) throw std::runtime_error("no queue (alfsh_init not called)");
    f();
    return 0;
  } catch (lp::exception &e) {
    if (e.info() <= 0)
      std::fprintf(stderr, "alfsh: %s: oneMKL LAPACK exception info=%lld: %s\n", what, (long long)e.info(), e.what());
    return e.info() != 0 ? (int)e.info() : -1;
  } catch (sycl::exception &e) {
    std::fprintf(stderr, "alfsh: %s: SYCL exception: %s\n", what, e.what());
    return -1;
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfsh: %s: exception: %s\n", what, e.what());
    return -1;
  }
}
sycl::queue &Q() { return *g_q; }
}  // namespace

extern "C" {

int alfsh_init(void *q) {
  g_q = static_cast<sycl::queue *>(q);
  return g_q ? 0 : -1;
}

void alfsh_finalize() { g_q = nullptr; }

// device name into buf; *is_l0 = 1 for a Level Zero GPU
int alfsh_device_info(char *buf, int n, int *is_l0) {
  return guarded("device_info", [&] {
    auto d = Q().get_device();
    std::snprintf(buf, (size_t)n, "%s", d.get_info<sycl::info::device::name>().c_str());
    *is_l0 = (d.get_backend() == sycl::backend::ext_oneapi_level_zero && d.is_gpu()) ? 1 : 0;
  });
}

// 1 host USM, 2 device USM of the queue's own device, 3 shared, 4 device USM of another device, 0 unknown, -1 error
int alfsh_ptr_kind(const void *p) {
  int r = -1;
  guarded("ptr_kind", [&] {
    switch (sycl::get_pointer_type(p, Q().get_context())) {
      case sycl::usm::alloc::host: r = 1; break;
      case sycl::usm::alloc::device:
        r = (sycl::get_pointer_device(p, Q().get_context()) == Q().get_device()) ? 2 : 4;
        break;
      case sycl::usm::alloc::shared: r = 3; break;
      default: r = 0;
    }
  });
  return r;
}

int alfsh_wait() {
  return guarded("wait", [&] { Q().wait_and_throw(); });
}

int alfsh_zgeqrf_ws(int64_t m, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("geqrf_ws", [&] { *ws = lp::geqrf_scratchpad_size<cd>(Q(), m, n, lda); });
}
int alfsh_zgeqrf(int64_t m, int64_t n, void *a, int64_t lda, void *tau, void *work, int64_t lwork) {
  return guarded("geqrf", [&] {
    lp::geqrf(Q(), m, n, (cd *)a, lda, (cd *)tau, (cd *)work, lwork).wait_and_throw();
  });
}

int alfsh_zunmqr_ws(int side, int trans, int64_t m, int64_t n, int64_t k, int64_t lda, int64_t ldc, int64_t *ws) {
  return guarded("unmqr_ws", [&] { *ws = lp::unmqr_scratchpad_size<cd>(Q(), sd(side), tr(trans), m, n, k, lda, ldc); });
}
int alfsh_zunmqr(int side, int trans, int64_t m, int64_t n, int64_t k, void *a, int64_t lda, void *tau, void *c,
                 int64_t ldc, void *work, int64_t lwork) {
  return guarded("unmqr", [&] {
    lp::unmqr(Q(), sd(side), tr(trans), m, n, k, (cd *)a, lda, (cd *)tau, (cd *)c, ldc, (cd *)work, lwork)
        .wait_and_throw();
  });
}

int alfsh_zgetrf_ws(int64_t m, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("getrf_ws", [&] { *ws = lp::getrf_scratchpad_size<cd>(Q(), m, n, lda); });
}
int alfsh_zgetrf(int64_t m, int64_t n, void *a, int64_t lda, int64_t *ipiv, void *work, int64_t lwork) {
  return guarded("getrf", [&] { lp::getrf(Q(), m, n, (cd *)a, lda, ipiv, (cd *)work, lwork).wait_and_throw(); });
}

int alfsh_zgetri_ws(int64_t n, int64_t lda, int64_t *ws) {
  return guarded("getri_ws", [&] { *ws = lp::getri_scratchpad_size<cd>(Q(), n, lda); });
}
int alfsh_zgetri(int64_t n, void *a, int64_t lda, int64_t *ipiv, void *work, int64_t lwork) {
  return guarded("getri", [&] { lp::getri(Q(), n, (cd *)a, lda, ipiv, (cd *)work, lwork).wait_and_throw(); });
}

}  // extern "C"

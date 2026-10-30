// alf_la SYCL/oneMKL shim (ticket 302a): the LAPACK routines of alf_la's ONEMKL backend (GPU=intel) that oneMKL's
// Fortran "!$omp dispatch" would run (partly) on the host, called from Fortran through bind(C).
// Derived from the 103 prototype (gpu_bench/src/bench_sycl_shim.cpp).
//
// Conventions
//  * alfsh_init(q): q is the opaque sycl::queue* of the OpenMP offload runtime (OpenMP interop, prefer_type("sycl")).
//    The queue is used by every later call and is NOT owned by the shim (a null q is refused: alf_la always passes the
//    interop queue, so that the shim calls are ordered with the OpenMP target regions).
//  * All matrix/vector arguments are USM pointers (device addresses from !$omp target data use_device_addr).
//  * Each call waits for completion before it returns (host-synchronous): the caller's OpenMP target regions and the
//    dispatched BLAS calls then see the result without any further ordering argument. oneMKL's USM API reports
//    failures by exception; it is translated into the return value (0 = ok, > 0 = LAPACK info, -1 = other error; a negative LAPACK info is passed through unchanged).
//  * `_ws` functions return the scratchpad size in complex<double> elements; the caller owns the scratchpad
//    (device-resident). ipiv is int64_t here (oneMKL USM API); the Fortran side converts from/to its int32 ipiv.
//  * Test hook (testsuite/alf_la.tests, ALF_LA_TEST_BREAK_ORDER): alfsh_zgemm_unordered runs a GEMM on a SEPARATE
//    queue (same context and device) and returns without waiting, i.e. it deliberately breaks the ordering.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <oneapi/mkl/lapack.hpp>
#include <atomic>
#include <chrono>
#include <complex>
#include <thread>
#include <optional>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using cd = std::complex<double>;
namespace lp = oneapi::mkl::lapack;

sycl::queue *g_q = nullptr;             // the queue every call uses
std::unique_ptr<sycl::queue> g_other;   // separate queue of the ordering-break test hook
std::atomic<bool> g_async_err{false};   // set by the async handler of the shim-owned queues
std::atomic<bool> g_release{false};     // test hook: releases the gate that holds back the unordered GEMMs
std::optional<sycl::event> g_gate;      // test hook: host task that every unordered GEMM depends on

oneapi::mkl::transpose tr(int t) {  // Fortran side: 0=N, 1=T, 2=C (the LA_* values)
  return t == 0 ? oneapi::mkl::transpose::nontrans
                : (t == 1 ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans);
}
oneapi::mkl::side sd(int s) { return s == 0 ? oneapi::mkl::side::left : oneapi::mkl::side::right; }

// Asynchronous errors of the shim-owned queue (the test-hook queue): recorded, reported on stderr, and turned into the
// return code -1 of the routine that waited (no handler would mean std::terminate). The interop queue of the OpenMP
// runtime keeps the runtime's own handling.
void async_handler(sycl::exception_list l) {
  for (auto &p : l) {
    g_async_err.store(true);
    try {
      std::rethrow_exception(p);
    } catch (std::exception &e) {
      std::fprintf(stderr, "alfsh: asynchronous SYCL error: %s\n", e.what());
    } catch (...) {
      std::fprintf(stderr, "alfsh: asynchronous error of unknown type\n");
    }
  }
}

template <class F> int guarded(const char *what, F &&f) {
  g_async_err.store(false);
  try {
    f();
    return g_async_err.load() ? -1 : 0;
  } catch (lp::exception &e) {
    // positive info = LAPACK computation error (e.g. exactly singular U in getrf); the data are still complete
    if (e.info() <= 0)
      std::fprintf(stderr, "alfsh: %s: oneMKL LAPACK exception info=%lld: %s\n", what, (long long)e.info(), e.what());
    return e.info() != 0 ? (int)e.info() : -1;   // the oneMKL info value is passed through unchanged (<0: argument -info illegal)
  } catch (sycl::exception &e) {
    std::fprintf(stderr, "alfsh: %s: SYCL exception: %s\n", what, e.what());
    return -1;
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfsh: %s: exception: %s\n", what, e.what());
    return -1;
  } catch (...) {
    std::fprintf(stderr, "alfsh: %s: exception of unknown type\n", what);
    return -1;
  }
}
sycl::queue &Q() { return *g_q; }
}  // namespace

extern "C" {

// Select the queue (see header). Returns 0, or -1 if no queue could be made.
int alfsh_init(void *q) {
  return guarded("init", [&] {
    if (!q) throw std::runtime_error("alfsh_init: no queue");
    g_q = static_cast<sycl::queue *>(q);
  });
}

// Opens the gate of the ordering-break test hook (idempotent, async-signal-safe enough: one atomic store). Called by
// alfsh_wait, alfsh_finalize, the Fortran failure path (om_fail) and at process exit, so that an error exit while the gate
// is closed cannot leave the host task spinning (the runtime shutdown would wait for it).
void alfsh_release() { g_release.store(true); }

void alfsh_finalize() {
  g_release.store(true);
  g_gate.reset();
  g_other.reset();     // the queue destructor waits for its (now released) tasks
  g_q = nullptr;
}

// Name of the device behind the queue (NUL-terminated into buf[n]); *is_l0 = 1 if its backend is Level Zero.
int alfsh_device_info(char *buf, int n, int *is_l0) {
  return guarded("device_info", [&] {
    auto &qq = Q();
    auto d = qq.get_device();
    std::snprintf(buf, (size_t)n, "%s", d.get_info<sycl::info::device::name>().c_str());
    *is_l0 = (d.get_backend() == sycl::backend::ext_oneapi_level_zero && d.is_gpu()) ? 1 : 0;
  });
}

// USM kind of p in the shim queue's context: 0 unknown/host-malloc, 1 host, 2 device USM of THE QUEUE'S DEVICE,
// 3 shared, 4 device USM of another device of the same context (e.g. another PVC tile); -1 on error.
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

// Wait for everything enqueued on the shim queue (and on the test-hook queue).
int alfsh_wait() {
  return guarded("wait", [&] {
    if (g_gate) g_release.store(true);     // test hook: release the held-back unordered GEMMs, then wait for them
    Q().wait_and_throw();
    if (g_other) g_other->wait_and_throw();    // its async_handler records errors (g_async_err)
    g_gate.reset();
  });
}

// ---- QR family -------------------------------------------------------------------------------------------
int alfsh_zgeqrf_ws(int64_t m, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("geqrf_ws", [&] { *ws = lp::geqrf_scratchpad_size<cd>(Q(), m, n, lda); });
}
int alfsh_zgeqrf(int64_t m, int64_t n, void *a, int64_t lda, void *tau, void *work, int64_t lwork) {
  return guarded("geqrf", [&] { lp::geqrf(Q(), m, n, (cd *)a, lda, (cd *)tau, (cd *)work, lwork).wait_and_throw(); });
}
int alfsh_zungqr_ws(int64_t m, int64_t n, int64_t k, int64_t lda, int64_t *ws) {
  return guarded("ungqr_ws", [&] { *ws = lp::ungqr_scratchpad_size<cd>(Q(), m, n, k, lda); });
}
int alfsh_zungqr(int64_t m, int64_t n, int64_t k, void *a, int64_t lda, void *tau, void *work, int64_t lwork) {
  return guarded("ungqr", [&] {
    lp::ungqr(Q(), m, n, k, (cd *)a, lda, (const cd *)tau, (cd *)work, lwork).wait_and_throw();
  });
}
int alfsh_zunmqr_ws(int side, int trans, int64_t m, int64_t n, int64_t k, int64_t lda, int64_t ldc, int64_t *ws) {
  return guarded("unmqr_ws", [&] { *ws = lp::unmqr_scratchpad_size<cd>(Q(), sd(side), tr(trans), m, n, k, lda, ldc); });
}
int alfsh_zunmqr(int side, int trans, int64_t m, int64_t n, int64_t k, void *a, int64_t lda, void *tau, void *c,
                 int64_t ldc, void *work, int64_t lwork) {
  return guarded("unmqr", [&] {
    lp::unmqr(Q(), sd(side), tr(trans), m, n, k, (const cd *)a, lda, (const cd *)tau, (cd *)c, ldc, (cd *)work, lwork)
        .wait_and_throw();
  });
}

// ---- LU family ---------------------------------------------------------------------------------------------
int alfsh_zgetrf_ws(int64_t m, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("getrf_ws", [&] { *ws = lp::getrf_scratchpad_size<cd>(Q(), m, n, lda); });
}
int alfsh_zgetrf(int64_t m, int64_t n, void *a, int64_t lda, int64_t *ipiv, void *work, int64_t lwork) {
  return guarded("getrf", [&] { lp::getrf(Q(), m, n, (cd *)a, lda, ipiv, (cd *)work, lwork).wait_and_throw(); });
}
int alfsh_zgetrs_ws(int trans, int64_t n, int64_t nrhs, int64_t lda, int64_t ldb, int64_t *ws) {
  return guarded("getrs_ws", [&] { *ws = lp::getrs_scratchpad_size<cd>(Q(), tr(trans), n, nrhs, lda, ldb); });
}
int alfsh_zgetrs(int trans, int64_t n, int64_t nrhs, void *a, int64_t lda, int64_t *ipiv, void *b, int64_t ldb,
                 void *work, int64_t lwork) {
  return guarded("getrs", [&] {
    lp::getrs(Q(), tr(trans), n, nrhs, (cd *)a, lda, ipiv, (cd *)b, ldb, (cd *)work, lwork).wait_and_throw();
  });
}

// ---- test hook: a GEMM that is NOT ordered against the caller's regions ------------------------------------------
// It goes to a SEPARATE queue (same context and device) and the call returns without waiting. Two in-order queues of
// one device may execute in submission order on the same hardware engine (observed on the B580: 0 of 5 t_order runs
// failed), which would hide the missing ordering, so the GEMMs also depend on a host task (the gate) that is released
// only by alfsh_wait (la_sync / la_finalize): they run LATE, after the caller's following target regions.
int alfsh_zgemm_unordered(int ta, int tb, int64_t m, int64_t n, int64_t k, double ar, double ai, void *a, int64_t lda,
                          void *b, int64_t ldb, double br, double bi, void *c, int64_t ldc) {
  return guarded("zgemm_unordered", [&] {
    if (!g_other) g_other.reset(new sycl::queue{Q().get_context(), Q().get_device(), async_handler});   // separate, out-of-order
    if (!g_gate) {
      static bool at_exit = false;
      if (!at_exit) {     // an error exit with the gate closed must not hang the runtime shutdown
        std::atexit([] { g_release.store(true); });
        at_exit = true;
      }
      g_release.store(false);
      g_gate = g_other->submit([](sycl::handler &h) {
        h.host_task([] {
          while (!g_release.load()) std::this_thread::sleep_for(std::chrono::microseconds(50));
        });
      });
    }
    oneapi::mkl::blas::column_major::gemm(*g_other, tr(ta), tr(tb), m, n, k, cd(ar, ai), (const cd *)a, lda,
                                          (const cd *)b, ldb, cd(br, bi), (cd *)c, ldc, {*g_gate});   // no wait
  });
}

}  // extern "C"

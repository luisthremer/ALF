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
#include <dlfcn.h>
#include <unistd.h>
#include <atomic>
#include <chrono>
#include <complex>
#include <thread>
#include <optional>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
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

// Profiling aid (team session 2026-11-11): ALFSH_TIMING=1 accumulates call count and wall time per shim routine and prints
// them at alfsh_finalize. Off by default (one getenv at the first call); no effect on results.
struct TimAcc { long long n = 0; double s = 0.0; };
std::map<std::string, TimAcc> g_tim;
int g_tim_on = -1;
template <class F> int guarded(const char *what, F &&f) {
  g_async_err.store(false);
  try {
    if (g_tim_on < 0) { const char *e = std::getenv("ALFSH_TIMING"); g_tim_on = (e && e[0] == '1') ? 1 : 0; }
    if (g_tim_on == 1) {
      auto t0 = std::chrono::steady_clock::now();
      f();
      auto &a = g_tim[what];
      a.n++; a.s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    } else {
      f();
    }
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

// Device-only Householder QR / ZUNMQR (alf_la_sycl_devqr.hpp, alf_la_sycl_qr2.hpp; in libalfdevqr.so because they contain SYCL device code, which a plain
// host object cannot carry) instead of oneMKL's hybrid ZGEQRF (panels on the host CPU) and ZUNMQR. Run-time switches (read once; default off = oneMKL,
// the library is then never loaded):
//   ALFSH_DEVQR=1  phase-4 ZGEQRF (panel + ZLARFT in one kernel, ZLARFB by 3 GEMMs)     ALFSH_DEVQR=2  phase-5 ZGEQRF (T-free panel + apply kernel)
//   ALFSH_DEVUNMQR=1  phase-5 device ZUNMQR (apply kernel; side L/R, trans N/C)
// block size ALFSH_QRNB (default 32; 16 for level 2 unless set), library path ALFSH_DEVQR_LIB (default: ../Libraries/Modules/libalfdevqr.so relative to
// the directory of the executable, else libalfdevqr.so on the loader path).
int g_devqr = -1, g_qrnb = 32, g_devunmqr = 0;
bool g_qrnb_set = false;
int64_t (*g_dq_scratch)(int64_t, int64_t, int) = nullptr;
int (*g_dq_geqrf)(void *, int64_t, int64_t, void *, int64_t, void *, void *, int64_t, int) = nullptr;
int (*g_dq_geqrf2)(void *, int64_t, int64_t, void *, int64_t, void *, int) = nullptr;
int (*g_dq_unmqr)(void *, int, int, int64_t, int64_t, int64_t, const void *, int64_t, const void *, void *, int64_t) = nullptr;
void (*g_dq_bench)(void *, void *) = nullptr;
void devqr_init() {
  if (g_devqr >= 0) return;
  const char *e = std::getenv("ALFSH_DEVQR");
  g_devqr = (e && e[0] == '1') ? 1 : ((e && e[0] == '2') ? 2 : 0);
  const char *u = std::getenv("ALFSH_DEVUNMQR");
  g_devunmqr = (u && u[0] == '1') ? 1 : 0;
  const char *b = std::getenv("ALFSH_QRNB");
  if (b && std::atoi(b) > 0) { g_qrnb = std::atoi(b); g_qrnb_set = true; }
  if (g_devqr == 2 && !g_qrnb_set) g_qrnb = 16;
  if (g_devqr > 0 || g_devunmqr) {
    std::string path;
    if (const char *l = std::getenv("ALFSH_DEVQR_LIB")) path = l;
    else {
      char buf[4096];
      ssize_t k = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
      if (k > 0) { buf[k] = 0; std::string d(buf); d = d.substr(0, d.rfind('/')); path = d + "/../Libraries/Modules/libalfdevqr.so"; }
      if (path.empty() || access(path.c_str(), R_OK) != 0) path = "libalfdevqr.so";
    }
    void *h = dlopen(path.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (h) {
      g_dq_scratch = (int64_t(*)(int64_t, int64_t, int))dlsym(h, "alfdevqr_scratch");
      g_dq_geqrf = (int (*)(void *, int64_t, int64_t, void *, int64_t, void *, void *, int64_t, int))dlsym(h, "alfdevqr_zgeqrf");
      g_dq_geqrf2 = (int (*)(void *, int64_t, int64_t, void *, int64_t, void *, int))dlsym(h, "alfdevqr_zgeqrf2");
      g_dq_unmqr = (int (*)(void *, int, int, int64_t, int64_t, int64_t, const void *, int64_t, const void *, void *, int64_t))dlsym(h, "alfdevqr_zunmqr");
      g_dq_bench = (void (*)(void *, void *))dlsym(h, "alfdevqr_launch_bench");
    }
    if (!g_dq_scratch || !g_dq_geqrf || (g_devqr == 2 && !g_dq_geqrf2) || (g_devunmqr && !g_dq_unmqr))
      throw std::runtime_error(std::string("ALFSH_DEVQR/ALFSH_DEVUNMQR set but cannot load the device QR library '") + path + "': " + (dlerror() ? dlerror() : "symbol missing"));
  }
}
bool devqr_on() { devqr_init(); return g_devqr > 0; }
bool devunmqr_on() { devqr_init(); return g_devunmqr == 1; }
}  // namespace

extern "C" {

// Select the queue (see header). Returns 0, or -1 if no queue could be made.
int alfsh_init(void *q) {
  return guarded("init", [&] {
    if (!q) throw std::runtime_error("alfsh_init: no queue");
    g_q = static_cast<sycl::queue *>(q);
  });
}

// The queue handed to alfsh_init (for the asynchronous block path, libalfsb: same in-order queue). nullptr before alfsh_init.
void *alfsh_get_queue() { return static_cast<void *>(g_q); }

// Opens the gate of the ordering-break test hook (idempotent, async-signal-safe enough: one atomic store). Called by
// alfsh_wait, alfsh_finalize, the Fortran failure path (om_fail) and at process exit, so that an error exit while the gate
// is closed cannot leave the host task spinning (the runtime shutdown would wait for it).
void alfsh_release() { g_release.store(true); }

void alfsh_finalize() {
  if (g_tim_on == 1)
    for (auto &kv : g_tim)
      std::fprintf(stderr, "ALFSH_TIMING %-12s calls %10lld  total %10.3f s  per call %9.3f ms\n", kv.first.c_str(), kv.second.n,
                   kv.second.s, kv.second.n ? 1e3 * kv.second.s / kv.second.n : 0.0);
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
  return guarded("geqrf_ws", [&] {
    *ws = lp::geqrf_scratchpad_size<cd>(Q(), m, n, lda);
    if (devqr_on()) *ws = std::max<int64_t>(*ws, g_dq_scratch(m, n, g_qrnb));   // the device QR uses the caller's scratch too
  });
}
int alfsh_zgeqrf(int64_t m, int64_t n, void *a, int64_t lda, void *tau, void *work, int64_t lwork) {
  return guarded(devqr_on() ? "geqrf_dev" : "geqrf", [&] {
    if (devqr_on()) {
      static bool lb = false;   // diagnostic: launch cost on this queue (once)
      if (!lb && (lb = true) && g_dq_bench && std::getenv("ALFSH_LAUNCHBENCH")) g_dq_bench(&Q(), work);
      int rc = g_devqr == 2 ? g_dq_geqrf2(&Q(), m, n, a, lda, tau, g_qrnb) : g_dq_geqrf(&Q(), m, n, a, lda, tau, work, lwork, g_qrnb);
      if (rc == 0) { Q().wait_and_throw(); return; }
      if (rc < -1) throw std::runtime_error("device QR failed");
    }
    lp::geqrf(Q(), m, n, (cd *)a, lda, (cd *)tau, (cd *)work, lwork).wait_and_throw();   // oneMKL (default, or the device QR declined)
  });
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
  return guarded(devunmqr_on() ? "unmqr_dev" : "unmqr", [&] {
    if (devunmqr_on()) {
      int rc = g_dq_unmqr(&Q(), side, trans, m, n, k, a, lda, tau, c, ldc);
      if (rc == 0) { Q().wait_and_throw(); return; }
      if (rc < -1) throw std::runtime_error("device ZUNMQR failed");
    }
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

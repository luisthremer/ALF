// SYCL/oneMKL shim (ticket 103): the LAPACK routines that oneMKL's Fortran "!$omp dispatch" does not run
// fully on the device (R3 / HARD.md item 2), called from Fortran through bind(C).
//
// Conventions
//  * Every function takes the SYCL queue as an opaque pointer (sycl::queue*).  The Fortran side gets it from
//    OpenMP interop (omp_get_interop_ptr(obj, omp_ipr_targetsync)).  A null pointer selects a private
//    fallback queue on the SYCL *CPU* device (sycl::cpu_selector_v; ALFSH_FALLBACK=default uses the default
//    selector instead): only meant for host fallback, where the "device pointers" are host pointers.
//    ONEAPI_DEVICE_SELECTOR must NOT be set for the process: the OpenMP offload runtime honours it too and
//    then tries to use the OpenCL CPU device as an offload device (fails/hangs on nividia).
//  * All matrix/vector arguments are USM pointers (device addresses obtained with use_device_addr).
//  * No device write of `info`: oneMKL's USM API reports errors by exception, which is translated into the
//    return value (0 = ok, >0 = LAPACK info, -1 = other error).
//  * `_ws` functions return the scratchpad size in complex<double> elements; the caller owns the scratchpad
//    (device-resident, same lifetime rules as the LAPACK `work` array).
//  * Synchronous: each call waits for completion before returning (matches the timing/sync model of the
//    benchmark).  ALF would use the returned event with an in-order queue instead.
#include <sycl/sycl.hpp>
#include <oneapi/mkl/lapack.hpp>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

namespace {
using cd = std::complex<double>;
namespace lp = oneapi::mkl::lapack;

sycl::queue &fallback_queue() {
  static sycl::queue q = [] {
    const char *e = std::getenv("ALFSH_FALLBACK");
    if (e && std::string(e) == "default")
      return sycl::queue{sycl::default_selector_v, sycl::property::queue::in_order()};
    return sycl::queue{sycl::cpu_selector_v, sycl::property::queue::in_order()};
  }();
  return q;
}
sycl::queue &Q(void *q) { return q ? *static_cast<sycl::queue *>(q) : fallback_queue(); }

oneapi::mkl::transpose tr(int t) {  // Fortran side uses the cuBLAS/LA_* values 0=N, 1=T, 2=C
  return t == 0 ? oneapi::mkl::transpose::nontrans
                : (t == 1 ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans);
}
oneapi::mkl::side sd(int s) { return s == 0 ? oneapi::mkl::side::left : oneapi::mkl::side::right; }

template <class F> int guarded(const char *what, F &&f) {
  try {
    f();
    return 0;
  } catch (lp::exception &e) {
    std::fprintf(stderr, "alfsh: %s: oneMKL LAPACK exception info=%lld detail=%lld: %s\n", what,
                 (long long)e.info(), (long long)e.detail(), e.what());
    return e.info() > 0 ? (int)e.info() : -1;
  } catch (sycl::exception &e) {
    std::fprintf(stderr, "alfsh: %s: SYCL exception: %s\n", what, e.what());
    return -1;
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfsh: %s: exception: %s\n", what, e.what());
    return -1;
  }
}
}  // namespace

extern "C" {

// Print the device behind a queue (null = fallback queue).  Returns 0.
int alfsh_describe(void *q) {
  return guarded("describe", [&] {
    auto &qq = Q(q);
    auto d = qq.get_device();
    std::printf("alfsh: queue %s, device '%s' (%s), in_order=%d\n", q ? "from OpenMP interop" : "private fallback",
                d.get_info<sycl::info::device::name>().c_str(),
                d.is_gpu() ? "gpu" : (d.is_cpu() ? "cpu" : "other"), (int)qq.is_in_order());
  });
}

// ---- QR family ------------------------------------------------------------------------------------------
int alfsh_zgeqrf_ws(void *q, int64_t m, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("geqrf_ws", [&] { *ws = lp::geqrf_scratchpad_size<cd>(Q(q), m, n, lda); });
}
int alfsh_zgeqrf(void *q, int64_t m, int64_t n, void *a, int64_t lda, void *tau, void *work, int64_t lwork) {
  return guarded("geqrf", [&] {
    lp::geqrf(Q(q), m, n, (cd *)a, lda, (cd *)tau, (cd *)work, lwork).wait_and_throw();
  });
}
int alfsh_zungqr_ws(void *q, int64_t m, int64_t n, int64_t k, int64_t lda, int64_t *ws) {
  return guarded("ungqr_ws", [&] { *ws = lp::ungqr_scratchpad_size<cd>(Q(q), m, n, k, lda); });
}
int alfsh_zungqr(void *q, int64_t m, int64_t n, int64_t k, void *a, int64_t lda, void *tau, void *work,
                 int64_t lwork) {
  return guarded("ungqr", [&] {
    lp::ungqr(Q(q), m, n, k, (cd *)a, lda, (const cd *)tau, (cd *)work, lwork).wait_and_throw();
  });
}
int alfsh_zunmqr_ws(void *q, int side, int trans, int64_t m, int64_t n, int64_t k, int64_t lda, int64_t ldc,
                    int64_t *ws) {
  return guarded("unmqr_ws", [&] {
    *ws = lp::unmqr_scratchpad_size<cd>(Q(q), sd(side), tr(trans), m, n, k, lda, ldc);
  });
}
int alfsh_zunmqr(void *q, int side, int trans, int64_t m, int64_t n, int64_t k, void *a, int64_t lda, void *tau,
                 void *c, int64_t ldc, void *work, int64_t lwork) {
  return guarded("unmqr", [&] {
    lp::unmqr(Q(q), sd(side), tr(trans), m, n, k, (const cd *)a, lda, (const cd *)tau, (cd *)c, ldc, (cd *)work,
              lwork)
        .wait_and_throw();
  });
}

// ---- LU family (not used by the benchmark sections; exercised by shim_test) ---------------------------------
// ipiv is int64_t (oneMKL USM API), unlike LAPACK's int32 ipiv.
int alfsh_zgetrf_ws(void *q, int64_t m, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("getrf_ws", [&] { *ws = lp::getrf_scratchpad_size<cd>(Q(q), m, n, lda); });
}
int alfsh_zgetrf(void *q, int64_t m, int64_t n, void *a, int64_t lda, int64_t *ipiv, void *work, int64_t lwork) {
  return guarded("getrf", [&] {
    lp::getrf(Q(q), m, n, (cd *)a, lda, ipiv, (cd *)work, lwork).wait_and_throw();
  });
}
int alfsh_zgetri_ws(void *q, int64_t n, int64_t lda, int64_t *ws) {
  return guarded("getri_ws", [&] { *ws = lp::getri_scratchpad_size<cd>(Q(q), n, lda); });
}
int alfsh_zgetri(void *q, int64_t n, void *a, int64_t lda, const int64_t *ipiv, void *work, int64_t lwork) {
  return guarded("getri", [&] {
    lp::getri(Q(q), n, (cd *)a, lda, ipiv, (cd *)work, lwork).wait_and_throw();
  });
}

}  // extern "C"

// libalfdevqr.so: the device-only QR of alf_la_sycl_devqr.hpp as a SYCL shared library. The SYCL device code needs a SYCL-aware
// link (icpx -fsycl); the main shim object stays a plain host object that ifx can link, and loads this library by dlopen
// when ALFSH_DEVQR=1 (see alf_la_sycl_shim.cpp). Build: see Libraries/Modules/Makefile (GPU=intel).
#include "alf_la_sycl_devqr.hpp"
#include "alf_la_sycl_qr2.hpp"   // phase 5: T-free panel + apply kernels (ZGEQRF level 2, ZUNMQR)
#include <chrono>
#include <cstdio>

extern "C" {
// scratch size (complex<double> elements) the caller must provide for an m x n matrix and block size nb
int64_t alfdevqr_scratch(int64_t m, int64_t n, int nb) {
  if (nb < 1) nb = 1;
  if (nb > alfdevqr::MAX_NB) nb = alfdevqr::MAX_NB;
  return alfdevqr::scratch_size(m, n, nb);
}
// Diagnostic (ALFSH_LAUNCHBENCH=1, first geqrf call): cost per submission of n in-order kernels / GEMMs on the shim queue q (one wait at the end)
void alfdevqr_launch_bench(void *qq, void *scratch /* >= 64 complex<double> */) {
  sycl::queue &q = *static_cast<sycl::queue *>(qq);
  double *d = (double *)scratch;
  auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  for (int rep = 0; rep < 2; ++rep) {
    for (int n : {1, 20, 100}) {
      double t0 = now();
      for (int i = 0; i < n; ++i) q.parallel_for(sycl::nd_range<1>(256, 256), [=](sycl::nd_item<1> it) { if (it.get_global_id(0) == 0) d[0] += 1.0; });
      q.wait();
      double t1 = now();
      for (int i = 0; i < n; ++i) oneapi::mkl::blas::column_major::gemm(q, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans, 4, 4, 4, alfdevqr::cd(1.0), (alfdevqr::cd *)d, 4, (alfdevqr::cd *)d + 16, 4, alfdevqr::cd(0.0), (alfdevqr::cd *)d + 32, 4);
      q.wait();
      double t2 = now();
      if (rep) std::fprintf(stderr, "ALFSH_LAUNCHBENCH n=%3d: %8.1f us per trivial kernel   %8.1f us per tiny GEMM (incl. one wait)\n", n, 1e6 * (t1 - t0) / n, 1e6 * (t2 - t1) / n);
    }
  }
}
// Diagnostic (ALFSH_LAUNCHBENCH=1): time of one panel kernel (first panel, M = m) on a COPY of the real matrix
void alfdevqr_probe(void *qq, void *A, int64_t lda, int64_t mfull, int nb) {
  sycl::queue &q = *static_cast<sycl::queue *>(qq);
  auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  {
    auto d = q.get_device();
    std::fprintf(stderr, "ALFSH_DEVINFO max_wg=%zu local_mem=%zu max_cu=%u max_clock=%u sub_group_sizes:", d.get_info<sycl::info::device::max_work_group_size>(), (size_t)d.get_info<sycl::info::device::local_mem_size>(),
                 d.get_info<sycl::info::device::max_compute_units>(), d.get_info<sycl::info::device::max_clock_frequency>());
    for (auto z : d.get_info<sycl::info::device::sub_group_sizes>()) std::fprintf(stderr, " %zu", (size_t)z);
    std::fprintf(stderr, "  name=%s\n", d.get_info<sycl::info::device::name>().c_str());
  }
  double *P0 = (double *)A;
  for (int64_t m : {128, 256, 384, 512, 576}) {
    if (m > mfull) continue;
    for (int jb : {16, 32}) {
      double *P = sycl::malloc_device<double>(2 * (size_t)m * jb, q), *tau = sycl::malloc_device<double>(2 * (size_t)jb, q), *Vx = sycl::malloc_device<double>(2 * (size_t)m * jb, q), *Tm = sycl::malloc_device<double>(2 * (size_t)jb * jb, q);
      double best = 1e9, first = 0;
      for (int rep = 0; rep < 6; ++rep) {
        for (int c = 0; c < jb; ++c) q.memcpy(P + 2 * (size_t)c * m, P0 + 2 * (size_t)c * lda, 2 * (size_t)m * 8);
        q.wait();
        double t0 = now();
        alfdevqr::panel_kernel(q, (int)m, jb, P, m, tau, Vx, Tm, true);
        q.wait();
        double t = now() - t0;
        if (rep == 0) first = t; else best = std::min(best, t);
      }
      std::fprintf(stderr, "ALFSH_PROBE panel M=%lld jb=%d: first call %.1f us, then %.1f us  (%.2f us per column)\n", (long long)m, jb, 1e6 * first, 1e6 * best, 1e6 * best / jb);
      sycl::free(P, q); sycl::free(tau, q); sycl::free(Vx, q); sycl::free(Tm, q);
    }
  }
}
// 0 = done (asynchronously on the queue: the caller waits), -1 = declined before any change (m too large, scratch too small),
// -2 = SYCL error (A may be partially modified)
int alfdevqr_zgeqrf(void *q, int64_t m, int64_t n, void *a, int64_t lda, void *tau, void *work, int64_t lwork, int nb) {
  try {
    static bool once = false;
    if (!once && m == 576 && std::getenv("ALFSH_LAUNCHBENCH")) { once = true; alfdevqr_probe(q, a, lda, m, nb); }
    return alfdevqr::geqrf(*static_cast<sycl::queue *>(q), m, n, (alfdevqr::cd *)a, lda, (alfdevqr::cd *)tau, (alfdevqr::cd *)work, lwork, nb, nullptr);
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfdevqr: exception: %s\n", e.what());
    return -2;
  }
}
// Phase 5. ZGEQRF with the T-free panel and the apply kernel (alf_la_sycl_qr2.hpp): no scratch. Same return codes as alfdevqr_zgeqrf.
int alfdevqr_zgeqrf2(void *q, int64_t m, int64_t n, void *a, int64_t lda, void *tau, int nb) {
  try {
    return alfdevqr::geqrf2(*static_cast<sycl::queue *>(q), m, n, (alfdevqr::cd *)a, lda, (alfdevqr::cd *)tau, nb);
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfdevqr: exception: %s\n", e.what());
    return -2;
  }
}
// Phase 5. ZUNMQR (side 0 = L, 1 = R; trans 0 = N, 2 = C). 0 = done (asynchronously on the queue), -1 = declined before any change, -2 = SYCL error.
int alfdevqr_zunmqr(void *q, int side, int trans, int64_t m, int64_t n, int64_t k, const void *a, int64_t lda, const void *tau, void *c, int64_t ldc) {
  try {
    return alfdevqr::unmqr2(*static_cast<sycl::queue *>(q), side, trans, m, n, k, (const alfdevqr::cd *)a, lda, (const alfdevqr::cd *)tau, (alfdevqr::cd *)c, ldc);
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfdevqr: exception: %s\n", e.what());
    return -2;
  }
}
}

// Team session microbenchmark: wall time per call of the oneMKL SYCL LAPACK/BLAS routines ALF's ONEMKL backend uses,
// complex double, square N x N, device USM on one tile (as alf_la_sycl_shim.cpp calls them: one call + wait).
// Also host LAPACK zgeqrf on the same matrix for comparison. Build: icpx -fsycl -qmkl -O2 lapack_bench.cpp -o lapack_bench
// usage: lapack_bench N [reps]
#include <sycl/sycl.hpp>
#include <oneapi/mkl.hpp>
#include <mkl.h>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
using cd = std::complex<double>;
namespace lp = oneapi::mkl::lapack;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc, char **argv) {
  int64_t n = argc > 1 ? atoll(argv[1]) : 576;
  int reps = argc > 2 ? atoi(argv[2]) : 20;
  sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order()};
  std::printf("device %s  N=%lld reps=%d omp_threads=%d\n", q.get_device().get_info<sycl::info::device::name>().c_str(),
              (long long)n, reps, mkl_get_max_threads());
  std::vector<cd> h(n * n);
  std::mt19937_64 g(1);
  std::uniform_real_distribution<double> u(-1, 1);
  for (auto &x : h) x = cd(u(g), u(g));
  cd *A = sycl::malloc_device<cd>(n * n, q), *A0 = sycl::malloc_device<cd>(n * n, q), *C = sycl::malloc_device<cd>(n * n, q);
  cd *tau = sycl::malloc_device<cd>(n, q);
  int64_t *piv = sycl::malloc_device<int64_t>(n, q);
  q.memcpy(A0, h.data(), n * n * sizeof(cd)).wait();
  int64_t ws = std::max({lp::geqrf_scratchpad_size<cd>(q, n, n, n),
                         lp::unmqr_scratchpad_size<cd>(q, oneapi::mkl::side::left, oneapi::mkl::transpose::nontrans, n, n, n, n, n),
                         lp::getrf_scratchpad_size<cd>(q, n, n, n),
                         lp::getrs_scratchpad_size<cd>(q, oneapi::mkl::transpose::nontrans, n, n, n, n)});
  cd *w = sycl::malloc_device<cd>(ws, q);
  auto reset = [&] { q.memcpy(A, A0, n * n * sizeof(cd)); q.wait(); };
  auto bench = [&](const char *name, auto f, bool rs) {
    double t = 0, tmin = 1e9;
    for (int r = 0; r < reps + 2; ++r) {
      if (rs) reset();
      double t0 = now();
      f();
      double dt = now() - t0;
      if (r >= 2) { t += dt; tmin = std::min(tmin, dt); }
    }
    std::printf("RESULT N=%lld %-22s mean %8.3f ms  min %8.3f ms\n", (long long)n, name, 1e3 * t / reps, 1e3 * tmin);
  };
  bench("dev zgeqrf", [&] { lp::geqrf(q, n, n, A, n, tau, w, ws); q.wait(); }, true);
  reset(); lp::geqrf(q, n, n, A, n, tau, w, ws).wait();
  bench("dev zunmqr L N", [&] {
    lp::unmqr(q, oneapi::mkl::side::left, oneapi::mkl::transpose::nontrans, n, n, n, A, n, tau, C, n, w, ws); q.wait(); }, false);
  bench("dev zgetrf", [&] { lp::getrf(q, n, n, A, n, piv, w, ws); q.wait(); }, true);
  reset(); lp::getrf(q, n, n, A, n, piv, w, ws).wait();
  bench("dev zgetrs", [&] { lp::getrs(q, oneapi::mkl::transpose::nontrans, n, n, A, n, piv, C, n, w, ws); q.wait(); }, false);
  // batched variants (strided, batch 1 and 2) and their deviation from the plain geqrf/getrf result
  {
    std::vector<cd> r1(n * n), r2(n * n);
    reset(); lp::geqrf(q, n, n, A, n, tau, w, ws); q.wait(); q.memcpy(r1.data(), A, n * n * sizeof(cd)).wait();
    for (int64_t bt : {1, 2}) {
      cd *Ab = sycl::malloc_device<cd>(bt * n * n, q); cd *tb = sycl::malloc_device<cd>(bt * n, q);
      int64_t *pb = sycl::malloc_device<int64_t>(bt * n, q);
      int64_t wq = std::max(lp::geqrf_batch_scratchpad_size<cd>(q, n, n, n, n * n, n, bt),
                            lp::getrf_batch_scratchpad_size<cd>(q, n, n, n, n * n, n, bt));
      cd *wb = sycl::malloc_device<cd>(wq, q);
      auto rs = [&] { for (int64_t i = 0; i < bt; ++i) q.memcpy(Ab + i * n * n, A0, n * n * sizeof(cd)); q.wait(); };
      char nm[64];
      std::snprintf(nm, 64, "dev zgeqrf_batch b=%lld", (long long)bt);
      double t = 0; for (int r = 0; r < reps + 2; ++r) { rs(); double t0 = now(); lp::geqrf_batch(q, n, n, Ab, n, n * n, tb, n, bt, wb, wq); q.wait(); if (r >= 2) t += now() - t0; }
      q.memcpy(r2.data(), Ab, n * n * sizeof(cd)).wait();
      double md = 0, mx = 0; for (int64_t i = 0; i < n * n; ++i) { md = std::max(md, std::abs(r1[i] - r2[i])); mx = std::max(mx, std::abs(r1[i])); }
      std::printf("RESULT N=%lld %-22s mean %8.3f ms  (max rel dev vs geqrf %.2e)\n", (long long)n, nm, 1e3 * t / reps, md / mx);
      std::snprintf(nm, 64, "dev zgetrf_batch b=%lld", (long long)bt);
      t = 0; for (int r = 0; r < reps + 2; ++r) { rs(); double t0 = now(); lp::getrf_batch(q, n, n, Ab, n, n * n, pb, n, bt, wb, wq); q.wait(); if (r >= 2) t += now() - t0; }
      std::printf("RESULT N=%lld %-22s mean %8.3f ms\n", (long long)n, nm, 1e3 * t / reps);
      sycl::free(Ab, q); sycl::free(tb, q); sycl::free(pb, q); sycl::free(wb, q);
    }
  }
  cd one(1, 0), zero(0, 0);
  bench("dev zgemm NN", [&] {
    oneapi::mkl::blas::column_major::gemm(q, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans, n, n, n, one,
                                          A0, n, A, n, zero, C, n); q.wait(); }, false);
  // host LAPACK zgeqrf (MKL, current thread count)
  std::vector<cd> hA(n * n), htau(n);
  MKL_Complex16 wq; lapack_int info, lw = -1;
  zgeqrf((lapack_int *)&n, (lapack_int *)&n, (MKL_Complex16 *)hA.data(), (lapack_int *)&n, (MKL_Complex16 *)htau.data(), &wq, &lw, &info);
  lw = (lapack_int)wq.real; std::vector<cd> hw(lw);
  bench("host zgeqrf", [&] {
    hA = h; zgeqrf((lapack_int *)&n, (lapack_int *)&n, (MKL_Complex16 *)hA.data(), (lapack_int *)&n, (MKL_Complex16 *)htau.data(),
                   (MKL_Complex16 *)hw.data(), &lw, &info); }, false);
  sycl::free(A, q); sycl::free(A0, q); sycl::free(C, q); sycl::free(tau, q); sycl::free(piv, q); sycl::free(w, q);
}

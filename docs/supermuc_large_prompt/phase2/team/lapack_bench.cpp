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

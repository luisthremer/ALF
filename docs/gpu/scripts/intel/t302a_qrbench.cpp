// development only: what does one tiny oneMKL zgeqrf cost on the B580?  icpx -fsycl -O2 -std=c++17 -I$MKLROOT/include ... -lmkl_sycl_lapack -lmkl_sycl_blas
#include <sycl/sycl.hpp>
#include <oneapi/mkl/lapack.hpp>
#include <chrono>
#include <complex>
#include <cstdio>
#include <vector>
extern "C" void zgeqrf_(int *, int *, void *, int *, void *, void *, int *, int *);
using cd = std::complex<double>;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main() {
  setvbuf(stdout, nullptr, _IOLBF, 0);
  
  for (int mode = 0; mode < 2; ++mode) {
    sycl::queue q = mode == 0 ? sycl::queue{sycl::gpu_selector_v, sycl::property::queue::in_order()}
                              : sycl::queue{sycl::gpu_selector_v};
    for (int n : {2, 8, 32, 64, 144, 400}) {
      const int NR = n >= 144 ? 50 : 300;
      std::vector<cd> h(n * n);
      for (int i = 0; i < n * n; ++i) h[i] = cd(0.1 * ((i * 7919) % 101) - 5, 0.05 * ((i * 104729) % 97));
      cd *a = sycl::malloc_device<cd>(n * n, q), *tau = sycl::malloc_device<cd>(n, q);
      int64_t lw = oneapi::mkl::lapack::geqrf_scratchpad_size<cd>(q, n, n, n);
      cd *w = sycl::malloc_device<cd>(lw, q);
      q.memcpy(a, h.data(), n * n * sizeof(cd)).wait();
      oneapi::mkl::lapack::geqrf(q, n, n, a, n, tau, w, lw).wait();
      double t0 = now();
      for (int i = 0; i < NR; ++i) oneapi::mkl::lapack::geqrf(q, n, n, a, n, tau, w, lw).wait_and_throw();
      double t1 = now();
      for (int i = 0; i < NR; ++i) oneapi::mkl::lapack::geqrf(q, n, n, a, n, tau, w, lw);
      q.wait();
      double t2 = now();
      {  // host fallback: copy to the host, LAPACK zgeqrf on the host, copy back
        std::vector<cd> hh(n * n), ht(n), hw(64 * n + 64);
        int info, lwh = 64 * n + 64, nn = n;
        double t3 = now();
        for (int i = 0; i < NR; ++i) {
          q.memcpy(hh.data(), a, (size_t)n * n * sizeof(cd)).wait();
          zgeqrf_(&nn, &nn, hh.data(), &nn, ht.data(), hw.data(), &lwh, &info);
          q.memcpy(a, hh.data(), (size_t)n * n * sizeof(cd));
          q.memcpy(tau, ht.data(), (size_t)n * sizeof(cd)).wait();
        }
        double t4 = now();
        std::printf("   host fallback n=%d: %.1f us per call\n", n, (t4 - t3) / NR * 1e6);
      }
      std::printf("%s n=%d: wait each %.1f us, no wait %.1f us per call (lwork %lld)\n", mode ? "out-of-order" : "in-order", n,
                  (t1 - t0) / NR * 1e6, (t2 - t1) / NR * 1e6, (long long)lw);
      sycl::free(a, q); sycl::free(tau, q); sycl::free(w, q);
    }
  }
}

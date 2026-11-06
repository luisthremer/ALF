// B580: zgemm and zgeqrf throughput with oneMKL SYCL (device) for N = 144, 400, 1024
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <oneapi/mkl/lapack.hpp>
#include <chrono>
#include <complex>
#include <cstdio>
#include <vector>
using cd = std::complex<double>;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main() {
  sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order()};
  printf("device: %s\n", q.get_device().get_info<sycl::info::device::name>().c_str());
  for (int n : {144, 400, 1024}) {
    std::vector<cd> h((size_t)n * n);
    for (size_t i = 0; i < h.size(); ++i) h[i] = cd(0.001 * ((i * 7919) % 101) - 0.05, 0.0005 * ((i * 104729) % 97));
    cd *a = sycl::malloc_device<cd>((size_t)n * n, q), *b = sycl::malloc_device<cd>((size_t)n * n, q), *c = sycl::malloc_device<cd>((size_t)n * n, q), *tau = sycl::malloc_device<cd>(n, q);
    q.memcpy(a, h.data(), h.size() * sizeof(cd)); q.memcpy(b, h.data(), h.size() * sizeof(cd)); q.wait();
    const int NR = n >= 1024 ? 5 : 30;
    oneapi::mkl::blas::column_major::gemm(q, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans, n, n, n, cd(1, 0), a, n, b, n, cd(0, 0), c, n).wait();
    double t0 = now();
    for (int i = 0; i < NR; ++i) oneapi::mkl::blas::column_major::gemm(q, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::nontrans, n, n, n, cd(1, 0), a, n, b, n, cd(0, 0), c, n);
    q.wait(); double t1 = now();
    double tg = (t1 - t0) / NR;
    printf("N=%4d zgemm %.3f ms  %.1f GFlop/s (8 N^3 flops)\n", n, tg * 1e3, 8.0 * n * n * n / tg / 1e9);
    int64_t lw = oneapi::mkl::lapack::geqrf_scratchpad_size<cd>(q, n, n, n);
    cd *w = sycl::malloc_device<cd>(lw, q);
    q.memcpy(a, h.data(), h.size() * sizeof(cd)).wait();
    oneapi::mkl::lapack::geqrf(q, n, n, a, n, tau, w, lw).wait();
    t0 = now();
    for (int i = 0; i < NR; ++i) { q.memcpy(a, h.data(), h.size() * sizeof(cd)); oneapi::mkl::lapack::geqrf(q, n, n, a, n, tau, w, lw); }
    q.wait(); t1 = now();
    double tq = (t1 - t0) / NR;
    printf("N=%4d zgeqrf %.3f ms (incl. H->D copy of the input)\n", n, tq * 1e3);
    sycl::free(a, q); sycl::free(b, q); sycl::free(c, q); sycl::free(tau, q); sycl::free(w, q);
  }
}

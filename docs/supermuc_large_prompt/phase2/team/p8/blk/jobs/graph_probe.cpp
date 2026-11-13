// probe: can a block-like op sequence (own kernels + oneMKL zgemm, USM) be recorded into a SYCL command graph and replayed? (in-order queue)
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <complex>
#include <chrono>
#include <cstdio>
namespace sx = sycl::ext::oneapi::experimental;
using cd = std::complex<double>;
int main() {
  sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order{}};
  std::printf("device: %s\n", q.get_device().get_info<sycl::info::device::name>().c_str());
  std::printf("graph support: %d\n", (int)q.get_device().has(sycl::aspect::ext_oneapi_graph));
  const int n = 576, k = 12;
  cd *A = sycl::malloc_device<cd>(n * k, q), *B = sycl::malloc_device<cd>(n * k, q), *C = sycl::malloc_device<cd>((size_t)n * n, q);
  q.fill(A, cd(1, 0), n * k).wait(); q.fill(B, cd(1, 0), n * k).wait(); q.fill(C, cd(0, 0), (size_t)n * n).wait();
  auto seq = [&](sycl::queue &qq) {
    qq.parallel_for(sycl::range<1>(n * k), [=](sycl::id<1> i) { A[i] = A[i] * 1.0; });
    oneapi::mkl::blas::column_major::gemm(qq, oneapi::mkl::transpose::nontrans, oneapi::mkl::transpose::trans, n, n, k, cd(-1, 0), A, n, B, n, cd(1, 0), C, n);
  };
  auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < 1000; ++i) seq(q);
  q.wait();
  double tdir = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("direct: 1000 x (kernel + zgemm) %.3f ms total; host submit+run %.2f us per sequence\n", 1e3 * tdir, 1e6 * tdir / 1000);
  try {
    sx::command_graph<sx::graph_state::modifiable> g(q.get_context(), q.get_device());
    g.begin_recording(q);
    seq(q);
    g.end_recording(q);
    auto ge = g.finalize();
    std::printf("graph recorded and finalized\n");
    auto t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000; ++i) q.ext_oneapi_graph(ge);
    auto t2 = std::chrono::steady_clock::now();
    q.wait();
    auto t3 = std::chrono::steady_clock::now();
    std::printf("graph: 1000 replays: submit %.2f us each, total incl. wait %.2f us each\n", 1e6 * std::chrono::duration<double>(t2 - t1).count() / 1000,
                1e6 * std::chrono::duration<double>(t3 - t1).count() / 1000);
  } catch (sycl::exception &e) {
    std::printf("graph: SYCL exception: %s (code %d)\n", e.what(), e.code().value());
  } catch (std::exception &e) {
    std::printf("graph: exception: %s\n", e.what());
  }
  return 0;
}

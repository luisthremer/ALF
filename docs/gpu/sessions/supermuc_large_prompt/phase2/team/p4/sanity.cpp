// harness sanity: upload -> kernel(s) -> readback; counts wrong results. Variants of copy path / waiting.
#include <sycl/sycl.hpp>
#include <chrono>
#include <thread>
#include <cstdio>
#include <vector>
#include <cstdlib>
int main(int argc, char **argv) {
  sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order()};
  const size_t n = 576 * 128 * 2;
  double *d = sycl::malloc_device<double>(n, q);
  const char *names[] = {"memcpy up/down", "memcpy + sleep 20ms before check", "kernel-copy readback", "kernel up + kernel down", "memcpy, 2nd readback", "shared USM no copy"};
  for (int mode = 0; mode < 6; ++mode) {
    int bad = 0, badruns = 0;
    for (int rep = 0; rep < 100; ++rep) {
      double *hp = sycl::malloc_host<double>(n, q), *op = sycl::malloc_host<double>(n, q);
      double *sh = mode == 5 ? sycl::malloc_shared<double>(n, q) : nullptr;
      for (size_t i = 0; i < n; ++i) { hp[i] = 1.0 + i; op[i] = 0.0; }
      if (mode == 5) { for (size_t i = 0; i < n; ++i) sh[i] = 1.0 + i; q.parallel_for(n, [=](sycl::id<1> i) { sh[i] += 1.0; }).wait(); for (size_t i = 0; i < n; ++i) op[i] = sh[i]; }
      else {
        if (mode == 3) q.parallel_for(n, [=](sycl::id<1> i) { d[i] = hp[i]; }).wait(); else q.memcpy(d, hp, n * 8).wait();
        q.parallel_for(n, [=](sycl::id<1> i) { d[i] = d[i] + 1.0; }).wait();
        if (mode == 2 || mode == 3) q.parallel_for(n, [=](sycl::id<1> i) { op[i] = d[i]; }).wait(); else q.memcpy(op, d, n * 8).wait();
        if (mode == 1) std::this_thread::sleep_for(std::chrono::milliseconds(20));
        if (mode == 4) { q.memcpy(op, d, n * 8).wait(); }
      }
      int b = 0;
      for (size_t i = 0; i < n; ++i) if (op[i] != 2.0 + i) ++b;
      if (b) { ++badruns; bad += b; }
      sycl::free(hp, q); sycl::free(op, q); if (sh) sycl::free(sh, q);
    }
    std::printf("mode %d (%s): bad runs %d / 100, bad elements %d\n", mode, names[mode], badruns, bad);
  }
}

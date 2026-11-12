// microbenchmark: cost of a work-group barrier and of a flag handoff on PVC (true time via host-USM flag)
#include <sycl/sycl.hpp>
#include <sycl/ext/intel/experimental/grf_size_properties.hpp>
#include <chrono>
#include <cstdio>
using namespace sycl;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
template <class F> static double run(queue &q, F &&f) {
  int *fl = malloc_host<int>(1, q); *fl = 0; volatile int *vf = fl;
  double t0 = now(); f(); q.single_task([=] { *fl = 1; }); while (*vf == 0) {} double t = now() - t0; free(fl, q); return t;
}
int main() {
  queue q{gpu_selector_v, property::queue::in_order()};
  const int N = 2000;
  for (int wi : {64, 256, 512}) {
    for (int mode = 0; mode < 5; ++mode) {
      auto launch = [&] {
        q.submit([&](handler &h) {
          local_accessor<int, 1> fl(range<1>(64), h);
          local_accessor<double, 1> dat(range<1>(64), h);
          ext::oneapi::experimental::properties props{ext::intel::experimental::grf_size<256>};
          h.parallel_for(nd_range<1>(wi, wi), props, [=](nd_item<1> it) [[sycl::reqd_sub_group_size(16)]] {
            auto sg = it.get_sub_group();
            int sid = sg.get_group_id()[0], nsg = wi / 16;
            if (it.get_local_linear_id() < 64) fl[it.get_local_linear_id()] = 0;
            it.barrier(access::fence_space::local_space);
            double acc = 0;
            for (int i = 0; i < N; ++i) {
              if (mode == 0) it.barrier(access::fence_space::local_space);
              else if (mode == 1) group_barrier(it.get_group());
              else if (mode == 2) { sycl::group_barrier(sg); }
              else if (mode == 3) {   // relay: sub-group (i % nsg) writes a flag, all others wait for it (spin on SLM)
                atomic_ref<int, memory_order::relaxed, memory_scope::work_group, access::address_space::local_space> f(fl[0]);
                if (sid == i % nsg) { if (sg.get_local_id()[0] == 0) { dat[0] = i; atomic_fence(memory_order::release, memory_scope::work_group); f.store(i + 1); } }
                else { if (sg.get_local_id()[0] == 0) { while (f.load() < i + 1) {} atomic_fence(memory_order::acquire, memory_scope::work_group); acc += dat[0]; } }
                // consumers do not signal back (no WAR protection): measures the handoff latency only
              } else if (mode == 4) {   // sub-group reduction chain (for reference)
                acc = reduce_over_group(sg, acc + i, plus<double>());
              }
            }
            if (acc == 12345.678) fl[1] = 1;
          });
        });
      };
      run(q, launch);
      double t = run(q, launch);
      const char *nm[] = {"it.barrier(local)", "group_barrier(wg) [local+global fence]", "sub-group barrier", "SLM flag relay (1-way handoff)", "reduce_over_group(double)"};
      std::printf("WI=%3d (%2d sub-groups) %-42s %8.3f us per iteration\n", wi, wi / 16, nm[mode], 1e6 * t / N);
    }
  }
}

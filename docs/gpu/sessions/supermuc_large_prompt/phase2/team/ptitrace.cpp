// Team session: minimal GPU kernel/copy tracer on Intel PTI (pti_view), LD_PRELOADed into ALF (VTune collects no GPU
// data on this system). At exit prints: number of kernels, summed device time, top kernels by device time, copies.
// Build: icpx -O2 -shared -fPIC ptitrace.cpp -I$PTI/include -L$PTI/lib -lpti_view -o libptitrace.so
// Run:   ZE_ENABLE_TRACING_LAYER=1 LD_PRELOAD=libptitrace.so:$PTI/lib/libpti_view.so <app>
#include <pti/pti_view.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace {
struct Acc { long long n = 0; double dev = 0, lat = 0; };
std::map<std::string, Acc> g_k;
Acc g_copy;
std::mutex g_m;
uint64_t g_first = ~0ull, g_last = 0;

void req(unsigned char **buf, size_t *size) {
  *size = 64 * 1024 * 1024;
  *buf = static_cast<unsigned char *>(std::aligned_alloc(8, *size));
}
void done(unsigned char *buf, size_t, size_t used) {
  pti_view_record_base *r = nullptr;
  std::lock_guard<std::mutex> l(g_m);
  while (ptiViewGetNextRecord(buf, used, &r) == PTI_SUCCESS) {
    if (r->_view_kind == PTI_VIEW_DEVICE_GPU_KERNEL) {
      auto *k = reinterpret_cast<pti_view_record_kernel *>(r);
      auto &a = g_k[k->_name ? k->_name : "?"];
      a.n++; a.dev += 1e-9 * (k->_end_timestamp - k->_start_timestamp);
      a.lat += 1e-9 * (k->_end_timestamp - k->_append_timestamp);   // append (host) -> device end
      g_first = std::min(g_first, (uint64_t)k->_append_timestamp); g_last = std::max(g_last, (uint64_t)k->_end_timestamp);
    } else if (r->_view_kind == PTI_VIEW_DEVICE_GPU_MEM_COPY) {
      auto *c = reinterpret_cast<pti_view_record_memory_copy *>(r);
      g_copy.n++; g_copy.dev += 1e-9 * (c->_end_timestamp - c->_start_timestamp);
    }
  }
  std::free(buf);
}
struct Init {
  Init() {
    ptiViewSetCallbacks(req, done);
    ptiViewEnable(PTI_VIEW_DEVICE_GPU_KERNEL);
    ptiViewEnable(PTI_VIEW_DEVICE_GPU_MEM_COPY);
  }
  ~Init() {
    ptiViewDisable(PTI_VIEW_DEVICE_GPU_KERNEL);
    ptiViewDisable(PTI_VIEW_DEVICE_GPU_MEM_COPY);
    ptiFlushAllViews();
    std::vector<std::pair<std::string, Acc>> v(g_k.begin(), g_k.end());
    std::sort(v.begin(), v.end(), [](auto &a, auto &b) { return a.second.dev > b.second.dev; });
    long long n = 0; double dev = 0, lat = 0;
    for (auto &e : v) { n += e.second.n; dev += e.second.dev; lat += e.second.lat; }
    std::fprintf(stderr, "PTI TOTAL kernels %lld  device %.3f s  append->end %.3f s  window %.3f s  copies %lld (%.3f s)\n", n, dev,
                 lat, g_last > g_first ? 1e-9 * (g_last - g_first) : 0.0, g_copy.n, g_copy.dev);
    int i = 0;
    for (auto &e : v) {
      if (i++ >= 40) break;
      std::fprintf(stderr, "PTI K %8lld  dev %8.3f s  mean %8.2f us  %.150s\n", e.second.n, e.second.dev,
                   1e6 * e.second.dev / e.second.n, e.first.c_str());
    }
  }
} g_init;
}  // namespace
